!===============================================================================
! ETAPE 6c : factorisation LU par blocs DISTRIBUEE (bloc-cyclique 2D),
! SANS PIVOT -- premiere brique de la factorisation elle-meme, volontairement
! simplifiee (matrices a diagonale dominante uniquement dans les tests) pour
! valider la MECANIQUE de communication (diffusion du panneau L le long
! d'une ligne de grille, du panneau U le long d'une colonne, mise a jour
! locale du complement de Schur) avant d'ajouter le pivot partiel distribue
! (etape 6d, la partie la plus delicate -- recherche de pivot inter-
! processus dans une colonne de processus).
!
! ALGORITHME (analogue non pivote de PxGETRF, schema "right-looking") --
! pour chaque bloc-colonne/ligne global K (KSTART:KEND, largeur KW<=NB) :
!
!   1) FACTORISATION DIAGONALE : le processus PROPRIETAIRE du bloc
!      diagonal A_kk (le seul, en bloc-cyclique un bloc appartient a UN
!      processus) le factorise localement, SANS PIVOT (LU_NOPIV_UNBLOCKED,
!      routine maison -- LAPACK n'offre pas de GETRF sans pivot, donc on
!      l'ecrit a la main pour ce petit bloc KWxKW, cf. commentaire de la
!      routine plus bas).
!
!   2) DIFFUSION DE A_kk FACTORISE :
!        - le long de la ligne de grille (ROW_COMM, racine = colonne du
!          proprietaire) -- vers les processus qui vont calculer le
!          panneau U (meme ligne de grille que le bloc diagonal)
!        - le long de la colonne de grille (COL_COMM, racine = ligne du
!          proprietaire) -- vers les processus qui vont calculer le
!          panneau L (meme colonne de grille)
!
!   3) PANNEAU U (processus de la ligne PROW_K) : U_kj = L_kk^-1 * A_kj
!      par CTRSM (triangle inferieur unitaire), pour toutes les colonnes
!      locales "de queue" (indice global > KEND).
!
!   4) PANNEAU L (processus de la colonne PCOL_K) : L_ik = A_ik * U_kk^-1
!      par CTRSM (triangle superieur), pour toutes les lignes locales de
!      queue.
!
!   5) DIFFUSION DES PANNEAUX vers le reste de la grille : U (calcule par
!      la ligne PROW_K) diffuse le long de COL_COMM (racine=PROW_K) pour
!      atteindre toutes les lignes de la grille ; L (calcule par la
!      colonne PCOL_K) diffuse le long de ROW_COMM (racine=PCOL_K) pour
!      atteindre toutes les colonnes.
!
!   6) MISE A JOUR LOCALE DU COMPLEMENT DE SCHUR : CHAQUE processus (partie
!      de queue non vide) calcule A_ij -= L_ik * U_kj (CGEMM) avec SA part
!      locale des panneaux L et U recus.
!
! Repartition "de queue" (indices globaux > KEND) SANS BOUCLE : grace a la
! propriete de prefixe du bloc-cyclique (le decoupage local des KEND
! premiers elements d'une matrice de taille N est identique a celui d'une
! matrice de taille KEND), le nombre d'elements locaux de tete est donne
! directement par NUMROC(KEND,NB,MYROW/MYCOL,PR/PC) -- pas de recherche a
! faire, cf. RTRAIL0/CTRAIL0 ci-dessous.
!===============================================================================
MODULE MPI_LU_NOPIV_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  IMPLICIT NONE

CONTAINS

  SUBROUTINE MPI_LU_NOPIV_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, INFO)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)     :: GRID
    INTEGER, INTENT(IN)          :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(INOUT)   :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    INTEGER, INTENT(OUT)         :: INFO

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), UBUF(:,:), LBUF(:,:)
    INTEGER :: NBLK, K, KSTART, KEND, KW
    INTEGER :: PROW_K, PCOL_K, OWNER_RANK, RS, CS
    INTEGER :: RTRAIL0, CTRAIL0, NROWTRAIL, NCOLTRAIL
    INTEGER :: INFO_STEP, IERR

    INFO = 0
    NBLK = (N + NB - 1) / NB

    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1

      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)
      OWNER_RANK = PROW_K*GRID%PC + PCOL_K

      RTRAIL0 = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      CTRAIL0 = NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      NROWTRAIL = MLOC - RTRAIL0
      NCOLTRAIL = NLOC - CTRAIL0

      ! --- 1) Factorisation locale du bloc diagonal (proprietaire seul) ---
      INFO_STEP = 0
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        CS = INDXG2L(KSTART, NB, GRID%PC)
        CALL LU_NOPIV_UNBLOCKED(KW, A_LOCAL(RS,CS), MLOC, INFO_STEP)
        IF (INFO_STEP .NE. 0) INFO_STEP = INFO_STEP + KSTART - 1
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        RETURN
      END IF

      ! --- 2) Diffusion du bloc diagonal factorise ---
      ALLOCATE(DKK(KW,KW))
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        CS = INDXG2L(KSTART, NB, GRID%PC)
        DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
      END IF
      IF (GRID%MYROW .EQ. PROW_K) THEN
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      END IF
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      END IF

      ! --- 3) Panneau U : U_kj = L_kk^-1 * A_kj (ligne PROW_K, colonnes de queue) ---
      IF (GRID%MYROW.EQ.PROW_K .AND. NCOLTRAIL.GT.0) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NCOLTRAIL, &
                    CONE, DKK, KW, A_LOCAL(RS,CTRAIL0+1), MLOC)
      END IF

      ! --- 4) Panneau L : L_ik = A_ik * U_kk^-1 (colonne PCOL_K, lignes de queue) ---
      IF (GRID%MYCOL.EQ.PCOL_K .AND. NROWTRAIL.GT.0) THEN
        CS = INDXG2L(KSTART, NB, GRID%PC)
        CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', NROWTRAIL, KW, &
                    CONE, DKK, KW, A_LOCAL(RTRAIL0+1,CS), MLOC)
      END IF
      DEALLOCATE(DKK)

      ! --- 5) Diffusion des panneaux vers le reste de la grille ---
      IF (NCOLTRAIL .GT. 0) THEN
        ALLOCATE(UBUF(KW,NCOLTRAIL))
        IF (GRID%MYROW .EQ. PROW_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          UBUF = A_LOCAL(RS:RS+KW-1, CTRAIL0+1:NLOC)
        END IF
        CALL MPI_BCAST(UBUF, KW*NCOLTRAIL, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      END IF

      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(LBUF(NROWTRAIL,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          LBUF = A_LOCAL(RTRAIL0+1:MLOC, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(LBUF, NROWTRAIL*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      END IF

      ! --- 6) Mise a jour locale du complement de Schur ---
      IF (NROWTRAIL.GT.0 .AND. NCOLTRAIL.GT.0) THEN
        CALL CGEMM('No transpose', 'No transpose', NROWTRAIL, NCOLTRAIL, KW, &
                    -CONE, LBUF, NROWTRAIL, UBUF, KW, &
                    CONE, A_LOCAL(RTRAIL0+1,CTRAIL0+1), MLOC)
      END IF

      IF (NCOLTRAIL .GT. 0) DEALLOCATE(UBUF)
      IF (NROWTRAIL .GT. 0) DEALLOCATE(LBUF)
    END DO
  END SUBROUTINE MPI_LU_NOPIV_FACTOR

  !-----------------------------------------------------------------------
  ! LU SANS PIVOT, non blocuee, sur un petit bloc dense local KBxKB
  ! (LAPACK n'a pas de "GETRF sans pivot" -- routine maison, Doolittle
  ! classique par rang-1 : A(J,J) devient le pivot U(J,J), la colonne
  ! sous la diagonale est divisee par lui (= colonne de L, diagonale
  ! unitaire implicite), puis mise a jour rang-1 du reste du bloc). En
  ! sortie A contient L (sous la diagonale, unitaire implicite) et U (sur
  ! et au-dessus), exactement le format d'empaquetage de CGETRF.
  !-----------------------------------------------------------------------
  SUBROUTINE LU_NOPIV_UNBLOCKED(KB, A, LDA, INFO)
    INTEGER, INTENT(IN)          :: KB, LDA
    COMPLEX(WP), INTENT(INOUT)   :: A(LDA,KB)
    INTEGER, INTENT(OUT)         :: INFO

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    REAL(WP), PARAMETER :: PIVTOL = 1.0E-6_WP
    INTEGER :: J

    INFO = 0
    DO J = 1, KB
      IF (ABS(A(J,J)) .LT. PIVTOL) THEN
        INFO = J
        RETURN
      END IF
      IF (J .LT. KB) THEN
        CALL CSCAL(KB-J, CONE/A(J,J), A(J+1,J), 1)
        CALL CGERU(KB-J, KB-J, -CONE, A(J+1,J), 1, A(J,J+1), LDA, A(J+1,J+1), LDA)
      END IF
    END DO
  END SUBROUTINE LU_NOPIV_UNBLOCKED

END MODULE MPI_LU_NOPIV_MOD
