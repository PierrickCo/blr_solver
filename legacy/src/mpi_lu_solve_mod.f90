!===============================================================================
! RESOLUTION MPI DISTRIBUEE : ferme le manque mis en evidence par le
! benchmark N=10000 (bench_mpi.f90) -- resoudre sans jamais RASSEMBLER
! (GATHER_2D) la matrice factorisee complete sur un seul processus (pic
! memoire O(N^2) mesure a ~980 Mo dans le benchmark).
!
! SECOND MEMBRE B : DISTRIBUE PAR LIGNES, comme A (meme decoupage NB, PR --
! B_LOCAL(MLOC,NRHS) avec MLOC=NUMROC(N,NB,MYROW,PR), EXACTEMENT le meme
! que celui de A_LOCAL), REPLIQUE seulement AU SEIN d'une meme ligne de
! grille (meme processus MYROW, tous les MYCOL) -- PAS replique sur toute
! la grille. C'est EXACTEMENT la convention deja validee par MPI_MATVEC_2D
! (etape 6b, Y_LOCAL) et GATHER_ROW_VECTOR (mpi_grid_mod.f90) -- reprise
! ici telle quelle plutot que reinventee.
!
! POURQUOI CE CHOIX (evolution depuis la premiere version, qui repliquait
! B INTEGRALEMENT sur toute la grille) : l'utilisateur a signale que NRHS
! peut devenir tres grand -- une replication totale (N x NRHS sur CHAQUE
! processus) ne passe pas a l'echelle. La distribuer PAR LIGNES divise la
! memoire de B par PR (nombre de lignes de la grille) -- pas par P=PR*PC
! (une vraie distribution 2D exigerait de repositionner B "en lignes" vers
! "en colonnes" a chaque etape, la vraie complexite d'un PxGETRS a la
! ScaLAPACK qu'on continue a eviter), mais un facteur PR est deja un vrai
! gain qui grandit avec la taille de la grille, et surtout : cette
! distribution reutilise EXACTEMENT les memes diffusions ROW_COMM/COL_COMM
! deja eprouvees pendant la factorisation elle-meme (voir algorithme
! ci-dessous) -- bien plus sur qu'inventer un nouveau motif de
! communication pour une distribution 2D complete du second membre.
!
! ALGORITHME (substitution par blocs, MEME schema de communication que la
! factorisation -- litteralement les memes appels ROW_COMM/COL_COMM que
! mpi_lu_pivot_mod.f90/mpi_blr_lu_pivot_mod.f90, juste appliques a B) :
!
!   0) Appliquer P*B : B etant distribue PAR LIGNES (pas replique), un
!      echange de lignes GLOBALES peut desormais impliquer DEUX PROCESSUS-
!      LIGNES differents -- reutilise SWAP_ROWS_RANGE (mpi_lu_pivot_mod.f90,
!      deja valide pour A) directement sur B_LOCAL, en traitant les NRHS
!      colonnes comme "toute la plage de colonnes locale" (B n'est pas
!      distribue par colonnes).
!
!   1) L*Y = P*B (substitution AVANT, K=1..NBLK) :
!      a) diffuser le bloc diagonal DKK le long de ROW_COMM (racine=
!         PCOL_K) -- atteint la ligne PROW_K -- EXACTEMENT le meme appel
!         que pour le calcul du panneau U pendant la factorisation.
!      b) les processus de la ligne PROW_K resolvent (CTRSM) sur LEUR
!         copie (repliquee au sein de la ligne) de B(K,:) -- redondant PC
!         fois mais bon marche, garde la replication INTRA-LIGNE coherente.
!      c) (PROW_K,PCOL_K), qui vient de calculer X(K), le diffuse le long
!         de COL_COMM (racine=PROW_K) -- atteint tous les processus de la
!         colonne PCOL_K, EXACTEMENT comme la diffusion du panneau U
!         pendant la factorisation.
!      d) les processus de la colonne PCOL_K (qui ont maintenant X(K))
!         calculent leur contribution locale A(lignes trainantes locales,
!         K)*X(K) (CGEMM), la diffusent le long de LEUR PROPRE ROW_COMM
!         (racine=PCOL_K) pour la repliquer au sein de leur ligne --
!         EXACTEMENT la diffusion de LBUF pendant la factorisation.
!      e) chaque processus met a jour SA copie de B aux lignes GLOBALES
!         trainantes correspondantes -- local.
!
!   2) U*X = Y (substitution ARRIERE, K=NBLK..1) : meme schema, lignes "en
!      tete" (indice global < KSTART) au lieu de "trainantes".
!
! A LA FIN, B_LOCAL contient X, distribue PAR LIGNES et REPLIQUE au sein de
! chaque ligne de grille -- forme directement reutilisable par
! GATHER_ROW_VECTOR (mpi_grid_mod.f90) si une forme rassemblee est desiree.
!===============================================================================
MODULE MPI_LU_SOLVE_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE MPI_LU_PIVOT_MOD, ONLY: SWAP_ROWS_RANGE
  IMPLICIT NONE

CONTAINS

  SUBROUTINE MPI_LU_SOLVE_DISTRIBUTED(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                                       NRHS, B_LOCAL, LDB)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)     :: GRID
    INTEGER, INTENT(IN)          :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(IN)      :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    INTEGER, INTENT(IN)          :: IPIV_GLOBAL(N)
    INTEGER, INTENT(IN)          :: NRHS, LDB
    ! B_LOCAL : distribue PAR LIGNES comme A (MLOC lignes locales), replique
    ! au sein d'une meme ligne de grille (meme MYROW, tous les MYCOL).
    COMPLEX(WP), INTENT(INOUT)   :: B_LOCAL(MAX(MLOC,1),NRHS)

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK(:,:), CONTRIB(:,:)
    INTEGER :: NBLK, K, KSTART, KEND, KW
    INTEGER :: PROW_K, PCOL_K, RS, CS
    INTEGER :: RTRAIL0, NROWTRAIL, RHEAD, I, R1, R2, PROW2
    INTEGER :: IERR

    NBLK = (N + NB - 1) / NB

    ! === 0) P*B : echange de lignes, potentiellement inter-processus-ligne
    !     (SWAP_ROWS_RANGE, deja valide pour A) ===
    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      DO I = 1, KW
        R1 = KSTART - 1 + I
        R2 = KSTART - 1 + IPIV_GLOBAL(KSTART-1+I)
        IF (R1 .NE. R2) THEN
          PROW2 = INDXG2P(R2, NB, GRID%PR)
          CALL SWAP_ROWS_RANGE(GRID, NB, B_LOCAL, MLOC, NRHS, R1, R2, PROW_K, PROW2, 1, NRHS)
        END IF
      END DO
    END DO

    ! === 1) L*Y = P*B (substitution avant) ===
    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = INDXG2L(KSTART, NB, GRID%PR)
        CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS, &
                    CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        DEALLOCATE(DKK)
      END IF

      ! X(K) : diffuse le long de COL_COMM (racine=PROW_K) pour atteindre
      ! la colonne de grille PCOL_K (qui seule a besoin de A(:,K)).
      ALLOCATE(XK(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        XK = B_LOCAL(RS:RS+KW-1,1:NRHS)
      END IF
      CALL MPI_BCAST(XK, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      RTRAIL0   = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      NROWTRAIL = MLOC - RTRAIL0
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          CALL CGEMM('No transpose', 'No transpose', NROWTRAIL, NRHS, KW, &
                      CONE, A_LOCAL(RTRAIL0+1,CS), MLOC, XK, KW, &
                      (0.0_WP,0.0_WP), CONTRIB, NROWTRAIL)
        END IF
        CALL MPI_BCAST(CONTRIB, NROWTRAIL*NRHS, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        DO I = 1, NROWTRAIL
          B_LOCAL(RTRAIL0+I,1:NRHS) = B_LOCAL(RTRAIL0+I,1:NRHS) - CONTRIB(I,1:NRHS)
        END DO
        DEALLOCATE(CONTRIB)
      END IF
      DEALLOCATE(XK)
    END DO

    ! === 2) U*X = Y (substitution arriere) ===
    DO K = NBLK, 1, -1
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = INDXG2L(KSTART, NB, GRID%PR)
        CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS, &
                    CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        DEALLOCATE(DKK)
      END IF

      ALLOCATE(XK(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        XK = B_LOCAL(RS:RS+KW-1,1:NRHS)
      END IF
      CALL MPI_BCAST(XK, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      RHEAD = NUMROC(KSTART-1, NB, GRID%MYROW, GRID%PR)
      IF (RHEAD .GT. 0) THEN
        ALLOCATE(CONTRIB(RHEAD,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          CALL CGEMM('No transpose', 'No transpose', RHEAD, NRHS, KW, &
                      CONE, A_LOCAL(1,CS), MLOC, XK, KW, &
                      (0.0_WP,0.0_WP), CONTRIB, RHEAD)
        END IF
        CALL MPI_BCAST(CONTRIB, RHEAD*NRHS, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        DO I = 1, RHEAD
          B_LOCAL(I,1:NRHS) = B_LOCAL(I,1:NRHS) - CONTRIB(I,1:NRHS)
        END DO
        DEALLOCATE(CONTRIB)
      END IF
      DEALLOCATE(XK)
    END DO
  END SUBROUTINE MPI_LU_SOLVE_DISTRIBUTED

END MODULE MPI_LU_SOLVE_MOD
