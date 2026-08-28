!===============================================================================
! ETAPE 6d : factorisation LU par blocs DISTRIBUEE (bloc-cyclique 2D), AVEC
! PIVOT PARTIEL -- la partie la plus delicate du chantier MPI (voir
! discussion avec l'utilisateur avant cette etape). Construite sur le
! squelette de 6c (mpi_lu_nopiv_mod.f90), avec le PANNEAU DIAGONAL
! remplace par une vraie factorisation pivotee.
!
! CHOIX D'IMPLEMENTATION CENTRAL (documente car c'est une simplification
! deliberee, dans le meme esprit que "SVD plutot que QR-pivote" a l'etape
! 2) : plutot que de reimplementer a la main la recherche de pivot ET
! l'elimination scalaire colonne-par-colonne de facon distribuee (tres
! haut risque de bug, cf. la discussion sur la complexite de PxGETF2),
! on RASSEMBLE le panneau-colonne complet [KSTART:N,KSTART:KEND] (toutes
! les lignes restantes, largeur KW<=NB) sur SON processus proprietaire
! (le meme "owner" que dans 6c), on le factorise LOCALEMENT avec CGETRF
! -- EXACTEMENT le meme appel, deja valide, que dans block_lu_mod.f90
! (etape 1) -- puis on RENVOIE le panneau factorise a ses proprietaires
! d'origine. La recherche de pivot et l'elimination scalaire restent donc
! un algorithme sequentiel bien rode (LAPACK), pas quelque chose que l'on
! reecrit soi-meme de facon distribuee. Le cout reel de la distribution
! (et la vraie complexite de cette etape) se trouve ENTIEREMENT dans :
!   (a) le rassemblement/renvoi du panneau (communication, pas de calcul),
!   (b) la PROPAGATION des echanges de lignes trouves vers le reste de la
!       grille -- colonnes DEJA factorisees a gauche ET colonnes pas
!       encore touchees a droite, sur TOUS les processus, quelle que soit
!       leur colonne de grille (contrairement au panneau lui-meme qui
!       n'existe que sur UNE colonne de processus).
!
! ALGORITHME, pour chaque bloc-colonne/ligne global K :
!   A) RASSEMBLER le panneau [KSTART:N,KSTART:KEND] sur son proprietaire
!      (processus (PROW_K,PCOL_K)) -- seuls les processus de la colonne de
!      grille PCOL_K participent (envoi de leur part locale, via
!      COL_COMM).
!   B) Le proprietaire factorise LOCALEMENT ce panneau via CGETRF (pivot
!      partiel standard, sur TOUTES les lignes du panneau -- exactement le
!      pivot "confine au panneau" de l'etape 1/3, transpose ici a un
!      panneau qui, cette fois, est physiquement reparti avant/apres le
!      rassemblement). IPIV_PANEL_GLOBAL(i) reste au format LOCAL-AU-
!      PANNEAU exact que renvoie CGETRF (PAS converti en indice global) --
!      c'est le format que BLOCK_LU_SOLVE (etape 1) attend directement
!      dans IPIV_GLOBAL(KSTART:KEND), meme convention que block_lu_mod.f90.
!      La conversion en indice GLOBAL (KSTART-1+valeur) n'est faite que
!      LOCALEMENT, au moment de calculer R2 dans la boucle de propagation
!      (etape E) -- necessaire LA pour determiner quel processus possede
!      la ligne gagnante, mais ne doit surtout pas contaminer la valeur
!      stockee/diffusee, sous peine de casser BLOCK_LU_SOLVE en aval (bug
!      reellement rencontre et corrige : residu ~0.5-1e7 avant correction,
!      y compris en grille 1x1 -- preuve que ce n'etait pas un bug de
!      communication distribuee mais une confusion de convention).
!   C) RENVOYER le panneau factorise a ses proprietaires d'origine (meme
!      schema de communication que A, en sens inverse).
!   D) DIFFUSER IPIV_PANEL_GLOBAL(1:KW) a TOUTE LA GRILLE (MPI_BCAST sur
!      GRID%COMM, pas seulement COL_COMM -- indispensable puisque l'etape
!      suivante doit s'appliquer a des processus de N'IMPORTE QUELLE
!      colonne de grille).
!   E) PROPAGER les KW echanges de lignes trouves, DANS L'ORDRE (comme le
!      ferait un CLASWP unique avec ce IPIV), aux colonnes locales A
!      GAUCHE (indice global < KSTART) et A DROITE (indice global > KEND)
!      -- sur TOUS les processus (SWAP_ROWS_RANGE, ci-dessous : swap local
!      si les deux lignes globales impliquees vivent sur le meme MYROW,
!      MPI_SENDRECV entre les deux processus-lignes concernes sinon).
!   F) Le reste est IDENTIQUE a 6c (le panneau, une fois factorise et
!      propage, contient DEJA tout L21 -- CGETRF l'a calcule sur toute la
!      hauteur du panneau, pas seulement le bloc diagonal) : diffuser le
!      bloc diagonal le long de ROW_COMM pour calculer le panneau U par
!      CTRSM, diffuser L (deja pret) le long de ROW_COMM et U le long de
!      COL_COMM, mise a jour du complement de Schur par CGEMM.
!===============================================================================
MODULE MPI_LU_PIVOT_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_MILESTONE
  IMPLICIT NONE

CONTAINS

  ! SHOW_PROGRESS (optionnel) : affichage par PALIERS (PRINT_PROGRESS_
  ! MILESTONE, pas la barre a rafraichissement en place) -- ce module
  ! tourne sous MPI/mpirun, qui relaie le stdout de chaque processus via
  ! un pipe et casse le retour chariot (voir progress_bar_mod.f90 et
  ! mpi_blr_matrix_mod.f90 pour le detail du bug rencontre/corrige).
  ! K/NBLK est un pourcentage EXACT (boucle globalement synchronisee par
  ! les diffusions a chaque etape K) -- affiche uniquement par le
  ! processus 0.
  SUBROUTINE MPI_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                                  INFO, N_SWAPS, MIN_PIVOT, SHOW_PROGRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(INOUT)      :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    INTEGER, INTENT(OUT)            :: IPIV_GLOBAL(N)
    INTEGER, INTENT(OUT)            :: INFO
    INTEGER, INTENT(OUT), OPTIONAL  :: N_SWAPS
    REAL(WP), INTENT(OUT), OPTIONAL :: MIN_PIVOT
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: PANEL(:,:), RBUF(:,:), SBUF(:,:), DKK(:,:), LBUF(:,:), UBUF(:,:)
    INTEGER :: IPIV_LOCAL(NB), IPIV_PANEL_GLOBAL(NB)
    INTEGER :: NBLK, K, KSTART, KEND, KW, MPANEL
    INTEGER :: PROW_K, PCOL_K, OWNER_RANK, RS, CS, RSTART0, NROWS_MINE
    INTEGER :: R, RMLOC_R, RSTART0_R, RNROWS_R, I, J, GI
    INTEGER :: LEFT_COUNT, RIGHT_COUNT, RTRAIL0, CTRAIL0, NROWTRAIL, NCOLTRAIL
    INTEGER :: R1, R2, PROW2, IINFO, INFO_STEP, IERR, NS, LAST_PCT
    LOGICAL :: PROGRESS
    INTEGER :: STATUS(MPI_STATUS_SIZE)
    REAL(WP) :: MP

    INFO = 0
    NS = 0
    MP = HUGE(1.0_WP)
    NBLK = (N + NB - 1) / NB
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    LAST_PCT = -5

    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      MPANEL = N - KSTART + 1

      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)
      OWNER_RANK = PROW_K*GRID%PC + PCOL_K

      INFO_STEP = 0

      ! === A+B+C : rassembler le panneau, factoriser (CGETRF), renvoyer ===
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        RSTART0    = NUMROC(KSTART-1, NB, GRID%MYROW, GRID%PR)
        NROWS_MINE = MLOC - RSTART0
        CS         = INDXG2L(KSTART, NB, GRID%PC)

        IF (GRID%MYROW .EQ. PROW_K) THEN
          ! --- proprietaire : rassemble le panneau complet ---
          ALLOCATE(PANEL(MPANEL,KW))
          IF (NROWS_MINE .GT. 0) THEN
            DO J = 1, KW
              DO I = 1, NROWS_MINE
                GI = INDXL2G(RSTART0+I, NB, GRID%MYROW, GRID%PR)
                PANEL(GI-KSTART+1,J) = A_LOCAL(RSTART0+I, CS+J-1)
              END DO
            END DO
          END IF
          DO R = 0, GRID%PR-1
            IF (R .EQ. PROW_K) CYCLE
            RMLOC_R   = NUMROC(N, NB, R, GRID%PR)
            RSTART0_R = NUMROC(KSTART-1, NB, R, GRID%PR)
            RNROWS_R  = RMLOC_R - RSTART0_R
            IF (RNROWS_R .LE. 0) CYCLE
            ALLOCATE(RBUF(RNROWS_R,KW))
            CALL MPI_RECV(RBUF, RNROWS_R*KW, MPI_COMPLEX, R, 10, GRID%COL_COMM, STATUS, IERR)
            DO J = 1, KW
              DO I = 1, RNROWS_R
                GI = INDXL2G(RSTART0_R+I, NB, R, GRID%PR)
                PANEL(GI-KSTART+1,J) = RBUF(I,J)
              END DO
            END DO
            DEALLOCATE(RBUF)
          END DO

          ! --- factorisation locale, pivot partiel standard (deja valide,
          !     etape 1) ---
          CALL CGETRF(MPANEL, KW, PANEL, MPANEL, IPIV_LOCAL, IINFO)
          IF (IINFO .NE. 0) INFO_STEP = KSTART - 1 + IINFO
          ! IPIV_PANEL_GLOBAL reste au format LOCAL AU PANNEAU (comme le
          ! renvoie directement CGETRF -- valeur v = ligne v DU SOUS-TABLEAU
          ! debutant a KSTART, PAS un indice de ligne global). C'est
          ! EXACTEMENT le format que BLOCK_LU_SOLVE (etape 1) attend dans
          ! IPIV_GLOBAL(KSTART:KEND) puisqu'il rejoue les echanges via
          ! CLASWP(...,B(KSTART,1),...,IPIV(KSTART),1) -- un CLASWP relatif
          ! au meme sous-tableau. La conversion en indice GLOBAL (KSTART-1+v,
          ! necessaire pour NOTRE propagation a nous, cf. R2 plus bas) se
          ! fait SEULEMENT localement dans la boucle de propagation, jamais
          ! dans la valeur stockee/diffusee ici.
          DO J = 1, KW
            IPIV_PANEL_GLOBAL(J) = IPIV_LOCAL(J)
            MP = MIN(MP, ABS(PANEL(J,J)))
          END DO

          ! --- renvoi du panneau factorise a ses proprietaires d'origine ---
          IF (NROWS_MINE .GT. 0) THEN
            DO J = 1, KW
              DO I = 1, NROWS_MINE
                GI = INDXL2G(RSTART0+I, NB, GRID%MYROW, GRID%PR)
                A_LOCAL(RSTART0+I, CS+J-1) = PANEL(GI-KSTART+1,J)
              END DO
            END DO
          END IF
          DO R = 0, GRID%PR-1
            IF (R .EQ. PROW_K) CYCLE
            RMLOC_R   = NUMROC(N, NB, R, GRID%PR)
            RSTART0_R = NUMROC(KSTART-1, NB, R, GRID%PR)
            RNROWS_R  = RMLOC_R - RSTART0_R
            IF (RNROWS_R .LE. 0) CYCLE
            ALLOCATE(SBUF(RNROWS_R,KW))
            DO J = 1, KW
              DO I = 1, RNROWS_R
                GI = INDXL2G(RSTART0_R+I, NB, R, GRID%PR)
                SBUF(I,J) = PANEL(GI-KSTART+1,J)
              END DO
            END DO
            CALL MPI_SEND(SBUF, RNROWS_R*KW, MPI_COMPLEX, R, 11, GRID%COL_COMM, IERR)
            DEALLOCATE(SBUF)
          END DO
          DEALLOCATE(PANEL)
        ELSE
          ! --- autre membre de la colonne de grille PCOL_K : envoie sa
          !     part, recoit la part factorisee en retour ---
          IF (NROWS_MINE .GT. 0) THEN
            ALLOCATE(SBUF(NROWS_MINE,KW))
            SBUF = A_LOCAL(RSTART0+1:MLOC, CS:CS+KW-1)
            CALL MPI_SEND(SBUF, NROWS_MINE*KW, MPI_COMPLEX, PROW_K, 10, GRID%COL_COMM, IERR)
            DEALLOCATE(SBUF)

            ALLOCATE(RBUF(NROWS_MINE,KW))
            CALL MPI_RECV(RBUF, NROWS_MINE*KW, MPI_COMPLEX, PROW_K, 11, GRID%COL_COMM, STATUS, IERR)
            A_LOCAL(RSTART0+1:MLOC, CS:CS+KW-1) = RBUF
            DEALLOCATE(RBUF)
          END IF
        END IF
      END IF

      ! === D : diffuser INFO_STEP et IPIV_PANEL_GLOBAL a TOUTE la grille ===
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        RETURN
      END IF
      CALL MPI_BCAST(IPIV_PANEL_GLOBAL, KW, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IPIV_GLOBAL(KSTART:KEND) = IPIV_PANEL_GLOBAL(1:KW)

      ! === E : propager les echanges de lignes a gauche et a droite ===
      LEFT_COUNT  = NUMROC(KSTART-1, NB, GRID%MYCOL, GRID%PC)
      RIGHT_COUNT = NLOC - NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      DO I = 1, KW
        R1 = KSTART - 1 + I
        R2 = KSTART - 1 + IPIV_PANEL_GLOBAL(I)   ! conversion locale->globale, ICI seulement
        IF (R1 .NE. R2) THEN
          NS = NS + 1
          PROW2 = INDXG2P(R2, NB, GRID%PR)
          IF (LEFT_COUNT .GT. 0) THEN
            CALL SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW_K, PROW2, &
                                  1, LEFT_COUNT)
          END IF
          IF (RIGHT_COUNT .GT. 0) THEN
            CALL SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW_K, PROW2, &
                                  NLOC-RIGHT_COUNT+1, NLOC)
          END IF
        END IF
      END DO

      ! === F : identique a 6c (le panneau contient deja tout L21) ===
      RTRAIL0   = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      CTRAIL0   = NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      NROWTRAIL = MLOC - RTRAIL0
      NCOLTRAIL = NLOC - CTRAIL0

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NCOLTRAIL .GT. 0) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NCOLTRAIL, &
                      CONE, DKK, KW, A_LOCAL(RS,CTRAIL0+1), MLOC)
        END IF
        DEALLOCATE(DKK)
      END IF

      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(LBUF(NROWTRAIL,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          LBUF = A_LOCAL(RTRAIL0+1:MLOC, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(LBUF, NROWTRAIL*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      END IF

      IF (NCOLTRAIL .GT. 0) THEN
        ALLOCATE(UBUF(KW,NCOLTRAIL))
        IF (GRID%MYROW .EQ. PROW_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          UBUF = A_LOCAL(RS:RS+KW-1, CTRAIL0+1:NLOC)
        END IF
        CALL MPI_BCAST(UBUF, KW*NCOLTRAIL, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      END IF

      IF (NROWTRAIL.GT.0 .AND. NCOLTRAIL.GT.0) THEN
        CALL CGEMM('No transpose', 'No transpose', NROWTRAIL, NCOLTRAIL, KW, &
                    -CONE, LBUF, NROWTRAIL, UBUF, KW, &
                    CONE, A_LOCAL(RTRAIL0+1,CTRAIL0+1), MLOC)
      END IF

      IF (NROWTRAIL .GT. 0) DEALLOCATE(LBUF)
      IF (NCOLTRAIL .GT. 0) DEALLOCATE(UBUF)

      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('Factorisation LU MPI (pivot)', K, NBLK, LAST_PCT)
    END DO

    IF (PRESENT(N_SWAPS)) N_SWAPS = NS
    IF (PRESENT(MIN_PIVOT)) THEN
      CALL MPI_ALLREDUCE(MP, MIN_PIVOT, 1, MPI_REAL, MPI_MIN, GRID%COMM, IERR)
    END IF
  END SUBROUTINE MPI_LU_PIVOT_FACTOR

  !-----------------------------------------------------------------------
  ! Echange les lignes GLOBALES R1 et R2 (propriete respective PROW1/PROW2,
  ! 0-based) sur l'intervalle de COLONNES LOCALES [COL_LO:COL_HI] -- swap
  ! purement local si PROW1==PROW2 (les deux lignes vivent sur le meme
  ! MYROW), sinon MPI_SENDRECV entre les deux processus-lignes concernes
  ! (chaque processus determine independamment son role via son propre
  ! MYROW -- aucune coordination supplementaire necessaire, cf. discussion
  ! avec l'utilisateur sur pourquoi c'est la partie delicate de 6d).
  !-----------------------------------------------------------------------
  SUBROUTINE SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW1, PROW2, &
                              COL_LO, COL_HI)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)   :: GRID
    INTEGER, INTENT(IN)        :: NB, MLOC, NLOC, R1, R2, PROW1, PROW2, COL_LO, COL_HI
    COMPLEX(WP), INTENT(INOUT) :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))

    COMPLEX(WP), ALLOCATABLE :: SENDTMP(:), RECVTMP(:)
    INTEGER :: NCOLS, RL1, RL2, PARTNER, IERR
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    IF (R1 .EQ. R2 .OR. COL_HI .LT. COL_LO) RETURN
    NCOLS = COL_HI - COL_LO + 1

    IF (PROW1 .EQ. PROW2) THEN
      IF (GRID%MYROW .EQ. PROW1) THEN
        RL1 = INDXG2L(R1, NB, GRID%PR)
        RL2 = INDXG2L(R2, NB, GRID%PR)
        ALLOCATE(SENDTMP(NCOLS))
        SENDTMP = A_LOCAL(RL1,COL_LO:COL_HI)
        A_LOCAL(RL1,COL_LO:COL_HI) = A_LOCAL(RL2,COL_LO:COL_HI)
        A_LOCAL(RL2,COL_LO:COL_HI) = SENDTMP
        DEALLOCATE(SENDTMP)
      END IF
    ELSE IF (GRID%MYROW .EQ. PROW1) THEN
      RL1 = INDXG2L(R1, NB, GRID%PR)
      PARTNER = PROW2*GRID%PC + GRID%MYCOL
      ALLOCATE(SENDTMP(NCOLS), RECVTMP(NCOLS))
      SENDTMP = A_LOCAL(RL1,COL_LO:COL_HI)
      CALL MPI_SENDRECV(SENDTMP, NCOLS, MPI_COMPLEX, PARTNER, 20, &
                         RECVTMP, NCOLS, MPI_COMPLEX, PARTNER, 20, &
                         GRID%COMM, STATUS, IERR)
      A_LOCAL(RL1,COL_LO:COL_HI) = RECVTMP
      DEALLOCATE(SENDTMP, RECVTMP)
    ELSE IF (GRID%MYROW .EQ. PROW2) THEN
      RL2 = INDXG2L(R2, NB, GRID%PR)
      PARTNER = PROW1*GRID%PC + GRID%MYCOL
      ALLOCATE(SENDTMP(NCOLS), RECVTMP(NCOLS))
      SENDTMP = A_LOCAL(RL2,COL_LO:COL_HI)
      CALL MPI_SENDRECV(SENDTMP, NCOLS, MPI_COMPLEX, PARTNER, 20, &
                         RECVTMP, NCOLS, MPI_COMPLEX, PARTNER, 20, &
                         GRID%COMM, STATUS, IERR)
      A_LOCAL(RL2,COL_LO:COL_HI) = RECVTMP
      DEALLOCATE(SENDTMP, RECVTMP)
    END IF
  END SUBROUTINE SWAP_ROWS_RANGE

END MODULE MPI_LU_PIVOT_MOD
