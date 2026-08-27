!===============================================================================
! RESOLUTION MPI DISTRIBUEE, VRAIE DISTRIBUTION 2D DU SECOND MEMBRE (demande
! explicite : completer le facteur PR deja obtenu par mpi_lu_solve_mod.f90
! avec le facteur PC manquant, quitte a introduire un mecanisme de
! communication reellement nouveau -- ALLGATHERV, jamais utilise ailleurs
! dans ce projet).
!
! DISTRIBUTION DE B : B_LOCAL(MLOC,NRHS_LOC) -- MLOC exactement comme A
! (lignes, NUMROC sur PR, INCHANGE par rapport a mpi_lu_solve_mod.f90),
! NRHS_LOC = un DECOUPAGE SIMPLE PAR BLOCS (pas bloc-cyclique -- NRHS est un
! lot de seconds membres independants, pas un axe spatial de A, donc pas
! besoin du decoupage a la NB qui a du sens pour N) de NRHS sur PC
! (RHS_LOCAL_RANGE, ci-dessous) -- AUCUNE replication residuelle : chaque
! colonne de grille ne stocke NI ne calcule que SA PROPRE tranche de
! seconds membres.
!
! LE PROBLEME REEL (identifie AVANT d'ecrire ce module, avec l'utilisateur) :
! une fois X(K) resolu, chaque colonne de la ligne PROW_K n'en detient QUE
! SA PROPRE tranche (KW x NRHS_LOC) -- mais c'est le processus de colonne
! PCOL_K (proprietaire de L(I,K)) qui doit calculer L(I,K)*X(K), et il n'a
! QUE SA PROPRE tranche RHS, pas les PC-1 autres. Il faut donc RECONSTRUIRE
! X(K) en entier (KW x NRHS) avant de pouvoir l'utiliser -- exactement le
! "mismatch ligne/colonne" qu'une distribution 2D complete de B soulevait
! deja en theorie.
!
! ALGORITHME (par etape de bloc K, substitution avant -- l'arriere est
! symetrique, tetes au lieu de trainantes) :
!   a) DKK diffuse le long de ROW_COMM (racine=PCOL_K) -- IDENTIQUE a la
!      version 1D (mpi_lu_solve_mod.f90) et a la factorisation elle-meme.
!   b) CTRSM LOCAL sur NOTRE PROPRE tranche NRHS_LOC -- PLUS de redondance
!      inter-colonnes ici (contrairement a la version 1D ou B etait
!      replique dans la ligne et CHAQUE colonne refaisait le MEME CTRSM) :
!      chaque colonne calcule un morceau DIFFERENT, un vrai partage de
!      travail en plus du partage memoire.
!   c) MPI_ALLGATHERV DANS ROW_COMM (ligne PROW_K) : reconstruit X(K)
!      COMPLET (KW x NRHS) a partir des PC tranches NRHS_LOC -- LE SEUL
!      mecanisme reellement NOUVEAU de cette version.
!   d) diffusion de ce X(K) complet le long de COL_COMM (racine=PROW_K) --
!      IDENTIQUE a la version 1D (deja full-NRHS-width la-bas aussi).
!   e) le proprietaire de L(:,K) (colonne PCOL_K) calcule la contribution
!      COMPLETE (NROWTRAIL x NRHS) via CGEMM -- IDENTIQUE a la version 1D.
!   f) MPI_SCATTERV le long de ROW_COMM (racine=PCOL_K), PAS une diffusion :
!      chaque colonne ne RECOIT QUE sa propre tranche NRHS_LOC -- evite le
!      gaspillage de bande passante d'une diffusion complete suivie d'un
!      rejet local (contrairement a ce qu'on aurait fait par simplicite).
!
! COMPROMIS ATTENDU (signale a l'utilisateur avant d'ecrire ce module) :
! gain memoire reel et complet (facteur PR x PC au lieu de PR seul), au
! prix d'un ALLGATHERV + 2 diffusions par etape au lieu de 2 diffusions
! (mpi_lu_solve_mod.f90) -- strictement PLUS d'appels collectifs par etape.
! Sur ce banc d'essai (WSL2, MPI en boucle locale), l'experience recente
! avec la proposition 3 MPI (BCAST_BLR_BLOCK) suggere que ce surcout de
! latence pourrait depasser le gain de volume -- A MESURER, pas supposer.
!===============================================================================
MODULE MPI_LU_SOLVE_2D_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE MPI_LU_PIVOT_MOD, ONLY: SWAP_ROWS_RANGE
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! Decoupage de [1:NRHS] en PC tranches CONTIGUES aussi equilibrees que
  ! possible (PAS bloc-cyclique, contrairement a BLOCK_PARTITION/NUMROC --
  ! NRHS n'a pas la structure de "nombreux petits panneaux" qui justifie le
  ! bloc-cyclique pour N). Calcul LOCAL, sans communication, pour N'IMPORTE
  ! QUELLE colonne C (pas seulement la sienne) -- necessaire pour que
  ! chaque processus construise seul les tableaux RECVCOUNTS/DISPLS d'un
  ! ALLGATHERV/SCATTERV, comme NUMROC le permet deja pour les lignes de A.
  !-----------------------------------------------------------------------
  SUBROUTINE RHS_LOCAL_RANGE(NRHS, PC, C, RSTART, RLOC)
    INTEGER, INTENT(IN)  :: NRHS, PC, C
    INTEGER, INTENT(OUT) :: RSTART, RLOC
    INTEGER :: BASE, REM
    BASE = NRHS / PC
    REM  = MOD(NRHS, PC)
    IF (C .LT. REM) THEN
      RLOC   = BASE + 1
      RSTART = C*(BASE+1) + 1
    ELSE
      RLOC   = BASE
      RSTART = REM*(BASE+1) + (C-REM)*BASE + 1
    END IF
  END SUBROUTINE RHS_LOCAL_RANGE

  !-----------------------------------------------------------------------
  ! Nombre de seconds membres locaux pour CE processus (a appeler par
  ! l'appelant AVANT d'allouer B_LOCAL, meme convention que NUMROC pour A).
  !-----------------------------------------------------------------------
  FUNCTION RHS_NUMROC(NRHS, PC, MYCOL) RESULT(NRHS_LOC)
    INTEGER, INTENT(IN) :: NRHS, PC, MYCOL
    INTEGER :: NRHS_LOC, DUMMY_START
    CALL RHS_LOCAL_RANGE(NRHS, PC, MYCOL, DUMMY_START, NRHS_LOC)
  END FUNCTION RHS_NUMROC

  ! B_LOCAL : distribue PAR LIGNES comme A (MLOC, meme NUMROC/PR que A) ET
  ! PAR COLONNES sur les seconds membres (NRHS_LOC, RHS_NUMROC/PC ci-dessus)
  ! -- AUCUNE replication residuelle, contrairement a mpi_lu_solve_mod.f90.
  SUBROUTINE MPI_LU_SOLVE_DISTRIBUTED_2D(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                                          NRHS, NRHS_LOC, B_LOCAL, LDB)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)     :: GRID
    INTEGER, INTENT(IN)          :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(IN)      :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    INTEGER, INTENT(IN)          :: IPIV_GLOBAL(N)
    INTEGER, INTENT(IN)          :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)   :: B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1))

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:)
    COMPLEX(WP), ALLOCATABLE :: CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, KSTART, KEND, KW
    INTEGER :: PROW_K, PCOL_K, RS, CS
    INTEGER :: RTRAIL0, NROWTRAIL, RHEAD, I, R1, R2, PROW2, C
    INTEGER :: IERR

    NBLK = (N + NB - 1) / NB

    ! Tranches RHS de TOUTES les colonnes de grille -- calcul local, sans
    ! communication, reutilise a chaque etape pour les RECVCOUNTS/DISPLS.
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO
    ALLOCATE(RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))

    ! === 0) P*B : sur la tranche RHS locale uniquement ===
    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      DO I = 1, KEND-KSTART+1
        R1 = KSTART - 1 + I
        R2 = KSTART - 1 + IPIV_GLOBAL(KSTART-1+I)
        IF (R1 .NE. R2 .AND. NRHS_LOC .GT. 0) THEN
          PROW2 = INDXG2P(R2, NB, GRID%PR)
          CALL SWAP_ROWS_RANGE(GRID, NB, B_LOCAL, MLOC, NRHS_LOC, R1, R2, PROW_K, PROW2, 1, NRHS_LOC)
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

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = INDXG2L(KSTART, NB, GRID%PR)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)

        ! c) ALLGATHERV : reconstruit X(K) complet a partir des tranches
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C)
          RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF

      ! d) diffusion de X(K) complet le long de COL_COMM
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      ! e) contribution complete par le proprietaire de L(:,K)
      RTRAIL0   = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      NROWTRAIL = MLOC - RTRAIL0
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          CALL CGEMM('No transpose', 'No transpose', NROWTRAIL, NRHS, KW, &
                      CONE, A_LOCAL(RTRAIL0+1,CS), MLOC, XK_FULL, KW, &
                      CZERO, CONTRIB, NROWTRAIL)
        END IF

        ! f) SCATTERV : chaque colonne ne recoit QUE sa propre tranche
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C)
          RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          DO I = 1, NROWTRAIL
            B_LOCAL(RTRAIL0+I,1:NRHS_LOC) = B_LOCAL(RTRAIL0+I,1:NRHS_LOC) - CONTRIB_LOC(I,1:NRHS_LOC)
          END DO
        END IF
        DEALLOCATE(CONTRIB, CONTRIB_LOC)
      END IF
      DEALLOCATE(XK_FULL)
    END DO

    ! === 2) U*X = Y (substitution arriere) -- meme schema, tetes au lieu
    !     de trainantes ===
    DO K = NBLK, 1, -1
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = INDXG2L(KSTART, NB, GRID%PR)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)

        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C)
          RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF

      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      RHEAD = NUMROC(KSTART-1, NB, GRID%MYROW, GRID%PR)
      IF (RHEAD .GT. 0) THEN
        ALLOCATE(CONTRIB(RHEAD,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          CALL CGEMM('No transpose', 'No transpose', RHEAD, NRHS, KW, &
                      CONE, A_LOCAL(1,CS), MLOC, XK_FULL, KW, &
                      CZERO, CONTRIB, RHEAD)
        END IF

        ALLOCATE(CONTRIB_LOC(RHEAD,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = RHEAD*RLOCS(C)
          RDISPLS(C) = RHEAD*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, RHEAD*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          DO I = 1, RHEAD
            B_LOCAL(I,1:NRHS_LOC) = B_LOCAL(I,1:NRHS_LOC) - CONTRIB_LOC(I,1:NRHS_LOC)
          END DO
        END IF
        DEALLOCATE(CONTRIB, CONTRIB_LOC)
      END IF
      DEALLOCATE(XK_FULL)
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_LU_SOLVE_DISTRIBUTED_2D

  !-----------------------------------------------------------------------
  ! Rassemble sur ROOT un second membre distribue PAR LIGNES (NUMROC/PR,
  ! comme A) ET PAR COLONNES (RHS_LOCAL_RANGE/PC, PAS bloc-cyclique) --
  ! meme schema point-a-point que GATHER_2D (mpi_grid_mod.f90), colonne
  ! mappee via RHS_LOCAL_RANGE au lieu de INDXL2G/NB. Utilite reelle : test
  ! et benchmark (validation/mesure), pas le chemin chaud du solveur.
  !-----------------------------------------------------------------------
  SUBROUTINE GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MLOC, NRHS_LOC, ROOT, B_GLOBAL)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)   :: GRID
    INTEGER, INTENT(IN)        :: N, NB, NRHS, MLOC, NRHS_LOC, ROOT
    COMPLEX(WP), INTENT(IN)    :: B_LOCAL(MLOC,MAX(NRHS_LOC,1))
    COMPLEX(WP), INTENT(OUT)   :: B_GLOBAL(N,NRHS)

    COMPLEX(WP), ALLOCATABLE :: RECVBUF(:,:)
    INTEGER :: P, RMLOC, RNRHSLOC, RSTART, RSTART_MINE, NRHS_LOC_MINE
    INTEGER :: RMYROW, RMYCOL, II, JJ, GII, GJJ, IERR
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    IF (GRID%MYID .EQ. ROOT) THEN
      CALL RHS_LOCAL_RANGE(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC_MINE)
      DO JJ = 1, NRHS_LOC_MINE
        GJJ = RSTART_MINE + JJ - 1
        DO II = 1, MLOC
          GII = INDXL2G(II, NB, GRID%MYROW, GRID%PR)
          B_GLOBAL(GII,GJJ) = B_LOCAL(II,JJ)
        END DO
      END DO
      DO P = 0, GRID%NPROCS-1
        IF (P .EQ. ROOT) CYCLE
        RMYROW = P / GRID%PC
        RMYCOL = MOD(P, GRID%PC)
        RMLOC = NUMROC(N, NB, RMYROW, GRID%PR)
        CALL RHS_LOCAL_RANGE(NRHS, GRID%PC, RMYCOL, RSTART, RNRHSLOC)
        IF (RMLOC.EQ.0 .OR. RNRHSLOC.EQ.0) CYCLE
        ALLOCATE(RECVBUF(RMLOC,RNRHSLOC))
        CALL MPI_RECV(RECVBUF, RMLOC*RNRHSLOC, MPI_COMPLEX, P, 0, GRID%COMM, &
                       STATUS, IERR)
        DO JJ = 1, RNRHSLOC
          GJJ = RSTART + JJ - 1
          DO II = 1, RMLOC
            GII = INDXL2G(II, NB, RMYROW, GRID%PR)
            B_GLOBAL(GII,GJJ) = RECVBUF(II,JJ)
          END DO
        END DO
        DEALLOCATE(RECVBUF)
      END DO
    ELSE
      IF (MLOC.GT.0 .AND. NRHS_LOC.GT.0) THEN
        CALL MPI_SEND(B_LOCAL, MLOC*NRHS_LOC, MPI_COMPLEX, ROOT, 0, GRID%COMM, IERR)
      END IF
    END IF
  END SUBROUTINE GATHER_2D_RHS

END MODULE MPI_LU_SOLVE_2D_MOD
