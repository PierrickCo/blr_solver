!===============================================================================
! Resolution MPI distribuee, distribution 2D complete du second membre --
! B est reparti a la fois sur les lignes (facteur PR) et les colonnes
! (facteur PC), sans replication residuelle.
!
! DISTRIBUTION DE B : B_LOCAL(MLOC,NRHS_LOC) -- MLOC exactement comme A
! (lignes, NUMROC sur PR), NRHS_LOC = un decoupage simple par blocs
! contigus (pas bloc-cyclique : NRHS est un lot de seconds membres
! independants, pas un axe spatial de A, donc pas besoin du decoupage a la
! NB qui a du sens pour N) de NRHS sur PC (RHS_LOCAL_RANGE, ci-dessous) --
! chaque colonne de grille ne stocke ni ne calcule que sa propre tranche
! de seconds membres.
!
! LE PROBLEME A RESOUDRE : une fois X(K) resolu, chaque colonne de la ligne
! PROW_K n'en detient que sa propre tranche (KW x NRHS_LOC) -- mais c'est
! le processus de colonne PCOL_K (proprietaire de L(I,K)) qui doit calculer
! L(I,K)*X(K), et il n'a que sa propre tranche RHS, pas les PC-1 autres. Il
! faut donc reconstruire X(K) en entier (KW x NRHS) avant de pouvoir
! l'utiliser -- le "mismatch ligne/colonne" inherent a une distribution 2D
! complete de B.
!
! ALGORITHME (par etape de bloc K, substitution avant -- l'arriere est
! symetrique, tetes au lieu de trainantes) :
!   a) DKK diffuse le long de ROW_COMM (racine=PCOL_K), comme lors de la
!      factorisation.
!   b) CTRSM local sur la propre tranche NRHS_LOC de chaque processus --
!      chaque colonne calcule un morceau different, sans redondance de
!      calcul entre colonnes.
!   c) MPI_ALLGATHERV dans ROW_COMM (ligne PROW_K) : reconstruit X(K)
!      complet (KW x NRHS) a partir des PC tranches NRHS_LOC.
!   d) diffusion de ce X(K) complet le long de COL_COMM (racine=PROW_K).
!   e) le proprietaire de L(:,K) (colonne PCOL_K) calcule la contribution
!      complete (NROWTRAIL x NRHS) via CGEMM.
!   f) MPI_SCATTERV le long de ROW_COMM (racine=PCOL_K), pas une diffusion :
!      chaque colonne ne recoit que sa propre tranche NRHS_LOC, evitant le
!      gaspillage de bande passante d'une diffusion complete suivie d'un
!      rejet local.
!
! COMPROMIS : gain memoire complet (facteur PR x PC au lieu de PR seul), au
! prix d'un ALLGATHERV et de deux diffusions supplementaires par etape par
! rapport a une distribution 1D de B -- mesure sur ce projet plus rapide ET
! moins gourmande en memoire qu'une distribution 1D a toutes les echelles
! testees, malgre le surcout apparent en nombre d'appels collectifs.
!===============================================================================
MODULE MPI_LU_SOLVE_2D_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_MILESTONE
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
