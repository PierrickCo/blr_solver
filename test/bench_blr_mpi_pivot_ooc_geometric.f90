!===============================================================================
! BENCHMARK COMPARATIF -- solveur BLR+MPI+pivot, VERSION CLASSIQUE vs OUT-OF-
! CORE, MEME structure et MEME protocole de validation que bench_blr_mpi_
! pivot_ooc_compare.f90 -- mais avec une matrice MOINS SYNTHETIQUE :
!
!   Cote bench_blr_mpi_pivot_ooc_compare.f90, le noyau est un simple 1/(1+
!   |ecart d'INDICE|) : la compressibilite d'un bloc ne depend que de la
!   distance entre indices de ligne/colonne, une abstraction commode mais
!   sans rapport avec ce qui rend un VRAI bloc BLR compressible en pratique
!   (methode des elements de frontiere, BEM -- separation GEOMETRIQUE dans
!   l'espace physique).
!
!   ICI : N points repartis sur un CERCLE unite (P(i) = (cos theta_i, sin
!   theta_i), theta_i = 2*pi*(i-1)/N -- cas academique standard pour
!   illustrer BLR/H-matrices), noyau de Helmholtz en espace libre regularise
!   G(r) = exp(i*K*r) / (1+r) avec r = DISTANCE EUCLIDIENNE (et non plus un
!   ecart d'indice) entre P(i) et P(j) -- noyau complexe, OSCILLANT (K =
!   nombre d'onde), utilise tel quel dans la litterature BLR/H-matrices pour
!   les problemes de diffusion d'ondes (acoustique/electromagnetisme). Le
!   parametrage par indice croissant = angle croissant sur le cercle assure
!   que la partition BLR (blocs d'indices consecutifs) reste GEOMETRIQUEMENT
!   admissible (indices voisins = points voisins sur le cercle = champ
!   proche ; indices eloignes = points opposes sur le cercle = champ
!   lointain), exactement l'hypothese que suppose l'algorithme.
!
!   Anti-diagonale forte CONFINEE au bloc diagonal conservee A L'IDENTIQUE
!   (voir bench_blr_mpi_pivot_ooc_compare.f90) pour continuer a exercer un
!   vrai pivotage visible du pivot confine.
!
! Reste de l'infrastructure (comparaison classique/OOC, Xexact, mesure RSS/
! temps par etape, garde-fou "un seul MODE par processus") IDENTIQUE a
! bench_blr_mpi_pivot_ooc_compare.f90 -- voir son commentaire de tete pour
! le detail de ces choix, non repete ici.
!
! Usage : mpirun -np <PR*PC> ./bench_blr_mpi_pivot_ooc_geometric PR PC MODE [N] [NB] \
!         [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [WAVENUMBER] [USE_RSVD]
!         MODE = 'classical' ou 'ooc' (obligatoire).
!         (defauts : N=8000 NB=400 NRHS=6 TOL=1e-4 MAX_MB_PER_PROCESS=16
!          NRHS_BATCH=2 SCRATCH_DIR=./blr_ooc_scratch_bench_geo WAVENUMBER=8.0 USE_RSVD=1)
!         USE_RSVD=1 : compression a l'assemblage par SVD randomisee adaptative
!         (COMPRESS_BLOCK_RSVD_ADAPTIVE) au lieu de CGESVD -- dans les DEUX modes.
!         En mode classical, B=A*Xexact reste construit avec CGESVD (A_CL_FOR_B) :
!         l'erreur inclut donc l'ecart de compression RSVD/CGESVD ; en mode ooc,
!         B vient du noyau exact (FILL_RHS), l'erreur inclut toute la compression.
!===============================================================================
PROGRAM BENCH_BLR_MPI_PIVOT_OOC_GEOMETRIC
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0, USE_RSVD_INT
  LOGICAL :: USE_RSVD_FLAG
  CHARACTER(LEN=40) :: COMPRESS_LABEL
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS, WAVENUMBER
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  CHARACTER(LEN=16) :: MODE
  REAL(WP), PARAMETER :: PI = 3.14159265358979_WP

  ! --- chemin CLASSIQUE (tout en RAM) ---
  TYPE(MPI_BLR_MATRIX_T) :: A_CL
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), XEX_LOCAL(:,:), X_GLOBAL_CL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO_CL, N_SWAPS_CL
  INTEGER :: N_COMPRESSED_CL, N_OFFDIAG_CL, li, I
  INTEGER(8) :: RANK_SUM_CL
  DOUBLE PRECISION :: T_ASM_CL, T_FACTOR_CL, T_SOLVE_CL
  DOUBLE PRECISION :: RSS_ASM_CL, RSS_FACT_CL, RSS_SOLVE_CL
  REAL(WP) :: ERR_CL

  ! --- chemin OUT-OF-CORE ---
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_OOC
  TYPE(OOC_RHS_T) :: RHS
  COMPLEX(WP), ALLOCATABLE :: X_GLOBAL_OOC(:,:)
  INTEGER :: INFO_OOC, N_SWAPS_OOC, N_COMPRESSED_OOC, N_OFFDIAG_OOC
  INTEGER(8) :: RANK_SUM_OOC
  DOUBLE PRECISION :: T_ASM_A_OOC, T_ASM_B_OOC, T_FACTOR_OOC, T_SOLVE_OOC
  DOUBLE PRECISION :: RSS_ASMA_OOC, RSS_ASMB_OOC, RSS_FACT_OOC, RSS_SOLVE_OOC
  DOUBLE PRECISION :: T_CACHE_GET, T_CACHE_EVICT, T_COMPUTE, T_MPICALL
  DOUBLE PRECISION :: T_CACHE_GET_MIN, T_CACHE_EVICT_MIN, T_COMPUTE_MIN, T_MPICALL_MIN
  DOUBLE PRECISION :: T_GETBLK_ASM, T_COMPRESS_ASM, T_IO_ASM
  DOUBLE PRECISION :: T_GETBLK_ASM_MIN, T_COMPRESS_ASM_MIN, T_IO_ASM_MIN
  REAL(WP) :: ERR_OOC
  DOUBLE PRECISION :: T0B, T1B

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  N = 8000; NB = 400; NRHS = 6; TOL = 1.0E-4_WP
  MAX_MB_PER_PROCESS = 16.0_WP; NRHS_BATCH = 2; SCRATCH_DIR = './blr_ooc_scratch_bench_geo'
  WAVENUMBER = 8.0_WP; USE_RSVD_INT = 1

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./bench_blr_mpi_pivot_ooc_geometric PR PC MODE [N] [NB] [NRHS] [TOL] '// &
         '[MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [WAVENUMBER] [USE_RSVD]   (MODE = classical|ooc, obligatoire)'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, MODE)
  IF (MODE .NE. 'classical' .AND. MODE .NE. 'ooc') THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : MODE doit etre "classical" ou "ooc", recu : ', TRIM(MODE)
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN
    CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) MAX_MB_PER_PROCESS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) THEN
    CALL GET_COMMAND_ARGUMENT(9, ARG); READ(ARG,*) NRHS_BATCH
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN
    CALL GET_COMMAND_ARGUMENT(10, SCRATCH_DIR)
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN
    CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) WAVENUMBER
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 12) THEN
    CALL GET_COMMAND_ARGUMENT(12, ARG); READ(ARG,*) USE_RSVD_INT
  END IF
  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
  IF (USE_RSVD_FLAG) THEN
    COMPRESS_LABEL = 'SVD randomisee adaptative'
  ELSE
    COMPRESS_LABEL = 'SVD complete, CGESVD'
  END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1,A,F6.2)') &
       ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
       ' , MAX_MB_PER_PROCESS(OOC)=', MAX_MB_PER_PROCESS, ' , WAVENUMBER=', WAVENUMBER
  IF (MYID0.EQ.0) WRITE(6,'(A,A)') ' Compression a l''assemblage : ', TRIM(COMPRESS_LABEL)

  IF (MODE .EQ. 'classical') THEN
    !===========================================================================
    ! CHEMIN CLASSIQUE (tout en RAM)
    !===========================================================================
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_CL, T_ASM_CL, &
                                  N_COMPRESSED_CL, N_OFFDIAG_CL, RANK_SUM_CL, &
                                  T_GET_BLOCK=T_GETBLK_ASM, T_COMPRESS=T_COMPRESS_ASM, USE_RSVD=USE_RSVD_FLAG)
    RSS_ASM_CL = RSS_NOW_MB()

    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_CL, TOL, INFO_CL, &
                                         N_COMPRESSED_CL, N_OFFDIAG_CL, RANK_SUM_CL, N_SWAPS_CL, T_FACTOR_CL)
    IF (INFO_CL .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO_CL
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    RSS_FACT_CL = RSS_NOW_MB()

    MLOC = 0
    DO li = 1, A_CL%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + (A_CL%BEND(I)-A_CL%BSTART(I)+1)
    END DO
    NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
    CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
    ALLOCATE(XEX_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)), B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
    CALL FILL_XEX(A_CL, GRID, NRHS_LOC, RSTART_MINE, XEX_LOCAL, MLOC)

    ! B = A*Xexact -- construction hors mesure (meme convention que les autres
    ! bench_blr_mpi_pivot_*.f90).
    BLOCK
      TYPE(MPI_BLR_MATRIX_T) :: A_CL_FOR_B
      DOUBLE PRECISION :: T_ASM_DUMMY
      INTEGER :: NC_DUMMY, NO_DUMMY
      INTEGER(8) :: RS_DUMMY
      CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_CL_FOR_B, T_ASM_DUMMY, &
                                    NC_DUMMY, NO_DUMMY, RS_DUMMY)
      CALL APPLY_MPI_BLR_MATRIX(GRID, A_CL_FOR_B, NRHS, NRHS_LOC, XEX_LOCAL, B_LOCAL, MLOC)
    END BLOCK

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0B = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_CL, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
    T1B = MPI_WTIME()
    T_SOLVE_CL = T1B - T0B
    RSS_SOLVE_CL = RSS_NOW_MB()

    ALLOCATE(X_GLOBAL_CL(N,NRHS))
    CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL_CL)
    ERR_CL = COMPUTE_ERR(X_GLOBAL_CL)
    DEALLOCATE(B_LOCAL, XEX_LOCAL)

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' '
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A)') '  [CLASSIQUE, noyau geometrique] Temps (s) et RSS max/processus (Mo)'
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Assemblage    : temps(MAX) = ', T_ASM_CL,   ' s  , RSS = ', RSS_ASM_CL,   ' Mo'
      WRITE(6,'(A,F10.3,A)') '    dont GET_BLOCK (evaluation du noyau) = ', T_GETBLK_ASM, ' s'
      WRITE(6,'(A,F10.3,A)') '    dont compression (' // TRIM(COMPRESS_LABEL) // ') = ', T_COMPRESS_ASM, ' s'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Factorisation : temps(MAX) = ', T_FACTOR_CL,' s  , RSS = ', RSS_FACT_CL,  ' Mo'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Resolution    : temps = ', T_SOLVE_CL, ' s  , RSS = ', RSS_SOLVE_CL, ' Mo'
      WRITE(6,'(A,F10.3,A)') '  TOTAL         : temps = ', T_ASM_CL+T_FACTOR_CL+T_SOLVE_CL, ' s'
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Blocs compresses : ', N_COMPRESSED_CL, ' / ', N_OFFDIAG_CL, ' (', &
           100.0_WP*REAL(N_COMPRESSED_CL,WP)/REAL(MAX(N_OFFDIAG_CL,1),WP), ' %)'
      WRITE(6,'(A,F6.1,A,I0)') '  Rang moyen des blocs compresses = ', &
           REAL(RANK_SUM_CL - (N_OFFDIAG_CL-N_COMPRESSED_CL)*NB,WP)/REAL(MAX(N_COMPRESSED_CL,1),WP), ' / NB=', NB
      WRITE(6,'(A,F6.2,A)') '  Taux de compression REEL (stockage bas-rang / stockage dense equivalent) = ', &
           REAL_COMPRESSION_PCT(N_COMPRESSED_CL, N_OFFDIAG_CL, RANK_SUM_CL), ' %'
      WRITE(6,'(A,I0)') '  Swaps de pivot (doit etre > 0) = ', N_SWAPS_CL
      WRITE(6,'(A,ES12.4)') '  Erreur relative ||X-Xexact|| = ', ERR_CL
      IF (ERR_CL.LT.1.0E-3_WP .AND. N_SWAPS_CL.GT.0) THEN
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_GEOMETRIC [CLASSIQUE] : OK ==='
      ELSE
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_GEOMETRIC [CLASSIQUE] : ECHEC ==='
      END IF
    END IF
    DEALLOCATE(X_GLOBAL_CL)

  ELSE
    !===========================================================================
    ! CHEMIN OUT-OF-CORE
    !===========================================================================
    CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_OOC, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                      T_ASM_A_OOC, N_COMPRESSED_OOC, N_OFFDIAG_OOC, RANK_SUM_OOC, &
                                      T_GET_BLOCK=T_GETBLK_ASM, T_COMPRESS=T_COMPRESS_ASM, T_IO=T_IO_ASM, &
                                      T_GET_BLOCK_MIN=T_GETBLK_ASM_MIN, T_COMPRESS_MIN=T_COMPRESS_ASM_MIN, &
                                      T_IO_MIN=T_IO_ASM_MIN, USE_RSVD=USE_RSVD_FLAG)
    RSS_ASMA_OOC = RSS_NOW_MB()

    T0B = MPI_WTIME()
    CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_OOC, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
    T1B = MPI_WTIME()
    T_ASM_B_OOC = T1B - T0B
    RSS_ASMB_OOC = RSS_NOW_MB()

    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_OOC, TOL, INFO_OOC, &
                                             N_COMPRESSED_OOC, N_OFFDIAG_OOC, RANK_SUM_OOC, N_SWAPS_OOC, T_FACTOR_OOC, &
                                             T_CACHE_GET=T_CACHE_GET, T_CACHE_EVICT=T_CACHE_EVICT, &
                                             T_COMPUTE=T_COMPUTE, T_MPICALL=T_MPICALL, &
                                             T_CACHE_GET_MIN=T_CACHE_GET_MIN, T_CACHE_EVICT_MIN=T_CACHE_EVICT_MIN, &
                                             T_COMPUTE_MIN=T_COMPUTE_MIN, T_MPICALL_MIN=T_MPICALL_MIN)
    IF (INFO_OOC .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO_OOC
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    RSS_FACT_OOC = RSS_NOW_MB()

    ALLOCATE(X_GLOBAL_OOC(N,NRHS))
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0B = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_OOC, RHS, CONSUME_X_BATCH)
    T1B = MPI_WTIME()
    T_SOLVE_OOC = T1B - T0B
    RSS_SOLVE_OOC = RSS_NOW_MB()

    ERR_OOC = COMPUTE_ERR(X_GLOBAL_OOC)

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' '
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A)') '  [OOC, noyau geometrique] Temps (s) et RSS max/processus (Mo)'
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Assemblage A  : temps(MAX) = ', T_ASM_A_OOC, ' s  , RSS = ', RSS_ASMA_OOC, ' Mo'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont GET_BLOCK (evaluation du noyau)          : MAX=', T_GETBLK_ASM, &
           ' s , MIN=', T_GETBLK_ASM_MIN, ' s'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont compression (' // TRIM(COMPRESS_LABEL) // ') : MAX=', T_COMPRESS_ASM, &
           ' s , MIN=', T_COMPRESS_ASM_MIN, ' s'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont E/S disque (ecriture diag+cache+flush)   : MAX=', T_IO_ASM, &
           ' s , MIN=', T_IO_ASM_MIN, ' s'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Assemblage B  : temps = ', T_ASM_B_OOC, ' s  , RSS = ', RSS_ASMB_OOC, ' Mo'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Factorisation : temps(MAX) = ', T_FACTOR_OOC,' s  , RSS = ', RSS_FACT_OOC, ' Mo'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont cache GET (lookup+lecture disque)        : MAX=', T_CACHE_GET, &
           ' s , MIN=', T_CACHE_GET_MIN, ' s'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont eviction panneau (ecriture disque)       : MAX=', T_CACHE_EVICT, &
           ' s , MIN=', T_CACHE_EVICT_MIN, ' s'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont calcul pur (CGETRF/TRSM/Schur/swap)      : MAX=', T_COMPUTE, &
           ' s , MIN=', T_COMPUTE_MIN, ' s'
      WRITE(6,'(A,F10.3,A,F10.3,A)') '    dont communications MPI (BCAST)               : MAX=', T_MPICALL, &
           ' s , MIN=', T_MPICALL_MIN, ' s'
      WRITE(6,'(A,F10.3,A)')         '    somme des 4 categories (MAX) = ', T_CACHE_GET+T_CACHE_EVICT+T_COMPUTE+T_MPICALL, ' s'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Resolution    : temps = ', T_SOLVE_OOC, ' s  , RSS = ', RSS_SOLVE_OOC,' Mo'
      WRITE(6,'(A,F10.3,A)') '  TOTAL         : temps = ', T_ASM_A_OOC+T_ASM_B_OOC+T_FACTOR_OOC+T_SOLVE_OOC, ' s'
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Blocs compresses : ', N_COMPRESSED_OOC, ' / ', N_OFFDIAG_OOC, ' (', &
           100.0_WP*REAL(N_COMPRESSED_OOC,WP)/REAL(MAX(N_OFFDIAG_OOC,1),WP), ' %)'
      WRITE(6,'(A,F6.1,A,I0)') '  Rang moyen des blocs compresses = ', &
           REAL(RANK_SUM_OOC - (N_OFFDIAG_OOC-N_COMPRESSED_OOC)*NB,WP)/REAL(MAX(N_COMPRESSED_OOC,1),WP), ' / NB=', NB
      WRITE(6,'(A,F6.2,A)') '  Taux de compression REEL (stockage bas-rang / stockage dense equivalent) = ', &
           REAL_COMPRESSION_PCT(N_COMPRESSED_OOC, N_OFFDIAG_OOC, RANK_SUM_OOC), ' %'
      WRITE(6,'(A,I0)') '  Swaps de pivot (doit etre > 0) = ', N_SWAPS_OOC
      WRITE(6,'(A,ES12.4)') '  Erreur relative ||X-Xexact|| = ', ERR_OOC
      IF (ERR_OOC.LT.1.0E-3_WP .AND. N_SWAPS_OOC.GT.0) THEN
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_GEOMETRIC [OOC] : OK ==='
      ELSE
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_GEOMETRIC [OOC] : ECHEC ==='
      END IF
    END IF
    DEALLOCATE(X_GLOBAL_OOC)
    CALL OOC_RHS_FREE(RHS)
    CALL MPI_BLR_MATRIX_OOC_FREE(A_OOC)
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Taux de compression REEL (stockage bas-rang effectif / stockage dense
  ! equivalent) -- voir examples/example_blr_mpi_pivot.f90 pour le detail
  ! du correctif : RANK_SUM melange rang REEL (blocs compresses) et un
  ! "rang virtuel"=NB par bloc en repli DENSE (convention de comptage
  ! interne a MPI_ASSEMBLE_BLR_MATRIX[_OOC], pas un vrai stockage UV) --
  ! une formule naive 2*NB*RANK_SUM/(NB^2*N_OFFDIAG) donnerait a tort 200%
  ! si tous les blocs tombaient dense (invisible ici tant que ce banc
  ! reste a 100% de blocs compresses, mais faux en general).
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  ! RSS reelle du processus, max sur toute la grille.
  DOUBLE PRECISION FUNCTION RSS_NOW_MB()
    DOUBLE PRECISION :: TMP_MB
    INTEGER :: IERR3
    TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
    CALL MPI_ALLREDUCE(TMP_MB, RSS_NOW_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR3)
  END FUNCTION RSS_NOW_MB

  REAL(WP) FUNCTION COMPUTE_ERR(XG)
    COMPLEX(WP), INTENT(IN) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: B_REF(:,:)
    INTEGER :: GI, R
    COMPUTE_ERR = -1.0_WP
    IF (MYID0 .NE. 0) RETURN
    ALLOCATE(B_REF(N,NRHS))
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = XEX_FORMULA(GI,R)
      END DO
    END DO
    COMPUTE_ERR = SQRT(SUM(ABS(XG-B_REF)**2)) / SQRT(SUM(ABS(B_REF)**2))
    DEALLOCATE(B_REF)
  END FUNCTION COMPUTE_ERR

  ! Position du point GI sur le cercle unite -- indice croissant = angle
  ! croissant, pour que la partition par blocs d'indices CONSECUTIFS reste
  ! geometriquement admissible (voir commentaire de tete).
  REAL(WP) FUNCTION PTX(GI)
    INTEGER, INTENT(IN) :: GI
    PTX = COS(2.0_WP*PI*REAL(GI-1,WP)/REAL(N,WP))
  END FUNCTION PTX
  REAL(WP) FUNCTION PTY(GI)
    INTEGER, INTENT(IN) :: GI
    PTY = SIN(2.0_WP*PI*REAL(GI-1,WP)/REAL(N,WP))
  END FUNCTION PTY

  ! Noyau de Helmholtz en espace libre, REGULARISE : exp(i*K*r)/(1+r), r =
  ! distance EUCLIDIENNE (et non plus un ecart d'indice) entre les points GI
  ! et GJ places sur le cercle -- oscillant (K = WAVENUMBER), decroissance
  ! reelle en 1/r pour r grand, meme famille que les noyaux utilises pour
  ! les problemes de diffusion d'ondes en BLR/H-matrices. R=0 impossible ici
  ! (GI!=GJ sur cette branche, deux points distincts du cercle).
  COMPLEX(WP) FUNCTION HELMHOLTZ_KERNEL(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: R
    R = SQRT((PTX(GI)-PTX(GJ))**2 + (PTY(GI)-PTY(GJ))**2)
    HELMHOLTZ_KERNEL = EXP(CMPLX(0.0_WP, WAVENUMBER*R, KIND=WP)) / CMPLX(1.0_WP+R, 0.0_WP, KIND=WP)
  END FUNCTION HELMHOLTZ_KERNEL

  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ, NN)
    INTEGER, INTENT(IN) :: GI, GJ, NN
    INTEGER :: BSTART0, BEND0, GJ_STRONG
    BSTART0 = ((GI-1)/NB)*NB + 1
    BEND0 = MIN(BSTART0+NB-1, NN)
    GJ_STRONG = BEND0 - (GI-BSTART0)
    IF (GJ .EQ. GJ_STRONG) THEN
      LOWRANK_A = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
    ELSE
      LOWRANK_A = HELMHOLTZ_KERNEL(GI, GJ)
    END IF
  END FUNCTION LOWRANK_A

  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        BLK(II-ISTART+1, JJ-JSTART+1) = LOWRANK_A(II, JJ, N)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  SUBROUTINE FILL_RHS(ISTART, IEND, RSTART, REND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RSTART, REND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, REND-RSTART+1)
    INTEGER :: II, JJ, RR
    COMPLEX(WP) :: ACC
    COMPLEX(WP), PARAMETER :: CZERO = (0.0_WP,0.0_WP)
    DO RR = RSTART, REND
      DO II = ISTART, IEND
        ACC = CZERO
        DO JJ = 1, N
          ACC = ACC + LOWRANK_A(II,JJ,N) * XEX_FORMULA(JJ,RR)
        END DO
        BLK(II-ISTART+1, RR-RSTART+1) = ACC
      END DO
    END DO
  END SUBROUTINE FILL_RHS

  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    COMPLEX(WP), ALLOCATABLE :: X_BATCH_GLOBAL(:,:)
    INTEGER :: WIDTH
    WIDTH = BATCH_REND - BATCH_RSTART + 1
    ALLOCATE(X_BATCH_GLOBAL(N,WIDTH))
    CALL GATHER_2D_RHS(GRID, N, NB, WIDTH, X_LOCAL, LDX, NRHS_LOC_B, 0, X_BATCH_GLOBAL)
    IF (MYID0 .EQ. 0) X_GLOBAL_OOC(:,BATCH_RSTART:BATCH_REND) = X_BATCH_GLOBAL
    DEALLOCATE(X_BATCH_GLOBAL)
  END SUBROUTINE CONSUME_X_BATCH

  SUBROUTINE FILL_XEX(AB, GD, NRL, RSTART, XL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: XL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, GRR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          GRR = RSTART + RR - 1
          XL(RS + (GII-AB%BSTART(II)), RR) = XEX_FORMULA(GII, GRR)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_XEX

  SUBROUTINE APPLY_MPI_BLR_MATRIX(GD, AB, NR, NRL, XL, YL, ML)
    TYPE(GRID_T), INTENT(IN)           :: GD
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    INTEGER, INTENT(IN)                :: NR, NRL, ML
    COMPLEX(WP), INTENT(IN)            :: XL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), INTENT(OUT)           :: YL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), PARAMETER :: CONE=(1.0_WP,0.0_WP), CZERO=(0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: XJ_LOC(:,:), XJ_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: J, lj, li, I, NJ, MI, RS, RE, C, IERR2, PROWJ, PCOLJ, MLOC_ALL

    IF (NRL .GT. 0) YL = CZERO

    ALLOCATE(RSTARTS(0:GD%PC-1), RLOCS(0:GD%PC-1), RCOUNTS(0:GD%PC-1), RDISPLS(0:GD%PC-1))
    DO C = 0, GD%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NR, GD%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    DO J = 1, AB%NBLK
      NJ = AB%BEND(J)-AB%BSTART(J)+1
      PROWJ = MOD(J-1, GD%PR); PCOLJ = MOD(J-1, GD%PC)

      ALLOCATE(XJ_FULL(NJ,NR))
      IF (GD%MYROW .EQ. PROWJ) THEN
        RS = LOCAL_ROW_OFFSET(AB, GD, J)
        ALLOCATE(XJ_LOC(NJ,MAX(NRL,1)))
        IF (NRL .GT. 0) XJ_LOC = XL(RS:RS+NJ-1,1:NRL)
        DO C = 0, GD%PC-1
          RCOUNTS(C) = NJ*RLOCS(C); RDISPLS(C) = NJ*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XJ_LOC, NJ*NRL, MPI_COMPLEX, XJ_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                             GD%ROW_COMM, IERR2)
        DEALLOCATE(XJ_LOC)
      END IF
      CALL MPI_BCAST(XJ_FULL, NJ*NR, MPI_COMPLEX, PROWJ, GD%COL_COMM, IERR2)

      MLOC_ALL = 0
      DO li = 1, AB%NBLKR_LOC
        I = (li-1)*GD%PR + GD%MYROW + 1
        MLOC_ALL = MLOC_ALL + (AB%BEND(I)-AB%BSTART(I)+1)
      END DO
      ALLOCATE(CONTRIB(MLOC_ALL,NR))
      IF (GD%MYCOL .EQ. PCOLJ) THEN
        lj = (J-1)/GD%PC + 1
        CONTRIB = CZERO
        RE = 0
        DO li = 1, AB%NBLKR_LOC
          I = (li-1)*GD%PR + GD%MYROW + 1
          MI = AB%BEND(I)-AB%BSTART(I)+1
          IF (I .EQ. J) THEN
            CALL CGEMM('No transpose','No transpose', MI, NR, NJ, CONE, AB%DIAGBLK(I)%D, MI, &
                        XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
          ELSE
            CALL BLR_BLOCK_APPLY(AB%OFFBLK(li,lj), NR, CONE, XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
          END IF
          RE = RE + MI
        END DO
      END IF
      ALLOCATE(CONTRIB_LOC(MLOC_ALL,MAX(NRL,1)))
      DO C = 0, GD%PC-1
        RCOUNTS(C) = MLOC_ALL*RLOCS(C); RDISPLS(C) = MLOC_ALL*(RSTARTS(C)-1)
      END DO
      CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, CONTRIB_LOC, MLOC_ALL*NRL, &
                         MPI_COMPLEX, PCOLJ, GD%ROW_COMM, IERR2)
      IF (NRL .GT. 0) YL(1:MLOC_ALL,1:NRL) = YL(1:MLOC_ALL,1:NRL) + CONTRIB_LOC(1:MLOC_ALL,1:NRL)
      DEALLOCATE(CONTRIB, CONTRIB_LOC)
      DEALLOCATE(XJ_FULL)
    END DO
    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE APPLY_MPI_BLR_MATRIX

END PROGRAM BENCH_BLR_MPI_PIVOT_OOC_GEOMETRIC
