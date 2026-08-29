!===============================================================================
! BENCHMARK -- impact de COMPRESS_WORTHWHILE_RATIO (blr_compress_mod.f90,
! desormais une variable de module modifiable a l'execution via
! SET_COMPRESS_WORTHWHILE_RATIO, PLUS un PARAMETER fige) sur le solveur
! BLR+MPI de reference AVEC PIVOT, DANS LE CAS MARGINAL identifie par
! bench_blr_mpi_pivot_compressibility.f90 : a EPS fixe proche du point de
! transition compressible/incompressible (ex. EPS=1e-6 sur la matrice de
! reference -- voir la memoire du projet, sweep_compressibility.sh), une
! grande partie des blocs hors-diagonale atteignent un rang final PROCHE du
! plafond RMAX = MAX(1, INT(RATIO*MIN(MI,NJ))-1) -- c'est PRECISEMENT le
! regime ou RATIO influe le plus (ailleurs -- tres compressible ou tres
! incompressible -- RATIO ne change presque rien, voir le sweep de
! compressibilite).
!
! MESURE, a EPS/N/NB/NRHS/TOL/grille FIXES, l'effet d'un SEUL parametre
! (RATIO) sur :
!   - RSS reelle (mesuree) et memoire theorique (BLR_BLOCK_T), + leur ratio
!   - precision numerique (||X_calcule - X_exact|| / ||X_exact||)
!   - temps par etape (assemblage, factorisation, resolution)
!
! MATRICE : IDENTIQUE a bench_blr_mpi_pivot_compressibility.f90 (anti-
! diagonale forte + noyau de Cauchy amorti + bruit plein-rang EPS*NOISE).
!
! Usage : mpirun -np <PR*PC> ./bench_blr_mpi_pivot_ratio PR PC N NB NRHS TOL EPS RATIO
!===============================================================================
PROGRAM BENCH_BLR_MPI_PIVOT_RATIO
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL, EPS, RATIO
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), XEX_LOCAL(:,:), X_GLOBAL(:,:)
  COMPLEX(WP), ALLOCATABLE :: B_REF(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I, GI, GJ, R
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  DOUBLE PRECISION :: RSS_FACT_MB, MEM_ANALYTIC_MB, TMP_MB
  REAL(WP) :: ERR, AVG_RANK

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 8) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./bench_blr_mpi_pivot_ratio PR PC N NB NRHS TOL EPS RATIO'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL
  CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) EPS
  CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) RATIO

  ! SEULE ligne qui differe de bench_blr_mpi_pivot_compressibility.f90 --
  ! surcharge COMPRESS_WORTHWHILE_RATIO AVANT tout assemblage/factorisation
  ! (voir blr_compress_mod.f90/SET_COMPRESS_WORTHWHILE_RATIO). Un run qui
  ! passerait RATIO = valeur par defaut reproduirait EXACTEMENT le
  ! comportement actuel.
  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ! === 1) assemblage BLR distribue ===
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM)

  ! === 2) factorisation LU, AVEC pivot confine au bloc diagonal ===
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                       N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, RSS_FACT_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
  TMP_MB = REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR),8)/1024.0D0/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, MEM_ANALYTIC_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)

  ! === 3) second membre B = A_BLR_compresse * Xexact, distribue lignes ET colonnes ===
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(XEX_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)), B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_XEX(A_BLR, GRID, NRHS_LOC, RSTART_MINE, XEX_LOCAL, MLOC)

  BLOCK
    TYPE(MPI_BLR_MATRIX_T) :: A_BLR_FOR_B
    DOUBLE PRECISION :: T_ASM_DUMMY
    INTEGER :: NC_DUMMY, NO_DUMMY
    INTEGER(8) :: RS_DUMMY
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR_FOR_B, T_ASM_DUMMY, &
                                  NC_DUMMY, NO_DUMMY, RS_DUMMY)
    CALL APPLY_MPI_BLR_MATRIX(GRID, A_BLR_FOR_B, NRHS, NRHS_LOC, XEX_LOCAL, B_LOCAL, MLOC)
  END BLOCK

  ! === 4) resolution 2D, AVEC pivot ===
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
  T1 = MPI_WTIME()
  T_SOLVE = T1 - T0

  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

  ERR = -1.0_WP
  IF (MYID0 .EQ. 0) THEN
    ALLOCATE(B_REF(N,NRHS))
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = XEX_FORMULA(GI,R)
      END DO
    END DO
    ERR = SQRT(SUM(ABS(X_GLOBAL-B_REF)**2)) / SQRT(SUM(ABS(B_REF)**2))
    AVG_RANK = REAL(RANK_SUM,WP) / REAL(MAX(N_OFFDIAG,1),WP)

    WRITE(6,'(A,ES9.2,A,F5.2,A,I0)') ' EPS = ', EPS, ' , RATIO = ', RATIO, &
         ' , RMAX = ', MAX(1, INT(RATIO*REAL(NB,WP)) - 1)
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, &
         ' , NRHS=', NRHS, ' , TOL=', TOL
    WRITE(6,'(A,F8.3,A,F8.3,A,F8.3,A)') '   Temps assemblage = ', T_ASM, ' s , factorisation = ', T_FACTOR, &
         ' s , resolution = ', T_SOLVE, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Compression : ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', &
         100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F8.2)') '   Rang moyen = ', AVG_RANK
    WRITE(6,'(A,I0)') '   Nombre de swaps de pivot effectues (total) = ', N_SWAPS
    WRITE(6,'(A,F10.1,A,F10.1,A)') '   Memoire ANALYTIQUE (max/processus) = ', MEM_ANALYTIC_MB, &
         ' Mo  ;  RSS apres factorisation (max/processus) = ', RSS_FACT_MB, ' Mo'
    WRITE(6,'(A,F6.2)') '   Ratio RSS/theorique = ', RSS_FACT_MB/MAX(MEM_ANALYTIC_MB,1.0D-6)
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_calcule - X_exact|| = ', ERR
    DEALLOCATE(B_REF)
  END IF

  DEALLOCATE(B_LOCAL, XEX_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Hash deterministe (I,J) -> [-1,1) -- pas de RNG a etat, recalculable
  ! independamment pour n'importe quel bloc par n'importe quel processus
  ! (necessaire : GET_BLOCK peut etre appelee dans n'importe quel ordre,
  ! par n'importe quel processus, pour n'importe quel sous-bloc).
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  ! A(GI,GJ) = 50 sur l'anti-diagonale CONFINEE AU BLOC DIAGONAL de GI (force
  ! le pivot), sinon noyau de Cauchy amorti + bruit plein-rang EPS*HASH_NOISE
  ! -- EPS reste FIXE pour tout ce benchmark (voir commentaire de tete),
  ! seul RATIO (COMPRESS_WORTHWHILE_RATIO) varie d'un run a l'autre.
  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ, NN)
    INTEGER, INTENT(IN) :: GI, GJ, NN
    REAL(WP) :: DF
    INTEGER :: BSTART0, BEND0, GJ_STRONG
    BSTART0 = ((GI-1)/NB)*NB + 1
    BEND0 = MIN(BSTART0+NB-1, NN)
    GJ_STRONG = BEND0 - (GI-BSTART0)
    IF (GJ .EQ. GJ_STRONG) THEN
      LOWRANK_A = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
    ELSE
      DF = 1.0_WP / (1.0_WP + ABS(REAL(GJ_STRONG-GJ,WP))) + EPS*HASH_NOISE(GI,GJ)
      LOWRANK_A = CMPLX(DF, 0.3_WP*DF, KIND=WP)
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

END PROGRAM BENCH_BLR_MPI_PIVOT_RATIO
