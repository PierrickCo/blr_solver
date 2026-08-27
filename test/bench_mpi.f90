!===============================================================================
! REFERENCE "MPI" : solveur MPI dense avec pivot (mpi_lu_pivot_mod.f90) +
! resolution 2D distribuee complete (mpi_lu_solve_2d_mod.f90 -- B distribue
! PAR LIGNES ET PAR COLONNES sur NRHS, aucune replication residuelle). Fait
! partie du jeu de reference standard a 5 solveurs de ce projet (A/B/MPI/
! C/BLR+MPI, voir le README).
!
! LIMITE IMPORTANTE, mesuree explicitement plutot que masquee : ce solveur
! n'a pas de stockage compresse -- A_LOCAL reste dense sur chaque processus
! (utile comme reference "MPI classique", a comparer a BLR+MPI qui, lui,
! reste compresse de bout en bout -- voir mpi_blr_matrix_mod.f90).
!
! (Les variantes BLR+pivot historiquement comparees ici -- 6e et proposition
! 3 -- ont ete deplacees vers legacy/, ou elles restent utilisables telles
! quelles : voir legacy/test/bench_mpi_full.f90, qui compare aussi
! l'ancienne resolution par rassemblement et la resolution 1D.)
!
! Usage : mpirun -np <PR*PC> ./bench_mpi PR PC [N] [NRHS] [NB] [TOL]
!         (defauts : N=10000 NRHS=100 NB=100 TOL=1e-4 -- TOL inutilise ici,
!         conserve pour la compatibilite de ligne de commande avec les
!         autres benchmarks du projet)
!===============================================================================
PROGRAM BENCH_MPI
  USE MPI_LU_PIVOT_MOD
  USE MPI_LU_SOLVE_2D_MOD
  USE BLOCK_LU_MOD, ONLY: BLOCK_LU_SOLVE
  USE BENCH_UTILS_MOD
  USE OMP_LIB, ONLY: OMP_GET_WTIME
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NRHS, NB, IERR, MYID0
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  INTEGER, PARAMETER :: RGLOBAL = 6

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./bench_mpi PR PC [N] [NRHS] [NB] [TOL]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  N = 10000; NRHS = 100; NB = 100; TOL = 1.0E-4_WP
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL
  END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL REPORT_THREAD_AFFINITY(MPI_COMM_WORLD)

  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES10.3)') ' === BENCHMARK MPI (reference MPI) : grille ', PR, 'x', PC, &
         ' , N=', N, ' , NRHS=', NRHS, ' , NB=', NB, ' , TOL=', TOL
    WRITE(6,*)
  END IF

  CALL RUN_DENSE_PIVOT(N, NRHS, NB)

  CALL MPI_FINALIZE(IERR)

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE BUILD_LOCAL_MATRIX(N, NB, A_LOCAL, MLOC, NLOC, T_ASM)
    INTEGER, INTENT(IN)        :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(OUT)   :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    DOUBLE PRECISION, INTENT(OUT) :: T_ASM
    INTEGER :: I, J, GI, GJ, R
    COMPLEX(WP) :: AVAL
    DOUBLE PRECISION :: T0, T1, T_LOCAL

    T0 = OMP_GET_WTIME()
    DO J = 1, NLOC
      GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
      DO I = 1, MLOC
        GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
        END DO
        IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        A_LOCAL(I,J) = AVAL
      END DO
    END DO
    T1 = OMP_GET_WTIME()
    T_LOCAL = T1 - T0
    CALL MPI_ALLREDUCE(T_LOCAL, T_ASM, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
  END SUBROUTINE BUILD_LOCAL_MATRIX

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

  FUNCTION MAX_RSS_MB() RESULT(MB)
    DOUBLE PRECISION :: MB, LOCAL_MB
    INTEGER(8) :: KB
    KB = GET_RSS_KB()
    LOCAL_MB = REAL(KB,8) / 1024.0D0
    CALL MPI_ALLREDUCE(LOCAL_MB, MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
  END FUNCTION MAX_RSS_MB

  SUBROUTINE RANDOM_COMPLEX(M, K, X)
    INTEGER, INTENT(IN) :: M, K
    COMPLEX(WP), INTENT(OUT) :: X(M,K)
    REAL(WP), ALLOCATABLE :: RP(:,:), IP(:,:)
    ALLOCATE(RP(M,K), IP(M,K))
    CALL RANDOM_NUMBER(RP)
    CALL RANDOM_NUMBER(IP)
    X = CMPLX(RP-0.5_WP, IP-0.5_WP, KIND=WP)
    DEALLOCATE(RP, IP)
  END SUBROUTINE RANDOM_COMPLEX

  FUNCTION FROB_ERR(N, NRHS, XCALC, XEX) RESULT(ERR)
    INTEGER, INTENT(IN) :: N, NRHS
    COMPLEX(WP), INTENT(IN) :: XCALC(N,NRHS), XEX(N,NRHS)
    REAL(WP) :: ERR
    ERR = SQRT(SUM(ABS(XCALC-XEX)**2)) / SQRT(SUM(ABS(XEX)**2))
  END FUNCTION FROB_ERR

  SUBROUTINE REPORT(LABEL, T_ASM, T_COMPRESS, T_BCAST, T_FACT_PURE, RSS_FACT_MB, T_GATHER, RSS_GATHER_MB, T_SOLVE, ERR)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    DOUBLE PRECISION, INTENT(IN) :: T_ASM, T_COMPRESS, T_BCAST, T_FACT_PURE, RSS_FACT_MB, T_GATHER, RSS_GATHER_MB, T_SOLVE
    REAL(WP), INTENT(IN) :: ERR
    IF (MYID0 .NE. 0) RETURN
    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'
    WRITE(6,'(A,F8.3,A)')  '   Temps assemblage (max/processus)      = ', T_ASM, ' s'
    WRITE(6,'(A,F8.3,A)')  '   Temps compression basse-rang (max/processus) = ', T_COMPRESS, ' s'
    WRITE(6,'(A,F8.3,A)')  '   Temps diffusion des panneaux L/U (max/processus) = ', T_BCAST, ' s'
    WRITE(6,'(A,F8.3,A)')  '   Temps factorisation (hors compression, distribuee) = ', T_FACT_PURE, ' s'
    WRITE(6,'(A,F10.1,A)') '   Memoire PAR PROCESSUS pendant factor. = ', RSS_FACT_MB, ' Mo (max sur tous)'
    WRITE(6,'(A)') '   [ANCIEN : resolution via rassemblement]'
    WRITE(6,'(A,F8.3,A)')  '     Temps rassemblement (GATHER_2D)       = ', T_GATHER, ' s'
    WRITE(6,'(A,F10.1,A)') '     Memoire RANG 0 apres rassemblement    = ', RSS_GATHER_MB, ' Mo'
    WRITE(6,'(A,F8.3,A)')  '     Temps resolution (sequentiel, rang 0) = ', T_SOLVE, ' s'
    WRITE(6,'(A,ES12.4)')  '     Erreur relative                       = ', ERR
  END SUBROUTINE REPORT

  !-----------------------------------------------------------------------
  ! Resolution DISTRIBUEE 2D (MPI_LU_SOLVE_DISTRIBUTED_2D) : B distribue
  ! PAR LIGNES (comme A) ET PAR COLONNES sur NRHS (RHS_NUMROC/PC) -- AUCUNE
  ! replication residuelle. C'est LA resolution de la reference "MPI".
  !-----------------------------------------------------------------------
  SUBROUTINE RUN_DISTRIBUTED_SOLVE_2D(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    INTEGER, INTENT(IN) :: N, NB, NRHS, MLOC, NLOC
    COMPLEX(WP), INTENT(IN) :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    INTEGER, INTENT(IN) :: IPIV_GLOBAL(N)
    COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:), A_REF(:,:), B_REF(:,:)
    INTEGER, ALLOCATABLE :: IPIV_REF(:)
    DOUBLE PRECISION :: T0, T1, T_SOLVE, RSS_MB, B_OLD_MB, B_NEW_MB
    REAL(WP) :: ERR, NORMB
    INTEGER :: GI, GJ, R, INFO_REF, NRHS_LOC, RSTART_MINE

    CALL RHS_LOCAL_RANGE(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
    ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
    DO R = 1, NRHS_LOC
      GJ = RSTART_MINE + R - 1
      DO GI = 1, MLOC
        BLOCK
          INTEGER :: GIROW
          GIROW = INDXL2G(GI, NB, GRID%MYROW, GRID%PR)
          B_LOCAL(GI,R) = CMPLX(REAL(GIROW,WP), -0.5_WP*REAL(GIROW,WP)+REAL(GJ,WP), KIND=WP)
        END BLOCK
      END DO
    END DO

    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = OMP_GET_WTIME()
    CALL MPI_LU_SOLVE_DISTRIBUTED_2D(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                                      NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
    T1 = OMP_GET_WTIME()
    T_SOLVE = T1 - T0
    RSS_MB = MAX_RSS_MB()

    ALLOCATE(X_GLOBAL(N,NRHS))
    CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

    ERR = 0.0_WP
    IF (MYID0 .EQ. 0) THEN
      ALLOCATE(A_REF(N,N), B_REF(N,NRHS), IPIV_REF(N))
      DO GJ = 1, N
        DO GI = 1, N
          A_REF(GI,GJ) = LOWRANK_A(GI,GJ,N)
        END DO
      END DO
      DO R = 1, NRHS
        DO GI = 1, N
          B_REF(GI,R) = CMPLX(REAL(GI,WP), -0.5_WP*REAL(GI,WP)+REAL(R,WP), KIND=WP)
        END DO
      END DO
      CALL CGESV(N, NRHS, A_REF, N, IPIV_REF, B_REF, N, INFO_REF)
      NORMB = SQRT(SUM(ABS(B_REF)**2))
      ERR = SQRT(SUM(ABS(X_GLOBAL-B_REF)**2)) / NORMB
      DEALLOCATE(A_REF, B_REF, IPIV_REF)
    END IF

    B_OLD_MB = 8.0D0*REAL(N,8)*REAL(NRHS,8)/1024.0D0/1024.0D0
    B_NEW_MB = 8.0D0*REAL(MAX(MLOC,1),8)*REAL(MAX(NRHS_LOC,1),8)/1024.0D0/1024.0D0

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') '   [2D : resolution DISTRIBUEE, B distribue lignes ET colonnes (PRxPC)]'
      WRITE(6,'(A,F8.3,A)')  '     Temps resolution (2D)                 = ', T_SOLVE, ' s'
      WRITE(6,'(A,F10.1,A)') '     Memoire PAR PROCESSUS pendant solve   = ', RSS_MB, ' Mo (max sur tous)'
      WRITE(6,'(A,ES12.4)')  '     Erreur relative                       = ', ERR
      WRITE(6,'(A,F8.1,A,F8.1,A)') '     B seul : ', B_OLD_MB, ' Mo (replication totale) vs ', &
           B_NEW_MB, ' Mo (distribution 2D, ce processus)'
      WRITE(6,*)
    END IF
    DEALLOCATE(B_LOCAL, X_GLOBAL)
  END SUBROUTINE RUN_DISTRIBUTED_SOLVE_2D

  FUNCTION LOWRANK_A(GI, GJ, N) RESULT(AVAL)
    INTEGER, INTENT(IN) :: GI, GJ, N
    COMPLEX(WP) :: AVAL
    INTEGER :: R
    AVAL = (0.0_WP,0.0_WP)
    DO R = 1, RGLOBAL
      AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
    END DO
    IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
  END FUNCTION LOWRANK_A

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_DENSE_PIVOT(N, NRHS, NB)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_FACT(:,:), XEX(:,:), B(:,:)
    INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: MLOC, NLOC, INFO
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT, T_GATHER, T_SOLVE, RSS_FACT, RSS_GATHER
    REAL(WP) :: ERR

    MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
    NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)
    ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
    CALL BUILD_LOCAL_MATRIX(N, NB, A_LOCAL, MLOC, NLOC, T_ASM)

    ALLOCATE(IPIV_GLOBAL(N))
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = OMP_GET_WTIME()
    CALL MPI_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, INFO)
    T1 = OMP_GET_WTIME()
    T_FACT = T1 - T0
    RSS_FACT = MAX_RSS_MB()

    ALLOCATE(A_FACT(N,N))
    T0 = OMP_GET_WTIME()
    CALL GATHER_2D(GRID, N, NB, A_LOCAL, MAX(MLOC,1), MAX(NLOC,1), 0, A_FACT)
    T1 = OMP_GET_WTIME()
    T_GATHER = T1 - T0
    RSS_GATHER = MAX_RSS_MB()

    T_SOLVE = 0.0D0; ERR = 0.0_WP
    IF (MYID0 .EQ. 0) THEN
      ALLOCATE(XEX(N,NRHS), B(N,NRHS))
      CALL RANDOM_COMPLEX(N, NRHS, XEX)
      ! B = A_original * Xexact -- A_original recalculee (A_FACT est deja
      ! factorisee) via la meme formule, comme dans tous les tests MPI
      CALL BUILD_REFERENCE_AND_RHS(N, NRHS, XEX, B)
      T0 = OMP_GET_WTIME()
      CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_FACT, N, IPIV_GLOBAL, B, N)
      T1 = OMP_GET_WTIME()
      T_SOLVE = T1 - T0
      ERR = FROB_ERR(N, NRHS, B, XEX)
      DEALLOCATE(XEX, B)
    END IF

    ! Pas de compression pour ce solveur dense : temps de compression nul
    ! par construction.
    CALL REPORT('MPI) Dense avec pivot + resolution 2D', T_ASM, 0.0D0, 0.0D0, T_FACT, RSS_FACT, T_GATHER, RSS_GATHER, T_SOLVE, ERR)
    ! A_FACT (rassemblee, ~N^2 octets sur le rang 0) liberee AVANT la
    ! mesure memoire de la resolution distribuee -- sinon celle-ci resterait
    ! contaminee par cette allocation encore residente, faussant la
    ! comparaison (voir note de tete de fichier).
    DEALLOCATE(A_FACT)
    CALL RUN_DISTRIBUTED_SOLVE_2D(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    DEALLOCATE(A_LOCAL, IPIV_GLOBAL)
  END SUBROUTINE RUN_DENSE_PIVOT

  !-----------------------------------------------------------------------
  ! Reconstruit la matrice ORIGINALE (meme formule) et B=A_originale*Xexact
  ! -- reference INDEPENDANTE du resultat factorise/gathered, meme
  ! methodologie que tous les tests MPI precedents.
  !-----------------------------------------------------------------------
  SUBROUTINE BUILD_REFERENCE_AND_RHS(N, NRHS, XEX, B)
    INTEGER, INTENT(IN) :: N, NRHS
    COMPLEX(WP), INTENT(IN)  :: XEX(N,NRHS)
    COMPLEX(WP), INTENT(OUT) :: B(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: A_REF(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: GI, GJ, R
    COMPLEX(WP) :: AVAL

    ALLOCATE(A_REF(N,N))
    DO GJ = 1, N
      DO GI = 1, N
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
        END DO
        IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        A_REF(GI,GJ) = AVAL
      END DO
    END DO
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_REF, N, XEX, N, CZERO, B, N)
    DEALLOCATE(A_REF)
  END SUBROUTINE BUILD_REFERENCE_AND_RHS

END PROGRAM BENCH_MPI
