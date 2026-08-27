!===============================================================================
! BENCHMARK MPI : meme systeme (N=10000, NRHS=100, meme matrice EXACTEMENT
! que bench_seq.f90) sur trois solveurs distribues :
!   1) Dense avec pivot (mpi_lu_pivot_mod.f90, 6d) -- reference distribuee.
!   2) BLR avec pivot, diffusion DENSE-puis-compression (mpi_blr_lu_pivot_
!      mod.f90, 6e).
!   3) BLR avec pivot, diffusion COMPRESSEE (mpi_blr_lu_pivot_bcast_mod.f90,
!      proposition 3) -- LA comparaison demandee contre le cas 2.
!
! LIMITE IMPORTANTE, mesuree explicitement plutot que masquee : aucune de
! ces trois versions n'a de RESOLUTION distribuee -- resoudre necessite
! aujourd'hui de RASSEMBLER (GATHER_2D) la matrice factorisee complete sur
! le rang 0 (cout memoire O(N^2) sur CE processus, cout de communication
! en etoile), puis d'appeler BLOCK_LU_SOLVE (etape 1) sequentiellement
! la-bas. Ce cout est mesure et rapporte SEPAREMENT du temps de
! factorisation, avec la memoire pic sur le rang 0 clairement etiquetee --
! c'est le cout REEL qu'un utilisateur paierait aujourd'hui, pas quelque
! chose a cacher.
!
! MEMOIRE PAR PROCESSUS PENDANT LA FACTORISATION (avant rassemblement) :
! max du RSS sur tous les processus (MPI_ALLREDUCE MAX) -- c'est LA
! comparaison utile pour la proposition 2 eventuellement etendue au
! distribue (pas encore fait, voir memoire du projet) : actuellement les
! 3 cas MPI devraient montrer une memoire par processus SIMILAIRE (aucun
! ne stocke de facteurs compresses de facon persistante) -- la
! proposition 3 ameliore la COMMUNICATION (BYTES_SENT, rapporte pour le
! cas 3), pas la memoire, contrairement a la proposition 2 (sequentielle
! uniquement, voir bench_seq.f90).
!
! Usage : mpirun -np <PR*PC> ./bench_mpi PR PC [N] [NRHS] [NB] [TOL]
!         (defauts : N=10000 NRHS=100 NB=100 TOL=1e-4)
!===============================================================================
PROGRAM BENCH_MPI
  USE MPI_LU_PIVOT_MOD
  USE MPI_BLR_LU_PIVOT_MOD
  USE MPI_BLR_LU_PIVOT_BCAST_MOD
  USE MPI_LU_SOLVE_MOD
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
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES10.3)') ' === BENCHMARK MPI : grille ', PR, 'x', PC, &
         ' , N=', N, ' , NRHS=', NRHS, ' , NB=', NB, ' , TOL=', TOL
    WRITE(6,*)
  END IF

  CALL RUN_DENSE_PIVOT(N, NRHS, NB)
  CALL RUN_BLR_PIVOT(N, NRHS, NB, TOL)
  CALL RUN_BLR_PIVOT_BCAST(N, NRHS, NB, TOL)

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
  ! Resolution DISTRIBUEE (MPI_LU_SOLVE_DISTRIBUTED) sur les MEMES facteurs
  ! -- aucun rassemblement. B DISTRIBUE PAR LIGNES (comme A), replique
  ! seulement au sein d'une ligne de grille -- PAS replique sur toute la
  ! grille (evolution demandee par l'utilisateur pour un NRHS pouvant
  ! devenir tres grand). Rapporte AUSSI la comparaison analytique contre
  ! l'ancienne replication totale, pour un NRHS realiste ou tres grand.
  !-----------------------------------------------------------------------
  SUBROUTINE RUN_DISTRIBUTED_SOLVE(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    INTEGER, INTENT(IN) :: N, NB, NRHS, MLOC, NLOC
    COMPLEX(WP), INTENT(IN) :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    INTEGER, INTENT(IN) :: IPIV_GLOBAL(N)
    COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:), A_REF(:,:), B_REF(:,:)
    INTEGER, ALLOCATABLE :: IPIV_REF(:)
    DOUBLE PRECISION :: T0, T1, T_SOLVE, RSS_MB, B_OLD_MB, B_NEW_MB
    REAL(WP) :: ERR, NORMB
    INTEGER :: GI, GJ, R, INFO_REF

    ALLOCATE(B_LOCAL(MAX(MLOC,1),NRHS))
    DO R = 1, NRHS
      DO GI = 1, MLOC
        GJ = INDXL2G(GI, NB, GRID%MYROW, GRID%PR)
        B_LOCAL(GI,R) = CMPLX(REAL(GJ,WP), -0.5_WP*REAL(GJ,WP)+REAL(R,WP), KIND=WP)
      END DO
    END DO

    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = OMP_GET_WTIME()
    CALL MPI_LU_SOLVE_DISTRIBUTED(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, NRHS, B_LOCAL, MAX(MLOC,1))
    T1 = OMP_GET_WTIME()
    T_SOLVE = T1 - T0
    RSS_MB = MAX_RSS_MB()

    ALLOCATE(X_GLOBAL(N,NRHS))
    DO R = 1, NRHS
      CALL GATHER_ROW_VECTOR(GRID, N, NB, B_LOCAL(:,R), MAX(MLOC,1), 0, X_GLOBAL(:,R))
    END DO

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

    ! B_OLD_MB : ce que l'ANCIENNE replication totale aurait coute PAR
    ! PROCESSUS (N x NRHS partout) ; B_NEW_MB : ce que la nouvelle
    ! distribution par lignes coute REELLEMENT sur CE processus.
    B_OLD_MB = 8.0D0*REAL(N,8)*REAL(NRHS,8)/1024.0D0/1024.0D0
    B_NEW_MB = 8.0D0*REAL(MAX(MLOC,1),8)*REAL(NRHS,8)/1024.0D0/1024.0D0

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') '   [NOUVEAU : resolution DISTRIBUEE, B distribue par lignes]'
      WRITE(6,'(A,F8.3,A)')  '     Temps resolution (distribuee)         = ', T_SOLVE, ' s'
      WRITE(6,'(A,F10.1,A)') '     Memoire PAR PROCESSUS pendant solve   = ', RSS_MB, ' Mo (max sur tous)'
      WRITE(6,'(A,ES12.4)')  '     Erreur relative                       = ', ERR
      WRITE(6,'(A,F8.1,A,F8.1,A)') '     B seul : ', B_OLD_MB, ' Mo (ancienne replication totale) vs ', &
           B_NEW_MB, ' Mo (nouvelle distribution par lignes, ce processus)'
      WRITE(6,*)
    END IF
    DEALLOCATE(B_LOCAL, X_GLOBAL)
  END SUBROUTINE RUN_DISTRIBUTED_SOLVE

  !-----------------------------------------------------------------------
  ! Resolution DISTRIBUEE 2D (MPI_LU_SOLVE_DISTRIBUTED_2D) : B distribue
  ! PAR LIGNES (comme A) ET PAR COLONNES sur NRHS (RHS_NUMROC/PC) --
  ! AUCUNE replication residuelle (contrairement a RUN_DISTRIBUTED_SOLVE
  ! ci-dessus, qui replique encore B au sein d'une ligne). Compare les 3
  ! approches (rassemblement, 1D-lignes, 2D complet) sur les MEMES facteurs.
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

    ! B_OLD_MB : replication totale (reference historique) ; B_NEW_MB : ce
    ! que la distribution 2D COMPLETE coute REELLEMENT sur CE processus.
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

    ! Pas de compression pour la version dense (6d) : temps de compression
    ! nul par construction.
    CALL REPORT('1) Dense avec pivot (6d)', T_ASM, 0.0D0, 0.0D0, T_FACT, RSS_FACT, T_GATHER, RSS_GATHER, T_SOLVE, ERR)
    ! A_FACT (rassemblee, ~N^2 octets sur le rang 0) liberee AVANT la
    ! mesure memoire de la resolution distribuee -- sinon celle-ci resterait
    ! contaminee par cette allocation encore residente, faussant la
    ! comparaison (voir note de tete de fichier).
    DEALLOCATE(A_FACT)
    CALL RUN_DISTRIBUTED_SOLVE(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    CALL RUN_DISTRIBUTED_SOLVE_2D(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    DEALLOCATE(A_LOCAL, IPIV_GLOBAL)
  END SUBROUTINE RUN_DENSE_PIVOT

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_BLR_PIVOT(N, NRHS, NB, TOL)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    REAL(WP), INTENT(IN) :: TOL
    COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_FACT(:,:), XEX(:,:), B(:,:)
    INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:)
    INTEGER :: MLOC, NLOC, INFO, N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8) :: RANK_SUM
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT_TOTAL, T_COMPRESS, T_BCAST, T_GATHER, T_SOLVE, RSS_FACT, RSS_GATHER
    REAL(WP) :: ERR

    MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
    NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)
    ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
    CALL BUILD_LOCAL_MATRIX(N, NB, A_LOCAL, MLOC, NLOC, T_ASM)

    ALLOCATE(IPIV_GLOBAL(N))
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = OMP_GET_WTIME()
    CALL MPI_BLR_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, TOL, IPIV_GLOBAL, INFO, &
                                  N_COMPRESSED, N_TOTAL_PANELS, RANK_SUM, T_COMPRESS=T_COMPRESS, T_BCAST=T_BCAST)
    T1 = OMP_GET_WTIME()
    T_FACT_TOTAL = T1 - T0
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
      CALL BUILD_REFERENCE_AND_RHS(N, NRHS, XEX, B)
      T0 = OMP_GET_WTIME()
      CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_FACT, N, IPIV_GLOBAL, B, N)
      T1 = OMP_GET_WTIME()
      T_SOLVE = T1 - T0
      ERR = FROB_ERR(N, NRHS, B, XEX)
      DEALLOCATE(XEX, B)
    END IF

    CALL REPORT('2) BLR avec pivot, diffusion DENSE (6e)', T_ASM, T_COMPRESS, T_BCAST, T_FACT_TOTAL-T_COMPRESS, &
                RSS_FACT, T_GATHER, RSS_GATHER, T_SOLVE, ERR)
    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   (compression : ', N_COMPRESSED, ' / ', &
           N_TOTAL_PANELS, ' panneaux, ', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_TOTAL_PANELS,1),WP), ' %)'
      WRITE(6,'(A,F8.2,A,I0,A)') '   (rang moyen = ', REAL(RANK_SUM,WP)/REAL(MAX(N_TOTAL_PANELS,1),WP), &
           ' -- panneaux locaux, taille variable selon le processus, NB=', NB, ')'
    END IF
    DEALLOCATE(A_FACT)
    CALL RUN_DISTRIBUTED_SOLVE(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    DEALLOCATE(A_LOCAL, IPIV_GLOBAL)
  END SUBROUTINE RUN_BLR_PIVOT

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_BLR_PIVOT_BCAST(N, NRHS, NB, TOL)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    REAL(WP), INTENT(IN) :: TOL
    COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_FACT(:,:), XEX(:,:), B(:,:)
    INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:)
    INTEGER :: MLOC, NLOC, INFO, N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8) :: RANK_SUM, BYTES_SENT, BYTES_DENSE_EQUIV
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT_TOTAL, T_COMPRESS, T_BCAST, T_GATHER, T_SOLVE, RSS_FACT, RSS_GATHER
    REAL(WP) :: ERR

    MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
    NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)
    ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
    CALL BUILD_LOCAL_MATRIX(N, NB, A_LOCAL, MLOC, NLOC, T_ASM)

    ALLOCATE(IPIV_GLOBAL(N))
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = OMP_GET_WTIME()
    CALL MPI_BLR_LU_PIVOT_BCAST_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, TOL, IPIV_GLOBAL, INFO, &
                                        N_COMPRESSED, N_TOTAL_PANELS, RANK_SUM, BYTES_SENT, BYTES_DENSE_EQUIV, &
                                        T_COMPRESS=T_COMPRESS, T_BCAST=T_BCAST)
    T1 = OMP_GET_WTIME()
    T_FACT_TOTAL = T1 - T0
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
      CALL BUILD_REFERENCE_AND_RHS(N, NRHS, XEX, B)
      T0 = OMP_GET_WTIME()
      CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_FACT, N, IPIV_GLOBAL, B, N)
      T1 = OMP_GET_WTIME()
      T_SOLVE = T1 - T0
      ERR = FROB_ERR(N, NRHS, B, XEX)
      DEALLOCATE(XEX, B)
    END IF

    CALL REPORT('3) BLR avec pivot, diffusion COMPRESSEE (proposition 3)', &
                T_ASM, T_COMPRESS, T_BCAST, T_FACT_TOTAL-T_COMPRESS, RSS_FACT, T_GATHER, RSS_GATHER, T_SOLVE, ERR)
    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   (compression : ', N_COMPRESSED, ' / ', &
           N_TOTAL_PANELS, ' panneaux, ', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_TOTAL_PANELS,1),WP), ' %)'
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   (octets diffuses : ', BYTES_SENT, ' vs equivalent dense (6e) ', &
           BYTES_DENSE_EQUIV, ' -> ', 100.0_WP*REAL(BYTES_SENT,WP)/REAL(MAX(BYTES_DENSE_EQUIV,1_8),WP), ' % du dense)'
      WRITE(6,'(A,F8.2,A,I0,A)') '   (rang moyen = ', REAL(RANK_SUM,WP)/REAL(MAX(N_TOTAL_PANELS,1),WP), &
           ' -- panneaux locaux, taille variable selon le processus, NB=', NB, ')'
    END IF
    DEALLOCATE(A_FACT)
    CALL RUN_DISTRIBUTED_SOLVE(N, NB, NRHS, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL)
    DEALLOCATE(A_LOCAL, IPIV_GLOBAL)
  END SUBROUTINE RUN_BLR_PIVOT_BCAST

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
