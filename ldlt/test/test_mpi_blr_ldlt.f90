!===============================================================================
! Test de validation du LDL^T BLR+MPI (Bunch-Kaufman confine) sur la
! matrice complexe symetrique de sym_testmat_mod.f90.
!
! Trois verifications INDEPENDANTES de ce que calcule le solveur :
!   1) Erreur ||X-Xexact||/||Xexact||, avec B = A*Xexact calcule par produit
!      BRUT a partir de la formule de A (jamais a partir des blocs
!      compresses) -- l'erreur inclut donc la compression.
!   2) (N <= 4000) Reference dense independante : CGESV (LU, autre
!      algorithme) sur A et B complets au rang 0 -- ||X-X_CGESV||/||X_CGESV||.
!   3) (grille 1x1 seulement, N <= 3000) Reconstruction directe du facteur :
!      ||Q*A*Q^T - L*D*L^T|| / ||A|| a partir de IPIV/DIAGP/OFFBLK/D --
!      isole la factorisation de la resolution.
!
! Usage : mpirun -np PR*PC ./test_mpi_blr_ldlt PR PC [N] [NB] [NRHS] [TOL] [STRONG] [DIAG]
!         (defauts 2000 200 8 1e-4 50 0.01 -- STRONG/DIAG : voir
!         SYM_TESTMAT_SETUP, choisissent le type de pivots 1x1/2x2)
!===============================================================================
PROGRAM TEST_MPI_BLR_LDLT
  USE MPI_BLR_LDLT_MOD
  USE SYM_TESTMAT_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  TYPE(MPI_BLR_SYM_T) :: A
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0, INFO, N_COMPRESSED, N_OFFDIAG, N_SWAPS, N_2X2
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE
  INTEGER(8) :: RANK_SUM
  REAL(WP) :: TOL, ERR, ERR_REF, ERR_FACT, ERR_REF_EXACT, STRONG, DIAG
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T0, T_SOLVE
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), XG(:,:)
  LOGICAL :: PASS

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)
  N = 2000; NB = 200; NRHS = 8; TOL = 1.0E-4_WP; STRONG = 50.0_WP; DIAG = 0.01_WP
  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_ldlt PR PC [N] [NB] [NRHS] [TOL]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN; CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N;    END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN; CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB;   END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN; CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN; CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL;  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN; CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) STRONG; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN; CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) DIAG; END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL SYM_TESTMAT_SETUP(N, NB, STRONG=STRONG, DIAG=DIAG)
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F6.1,A,F7.3)') ' LDLt BLR+MPI : grille ', PR, 'x', PC, &
       ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, ' , STRONG=', STRONG, ' , DIAG=', DIAG

  CALL MPI_ASSEMBLE_BLR_SYM(GRID, N, NB, TOL, SYM_GET_BLOCK, A, T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
  CALL MPI_BLR_SYM_LDLT_FACTOR(GRID, A, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, N_2X2, T_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_SYM_LDLT_FACTOR, INFO=', INFO
    IF (MYID0.EQ.0) WRITE(6,*) '=== TEST_MPI_BLR_LDLT : ECHEC ==='
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ERR_FACT = -1.0_WP
  IF (GRID%NPROCS .EQ. 1 .AND. N .LE. 3000) ERR_FACT = CHECK_FACTOR()

  MLOC = MPI_BLR_SYM_MLOC(A, GRID)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_B_EXACT()

  CALL MPI_BARRIER(GRID%COMM, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, A, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
  T_SOLVE = MPI_WTIME() - T0

  ALLOCATE(XG(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, XG)
  ERR_REF = -1.0_WP; ERR_REF_EXACT = -1.0_WP
  IF (MYID0 .EQ. 0) THEN
    ERR = REL_ERR_XEX(XG)
    IF (N .LE. 4000) CALL DENSE_REFERENCE(XG, ERR_REF, ERR_REF_EXACT)
    WRITE(6,'(A,F9.3,A,F9.3,A,F9.3,A)') '  Temps : assemblage ', T_ASM, ' s , factorisation ', T_FACTOR, &
         ' s , resolution ', T_SOLVE, ' s'
    WRITE(6,'(A,I0,A,I0,A,I0,A)') '  Blocs L compresses : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' , rang moyen ', INT(REAL(RANK_SUM,WP)/REAL(MAX(N_OFFDIAG,1),WP)), ' (inclut les blocs denses)'
    WRITE(6,'(A,I0,A,I0)') '  Echanges de pivot : ', N_SWAPS, ' , pivots 2x2 : ', N_2X2
    IF (ERR_FACT .GE. 0.0_WP) WRITE(6,'(A,ES12.4)') '  ||Q*A*Q^T - L*D*L^T|| / ||A||      = ', ERR_FACT
    WRITE(6,'(A,ES12.4)') '  ||X-Xexact|| / ||Xexact||           = ', ERR
    IF (ERR_REF .GE. 0.0_WP) THEN
      WRITE(6,'(A,ES12.4)') '  ||X-X_CGESV|| / ||X_CGESV|| (dense) = ', ERR_REF
      WRITE(6,'(A,ES12.4)') '  ||X_CGESV-Xexact|| / ||Xexact||     = ', ERR_REF_EXACT
    END IF
    PASS = (ERR .LT. 1.0E-3_WP)
    IF (ERR_FACT .GE. 0.0_WP) PASS = PASS .AND. (ERR_FACT .LT. 1.0E-3_WP)
    IF (ERR_REF .GE. 0.0_WP)  PASS = PASS .AND. (ERR_REF .LT. 1.0E-3_WP)
    IF (PASS) THEN
      WRITE(6,*) '=== TEST_MPI_BLR_LDLT : OK ==='
    ELSE
      WRITE(6,*) '=== TEST_MPI_BLR_LDLT : ECHEC ==='
    END IF
  END IF
  DEALLOCATE(B_LOCAL, XG)
  CALL MPI_BLR_SYM_FREE(A)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! B_LOCAL = A*Xexact (formule de A, produit brut par blocs de lignes).
  SUBROUTINE FILL_B_EXACT()
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: AROWS(:,:), XEX(:,:)
    INTEGER :: li, I, BSZ, RS, GI, R
    IF (NRHS_LOC .LE. 0) RETURN
    ALLOCATE(XEX(N,NRHS_LOC))
    DO R = 1, NRHS_LOC
      DO GI = 1, N
        XEX(GI,R) = SYM_XEX(GI, RSTART_MINE+R-1)
      END DO
    END DO
    RS = 1
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      BSZ = A%BEND(I) - A%BSTART(I) + 1
      ALLOCATE(AROWS(BSZ,N))
      CALL SYM_GET_BLOCK(A%BSTART(I), A%BEND(I), 1, N, AROWS)
      CALL CGEMM('N', 'N', BSZ, NRHS_LOC, N, CONE, AROWS, BSZ, XEX, N, CZERO, B_LOCAL(RS,1), MAX(MLOC,1))
      DEALLOCATE(AROWS)
      RS = RS + BSZ
    END DO
    DEALLOCATE(XEX)
  END SUBROUTINE FILL_B_EXACT

  REAL(WP) FUNCTION REL_ERR_XEX(X)
    COMPLEX(WP), INTENT(IN) :: X(N,NRHS)
    REAL(WP) :: NUM, DEN
    INTEGER :: GI, R
    NUM = 0.0_WP; DEN = 0.0_WP
    DO R = 1, NRHS
      DO GI = 1, N
        NUM = NUM + ABS(X(GI,R) - SYM_XEX(GI,R))**2
        DEN = DEN + ABS(SYM_XEX(GI,R))**2
      END DO
    END DO
    REL_ERR_XEX = SQRT(NUM/DEN)
  END FUNCTION REL_ERR_XEX

  ! Reference dense independante (rang 0) : CGESV sur A, B complets.
  SUBROUTINE DENSE_REFERENCE(X, ERR_X, ERR_EXACT)
    COMPLEX(WP), INTENT(IN) :: X(N,NRHS)
    REAL(WP), INTENT(OUT)   :: ERR_X, ERR_EXACT
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: AF(:,:), XEX(:,:), BF(:,:)
    INTEGER, ALLOCATABLE :: IP(:)
    INTEGER :: GI, R, INFO2
    ALLOCATE(AF(N,N), XEX(N,NRHS), BF(N,NRHS), IP(N))
    CALL SYM_GET_BLOCK(1, N, 1, N, AF)
    DO R = 1, NRHS
      DO GI = 1, N
        XEX(GI,R) = SYM_XEX(GI,R)
      END DO
    END DO
    CALL CGEMM('N', 'N', N, NRHS, N, CONE, AF, N, XEX, N, CZERO, BF, N)
    CALL CGESV(N, NRHS, AF, N, IP, BF, N, INFO2)
    ERR_X = SQRT(SUM(ABS(X-BF)**2) / SUM(ABS(BF)**2))
    ERR_EXACT = SQRT(SUM(ABS(BF-XEX)**2) / SUM(ABS(XEX)**2))
    DEALLOCATE(AF, XEX, BF, IP)
  END SUBROUTINE DENSE_REFERENCE

  ! Grille 1x1 : reconstruit L, D et Q en dense et mesure
  ! ||Q*A*Q^T - L*D*L^T||_F / ||A||_F (A = formule exacte, non compressee).
  REAL(WP) FUNCTION CHECK_FACTOR()
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: AQ(:,:), LF(:,:), DF(:,:), T(:,:), R(:,:), DK(:,:)
    INTEGER :: K, I, KW, P, II, JJ, B0
    ALLOCATE(AQ(N,N), LF(N,N), DF(N,N), T(N,N), R(N,N))
    CALL SYM_GET_BLOCK(1, N, 1, N, AQ)
    LF = CZERO; DF = CZERO
    DO K = 1, A%NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      B0 = A%BSTART(K) - 1
      ! Q*A*Q^T : echanges symetriques du bloc K, dans l'ordre
      DO I = 1, KW
        P = ABS(A%IPIV(I,K))
        IF (P .NE. I) THEN
          CALL CSWAP(N, AQ(B0+I,1), N, AQ(B0+P,1), N)
          CALL CSWAP(N, AQ(1,B0+I), 1, AQ(1,B0+P), 1)
        END IF
      END DO
      ! L_KK unitaire inferieure et D_K
      ALLOCATE(DK(KW,KW))
      DK = CZERO
      CALL SYM_UNPACK_LOWER(KW, A%DIAGP(K)%AP, DK, KW)
      DO JJ = 1, KW
        LF(B0+JJ,B0+JJ) = CONE
        DO II = JJ+1, KW
          LF(B0+II,B0+JJ) = DK(II,JJ)
        END DO
        DF(B0+JJ,B0+JJ) = A%DDIAG(JJ,K)
      END DO
      DEALLOCATE(DK)
      JJ = 1
      DO WHILE (JJ .LE. KW)
        IF (A%DSUB(JJ,K) .NE. CZERO .AND. JJ .LT. KW) THEN
          DF(B0+JJ+1,B0+JJ) = A%DSUB(JJ,K); DF(B0+JJ,B0+JJ+1) = A%DSUB(JJ,K)
          LF(B0+JJ+1,B0+JJ) = CZERO
          JJ = JJ + 2
        ELSE
          JJ = JJ + 1
        END IF
      END DO
      ! L(I,K), I>K (grille 1x1 : indices locaux = globaux)
      DO I = K+1, A%NBLK
        CALL EXPAND_BLOCK(A%OFFBLK(I,K), LF(A%BSTART(I),A%BSTART(K)), N)
      END DO
    END DO
    CALL CGEMM('N', 'N', N, N, N, CONE, LF, N, DF, N, CZERO, T, N)
    CALL CGEMM('N', 'T', N, N, N, CONE, T, N, LF, N, CZERO, R, N)
    CALL SYM_GET_BLOCK(1, N, 1, N, T)
    CHECK_FACTOR = SQRT(SUM(ABS(AQ-R)**2) / SUM(ABS(T)**2))
    DEALLOCATE(AQ, LF, DF, T, R)
  END FUNCTION CHECK_FACTOR

END PROGRAM TEST_MPI_BLR_LDLT
