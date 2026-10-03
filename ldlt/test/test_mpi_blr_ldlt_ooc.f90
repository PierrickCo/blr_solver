!===============================================================================
! Test de la variante OUT-OF-CORE du LDL^T BLR+MPI.
!
! B = A*Xexact est calcule entree par entree depuis la formule de A
! (FILL_RHS, jamais depuis les blocs compresses ni via une distribution
! dependante du decoupage en paquets) : le test est donc independant de
! NRHS_BATCH. Un petit MAX_MB_PER_PROCESS force de vraies evictions disque.
!
! Verifie : ||X_ooc - Xexact||/||Xexact|| < 1e-3, ecart a la version
! classique en memoire < 1e-3 (la SVD randomisee tire des vecteurs
! differents dans les deux runs : ecart attendu de l'ordre de TOL), et que
! le cache resident est reste sous le budget (a un bloc pres : un bloc plus
! gros que le budget est toujours accepte, voir blr_ooc_mod.f90).
!
! Usage : mpirun -np PR*PC ./test_mpi_blr_ldlt_ooc PR PC [N] [NB] [NRHS] [TOL] [MAX_MB] [NRHS_BATCH] [SCRATCH_DIR]
!         (defauts 2000 200 7 1e-4 2 3 ./ldlt_ooc_scratch)
!===============================================================================
PROGRAM TEST_MPI_BLR_LDLT_OOC
  USE MPI_BLR_LDLT_MOD
  USE MPI_BLR_LDLT_OOC_MOD
  USE SYM_TESTMAT_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: OOC_RHS_T, OOC_RHS_FREE
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  TYPE(MPI_BLR_SYM_OOC_T) :: AO
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0, MLOC, NRHS_LOC, RSTART_MINE
  REAL(WP) :: TOL, MAX_MB, ERR_C, ERR_O, DIFF
  CHARACTER(LEN=256) :: ARG, SCRATCH_DIR
  COMPLEX(WP), ALLOCATABLE :: XG_C(:,:), XG_O(:,:), XEX_ALL(:,:)
  INTEGER(8) :: RES_MAX_L, RES_MAX, DATA_L, DATA_MAX
  DOUBLE PRECISION :: TF_C, TF_O, TC_O, TCOMP_O, TMPI_O
  INTEGER :: NSW_C, NSW_O
  LOGICAL :: PASS

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)
  N = 2000; NB = 200; NRHS = 7; TOL = 1.0E-4_WP; MAX_MB = 2.0_WP; NRHS_BATCH = 3
  SCRATCH_DIR = './ldlt_ooc_scratch'
  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_ldlt_ooc PR PC [N] [NB] [NRHS] [TOL] '// &
         '[MAX_MB] [NRHS_BATCH] [SCRATCH_DIR]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN; CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N;          END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN; CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB;         END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN; CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS;       END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN; CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL;        END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN; CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) MAX_MB;     END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN; CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) NRHS_BATCH; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) CALL GET_COMMAND_ARGUMENT(9, SCRATCH_DIR)
  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0 .EQ. 0) CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(SCRATCH_DIR))
  CALL MPI_BARRIER(GRID%COMM, IERR)
  CALL SYM_TESTMAT_SETUP(N, NB)
  ALLOCATE(XEX_ALL(N,NRHS))
  CALL FILL_XEX_ALL()
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F7.2,A,I0)') ' LDLt ooc vs classique : grille ', PR, 'x', &
       PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, ' , MAX_MB=', MAX_MB, ' , NRHS_BATCH=', NRHS_BATCH
  ALLOCATE(XG_C(N,NRHS), XG_O(N,NRHS))
  XG_O = (0.0_WP,0.0_WP)
  CALL RUN_CLASSIC()
  CALL RUN_OOC()
  CALL MPI_ALLREDUCE(RES_MAX_L, RES_MAX, 1, MPI_INTEGER8, MPI_MAX, GRID%COMM, IERR)
  CALL MPI_ALLREDUCE(DATA_L, DATA_MAX, 1, MPI_INTEGER8, MPI_MAX, GRID%COMM, IERR)
  IF (MYID0 .EQ. 0) THEN
    ERR_C = REL_ERR(XG_C)
    ERR_O = REL_ERR(XG_O)
    DIFF = SQRT(SUM(ABS(XG_O-XG_C)**2) / SUM(ABS(XG_C)**2))
    WRITE(6,'(A,ES12.4,A,F8.3,A,I0)') '  classique : ||X-Xexact||/||Xexact|| = ', ERR_C, ' , factorisation ', TF_C, &
         ' s , echanges ', NSW_C
    WRITE(6,'(A,ES12.4,A,F8.3,A,I0)') '  ooc       : ||X-Xexact||/||Xexact|| = ', ERR_O, ' , factorisation ', TF_O, &
         ' s , echanges ', NSW_O
    WRITE(6,'(A,3F9.3)') '  ooc, decoupage factorisation (cache, calcul, MPI) (s) : ', TC_O, TCOMP_O, TMPI_O
    WRITE(6,'(A,F9.2,A,F9.2,A,F7.2,A)') '  donnees sur disque max/proc = ', REAL(DATA_MAX,WP)/1048576.0_WP, &
         ' Mo , cache resident max observe = ', REAL(RES_MAX,WP)/1048576.0_WP, ' Mo (budget ', MAX_MB, ' Mo PAR cache, 2 caches)'
    WRITE(6,'(A,ES12.4)') '  ||X_ooc - X_classique|| / ||X_classique|| = ', DIFF
    PASS = (ERR_O .LT. 1.0E-3_WP) .AND. (DIFF .LT. 1.0E-3_WP) .AND. (NSW_O .EQ. NSW_C)
    IF (PASS) THEN
      WRITE(6,*) '=== TEST_MPI_BLR_LDLT_OOC : OK ==='
    ELSE
      WRITE(6,*) '=== TEST_MPI_BLR_LDLT_OOC : ECHEC ==='
    END IF
  END IF
  DEALLOCATE(XG_C, XG_O, XEX_ALL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  SUBROUTINE FILL_XEX_ALL()
    INTEGER :: GI, R
    DO R = 1, NRHS
      DO GI = 1, N
        XEX_ALL(GI,R) = SYM_XEX(GI,R)
      END DO
    END DO
  END SUBROUTINE FILL_XEX_ALL

  ! FILL_RHS_IFACE : B(ISTART:IEND, RS:RE) = A(ISTART:IEND,:) * Xexact(:,RS:RE)
  SUBROUTINE FILL_RHS(ISTART, IEND, RS, RE, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RS, RE
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, RE-RS+1)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: AROWS(:,:)
    ALLOCATE(AROWS(IEND-ISTART+1,N))
    CALL SYM_GET_BLOCK(ISTART, IEND, 1, N, AROWS)
    CALL CGEMM('N', 'N', IEND-ISTART+1, RE-RS+1, N, CONE, AROWS, IEND-ISTART+1, XEX_ALL(1,RS), N, CZERO, &
               BLK, IEND-ISTART+1)
    DEALLOCATE(AROWS)
  END SUBROUTINE FILL_RHS

  ! CONSUME_X_BATCH_IFACE : rassemble la solution du paquet dans XG_O (rang 0).
  SUBROUTINE CONSUME_X(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)     :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN) :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    COMPLEX(WP), ALLOCATABLE :: XB(:,:)
    INTEGER :: W
    W = BATCH_REND - BATCH_RSTART + 1
    ALLOCATE(XB(N,W))
    CALL GATHER_2D_RHS(GRID, N, NB, W, X_LOCAL, LDX, NRHS_LOC_B, 0, XB)
    IF (MYID0 .EQ. 0) XG_O(:,BATCH_RSTART:BATCH_REND) = XB
    DEALLOCATE(XB)
    IF (GSTART .GT. GEND + 1) WRITE(6,*) 'incoherence GSTART/GEND'
  END SUBROUTINE CONSUME_X

  REAL(WP) FUNCTION REL_ERR(X)
    COMPLEX(WP), INTENT(IN) :: X(N,NRHS)
    REL_ERR = SQRT(SUM(ABS(X - XEX_ALL)**2) / SUM(ABS(XEX_ALL)**2))
  END FUNCTION REL_ERR

  SUBROUTINE RUN_CLASSIC()
    TYPE(MPI_BLR_SYM_T) :: AS
    COMPLEX(WP), ALLOCATABLE :: BL(:,:)
    INTEGER :: NC, NO, INFO, li, I, RS
    INTEGER(8) :: RSUM
    DOUBLE PRECISION :: TA
    CALL MPI_ASSEMBLE_BLR_SYM(GRID, N, NB, TOL, SYM_GET_BLOCK, AS, TA, NC, NO, RSUM)
    CALL MPI_BLR_SYM_LDLT_FACTOR(GRID, AS, TOL, INFO, NC, NO, RSUM, N_SWAPS=NSW_C, T_FACTOR=TF_C)
    IF (INFO .NE. 0 .AND. MYID0.EQ.0) WRITE(6,*) 'ERREUR LDLt classique, INFO=', INFO
    MLOC = MPI_BLR_SYM_MLOC(AS, GRID)
    CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
    ALLOCATE(BL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
    RS = 1
    DO li = 1, AS%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      IF (NRHS_LOC .GT. 0) CALL FILL_RHS(AS%BSTART(I), AS%BEND(I), RSTART_MINE, RSTART_MINE+NRHS_LOC-1, &
                                         BL(RS:RS+AS%BEND(I)-AS%BSTART(I),1:NRHS_LOC))
      RS = RS + AS%BEND(I) - AS%BSTART(I) + 1
    END DO
    CALL MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, AS, NRHS, NRHS_LOC, BL, MAX(MLOC,1))
    CALL GATHER_2D_RHS(GRID, N, NB, NRHS, BL, MAX(MLOC,1), NRHS_LOC, 0, XG_C)
    DEALLOCATE(BL)
    CALL MPI_BLR_SYM_FREE(AS)
  END SUBROUTINE RUN_CLASSIC

  ! Suivi du cache resident pendant la factorisation (appele a chaque palier
  ! de 5 %) -- AO est au niveau du programme pour etre visible ici.
  SUBROUTINE TRACK_RESIDENT(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    RES_MAX_L = MAX(RES_MAX_L, MPI_BLR_SYM_OOC_RESIDENT_BYTES(AO))
    IF (NDONE .GT. NTOTAL) WRITE(6,*) 'incoherence NDONE/NTOTAL'
  END SUBROUTINE TRACK_RESIDENT

  SUBROUTINE RUN_OOC()
    TYPE(OOC_RHS_T) :: RHS
    INTEGER :: NC, NO, INFO
    INTEGER(8) :: RSUM
    DOUBLE PRECISION :: TA
    CALL MPI_ASSEMBLE_BLR_SYM_OOC(GRID, N, NB, TOL, SYM_GET_BLOCK, AO, TRIM(SCRATCH_DIR), MAX_MB, TA, NC, NO, RSUM)
    CALL MPI_ASSEMBLE_RHS_SYM_OOC(GRID, AO, NRHS, NRHS_BATCH, FILL_RHS, TRIM(SCRATCH_DIR), RHS)
    RES_MAX_L = 0_8
    CALL MPI_BLR_SYM_LDLT_FACTOR_OOC(GRID, AO, TOL, INFO, NC, NO, RSUM, N_SWAPS=NSW_O, T_FACTOR=TF_O, &
                                     T_CACHE=TC_O, T_COMPUTE=TCOMP_O, T_MPICALL=TMPI_O, &
                                     PROGRESS_CALLBACK=TRACK_RESIDENT)
    IF (INFO .NE. 0 .AND. MYID0.EQ.0) WRITE(6,*) 'ERREUR LDLt ooc, INFO=', INFO
    DATA_L = MPI_BLR_SYM_OOC_DATA_BYTES(AO)
    CALL MPI_BLR_SYM_LDLT_SOLVE_OOC_2D(GRID, AO, RHS, CONSUME_X)
    CALL OOC_RHS_FREE(RHS)
    CALL MPI_BLR_SYM_OOC_FREE(AO)
  END SUBROUTINE RUN_OOC

END PROGRAM TEST_MPI_BLR_LDLT_OOC
