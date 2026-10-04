!===============================================================================
! Benchmark LDL^T BLR+MPI vs LU BLR+MPI (reference, pivot confine) sur la
! MEME matrice complexe symetrique (sym_testmat_mod.f90), meme grille,
! meme NB/TOL, meme second membre B = A*Xexact (produit brut sur la formule).
!
! Rapporte, pour chaque solveur : temps d'assemblage / factorisation /
! resolution (MAX sur les processus), memoire theorique de la matrice
! assemblee et des facteurs (MAX par processus ET somme sur la grille),
! erreur ||X-Xexact||/||Xexact||. Les deux solveurs tournent l'un apres
! l'autre dans le meme processus (le premier est libere avant le second).
!
! Usage : mpirun -np PR*PC ./bench_ldlt_vs_lu PR PC [N] [NB] [NRHS] [TOL] [USE_RSVD] [STRONG] [DIAG] [LOOKAHEAD] [DIAG_SPREAD]
!         (defauts 8000 400 8 1e-4 1 50 0.01 1 1)
!===============================================================================
PROGRAM BENCH_LDLT_VS_LU
  USE MPI_BLR_LDLT_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: MPI_BLR_MATRIX_T, MPI_ASSEMBLE_BLR_MATRIX, MPI_BLR_MATRIX_LU_FACTOR_PIVOT, &
                                MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D, MPI_BLR_MATRIX_MEMORY_BYTES
  USE SYM_TESTMAT_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0, USE_RSVD_INT, LA_INT, SPREAD_INT, MLOC, NRHS_LOC, RSTART_MINE
  REAL(WP) :: TOL, STRONG, DIAG
  CHARACTER(LEN=32) :: ARG
  ! resultats : 1 = LU, 2 = LDLt
  DOUBLE PRECISION :: T_ASM(2), T_FAC(2), T_SOL(2)
  INTEGER(8) :: MEM_A_MAX(2), MEM_A_SUM(2), MEM_F_MAX(2), MEM_F_SUM(2)
  REAL(WP) :: ERR(2)
  INTEGER :: NSW(2), N2X2
  DOUBLE PRECISION :: TBK(7,2)
  INTEGER :: JB
  CHARACTER(LEN=26), PARAMETER :: TBK_NAMES(7) = (/ 'facto. blocs diagonaux    ', 'echanges (pivot)          ', &
       'panneau (TRSM, D^-1, ...) ', 'MPI (attentes comprises)  ', 'Schur blocs diagonaux     ', &
       'Schur hors-diagonale      ', 'boucle complete           ' /)

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)
  N = 8000; NB = 400; NRHS = 8; TOL = 1.0E-4_WP; USE_RSVD_INT = 1; STRONG = 50.0_WP; DIAG = 0.01_WP; LA_INT = 1; SPREAD_INT = 1
  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./bench_ldlt_vs_lu PR PC [N] [NB] [NRHS] [TOL] [USE_RSVD] [STRONG] [DIAG]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN; CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN; CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN; CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN; CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN; CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) USE_RSVD_INT; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN; CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) STRONG; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) THEN; CALL GET_COMMAND_ARGUMENT(9, ARG); READ(ARG,*) DIAG; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN; CALL GET_COMMAND_ARGUMENT(10, ARG); READ(ARG,*) LA_INT; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN; CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) SPREAD_INT; END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL SYM_TESTMAT_SETUP(N, NB, STRONG=STRONG, DIAG=DIAG)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,I0,A,F6.1,A,F7.3)') ' Grille ', PR, 'x', PC, &
       ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, ' , USE_RSVD=', USE_RSVD_INT, &
       ' , STRONG=', STRONG, ' , DIAG=', DIAG
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0)') ' LOOKAHEAD=', LA_INT, ' , DIAG_SPREAD=', SPREAD_INT

  CALL RUN_LU()
  CALL RUN_LDLT()

  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') ' ======================================================================'
    WRITE(6,'(A)') '                         LU (pivot confine)     LDLt (Bunch-Kaufman)   LDLt/LU'
    WRITE(6,'(A)') ' ======================================================================'
    WRITE(6,'(A,F14.3,A,F14.3,A,F10.3)') '  Assemblage (s)        ', T_ASM(1), '      ', T_ASM(2), '     ', T_ASM(2)/T_ASM(1)
    WRITE(6,'(A,F14.3,A,F14.3,A,F10.3)') '  Factorisation (s)     ', T_FAC(1), '      ', T_FAC(2), '     ', T_FAC(2)/T_FAC(1)
    WRITE(6,'(A,F14.3,A,F14.3,A,F10.3)') '  Resolution (s)        ', T_SOL(1), '      ', T_SOL(2), '     ', T_SOL(2)/T_SOL(1)
    WRITE(6,'(A,F14.1,A,F14.1,A,F10.3)') '  Matrice A, max/proc (Mo)', MB(MEM_A_MAX(1)), '    ', MB(MEM_A_MAX(2)), '     ', &
         REAL(MEM_A_MAX(2),WP)/REAL(MEM_A_MAX(1),WP)
    WRITE(6,'(A,F14.1,A,F14.1,A,F10.3)') '  Matrice A, total (Mo)   ', MB(MEM_A_SUM(1)), '    ', MB(MEM_A_SUM(2)), '     ', &
         REAL(MEM_A_SUM(2),WP)/REAL(MEM_A_SUM(1),WP)
    WRITE(6,'(A,F14.1,A,F14.1,A,F10.3)') '  Facteurs, max/proc (Mo) ', MB(MEM_F_MAX(1)), '    ', MB(MEM_F_MAX(2)), '     ', &
         REAL(MEM_F_MAX(2),WP)/REAL(MEM_F_MAX(1),WP)
    WRITE(6,'(A,F14.1,A,F14.1,A,F10.3)') '  Facteurs, total (Mo)    ', MB(MEM_F_SUM(1)), '    ', MB(MEM_F_SUM(2)), '     ', &
         REAL(MEM_F_SUM(2),WP)/REAL(MEM_F_SUM(1),WP)
    WRITE(6,'(A,ES14.4,A,ES14.4)') '  ||X-Xexact||/||Xexact|| ', ERR(1), '      ', ERR(2)
    WRITE(6,'(A,I14,A,I14)')       '  Echanges de pivot       ', NSW(1), '      ', NSW(2)
    WRITE(6,'(A,I0)')              '  Pivots 2x2 (LDLt)       ', N2X2
    WRITE(6,'(A)') '  Decoupage de la factorisation (s, MAX sur les processus) :'
    DO JB = 1, 7
      WRITE(6,'(A,A,F12.4,A,F12.4)') '    ', TBK_NAMES(JB), TBK(JB,1), '      ', TBK(JB,2)
    END DO
    IF (ERR(1).LT.1.0E-3_WP .AND. ERR(2).LT.1.0E-3_WP) THEN
      WRITE(6,*) '=== BENCH_LDLT_VS_LU : OK ==='
    ELSE
      WRITE(6,*) '=== BENCH_LDLT_VS_LU : ECHEC (erreur > 1e-3) ==='
    END IF
  END IF
  CALL MPI_FINALIZE(IERR)

CONTAINS

  REAL(WP) FUNCTION MB(B)
    INTEGER(8), INTENT(IN) :: B
    MB = REAL(B,WP)/1024.0_WP/1024.0_WP
  END FUNCTION MB

  SUBROUTINE REDUCE_MEM(BYTES, BMAX, BSUM)
    INTEGER(8), INTENT(IN)  :: BYTES
    INTEGER(8), INTENT(OUT) :: BMAX, BSUM
    INTEGER :: IERR2
    CALL MPI_ALLREDUCE(BYTES, BMAX, 1, MPI_INTEGER8, MPI_MAX, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(BYTES, BSUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR2)
  END SUBROUTINE REDUCE_MEM

  ! B_LOCAL = A*Xexact sur les lignes locales (meme distribution pour les
  ! deux solveurs : blocs-lignes cycliques de NB lignes).
  SUBROUTINE FILL_B_EXACT(NBLKR_LOC, BSTART, BEND, BL)
    INTEGER, INTENT(IN)      :: NBLKR_LOC, BSTART(:), BEND(:)
    COMPLEX(WP), INTENT(OUT) :: BL(MAX(MLOC,1),MAX(NRHS_LOC,1))
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
    DO li = 1, NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      BSZ = BEND(I) - BSTART(I) + 1
      ALLOCATE(AROWS(BSZ,N))
      CALL SYM_GET_BLOCK(BSTART(I), BEND(I), 1, N, AROWS)
      CALL CGEMM('N', 'N', BSZ, NRHS_LOC, N, CONE, AROWS, BSZ, XEX, N, CZERO, BL(RS,1), MAX(MLOC,1))
      DEALLOCATE(AROWS)
      RS = RS + BSZ
    END DO
    DEALLOCATE(XEX)
  END SUBROUTINE FILL_B_EXACT

  REAL(WP) FUNCTION ERR_OF(BL)
    COMPLEX(WP), INTENT(IN) :: BL(MAX(MLOC,1),MAX(NRHS_LOC,1))
    COMPLEX(WP), ALLOCATABLE :: XG(:,:)
    REAL(WP) :: NUM, DEN
    INTEGER :: GI, R
    ALLOCATE(XG(N,NRHS))
    CALL GATHER_2D_RHS(GRID, N, NB, NRHS, BL, MAX(MLOC,1), NRHS_LOC, 0, XG)
    ERR_OF = 0.0_WP
    IF (MYID0 .EQ. 0) THEN
      NUM = 0.0_WP; DEN = 0.0_WP
      DO R = 1, NRHS
        DO GI = 1, N
          NUM = NUM + ABS(XG(GI,R) - SYM_XEX(GI,R))**2
          DEN = DEN + ABS(SYM_XEX(GI,R))**2
        END DO
      END DO
      ERR_OF = SQRT(NUM/DEN)
    END IF
    DEALLOCATE(XG)
  END FUNCTION ERR_OF

  SUBROUTINE RUN_LU()
    TYPE(MPI_BLR_MATRIX_T) :: AL
    COMPLEX(WP), ALLOCATABLE :: BL(:,:)
    INTEGER :: NC, NO, INFO, li, I
    INTEGER(8) :: RS
    DOUBLE PRECISION :: T0
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' --- LU BLR+MPI (reference, matrice entiere) ---'
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, SYM_GET_BLOCK, AL, T_ASM(1), NC, NO, RS, &
                                 SHOW_PROGRESS=.TRUE., USE_RSVD=(USE_RSVD_INT.NE.0))
    CALL REDUCE_MEM(MPI_BLR_MATRIX_MEMORY_BYTES(AL), MEM_A_MAX(1), MEM_A_SUM(1))
    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, AL, TOL, INFO, NC, NO, RS, NSW(1), T_FAC(1), SHOW_PROGRESS=.TRUE., &
                                        T_BREAKDOWN=TBK(:,1))
    IF (INFO .NE. 0 .AND. MYID0.EQ.0) WRITE(6,*) 'ERREUR LU, INFO=', INFO
    CALL REDUCE_MEM(MPI_BLR_MATRIX_MEMORY_BYTES(AL), MEM_F_MAX(1), MEM_F_SUM(1))
    MLOC = 0
    DO li = 1, AL%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + AL%BEND(I) - AL%BSTART(I) + 1
    END DO
    ALLOCATE(BL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
    CALL FILL_B_EXACT(AL%NBLKR_LOC, AL%BSTART, AL%BEND, BL)
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, AL, NRHS, NRHS_LOC, BL, MAX(MLOC,1))
    T_SOL(1) = MPI_WTIME() - T0
    ERR(1) = ERR_OF(BL)
    DEALLOCATE(BL)
  END SUBROUTINE RUN_LU

  SUBROUTINE RUN_LDLT()
    TYPE(MPI_BLR_SYM_T) :: AS
    DOUBLE PRECISION :: TPH_LDLT(3)
    COMPLEX(WP), ALLOCATABLE :: BL(:,:)
    INTEGER :: NC, NO, INFO
    INTEGER(8) :: RS
    DOUBLE PRECISION :: T0
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' --- LDLt BLR+MPI (triangle inferieur) ---'
    CALL MPI_ASSEMBLE_BLR_SYM(GRID, N, NB, TOL, SYM_GET_BLOCK, AS, T_ASM(2), NC, NO, RS, &
                              SHOW_PROGRESS=.TRUE., USE_RSVD=(USE_RSVD_INT.NE.0), DIAG_SPREAD=(SPREAD_INT.NE.0))
    CALL REDUCE_MEM(MPI_BLR_SYM_MEMORY_BYTES(AS), MEM_A_MAX(2), MEM_A_SUM(2))
    CALL MPI_BLR_SYM_LDLT_FACTOR(GRID, AS, TOL, INFO, NC, NO, RS, NSW(2), N2X2, T_FAC(2), SHOW_PROGRESS=.TRUE., &
                                 T_BREAKDOWN=TBK(:,2), LOOKAHEAD=(LA_INT.NE.0))
    IF (INFO .NE. 0 .AND. MYID0.EQ.0) WRITE(6,*) 'ERREUR LDLt, INFO=', INFO
    CALL REDUCE_MEM(MPI_BLR_SYM_MEMORY_BYTES(AS), MEM_F_MAX(2), MEM_F_SUM(2))
    MLOC = MPI_BLR_SYM_MLOC(AS, GRID)
    ALLOCATE(BL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
    CALL FILL_B_EXACT(AS%NBLKR_LOC, AS%BSTART, AS%BEND, BL)
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = MPI_WTIME()
    CALL MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, AS, NRHS, NRHS_LOC, BL, MAX(MLOC,1), T_PHASES=TPH_LDLT)
    T_SOL(2) = MPI_WTIME() - T0
    IF (MYID0.EQ.0) WRITE(6,'(A,3F10.4)') '  Resolution LDLt, phases (descente, D, remontee) (s) :', TPH_LDLT
    ERR(2) = ERR_OF(BL)
    DEALLOCATE(BL)
    CALL MPI_BLR_SYM_FREE(AS)
  END SUBROUTINE RUN_LDLT

END PROGRAM BENCH_LDLT_VS_LU
