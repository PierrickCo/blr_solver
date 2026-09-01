!===============================================================================
! Test du solveur LU MPI, BLR, AVEC PIVOT, VARIANTE OUT-OF-CORE (voir
! mpi_blr_matrix_mod.f90, type MPI_BLR_MATRIX_OOC_T, et blr_ooc_mod.f90) --
! MEME methodologie de validation que test_mpi_blr_matrix_pivot.f90/
! test_mpi_blr_matrix_pivot_arena.f90 (X_exact connu par formule
! deterministe, comparaison directe) -- MAIS le second membre B est calcule
! ICI par PRODUIT BRUT (boucle sur les N colonnes de A, PAS via le A_BLR
! compresse) car l'assemblage OOC de B (MPI_ASSEMBLE_RHS_OOC) decoupe les
! seconds membres par PAQUETS dont la distribution colonnes<->processus
! (RHS_LOCAL_RANGE_MBM sur la LARGEUR du paquet) peut DIFFERER de celle
! d'un appel a NRHS complet (RHS_LOCAL_RANGE_MBM sur NRHS total) -- MESURE
! et documentee dans la memoire du projet, PAS un bug (la solution restant
! mathematiquement correcte quelle que soit cette distribution). Un calcul
! brut, base UNIQUEMENT sur les indices GLOBAUX demandes par FILL_RHS, est
! donc la seule facon simple d'obtenir un B independant de CETTE
! distribution.
!
! Teste explicitement le PARAMETRE NRHS_BATCH (7e argument) : la meme
! matrice/Xexact doit rester correctement resolue quel que soit le
! decoupage en paquets, y compris NRHS_BATCH=1 (chaque second membre
! resolu separement).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_blr_matrix_pivot_ooc PR PC N NB NRHS TOL MAX_MB_PER_PROCESS NRHS_BATCH [SCRATCH_DIR]
!===============================================================================
PROGRAM TEST_MPI_BLR_MATRIX_PIVOT_OOC
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_BLR
  TYPE(OOC_RHS_T) :: RHS
  COMPLEX(WP), ALLOCATABLE :: X_GLOBAL(:,:), B_REF(:,:)
  INTEGER :: INFO, N_SWAPS, N_COMPRESSED, N_OFFDIAG, GI, R
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR
  REAL(WP) :: ERR

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 8) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./test_mpi_blr_matrix_pivot_ooc PR PC N NB NRHS TOL '// &
         'MAX_MB_PER_PROCESS NRHS_BATCH [SCRATCH_DIR]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL
  CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) MAX_MB_PER_PROCESS
  CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) NRHS_BATCH
  SCRATCH_DIR = './blr_ooc_scratch_test'
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) CALL GET_COMMAND_ARGUMENT(9, SCRATCH_DIR)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ! === 1) assemblage BLR distribue, out-of-core ===
  CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                    T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM)

  ! === 2) assemblage du second membre, out-of-core, par paquets de NRHS_BATCH ===
  CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)

  ! === 3) factorisation LU, AVEC pivot, out-of-core ===
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                           N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ! === 4) resolution, out-of-core, par paquets ===
  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH)

  ERR = -1.0_WP
  IF (MYID0 .EQ. 0) THEN
    ALLOCATE(B_REF(N,NRHS))
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = XEX_FORMULA(GI,R)
      END DO
    END DO
    ERR = SQRT(SUM(ABS(X_GLOBAL-B_REF)**2)) / SQRT(SUM(ABS(B_REF)**2))

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1,A,I0)') ' [OOC] Grille ', PR, 'x', PC, ' , N=', N, &
         ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, ' , MAX_MB_PER_PROCESS=', MAX_MB_PER_PROCESS, &
         ' , NRHS_BATCH=', NRHS_BATCH
    WRITE(6,'(A,F8.3,A,F8.3,A)') '   Temps assemblage = ', T_ASM, ' s , factorisation = ', T_FACTOR, ' s'
    WRITE(6,'(A,I0)') '   Nombre de swaps de pivot effectues (total) = ', N_SWAPS
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Compression : ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', &
         100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_calcule - X_exact|| = ', ERR
    IF (ERR.LT.1.0E-3_WP .AND. INFO.EQ.0 .AND. N_SWAPS.GT.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT_OOC : OK ==='
    ELSE IF (N_SWAPS.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT_OOC : SUSPECT (aucun pivot exerce, matrice pas assez agressive) ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT_OOC : ECHEC ==='
    END IF
    DEALLOCATE(B_REF)
  END IF

  DEALLOCATE(X_GLOBAL)
  CALL OOC_RHS_FREE(RHS)
  CALL MPI_BLR_MATRIX_OOC_FREE(A_BLR)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! A(GI,GJ) = 50 sur une anti-diagonale CONFINEE AU BLOC DIAGONAL de GI --
  ! IDENTIQUE a test_mpi_blr_matrix_pivot.f90/_arena.f90 (force un vrai
  ! pivotage QUE LE PIVOT CONFINE PEUT VOIR).
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
      DF = 1.0_WP / (1.0_WP + ABS(REAL(GJ_STRONG-GJ,WP)))
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

  ! B = A*Xexact, PAR PRODUIT BRUT (boucle sur les N colonnes de A) --
  ! INDEPENDANT de toute distribution colonnes<->processus (voir
  ! commentaire de tete) : correct pour N'IMPORTE QUELLE fenetre
  ! (ISTART,IEND,RSTART,REND) demandee par MPI_ASSEMBLE_RHS_OOC, quel que
  ! soit NRHS_BATCH.
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

  ! Rassemble CE paquet resolu dans X_GLOBAL (verification uniquement).
  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    COMPLEX(WP), ALLOCATABLE :: X_BATCH_GLOBAL(:,:)
    INTEGER :: WIDTH
    WIDTH = BATCH_REND - BATCH_RSTART + 1
    ALLOCATE(X_BATCH_GLOBAL(N,WIDTH))
    CALL GATHER_2D_RHS(GRID, N, NB, WIDTH, X_LOCAL, LDX, NRHS_LOC_B, 0, X_BATCH_GLOBAL)
    IF (MYID0 .EQ. 0) X_GLOBAL(:,BATCH_RSTART:BATCH_REND) = X_BATCH_GLOBAL
    DEALLOCATE(X_BATCH_GLOBAL)
  END SUBROUTINE CONSUME_X_BATCH

END PROGRAM TEST_MPI_BLR_MATRIX_PIVOT_OOC
