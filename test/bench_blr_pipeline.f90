!===============================================================================
! PIPELINE "MATRICE COMPRESSEE DE BOUT EN BOUT" (blr_matrix_mod.f90),
! organise EXACTEMENT selon les 5 etapes demandees par l'utilisateur :
!
!   1) Assemblage de A bloc par bloc, CHAQUE bloc compresse en low-rank
!      IMMEDIATEMENT avant de passer au bloc suivant -- mesure du temps
!      combine assemblage+compression, et du taux de compression de A.
!   2) Assemblage du second membre B -- mesure separee.
!   3) Factorisation LU de la matrice DEJA COMPRESSEE -- mesure du temps,
!      et du taux de compression des facteurs L/U finaux.
!   4) Resolution du systeme.
!   5) Verification de la precision EN APPLIQUANT LA MATRICE COMPRESSEE
!      (pas une reference dense recalculee a part) : Y = A_BLR_original*X,
!      compare a B.
!
! Matrice/second membre : meme recette deterministe "diagonale dominante +
! correction basse-rang" (rang global 6) que bench_seq.f90/bench_mpi.f90,
! pour rester comparable aux benchmarks precedents. SANS PIVOT (voir la
! note de portee de blr_matrix_mod.f90).
!
! Usage : ./bench_blr_pipeline [N] [NRHS] [NB] [TOL]  (defauts : 10000 100 100 1e-4)
!===============================================================================
PROGRAM BENCH_BLR_PIPELINE
  USE BLR_MATRIX_MOD
  USE BENCH_UTILS_MOD
  IMPLICIT NONE

  INTEGER :: N, NRHS, NB
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  INTEGER, PARAMETER :: RGLOBAL = 6

  N = 10000; NRHS = 100; NB = 100; TOL = 1.0E-4_WP
  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) TOL
  END IF

  WRITE(6,'(A,I0,A,I0,A,I0,A,ES10.3)') ' === PIPELINE MATRICE COMPRESSEE : N=', N, &
       ' , NRHS=', NRHS, ' , NB=', NB, ' , TOL=', TOL
  WRITE(6,*)

  CALL RUN_PIPELINE(N, NRHS, NB, TOL)

CONTAINS

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

  DOUBLE PRECISION FUNCTION WTIME()
    USE OMP_LIB, ONLY: OMP_GET_WTIME
    WTIME = OMP_GET_WTIME()
  END FUNCTION WTIME

  FUNCTION KB_TO_MB(KB) RESULT(MB)
    INTEGER(8), INTENT(IN) :: KB
    DOUBLE PRECISION :: MB
    MB = REAL(KB,8) / 1024.0D0
  END FUNCTION KB_TO_MB

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_PIPELINE(N, NRHS, NB, TOL)
    INTEGER, INTENT(IN)  :: N, NRHS, NB
    REAL(WP), INTENT(IN) :: TOL

    TYPE(BLR_MATRIX_T) :: A_BLR, A_BLR_ORIG
    COMPLEX(WP), ALLOCATABLE :: B(:,:), XEX(:,:), Y(:,:)
    INTEGER :: INFO, N_COMP_A, N_OFF_A, N_COMP_LU, N_OFF_LU
    INTEGER(8) :: RANK_SUM_A, RANK_SUM_LU
    DOUBLE PRECISION :: T0, T1
    DOUBLE PRECISION :: T_ASSEMBLE_COMPRESS_A, T_ASSEMBLE_B, T_FACTOR, T_SOLVE
    DOUBLE PRECISION :: MEM_A_MB, MEM_DENSE_MB, RSS_MB
    REAL(WP) :: ERR
    INTEGER :: GI, R

    ! === 1) Assemblage de A, bloc par bloc, compression immediate ===
    CALL ASSEMBLE_BLR_MATRIX(N, NB, TOL, GET_BLOCK, A_BLR, T_ASSEMBLE_COMPRESS_A, &
                              N_COMP_A, N_OFF_A, RANK_SUM_A)
    MEM_A_MB = KB_TO_MB(BLR_MATRIX_MEMORY_BYTES(A_BLR)/1024_8)
    MEM_DENSE_MB = 8.0D0*REAL(N,8)*REAL(N,8)/1024.0D0/1024.0D0
    RSS_MB = KB_TO_MB(GET_RSS_KB())

    WRITE(6,'(A)') ' --- 1) Assemblage de A (bloc par bloc, compression immediate) ---'
    WRITE(6,'(A,F8.3,A)')  '   Temps assemblage+compression = ', T_ASSEMBLE_COMPRESS_A, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Blocs compresses = ', N_COMP_A, ' / ', N_OFF_A, ' (', &
         100.0_WP*REAL(N_COMP_A,WP)/REAL(MAX(N_OFF_A,1),WP), ' %)'
    WRITE(6,'(A,F8.2)') '   Rang moyen = ', REAL(RANK_SUM_A,WP)/REAL(MAX(N_OFF_A,1),WP)
    WRITE(6,'(A,F10.1,A,F10.1,A,F6.1,A)') '   Memoire de A : ', MEM_A_MB, ' Mo (compressee) vs ', &
         MEM_DENSE_MB, ' Mo (dense) -> ', 100.0_WP*MEM_A_MB/MEM_DENSE_MB, ' % du dense'
    WRITE(6,'(A,F10.1,A)') '   RSS apres assemblage = ', RSS_MB, ' Mo'
    WRITE(6,*)

    ! garder une copie NON FACTORISEE pour la verification (etape 5) --
    ! peu couteux puisque A_BLR est deja compressee
    CALL COPY_BLR_MATRIX(A_BLR, A_BLR_ORIG)

    ! === 2) Assemblage du second membre B ===
    T0 = WTIME()
    ALLOCATE(B(N,NRHS), XEX(N,NRHS))
    CALL SET_FIXED_SEED()
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    T1 = WTIME()
    T_ASSEMBLE_B = T1 - T0
    ! B = A_original * Xexact (utilise la matrice COMPRESSEE, coherent
    ! avec l'etape 5 qui reutilise la meme routine d'application)
    ALLOCATE(Y(N,NRHS))
    CALL APPLY_BLR_MATRIX(A_BLR_ORIG, NRHS, XEX, N, Y, N)
    B = Y
    T1 = WTIME()

    WRITE(6,'(A)') ' --- 2) Assemblage du second membre B ---'
    WRITE(6,'(A,F8.3,A)') '   Temps assemblage (Xexact + B=A_BLR*Xexact) = ', T1-T0, ' s'
    WRITE(6,*)

    ! === 3) Factorisation LU de la matrice DEJA COMPRESSEE ===
    CALL BLR_MATRIX_LU_FACTOR(A_BLR, TOL, INFO, N_COMP_LU, N_OFF_LU, RANK_SUM_LU, T_FACTOR)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLR_MATRIX_LU_FACTOR, INFO=', INFO
      RETURN
    END IF
    RSS_MB = KB_TO_MB(GET_RSS_KB())

    WRITE(6,'(A)') ' --- 3) Factorisation LU (sur la matrice compressee) ---'
    WRITE(6,'(A,F8.3,A)') '   Temps factorisation = ', T_FACTOR, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Blocs L/U compresses = ', N_COMP_LU, ' / ', N_OFF_LU, ' (', &
         100.0_WP*REAL(N_COMP_LU,WP)/REAL(MAX(N_OFF_LU,1),WP), ' %)'
    WRITE(6,'(A,F8.2)') '   Rang moyen (L/U) = ', REAL(RANK_SUM_LU,WP)/REAL(MAX(N_OFF_LU,1),WP)
    WRITE(6,'(A,F10.1,A)') '   Memoire des facteurs L/U = ', &
         KB_TO_MB(BLR_MATRIX_MEMORY_BYTES(A_BLR)/1024_8), ' Mo'
    WRITE(6,'(A,F10.1,A)') '   RSS apres factorisation = ', RSS_MB, ' Mo'
    WRITE(6,*)

    ! === 4) Resolution ===
    T0 = WTIME()
    CALL BLR_MATRIX_LU_SOLVE(A_BLR, NRHS, B, N)
    T1 = WTIME()
    T_SOLVE = T1 - T0

    WRITE(6,'(A)') ' --- 4) Resolution ---'
    WRITE(6,'(A,F8.3,A)') '   Temps resolution = ', T_SOLVE, ' s'
    WRITE(6,*)

    ! === 5) Verification EN APPLIQUANT LA MATRICE COMPRESSEE ===
    ! Y = A_BLR_ORIG * Xcalc (B contient Xcalc apres l'etape 4), compare a
    ! B_original -- on a ecrase B, donc on recompare a XEX directement
    ! (equivalent, puisque B_original = A_BLR_ORIG*Xexact par construction).
    ERR = SQRT(SUM(ABS(B-XEX)**2)) / SQRT(SUM(ABS(XEX)**2))

    ! residu "physique" complementaire : Y=A_BLR_ORIG*Xcalc, compare a
    ! B_original (recalcule ici puisque B a ete ecrase par le solve)
    BLOCK
      COMPLEX(WP), ALLOCATABLE :: B_ORIG(:,:)
      REAL(WP) :: RESID
      ALLOCATE(B_ORIG(N,NRHS))
      CALL APPLY_BLR_MATRIX(A_BLR_ORIG, NRHS, XEX, N, B_ORIG, N)
      CALL APPLY_BLR_MATRIX(A_BLR_ORIG, NRHS, B, N, Y, N)
      RESID = SQRT(SUM(ABS(Y-B_ORIG)**2)) / SQRT(SUM(ABS(B_ORIG)**2))
      WRITE(6,'(A)') ' --- 5) Verification (via la matrice COMPRESSEE) ---'
      WRITE(6,'(A,ES12.4)') '   Erreur relative ||Xcalc-Xexact||/||Xexact||        = ', ERR
      WRITE(6,'(A,ES12.4)') '   Residu  ||A_BLR*Xcalc-A_BLR*Xexact||/||A_BLR*Xexact|| = ', RESID
      IF (ERR.LT.MAX(100.0_WP*TOL,1.0E-4_WP)) THEN
        WRITE(6,*) '=== PIPELINE : OK ==='
      ELSE
        WRITE(6,*) '=== PIPELINE : ECHEC ==='
      END IF
      DEALLOCATE(B_ORIG)
    END BLOCK

    DEALLOCATE(B, XEX, Y)
  END SUBROUTINE RUN_PIPELINE

  !-----------------------------------------------------------------------
  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ, R
    COMPLEX(WP) :: AVAL
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(II,R)*LOWRANK_V(R,JJ)
        END DO
        IF (II .EQ. JJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        BLK(II-ISTART+1, JJ-JSTART+1) = AVAL
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  !-----------------------------------------------------------------------
  SUBROUTINE COPY_BLR_MATRIX(SRC, DST)
    TYPE(BLR_MATRIX_T), INTENT(IN)  :: SRC
    TYPE(BLR_MATRIX_T), INTENT(OUT) :: DST
    INTEGER :: I, J, MI

    DST%N = SRC%N; DST%NB = SRC%NB; DST%NBLK = SRC%NBLK
    ALLOCATE(DST%BSTART(SRC%NBLK), DST%BEND(SRC%NBLK))
    DST%BSTART = SRC%BSTART; DST%BEND = SRC%BEND
    ALLOCATE(DST%DIAGBLK(SRC%NBLK))
    ALLOCATE(DST%OFFBLK(SRC%NBLK, SRC%NBLK))
    DO I = 1, SRC%NBLK
      MI = SRC%BEND(I) - SRC%BSTART(I) + 1
      ALLOCATE(DST%DIAGBLK(I)%D(MI,MI))
      DST%DIAGBLK(I)%D = SRC%DIAGBLK(I)%D
      DO J = 1, SRC%NBLK
        IF (I .NE. J) THEN
          DST%OFFBLK(I,J) = SRC%OFFBLK(I,J)   ! copie profonde (composantes allouables)
        END IF
      END DO
    END DO
  END SUBROUTINE COPY_BLR_MATRIX

  SUBROUTINE SET_FIXED_SEED()
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: NSEED
    CALL RANDOM_SEED(SIZE=NSEED)
    ALLOCATE(SEED(NSEED))
    SEED = 42
    CALL RANDOM_SEED(PUT=SEED)
    DEALLOCATE(SEED)
  END SUBROUTINE SET_FIXED_SEED

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

END PROGRAM BENCH_BLR_PIPELINE
