!===============================================================================
! REFERENCE "C" : solveur direct dense, sans decoupage en blocs -- appele
! avec NB=N, BLOCK_LU_FACTOR degenere alors en un seul CGETRF monolithique
! (aucune orchestration par blocs de ce projet, uniquement MKL). Fait partie
! du jeu de reference standard a 5 solveurs de ce projet (A/B/MPI/C/BLR+MPI,
! voir le README) -- mesure temps d'assemblage/factorisation/resolution,
! erreur relative, memoire (analytique + RSS reel via bench_utils_mod.f90).
!
! Matrice : recette deterministe "diagonale dominante + correction basse-
! rang" (rang global 6), la meme utilisee par tous les benchmarks du
! projet, pour rester directement comparable entre solveurs.
!
! (Les variantes BLR a stockage dense/compresse historiquement comparees
! ici -- etapes 2/5 et proposition 2 -- ont ete deplacees vers legacy/, ou
! elles restent utilisables telles quelles : voir legacy/test/bench_seq_full.f90.)
!
! Usage : ./bench_seq [N] [NRHS] [NB] [TOL]  (defauts : 10000 100 100 1e-4 ;
!         NB=N pour la reference C standard)
!===============================================================================
PROGRAM BENCH_SEQ
  USE BLOCK_LU_MOD
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

  WRITE(6,'(A,I0,A,I0,A,I0,A,ES10.3)') ' === BENCHMARK SEQUENTIEL (reference C) : N=', N, &
       ' , NRHS=', NRHS, ' , NB=', NB, ' , TOL=', TOL
  WRITE(6,*)

  CALL RUN_DENSE(N, NRHS, NB)

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE BUILD_MATRIX(N, A, T_ASSEMBLY)
    INTEGER, INTENT(IN)        :: N
    COMPLEX(WP), INTENT(OUT)   :: A(N,N)
    DOUBLE PRECISION, INTENT(OUT) :: T_ASSEMBLY
    INTEGER :: GI, GJ, R
    COMPLEX(WP) :: AVAL
    DOUBLE PRECISION :: T0, T1

    T0 = OMP_GET_WTIME_L()
    DO GJ = 1, N
      DO GI = 1, N
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
        END DO
        IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        A(GI,GJ) = AVAL
      END DO
    END DO
    T1 = OMP_GET_WTIME_L()
    T_ASSEMBLY = T1 - T0
  END SUBROUTINE BUILD_MATRIX

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

  DOUBLE PRECISION FUNCTION OMP_GET_WTIME_L()
    USE BLR_KINDS_MOD, ONLY: BLR_WALLTIME
    OMP_GET_WTIME_L = BLR_WALLTIME()
  END FUNCTION OMP_GET_WTIME_L

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

  SUBROUTINE SET_FIXED_SEED()
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: NSEED
    CALL RANDOM_SEED(SIZE=NSEED)
    ALLOCATE(SEED(NSEED))
    SEED = 42
    CALL RANDOM_SEED(PUT=SEED)
    DEALLOCATE(SEED)
  END SUBROUTINE SET_FIXED_SEED

  SUBROUTINE REPORT(LABEL, T_ASM, T_COMPRESS, T_FACT_PURE, T_SOLVE, ERR, &
                     MEM_ANALYTIC_MB, RSS_PEAK_MB, RSS_FINAL_MB)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    DOUBLE PRECISION, INTENT(IN) :: T_ASM, T_COMPRESS, T_FACT_PURE, T_SOLVE
    REAL(WP), INTENT(IN) :: ERR
    DOUBLE PRECISION, INTENT(IN) :: MEM_ANALYTIC_MB, RSS_PEAK_MB, RSS_FINAL_MB
    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'
    WRITE(6,'(A,F8.3,A)')   '   Temps assemblage              = ', T_ASM, ' s'
    WRITE(6,'(A,F8.3,A)')   '   Temps compression basse-rang  = ', T_COMPRESS, ' s'
    WRITE(6,'(A,F8.3,A)')   '   Temps factorisation (hors compression) = ', T_FACT_PURE, ' s'
    WRITE(6,'(A,F8.3,A)')   '   Temps resolution               = ', T_SOLVE, ' s'
    WRITE(6,'(A,ES12.4)')   '   Erreur relative                = ', ERR
    WRITE(6,'(A,F10.1,A)')  '   Memoire analytique (stockage final) = ', MEM_ANALYTIC_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)')  '   RSS pic (pendant factorisation)     = ', RSS_PEAK_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)')  '   RSS final (apres liberation scratch)= ', RSS_FINAL_MB, ' Mo'
    WRITE(6,*)
  END SUBROUTINE REPORT

  FUNCTION KB_TO_MB(KB) RESULT(MB)
    INTEGER(8), INTENT(IN) :: KB
    DOUBLE PRECISION :: MB
    MB = REAL(KB,8) / 1024.0D0
  END FUNCTION KB_TO_MB

  FUNCTION FROB_ERR(N, NRHS, XCALC, XEX) RESULT(ERR)
    INTEGER, INTENT(IN) :: N, NRHS
    COMPLEX(WP), INTENT(IN) :: XCALC(N,NRHS), XEX(N,NRHS)
    REAL(WP) :: ERR
    ERR = SQRT(SUM(ABS(XCALC-XEX)**2)) / SQRT(SUM(ABS(XEX)**2))
  END FUNCTION FROB_ERR

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_DENSE(N, NRHS, NB)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    COMPLEX(WP), ALLOCATABLE :: A(:,:), XEX(:,:), B(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT, T_SOLVE, RSS_PEAK, RSS_FINAL
    INTEGER :: INFO
    REAL(WP) :: ERR

    CALL SET_FIXED_SEED()
    ALLOCATE(A(N,N))
    CALL BUILD_MATRIX(N, A, T_ASM)

    ALLOCATE(XEX(N,NRHS), B(N,NRHS), IPIV(N))
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A, N, XEX, N, CZERO, B, N)

    T0 = OMP_GET_WTIME_L()
    CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO)
    T1 = OMP_GET_WTIME_L()
    T_FACT = T1 - T0
    RSS_PEAK = KB_TO_MB(GET_RSS_KB())

    T0 = OMP_GET_WTIME_L()
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, B, N)
    T1 = OMP_GET_WTIME_L()
    T_SOLVE = T1 - T0

    ERR = FROB_ERR(N, NRHS, B, XEX)
    RSS_FINAL = KB_TO_MB(GET_RSS_KB())

    ! Pas de compression dans ce cas (reference C pure) : temps de
    ! compression nul par construction, pas mesure separement.
    CALL REPORT('C) Direct dense, sans decoupage en blocs', T_ASM, 0.0D0, T_FACT, T_SOLVE, ERR, &
                8.0D0*REAL(N,8)*REAL(N,8)/1024.0D0/1024.0D0, RSS_PEAK, RSS_FINAL)

    DEALLOCATE(A, XEX, B, IPIV)
  END SUBROUTINE RUN_DENSE

END PROGRAM BENCH_SEQ
