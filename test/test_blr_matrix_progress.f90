!===============================================================================
! Test dedie a la barre de progression du pipeline "matrice compressee"
! (blr_matrix_mod.f90) : verifie que les 3 barres (assemblage+compression,
! factorisation, resolution -- SHOW_PROGRESS=.TRUE.) s'affichent chacune
! independamment SANS casser le resultat (compare a une reference dense
! CGESV independante), et UNIQUEMENT sur un petit cas test (N=200), comme
! demande -- ce n'est pas fait sur les gros benchmarks.
!===============================================================================
PROGRAM TEST_BLR_MATRIX_PROGRESS
  USE BLR_MATRIX_MOD
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 200, NB = 40, NRHS = 5, RGLOBAL = 4
  REAL(WP), PARAMETER :: TOL = 1.0E-4_WP

  TYPE(BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: ADENSE(:,:), B(:,:), BREF(:,:)
  COMPLEX(WP), ALLOCATABLE :: XEX(:,:)
  INTEGER, ALLOCATABLE :: IPIV(:)
  INTEGER :: INFO, N_COMP, N_OFF, I, J
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACT
  REAL(WP) :: ERR

  WRITE(6,'(A)') ' === Test barre de progression (petit cas, N=200) ==='
  WRITE(6,*)

  ! --- reference dense independante (CGESV) --------------------------
  ALLOCATE(ADENSE(N,N), XEX(N,NRHS), BREF(N,NRHS))
  CALL SET_FIXED_SEED()
  CALL BUILD_DENSE(N, ADENSE)
  CALL RANDOM_COMPLEX(N, NRHS, XEX)
  BREF = MATMUL(ADENSE, XEX)

  ALLOCATE(IPIV(N))
  BLOCK
    COMPLEX(WP), ALLOCATABLE :: ACOPY(:,:)
    ALLOCATE(ACOPY(N,N))
    ACOPY = ADENSE
    BLOCK
      COMPLEX(WP), ALLOCATABLE :: BCOPY(:,:)
      ALLOCATE(BCOPY(N,NRHS))
      BCOPY = BREF
      CALL CGESV(N, NRHS, ACOPY, N, IPIV, BCOPY, N, INFO)
      IF (INFO .NE. 0) THEN
        WRITE(6,'(A,I0)') '  ERREUR CGESV (reference), INFO=', INFO
        STOP 1
      END IF
      DEALLOCATE(BCOPY)
    END BLOCK
    DEALLOCATE(ACOPY)
  END BLOCK
  DEALLOCATE(IPIV)

  ! --- pipeline compresse, AVEC barres de progression -----------------
  WRITE(6,'(A)') ' -- 1) Assemblage + compression --'
  CALL ASSEMBLE_BLR_MATRIX(N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, N_COMP, N_OFF, RANK_SUM, &
                            SHOW_PROGRESS=.TRUE.)

  WRITE(6,'(A)') ' -- 3) Factorisation LU --'
  CALL BLR_MATRIX_LU_FACTOR(A_BLR, TOL, INFO, N_COMP, N_OFF, RANK_SUM, T_FACT, &
                             SHOW_PROGRESS=.TRUE.)
  IF (INFO .NE. 0) THEN
    WRITE(6,'(A,I0)') '  ERREUR BLR_MATRIX_LU_FACTOR, INFO=', INFO
    STOP 1
  END IF

  WRITE(6,'(A)') ' -- 4) Resolution --'
  ALLOCATE(B(N,NRHS))
  B = BREF
  CALL BLR_MATRIX_LU_SOLVE(A_BLR, NRHS, B, N, SHOW_PROGRESS=.TRUE.)

  ERR = SQRT(SUM(ABS(B-XEX)**2)) / SQRT(SUM(ABS(XEX)**2))
  WRITE(6,*)
  WRITE(6,'(A,ES12.4)') '  Erreur relative vs reference dense CGESV = ', ERR
  IF (ERR .LT. 1.0E-3_WP) THEN
    WRITE(6,*) '=== TEST BARRE DE PROGRESSION : OK ==='
  ELSE
    WRITE(6,*) '=== TEST BARRE DE PROGRESSION : ECHEC ==='
    STOP 1
  END IF

CONTAINS

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    BLK = ADENSE(ISTART:IEND, JSTART:JEND)
  END SUBROUTINE GET_BLOCK

  SUBROUTINE BUILD_DENSE(NN, A)
    INTEGER, INTENT(IN) :: NN
    COMPLEX(WP), INTENT(OUT) :: A(NN,NN)
    INTEGER :: II, JJ, R
    COMPLEX(WP) :: AVAL
    DO JJ = 1, NN
      DO II = 1, NN
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(II,R)*LOWRANK_V(R,JJ)
        END DO
        IF (II .EQ. JJ) AVAL = AVAL + CMPLX(REAL(NN,WP), 0.0_WP, KIND=WP)
        A(II,JJ) = AVAL
      END DO
    END DO
  END SUBROUTINE BUILD_DENSE

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

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

END PROGRAM TEST_BLR_MATRIX_PROGRESS
