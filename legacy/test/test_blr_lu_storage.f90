!===============================================================================
! Test de BLR_LU_FACTOR_COMPRESSED/BLR_LU_SOLVE_COMPRESSED (proposition 2 :
! stockage compresse des facteurs). Meme recette "diagonale + correction
! basse-rang" que test_blr_lu.f90 pour la correction, PLUS une mesure
! CONCRETE de l'empreinte memoire du stockage compresse contre l'equivalent
! dense N x N -- la promesse de cette proposition (reduction MEMOIRE, pas
! seulement temps de factorisation) doit etre MESUREE, pas seulement
! supposee (meme discipline que les mesures de temps de l'etape 5).
!
! Cas testes :
!   1) Matrice COMPRESSIBLE, TOL demande : residu de l'ordre de TOL, ET
!      empreinte memoire du stockage compresse < taille dense N x N.
!   2) REGRESSION -- meme matrice, TOL=0 : precision machine, ET empreinte
!      memoire >= taille dense (aucune compression, stockage plus gros
!      voire egal -- normal, chaque bloc reste dense).
!   3) REGRESSION -- matrice incompressible : compression ne doit
!      quasiment rien trouver (comme test_blr_lu.f90).
!
! Usage : ./test_blr_lu_storage N NB TOL
!===============================================================================
PROGRAM TEST_BLR_LU_STORAGE
  USE BLR_LU_STORAGE_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, INFO
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    WRITE(6,*) 'Usage: ./test_blr_lu_storage N NB TOL'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL

  CALL SET_FIXED_SEED()
  ALL_OK = .TRUE.

  WRITE(6,'(A,I0,A,I0,A,ES10.3)') ' N = ', N, ' , NB = ', NB, ' , TOL = ', TOL
  WRITE(6,*)

  CALL RUN_CASE('1) Compressible, TOL demande', TOL, &
                 COMPRESSIBLE=.TRUE., EXPECT_SMALLER=.TRUE., &
                 ERR_BOUND=MAX(TOL*100.0_WP, 1.0E-4_WP), OK=ALL_OK)
  CALL RUN_CASE('2) Regression : compressible, TOL=0', 0.0_WP, &
                 COMPRESSIBLE=.TRUE., EXPECT_SMALLER=.FALSE., &
                 ERR_BOUND=1.0E-4_WP, OK=ALL_OK)
  CALL RUN_CASE('3) Regression : incompressible, TOL demande', TOL, &
                 COMPRESSIBLE=.FALSE., EXPECT_SMALLER=.FALSE., &
                 ERR_BOUND=1.0E-4_WP, OK=ALL_OK)

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TOUS LES TESTS BLR_LU_STORAGE : OK ==='
  ELSE
    WRITE(6,*) '=== AU MOINS UN TEST BLR_LU_STORAGE A ECHOUE ==='
    STOP 1
  END IF

CONTAINS

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

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_CASE(LABEL, CASE_TOL, COMPRESSIBLE, EXPECT_SMALLER, ERR_BOUND, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    REAL(WP), INTENT(IN)         :: CASE_TOL, ERR_BOUND
    LOGICAL, INTENT(IN)          :: COMPRESSIBLE, EXPECT_SMALLER
    LOGICAL, INTENT(INOUT)       :: OK

    COMPLEX(WP), ALLOCATABLE :: A(:,:), A_ORIG(:,:), X(:,:), B(:,:), XCOMP(:,:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER       :: NRHS = 1, RGLOBAL = 6
    COMPLEX(WP), ALLOCATABLE :: UG(:,:), VG(:,:)
    TYPE(BLR_FACTORED_LU_T)  :: FACT
    INTEGER :: I, N_COMPRESSED, N_OFFDIAG
    INTEGER(8) :: RANK_SUM, MEM_COMPRESSED, MEM_DENSE
    REAL(WP) :: ERR, NORMB

    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'

    ALLOCATE(A(N,N), A_ORIG(N,N), X(N,NRHS), B(N,NRHS), XCOMP(N,NRHS))

    IF (COMPRESSIBLE) THEN
      ALLOCATE(UG(N,RGLOBAL), VG(RGLOBAL,N))
      CALL RANDOM_COMPLEX(N, RGLOBAL, UG)
      CALL RANDOM_COMPLEX(RGLOBAL, N, VG)
      CALL CGEMM('N','N', N, N, RGLOBAL, CONE, UG, N, VG, RGLOBAL, CZERO, A, N)
      DO I = 1, N
        A(I,I) = A(I,I) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
      END DO
      DEALLOCATE(UG, VG)
    ELSE
      CALL RANDOM_COMPLEX(N, N, A)
      DO I = 1, N
        A(I,I) = A(I,I) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
      END DO
    END IF
    A_ORIG = A

    CALL RANDOM_COMPLEX(N, NRHS, X)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, X, N, CZERO, B, N)
    NORMB = SQRT(SUM(ABS(B(:,1))**2))

    CALL BLR_LU_FACTOR_COMPRESSED(N, NB, A, N, CASE_TOL, FACT, INFO, &
                                   N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLR_LU_FACTOR_COMPRESSED, INFO = ', INFO
      OK = .FALSE.
      DEALLOCATE(A, A_ORIG, X, B, XCOMP)
      RETURN
    END IF

    XCOMP = B
    CALL BLR_LU_SOLVE_COMPRESSED(FACT, NRHS, XCOMP, N)

    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, XCOMP, N, CZERO, X, N)
    ERR = SQRT(SUM(ABS(X(:,1)-B(:,1))**2)) / NORMB

    MEM_COMPRESSED = BLR_FACTORED_LU_MEMORY_BYTES(FACT)
    MEM_DENSE = 8_8 * INT(N,8) * INT(N,8)

    WRITE(6,'(A,ES12.4)') '  Erreur relative ||A*Xcalc-B||/||B|| = ', ERR
    WRITE(6,'(A,I0,A,I0)') '  Blocs compresses : ', N_COMPRESSED, ' / ', N_OFFDIAG
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Memoire stockage compresse = ', MEM_COMPRESSED, &
         ' octets , dense N x N = ', MEM_DENSE, ' octets (', &
         100.0_WP*REAL(MEM_COMPRESSED,WP)/REAL(MEM_DENSE,WP), ' % du dense)'

    IF (ERR .GT. ERR_BOUND) THEN
      WRITE(6,*) '  -> ECHEC (residu trop grand)'
      OK = .FALSE.
    ELSE IF (EXPECT_SMALLER .AND. MEM_COMPRESSED .GE. MEM_DENSE) THEN
      WRITE(6,*) '  -> ECHEC (reduction memoire attendue mais absente)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A_ORIG, X, B, XCOMP)
  END SUBROUTINE RUN_CASE

END PROGRAM TEST_BLR_LU_STORAGE
