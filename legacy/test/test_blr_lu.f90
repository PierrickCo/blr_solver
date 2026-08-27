!===============================================================================
! Test de bout en bout de BLR_LU_FACTOR (etape 2) + BLOCK_LU_SOLVE
! (etape 1, reutilise tel quel puisque les facteurs restent stockes denses
! -- voir le commentaire en tete de blr_lu_mod.f90).
!
! Cas testes :
!   1) Matrice COMPRESSIBLE (A = c*I + U*V, rang global RGLOBAL petit,
!      meme principe que examples/cmumps_blr_lowrank.F90 : un correctif de
!      rang faible sur une diagonale dominante -- bien conditionnee par
!      construction, blocs hors-diagonale numeriquement de rang <=
!      RGLOBAL). Attendu : residu de l'ordre de TOL (PAS la precision
!      machine -- la factorisation BLR est approchee), ET une fraction
!      significative de blocs hors-diagonale effectivement compresses.
!   2) REGRESSION -- meme matrice compressible, TOL=0 (compression
!      desactivee) : doit redonner la precision machine, comme l'etape 1
!      pure (verifie qu'on n'a rien casse pour le cas non-compresse).
!   3) REGRESSION -- matrice ALEATOIRE INCOMPRESSIBLE (comme le test de
!      l'etape 1), TOL>0 : la compression doit se desactiver d'elle-meme
!      sur (quasi) tous les blocs (rien a compresser), residu proche de la
!      precision machine malgre TOL>0.
!
! Usage : ./test_blr_lu N NB TOL
!===============================================================================
PROGRAM TEST_BLR_LU
  USE BLOCK_LU_MOD
  USE BLR_LU_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, INFO
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    WRITE(6,*) 'Usage: ./test_blr_lu N NB TOL'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL

  CALL SET_FIXED_SEED()   ! reproductible d'une execution a l'autre
  ALL_OK = .TRUE.

  WRITE(6,'(A,I0,A,I0,A,ES10.3)') ' N = ', N, ' , NB = ', NB, ' , TOL = ', TOL
  WRITE(6,*)

  CALL RUN_CASE('1) Matrice compressible, TOL demande', N, NB, TOL, &
                 COMPRESSIBLE=.TRUE., EXPECT_COMPRESSION=.TRUE., &
                 ERR_BOUND=MAX(TOL*100.0_WP, 1.0E-4_WP), OK=ALL_OK)

  CALL RUN_CASE('2) Regression : compressible, TOL=0 (== etape 1)', N, NB, 0.0_WP, &
                 COMPRESSIBLE=.TRUE., EXPECT_COMPRESSION=.FALSE., &
                 ERR_BOUND=1.0E-4_WP, OK=ALL_OK)

  CALL RUN_CASE('3) Regression : matrice incompressible, TOL demande', N, NB, TOL, &
                 COMPRESSIBLE=.FALSE., EXPECT_COMPRESSION=.FALSE., &
                 ERR_BOUND=1.0E-4_WP, OK=ALL_OK)

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TOUS LES TESTS BLR_LU : OK ==='
  ELSE
    WRITE(6,*) '=== AU MOINS UN TEST BLR_LU A ECHOUE ==='
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

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_CASE(LABEL, N, NB, TOL, COMPRESSIBLE, EXPECT_COMPRESSION, &
                       ERR_BOUND, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(IN)          :: N, NB
    REAL(WP), INTENT(IN)         :: TOL, ERR_BOUND
    LOGICAL, INTENT(IN)          :: COMPRESSIBLE, EXPECT_COMPRESSION
    LOGICAL, INTENT(INOUT)       :: OK

    COMPLEX(WP), ALLOCATABLE :: A(:,:), A_ORIG(:,:), X(:,:), B(:,:), XCOMP(:,:)
    INTEGER, ALLOCATABLE     :: IPIV(:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER       :: NRHS = 1, RGLOBAL = 6
    COMPLEX(WP), ALLOCATABLE :: UG(:,:), VG(:,:)
    INTEGER :: I, N_COMPRESSED, N_OFFDIAG
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR, NORMB, AVG_RANK, COMPRESS_FRAC

    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'

    ALLOCATE(A(N,N), A_ORIG(N,N), X(N,NRHS), B(N,NRHS), XCOMP(N,NRHS), IPIV(N))

    IF (COMPRESSIBLE) THEN
      ! A = c*I + U*V, U/V aleatoires N x RGLOBAL / RGLOBAL x N -- rang
      ! global RGLOBAL, donc tout bloc hors-diagonale est numeriquement de
      ! rang <= RGLOBAL, quelle que soit sa taille.
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
    NORMB = 1.0_WP  ! recalcule ci-dessous a partir de B
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, X, N, CZERO, B, N)
    NORMB = SQRT(SUM(ABS(B(:,1))**2))

    CALL BLR_LU_FACTOR(N, NB, A, N, IPIV, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLR_LU_FACTOR, INFO = ', INFO
      OK = .FALSE.
      RETURN
    END IF

    XCOMP = B
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, XCOMP, N)

    ! Residu independant : ||A_ORIG*Xcalc - B|| / ||B|| (matrice NON
    ! compressee, meme methodologie que les exemples MUMPS/STRUMPACK)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, XCOMP, N, CZERO, X, N)
    ERR = SQRT(SUM(ABS(X(:,1)-B(:,1))**2)) / NORMB

    AVG_RANK = 0.0_WP
    COMPRESS_FRAC = 0.0_WP
    IF (N_OFFDIAG .GT. 0) THEN
      AVG_RANK = REAL(RANK_SUM,WP) / REAL(N_OFFDIAG,WP)
      COMPRESS_FRAC = 100.0_WP * REAL(N_COMPRESSED,WP) / REAL(N_OFFDIAG,WP)
    END IF

    WRITE(6,'(A,ES12.4)') '  Erreur relative ||A*Xcalc-B||/||B||   = ', ERR
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Blocs compresses : ', N_COMPRESSED, ' / ', &
         N_OFFDIAG, ' (', COMPRESS_FRAC, ' %)'
    WRITE(6,'(A,F8.2)') '  Rang moyen (blocs compresses inclus)  = ', AVG_RANK

    IF (ERR .GT. ERR_BOUND) THEN
      WRITE(6,*) '  -> ECHEC (residu trop grand)'
      OK = .FALSE.
    ELSE IF (EXPECT_COMPRESSION .AND. N_COMPRESSED .EQ. 0) THEN
      WRITE(6,*) '  -> ECHEC (compression attendue mais aucun bloc compresse)'
      OK = .FALSE.
    ELSE IF (.NOT.EXPECT_COMPRESSION .AND. N_COMPRESSED .GT. 0) THEN
      WRITE(6,*) '  -> ECHEC (compression inattendue)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A_ORIG, X, B, XCOMP, IPIV)
  END SUBROUTINE RUN_CASE

  !-----------------------------------------------------------------------
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

END PROGRAM TEST_BLR_LU
