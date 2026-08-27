!===============================================================================
! Test de BLR_LDLT_FACTOR (proposition 1 : LDLt + compression BLR) +
! BLOCK_LDLT_SOLVE (etape 4, REUTILISE TEL QUEL -- voir le commentaire en
! tete de blr_ldlt_mod.f90 sur pourquoi le format de sortie reste compatible).
!
! Cas testes :
!   1) Matrice COMPRESSIBLE et symetrique (A = c*I + (U*V + (U*V)^T), rang
!      global <= 2*RGLOBAL -- meme principe que test_blr_lu.f90, symetrise).
!      Attendu : residu de l'ordre de TOL (pas la precision machine), ET
!      une fraction significative de panneaux L21/DL21T compresses.
!   2) REGRESSION -- meme matrice, TOL=0 : doit redonner la precision
!      machine, comme l'etape 4 pure (BLOCK_LDLT_FACTOR).
!   3) Matrice SANS dominance diagonale + A(1,1)=0 force (exactement la
!      recette de test_ldlt.f90) : force de VRAIS pivots 2x2 de Bunch-
!      Kaufman -- verifie via le diagnostic N_2X2_PIVOTS>0 que la
!      reconstruction DL21T=D*L21^T pour le cas 2x2 (la partie la plus
!      delicate de ce module, ecrite a la main) est reellement exercisee,
!      pas seulement le cas 1x1.
!
! Usage : ./test_blr_ldlt N NB TOL
!===============================================================================
PROGRAM TEST_BLR_LDLT
  USE BLOCK_LDLT_MOD
  USE BLR_LDLT_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, INFO
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    WRITE(6,*) 'Usage: ./test_blr_ldlt N NB TOL'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL

  CALL SET_FIXED_SEED()
  ALL_OK = .TRUE.

  WRITE(6,'(A,I0,A,I0,A,ES10.3)') ' N = ', N, ' , NB = ', NB, ' , TOL = ', TOL
  WRITE(6,*)

  CALL RUN_COMPRESSIBLE('1) Compressible, TOL demande', TOL, .TRUE., &
                         MAX(TOL*100.0_WP, 1.0E-4_WP), ALL_OK)
  CALL RUN_COMPRESSIBLE('2) Regression : compressible, TOL=0 (== etape 4)', 0.0_WP, &
                         .FALSE., 1.0E-4_WP, ALL_OK)
  CALL RUN_2X2_FORCING('3) Sans dominance + A(1,1)=0 (force des pivots 2x2)', TOL, ALL_OK)

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TOUS LES TESTS BLR_LDLT : OK ==='
  ELSE
    WRITE(6,*) '=== AU MOINS UN TEST BLR_LDLT A ECHOUE ==='
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
  SUBROUTINE RUN_COMPRESSIBLE(LABEL, CASE_TOL, EXPECT_COMPRESSION, ERR_BOUND, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    REAL(WP), INTENT(IN)         :: CASE_TOL, ERR_BOUND
    LOGICAL, INTENT(IN)          :: EXPECT_COMPRESSION
    LOGICAL, INTENT(INOUT)       :: OK

    INTEGER, PARAMETER :: RGLOBAL = 6, NRHS = 1
    COMPLEX(WP), ALLOCATABLE :: UG(:,:), VG(:,:), A(:,:), A_ORIG(:,:)
    COMPLEX(WP), ALLOCATABLE :: X(:,:), B(:,:), XCOMP(:,:)
    INTEGER, ALLOCATABLE     :: IPIV(:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: I, N_COMPRESSED, N_PANELS, N2X2
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR, NORMX, AVG_RANK, COMPRESS_FRAC

    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'

    ALLOCATE(A(N,N), A_ORIG(N,N), X(N,NRHS), B(N,NRHS), XCOMP(N,NRHS), IPIV(N))

    ! A = c*I + (U*V + (U*V)^T) : symetrique (pas hermitienne), rang
    ! global <= 2*RGLOBAL, diagonale dominante par construction.
    ALLOCATE(UG(N,RGLOBAL), VG(RGLOBAL,N))
    CALL RANDOM_COMPLEX(N, RGLOBAL, UG)
    CALL RANDOM_COMPLEX(RGLOBAL, N, VG)
    CALL CGEMM('N','N', N, N, RGLOBAL, CONE, UG, N, VG, RGLOBAL, CZERO, A, N)
    A = A + TRANSPOSE(A)
    DO I = 1, N
      A(I,I) = A(I,I) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
    END DO
    DEALLOCATE(UG, VG)
    A_ORIG = A

    CALL RANDOM_COMPLEX(N, NRHS, X)
    NORMX = SQRT(SUM(ABS(X(:,1))**2))
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, X, N, CZERO, B, N)

    CALL BLR_LDLT_FACTOR(N, NB, A, N, IPIV, CASE_TOL, INFO, &
                          N_COMPRESSED, N_PANELS, RANK_SUM, N2X2)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLR_LDLT_FACTOR, INFO = ', INFO
      OK = .FALSE.
      DEALLOCATE(A, A_ORIG, X, B, XCOMP, IPIV)
      RETURN
    END IF

    XCOMP = B
    CALL BLOCK_LDLT_SOLVE(N, NRHS, A, N, IPIV, XCOMP, N, INFO)

    ERR = SQRT(SUM(ABS(XCOMP(:,1)-X(:,1))**2)) / NORMX

    AVG_RANK = 0.0_WP
    COMPRESS_FRAC = 0.0_WP
    IF (N_PANELS .GT. 0) THEN
      AVG_RANK = REAL(RANK_SUM,WP) / REAL(N_PANELS,WP)
      COMPRESS_FRAC = 100.0_WP * REAL(N_COMPRESSED,WP) / REAL(N_PANELS,WP)
    END IF

    WRITE(6,'(A,ES12.4)') '  Erreur relative ||Xcalc-Xexact||/||Xexact|| = ', ERR
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Panneaux compresses : ', N_COMPRESSED, ' / ', &
         N_PANELS, ' (', COMPRESS_FRAC, ' %)'
    WRITE(6,'(A,F8.2)') '  Rang moyen = ', AVG_RANK

    IF (ERR .GT. ERR_BOUND) THEN
      WRITE(6,*) '  -> ECHEC (residu trop grand)'
      OK = .FALSE.
    ELSE IF (EXPECT_COMPRESSION .AND. N_COMPRESSED .EQ. 0) THEN
      WRITE(6,*) '  -> ECHEC (compression attendue mais aucun panneau compresse)'
      OK = .FALSE.
    ELSE IF (.NOT.EXPECT_COMPRESSION .AND. N_COMPRESSED .GT. 0) THEN
      WRITE(6,*) '  -> ECHEC (compression inattendue)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A_ORIG, X, B, XCOMP, IPIV)
  END SUBROUTINE RUN_COMPRESSIBLE

  !-----------------------------------------------------------------------
  ! Meme recette que test_ldlt.f90 (symetrique sans dominance + A(1,1)=0)
  ! -- ici pour verifier que le cas 2x2 de BLR_LDLT_FACTOR est correct,
  ! pas pour re-tester le pivotage lui-meme (deja fait a l'etape 4).
  ! Seuil calibre sur l'erreur de la reference LAPACK elle-meme, meme
  ! raison qu'a l'etape 4 (matrice deliberement marginale).
  !-----------------------------------------------------------------------
  SUBROUTINE RUN_2X2_FORCING(LABEL, CASE_TOL, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    REAL(WP), INTENT(IN)         :: CASE_TOL
    LOGICAL, INTENT(INOUT)       :: OK

    INTEGER, PARAMETER :: NRHS = 1
    COMPLEX(WP), ALLOCATABLE :: R(:,:), A(:,:), A2(:,:), A_ORIG(:,:)
    COMPLEX(WP), ALLOCATABLE :: X(:,:), B(:,:), B2(:,:), CWORK(:)
    INTEGER, ALLOCATABLE     :: IPIV(:), IPIV2(:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: I, J, INFO2, LWORK, N_COMPRESSED, N_PANELS, N2X2
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR_OURS, ERR_LAPACK, NORMX

    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'

    ALLOCATE(R(N,N), A(N,N), A2(N,N), A_ORIG(N,N))
    ALLOCATE(X(N,NRHS), B(N,NRHS), B2(N,NRHS), IPIV(N), IPIV2(N))

    CALL RANDOM_COMPLEX(N, N, R)
    DO J = 1, N
      DO I = 1, N
        A(I,J) = R(I,J) + R(J,I)
      END DO
    END DO
    DEALLOCATE(R)
    A(1,1) = CMPLX(0.0_WP, 0.0_WP, KIND=WP)
    A_ORIG = A
    A2 = A

    CALL RANDOM_COMPLEX(N, NRHS, X)
    NORMX = SQRT(SUM(ABS(X(:,1))**2))
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, X, N, CZERO, B, N)
    B2 = B

    CALL BLR_LDLT_FACTOR(N, NB, A, N, IPIV, CASE_TOL, INFO, &
                          N_COMPRESSED, N_PANELS, RANK_SUM, N2X2)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLR_LDLT_FACTOR, INFO=', INFO
      OK = .FALSE.
      DEALLOCATE(A, A2, A_ORIG, X, B, B2, IPIV, IPIV2)
      RETURN
    END IF
    CALL BLOCK_LDLT_SOLVE(N, NRHS, A, N, IPIV, B, N, INFO)

    ALLOCATE(CWORK(1))
    LWORK = -1
    CALL CSYTRF('Lower', N, A2, N, IPIV2, CWORK, LWORK, INFO2)
    LWORK = INT(REAL(CWORK(1)))
    DEALLOCATE(CWORK)
    ALLOCATE(CWORK(MAX(LWORK,1)))
    CALL CSYTRF('Lower', N, A2, N, IPIV2, CWORK, LWORK, INFO2)
    DEALLOCATE(CWORK)
    CALL CSYTRS('Lower', N, NRHS, A2, N, IPIV2, B2, N, INFO2)

    ERR_OURS   = SQRT(SUM(ABS(B(:,1) -X(:,1))**2)) / NORMX
    ERR_LAPACK = SQRT(SUM(ABS(B2(:,1)-X(:,1))**2)) / NORMX

    WRITE(6,'(A,ES12.4)') '  ||X_nous - X_exact||/||X_exact||   = ', ERR_OURS
    WRITE(6,'(A,ES12.4)') '  ||X_LAPACK - X_exact||/||X_exact|| = ', ERR_LAPACK
    WRITE(6,'(A,I0)')     '  Pivots 2x2 rencontres = ', N2X2

    IF (ERR_OURS .GT. MAX(5.0_WP*ERR_LAPACK, 1.0E-4_WP)) THEN
      WRITE(6,*) '  -> ECHEC (residu trop grand)'
      OK = .FALSE.
    ELSE IF (N2X2 .EQ. 0) THEN
      WRITE(6,*) '  -> ECHEC (aucun pivot 2x2 rencontre, le cas delicat n''est pas teste)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A2, A_ORIG, X, B, B2, IPIV, IPIV2)
  END SUBROUTINE RUN_2X2_FORCING

END PROGRAM TEST_BLR_LDLT
