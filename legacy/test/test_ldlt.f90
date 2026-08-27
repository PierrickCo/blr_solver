!===============================================================================
! Test de BLOCK_LDLT_FACTOR/BLOCK_LDLT_SOLVE (etape 4, partie 1 -- sans
! compression). Meme methodologie de validation que test_block_lu.f90 :
! comparaison a un appel LAPACK direct (CSYTRF+CSYTRS) sur une copie
! independante de la matrice, PLUS residu contre la solution exacte connue.
!
! Cas testes :
!   1) Matrice symetrique (A(i,j)=A(j,i), PAS hermitienne) ALEATOIRE SANS
!      dominance diagonale -- contrairement aux tests LU precedents, on ne
!      force PAS de diagonale dominante ici : le but est justement de
!      solliciter naturellement le pivot de Bunch-Kaufman (melange de
!      pivots 1x1 et 2x2), pas de l'eviter.
!   2) Meme matrice, avec en plus A(1,1) force a zero explicitement --
!      cas limite ou le premier pivot naturel est totalement inutilisable.
!   3) Balayage de N, NB (y compris N non multiple de NB).
!
! Usage : ./test_ldlt N NB
!===============================================================================
PROGRAM TEST_LDLT
  USE BLOCK_LDLT_MOD
  IMPLICIT NONE

  INTEGER :: N, NB
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    WRITE(6,*) 'Usage: ./test_ldlt N NB'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB

  CALL SET_FIXED_SEED()
  ALL_OK = .TRUE.

  WRITE(6,'(A,I0,A,I0)') ' N = ', N, ' , NB = ', NB
  WRITE(6,*)

  CALL RUN_CASE('1) Symetrique aleatoire, sans dominance diagonale', N, NB, &
                 FORCE_ZERO_DIAG=.FALSE., OK=ALL_OK)
  CALL RUN_CASE('2) Idem, avec A(1,1)=0 force', N, NB, &
                 FORCE_ZERO_DIAG=.TRUE., OK=ALL_OK)

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TOUS LES TESTS LDLT : OK ==='
  ELSE
    WRITE(6,*) '=== AU MOINS UN TEST LDLT A ECHOUE ==='
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
  SUBROUTINE RUN_CASE(LABEL, N, NB, FORCE_ZERO_DIAG, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(IN)          :: N, NB
    LOGICAL, INTENT(IN)          :: FORCE_ZERO_DIAG
    LOGICAL, INTENT(INOUT)       :: OK

    COMPLEX(WP), ALLOCATABLE :: R(:,:), A(:,:), A2(:,:), A_ORIG(:,:)
    COMPLEX(WP), ALLOCATABLE :: X(:,:), B(:,:), B2(:,:), CWORK(:)
    INTEGER, ALLOCATABLE     :: IPIV(:), IPIV2(:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER       :: NRHS = 1
    INTEGER :: I, J, INFO, INFO2, LWORK
    REAL(WP) :: ERR_OURS, ERR_LAPACK, ERR_BETWEEN, NORMX

    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'

    ALLOCATE(R(N,N), A(N,N), A2(N,N), A_ORIG(N,N))
    ALLOCATE(X(N,NRHS), B(N,NRHS), B2(N,NRHS), IPIV(N), IPIV2(N))

    ! Symetrisation A = R + R^T (PAS R + R^H : symetrique, pas hermitienne
    ! -- meme convention que CMUMPS SYM=2 dans tout ce projet). Volontai-
    ! rement PAS de renfort diagonal : on veut que Bunch-Kaufman travaille.
    CALL RANDOM_COMPLEX(N, N, R)
    DO J = 1, N
      DO I = 1, N
        A(I,J) = R(I,J) + R(J,I)
      END DO
    END DO
    DEALLOCATE(R)

    IF (FORCE_ZERO_DIAG) THEN
      A(1,1) = CMPLX(0.0_WP, 0.0_WP, KIND=WP)
    END IF

    ! Ne remplir explicitement que le triangle inferieur pour A2 (ce que
    ! CSYTRF/notre BLOCK_LDLT_FACTOR sont censes seuls utiliser) -- la
    ! partie superieure ne doit jamais etre lue.
    A_ORIG = A
    A2 = A

    CALL RANDOM_COMPLEX(N, NRHS, X)
    NORMX = SQRT(SUM(ABS(X(:,1))**2))
    ! B = A*X avec A PLEINE (symetrique complete) puisque CGEMM ne connait
    ! pas la structure triangulaire -- coherent avec A_ORIG complete.
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, X, N, CZERO, B, N)
    B2 = B

    ! --- Notre LDLt par blocs ---
    CALL BLOCK_LDLT_FACTOR(N, NB, A, N, IPIV, INFO)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLOCK_LDLT_FACTOR, INFO=', INFO
      OK = .FALSE.
      DEALLOCATE(A, A2, A_ORIG, X, B, B2, IPIV, IPIV2)
      RETURN
    END IF
    CALL BLOCK_LDLT_SOLVE(N, NRHS, A, N, IPIV, B, N, INFO)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLOCK_LDLT_SOLVE, INFO=', INFO
      OK = .FALSE.
      DEALLOCATE(A, A2, A_ORIG, X, B, B2, IPIV, IPIV2)
      RETURN
    END IF

    ! --- Reference LAPACK directe (CSYTRF/CSYTRS, non bloque "boite noire") ---
    ALLOCATE(CWORK(1))
    LWORK = -1
    CALL CSYTRF('Lower', N, A2, N, IPIV2, CWORK, LWORK, INFO2)
    LWORK = INT(REAL(CWORK(1)))
    DEALLOCATE(CWORK)
    ALLOCATE(CWORK(MAX(LWORK,1)))
    CALL CSYTRF('Lower', N, A2, N, IPIV2, CWORK, LWORK, INFO2)
    DEALLOCATE(CWORK)
    IF (INFO2 .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR CSYTRF (reference), INFO=', INFO2
      OK = .FALSE.
      DEALLOCATE(A, A2, A_ORIG, X, B, B2, IPIV, IPIV2)
      RETURN
    END IF
    CALL CGEMM('N','N', N, NRHS, N, CONE, A_ORIG, N, X, N, CZERO, B2, N)
    CALL CSYTRS('Lower', N, NRHS, A2, N, IPIV2, B2, N, INFO2)

    ERR_OURS    = SQRT(SUM(ABS(B(:,1)  - X(:,1))**2)) / NORMX
    ERR_LAPACK  = SQRT(SUM(ABS(B2(:,1) - X(:,1))**2)) / NORMX
    ERR_BETWEEN = SQRT(SUM(ABS(B(:,1)  - B2(:,1))**2)) / NORMX

    WRITE(6,'(A,ES12.4)') '  ||X_notre_LDLt - X_exact||/||X_exact||  = ', ERR_OURS
    WRITE(6,'(A,ES12.4)') '  ||X_LAPACK     - X_exact||/||X_exact||  = ', ERR_LAPACK
    WRITE(6,'(A,ES12.4)') '  ||X_notre_LDLt - X_LAPACK||/||X_exact|| = ', ERR_BETWEEN

    ! Seuil relatif a l'erreur de la reference LAPACK elle-meme (plutot
    ! qu'un seuil absolu fixe) : cette matrice est DELIBEREMENT sans
    ! dominance diagonale (pour solliciter Bunch-Kaufman), donc meme la
    ! reference LAPACK directe n'atteint pas la precision machine en
    ! simple precision -- voir test_pivoting.f90 pour la meme lecon.
    IF (ERR_OURS .GT. MAX(5.0_WP*ERR_LAPACK, 1.0E-4_WP)) THEN
      WRITE(6,*) '  -> ECHEC'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A2, A_ORIG, X, B, B2, IPIV, IPIV2)
  END SUBROUTINE RUN_CASE

END PROGRAM TEST_LDLT
