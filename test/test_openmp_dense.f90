!===============================================================================
! Case "dense + OpenMP" de la grille dense/BLR x OpenMP/MPI-sans-pivot/
! MPI-avec-pivot : contrairement a blr_lu_mod.f90 (etape 5), block_lu_mod.f90
! (etape 1, SANS compression) n'a AUCUNE directive !$OMP a lui -- un seul
! gros CGEMM par panneau, deja parallelise EN INTERNE par la MKL threadee
! (voir l'explication en tete de block_lu_mod.f90). Ce test MESURE que
! cette acceleration existe reellement (via OMP_SET_NUM_THREADS, qui
! pilote directement le threading interne de MKL) plutot que de se
! contenter de l'affirmer -- meme discipline que test_openmp.f90 (etape 5,
! version BLR).
!
! Usage : ./test_openmp_dense N NB
!===============================================================================
PROGRAM TEST_OPENMP_DENSE
  USE BLOCK_LU_MOD
  USE OMP_LIB
  IMPLICIT NONE

  INTEGER :: N, NB, INFO, NCORES, NT, I
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  COMPLEX(WP), ALLOCATABLE :: A0(:,:), A(:,:)
  INTEGER, ALLOCATABLE     :: IPIV(:), THREAD_COUNTS(:)
  COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
  DOUBLE PRECISION :: T0, T1
  DOUBLE PRECISION, ALLOCATABLE :: TIMES(:)
  REAL(WP) :: ERR_REF, ERR

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    WRITE(6,*) 'Usage: ./test_openmp_dense N NB'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB

  CALL SET_FIXED_SEED()
  ALL_OK = .TRUE.

  NCORES = OMP_GET_NUM_PROCS()
  WRITE(6,'(A,I0,A,I0,A,I0,A)') ' N = ', N, ' , NB = ', NB, &
       ' , coeurs logiques disponibles = ', NCORES, ' (nproc)'

  ALLOCATE(THREAD_COUNTS(20))
  NT = 0
  I = 1
  DO WHILE (I .LE. NCORES)
    NT = NT + 1
    THREAD_COUNTS(NT) = I
    I = I * 2
  END DO
  IF (THREAD_COUNTS(NT) .NE. NCORES) THEN
    NT = NT + 1
    THREAD_COUNTS(NT) = NCORES
  END IF
  ALLOCATE(TIMES(NT))

  ! Matrice dense a diagonale dominante (pas besoin d'etre compressible
  ! ici -- pas de compression dans ce module).
  ALLOCATE(A0(N,N))
  CALL RANDOM_COMPLEX(N, N, A0)
  DO I = 1, N
    A0(I,I) = A0(I,I) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
  END DO

  ALLOCATE(A(N,N), IPIV(N))
  WRITE(6,*)
  WRITE(6,'(A)') ' Threads   Temps (s)   Speedup   Erreur relative (vs A0)'

  DO I = 1, NT
    CALL OMP_SET_NUM_THREADS(THREAD_COUNTS(I))
    A = A0
    T0 = OMP_GET_WTIME()
    CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO)
    T1 = OMP_GET_WTIME()
    TIMES(I) = T1 - T0

    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLOCK_LU_FACTOR, INFO=', INFO
      ALL_OK = .FALSE.
      CYCLE
    END IF

    CALL VERIFY_RESIDUAL(N, NB, A0, A, IPIV, ERR)
    IF (I .EQ. 1) ERR_REF = ERR

    WRITE(6,'(I6,A,F10.3,A,F7.2,A,ES12.4)') THREAD_COUNTS(I), '     ', &
         TIMES(I), '    ', TIMES(1)/TIMES(I), '    ', ERR
  END DO

  WRITE(6,*)
  DO I = 2, NT
    CALL OMP_SET_NUM_THREADS(THREAD_COUNTS(I))
    A = A0
    CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO)
    CALL VERIFY_RESIDUAL(N, NB, A0, A, IPIV, ERR)
    IF (ERR .GT. 10.0_WP*ERR_REF + 1.0E-5_WP) THEN
      WRITE(6,'(A,I0,A)') '  ECHEC : residu degrade avec ', THREAD_COUNTS(I), ' threads'
      ALL_OK = .FALSE.
    END IF
  END DO

  IF (TIMES(NT) .GE. TIMES(1)) THEN
    WRITE(6,*) 'ATTENTION : aucun gain de temps observe entre 1 et ', &
         THREAD_COUNTS(NT), ' threads (matrice trop petite, MKL threading', &
         ' pas declenche, ou machine trop chargee)'
  ELSE
    WRITE(6,'(A,I0,A,F5.2)') ' Speedup obtenu (1 -> ', THREAD_COUNTS(NT), &
         ' threads), SANS AUCUNE directive OMP a nous : x', TIMES(1)/TIMES(NT)
  END IF

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TEST OPENMP_DENSE : OK (correction preservee) ==='
  ELSE
    WRITE(6,*) '=== TEST OPENMP_DENSE : ECHEC ==='
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

  SUBROUTINE VERIFY_RESIDUAL(N, NB, A0, A_FACT, IPIV, ERR)
    INTEGER, INTENT(IN)        :: N, NB
    COMPLEX(WP), INTENT(IN)    :: A0(N,N), A_FACT(N,N)
    INTEGER, INTENT(IN)        :: IPIV(N)
    REAL(WP), INTENT(OUT)      :: ERR
    COMPLEX(WP), ALLOCATABLE :: XEX(:,:), B(:,:), XCALC(:,:)
    COMPLEX(WP), PARAMETER :: CONE_L = (1.0_WP,0.0_WP), CZERO_L = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER :: NRHS = 1

    ALLOCATE(XEX(N,NRHS), B(N,NRHS), XCALC(N,NRHS))
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    CALL CGEMM('N','N', N, NRHS, N, CONE_L, A0, N, XEX, N, CZERO_L, B, N)
    XCALC = B
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_FACT, N, IPIV, XCALC, N)
    ERR = SQRT(SUM(ABS(XCALC(:,1)-XEX(:,1))**2)) / SQRT(SUM(ABS(XEX(:,1))**2))
    DEALLOCATE(XEX, B, XCALC)
  END SUBROUTINE VERIFY_RESIDUAL

END PROGRAM TEST_OPENMP_DENSE
