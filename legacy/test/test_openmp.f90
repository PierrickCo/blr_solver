!===============================================================================
! ETAPE 5 : verifie que le parallelisme OpenMP ajoute a blr_lu_mod.f90 (1)
! ne casse pas la correction (residu identique, aux fluctuations de
! sommation flottante MKL pres, quel que soit le nombre de threads) et (2)
! apporte un VRAI gain de temps -- une directive OpenMP qui compile mais
! n'accelere rien serait un echec silencieux facile a manquer si on ne
! mesure pas explicitement.
!
! Methode : CALL OMP_SET_NUM_THREADS(n) pour n=1,2,4,...,NCORES avant
! chaque factorisation (pas besoin de relancer l'executable avec des
! variables d'environnement differentes), matrice assez grande et assez
! de blocs pour que la parallelisation ait quelque chose a repartir.
!
! Usage : ./test_openmp N NB TOL
!===============================================================================
PROGRAM TEST_OPENMP
  USE BLR_LU_MOD
  USE OMP_LIB
  IMPLICIT NONE

  INTEGER :: N, NB, INFO, NCORES, NT, I
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  COMPLEX(WP), ALLOCATABLE :: A0(:,:), A(:,:)
  INTEGER, ALLOCATABLE     :: IPIV(:), THREAD_COUNTS(:)
  COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
  INTEGER, PARAMETER       :: RGLOBAL = 6
  COMPLEX(WP), ALLOCATABLE :: UG(:,:), VG(:,:)
  INTEGER :: N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T0, T1
  DOUBLE PRECISION, ALLOCATABLE :: TIMES(:)
  REAL(WP) :: ERR_REF, ERR

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    WRITE(6,*) 'Usage: ./test_openmp N NB TOL'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL

  CALL SET_FIXED_SEED()
  ALL_OK = .TRUE.

  NCORES = OMP_GET_NUM_PROCS()
  WRITE(6,'(A,I0,A,I0,A,ES10.3,A,I0,A)') ' N = ', N, ' , NB = ', NB, &
       ' , TOL = ', TOL, ' , coeurs logiques disponibles = ', NCORES, &
       ' (nproc)'

  ! Liste des nombres de threads a essayer : 1, 2, 4, ... jusqu'a NCORES
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

  ! Matrice compressible (meme recette que test_blr_lu.f90) -- assez
  ! grande pour qu'il y ait un nombre de blocs significatif a repartir
  ! entre threads.
  ALLOCATE(A0(N,N), UG(N,RGLOBAL), VG(RGLOBAL,N))
  CALL RANDOM_COMPLEX(N, RGLOBAL, UG)
  CALL RANDOM_COMPLEX(RGLOBAL, N, VG)
  CALL CGEMM('N','N', N, N, RGLOBAL, CONE, UG, N, VG, RGLOBAL, CZERO, A0, N)
  DO I = 1, N
    A0(I,I) = A0(I,I) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
  END DO
  DEALLOCATE(UG, VG)

  ALLOCATE(A(N,N), IPIV(N))
  WRITE(6,*)
  WRITE(6,'(A)') ' Threads   Temps (s)   Speedup   Erreur relative (vs A0)'

  DO I = 1, NT
    CALL OMP_SET_NUM_THREADS(THREAD_COUNTS(I))
    A = A0
    T0 = OMP_GET_WTIME()
    CALL BLR_LU_FACTOR(N, NB, A, N, IPIV, TOL, INFO, N_COMPRESSED, &
                        N_OFFDIAG, RANK_SUM)
    T1 = OMP_GET_WTIME()
    TIMES(I) = T1 - T0

    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLR_LU_FACTOR, INFO=', INFO
      ALL_OK = .FALSE.
      CYCLE
    END IF

    ! Verification de coherence : residu de factorisation directe
    ! ||A0 - L*U (recompose)||... simplifie ici via un solve/verif rapide
    ! a la place (methodologie identique aux autres tests : residu contre
    ! A0 non compressee, calcule une seule fois avec un second membre
    ! connu -- voir CALL VERIFY_RESIDUAL ci-dessous).
    CALL VERIFY_RESIDUAL(N, NB, A0, A, IPIV, ERR)
    IF (I .EQ. 1) ERR_REF = ERR

    WRITE(6,'(I6,A,F10.3,A,F7.2,A,ES12.4)') THREAD_COUNTS(I), '     ', &
         TIMES(I), '    ', TIMES(1)/TIMES(I), '    ', ERR
  END DO

  WRITE(6,*)
  ! Verification de coherence numerique : le residu ne doit pas se
  ! degrader significativement quand le nombre de threads change (aucune
  ! section critique manquante ne devrait produire un resultat different).
  DO I = 2, NT
    CALL OMP_SET_NUM_THREADS(THREAD_COUNTS(I))
    A = A0
    CALL BLR_LU_FACTOR(N, NB, A, N, IPIV, TOL, INFO, N_COMPRESSED, &
                        N_OFFDIAG, RANK_SUM)
    CALL VERIFY_RESIDUAL(N, NB, A0, A, IPIV, ERR)
    IF (ERR .GT. 10.0_WP*ERR_REF + 1.0E-5_WP) THEN
      WRITE(6,'(A,I0,A)') '  ECHEC : residu degrade avec ', THREAD_COUNTS(I), ' threads'
      ALL_OK = .FALSE.
    END IF
  END DO

  IF (TIMES(NT) .GE. TIMES(1)) THEN
    WRITE(6,*) 'ATTENTION : aucun gain de temps observe entre 1 et ', &
         THREAD_COUNTS(NT), ' threads (voir discussion : matrice/NB trop', &
         ' petits, overhead domine, ou machine trop chargee)'
  ELSE
    WRITE(6,'(A,I0,A,F5.2)') ' Speedup obtenu (1 -> ', THREAD_COUNTS(NT), &
         ' threads) : x', TIMES(1)/TIMES(NT)
  END IF

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TEST OPENMP : OK (correction preservee) ==='
  ELSE
    WRITE(6,*) '=== TEST OPENMP : ECHEC ==='
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
  ! Residu ||A0*Xcalc - B||/||B|| pour un second membre B=A0*Xexact connu,
  ! Xcalc obtenu en resolvant avec les facteurs (A_FACT,IPIV) -- meme
  ! methodologie que test_blr_lu.f90, factorisee ici pour etre appelee a
  ! chaque nombre de threads sans dupliquer le code.
  !-----------------------------------------------------------------------
  SUBROUTINE VERIFY_RESIDUAL(N, NB, A0, A_FACT, IPIV, ERR)
    USE BLOCK_LU_MOD, ONLY: BLOCK_LU_SOLVE
    INTEGER, INTENT(IN)        :: N, NB
    COMPLEX(WP), INTENT(IN)    :: A0(N,N), A_FACT(N,N)
    INTEGER, INTENT(IN)        :: IPIV(N)
    REAL(WP), INTENT(OUT)      :: ERR
    COMPLEX(WP), ALLOCATABLE :: XEX(:,:), B(:,:), XCALC(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER :: NRHS = 1

    ALLOCATE(XEX(N,NRHS), B(N,NRHS), XCALC(N,NRHS))
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A0, N, XEX, N, CZERO, B, N)
    XCALC = B
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_FACT, N, IPIV, XCALC, N)
    ERR = SQRT(SUM(ABS(XCALC(:,1)-XEX(:,1))**2)) / SQRT(SUM(ABS(XEX(:,1))**2))
    DEALLOCATE(XEX, B, XCALC)
  END SUBROUTINE VERIFY_RESIDUAL

END PROGRAM TEST_OPENMP
