!===============================================================================
! Test de BLOCK_LU_FACTOR/BLOCK_LU_SOLVE (etape 1, sans compression) :
!
!   1) matrice A dense, complexe simple precision, N~500-2000, a diagonale
!      dominante (aleatoire) -- QUE de la correction d'algorithme visee ici,
!      pas encore de compressibilite : c'est pour ca qu'on utilise une
!      matrice aleatoire generique plutot qu'un noyau geometrique lisse
!      (voir la memoire du projet STRUMPACK sur les problemes de
!      conditionnement des noyaux lisses -- non pertinent a ce stade).
!   2) solution exacte X connue, second membre B = A*X (CGEMM).
!   3) copie A2 = A, on factorise/resout DEUX FOIS EN PARALLELE :
!        - avec NOTRE BLOCK_LU_FACTOR/BLOCK_LU_SOLVE (sur A)
!        - avec CGETRF/CGETRS de LAPACK directement (sur A2), comme
!          reference independante
!   4) on compare les deux solutions a X (residu relatif) ET entre elles.
!
! Usage : ./test_block_lu N NB
!===============================================================================
PROGRAM TEST_BLOCK_LU
  USE BLOCK_LU_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, I, J, INFO
  CHARACTER(LEN=32) :: ARG

  COMPLEX(WP), ALLOCATABLE :: A(:,:), A2(:,:), X(:,:), B(:,:), B2(:,:)
  INTEGER, ALLOCATABLE     :: IPIV(:), IPIV2(:)
  REAL(WP), ALLOCATABLE    :: RPART(:,:), IPART(:,:)

  COMPLEX(WP), PARAMETER :: CONE = (1.0_WP, 0.0_WP), CZERO = (0.0_WP, 0.0_WP)
  INTEGER, PARAMETER :: NRHS = 1

  REAL(WP) :: ERR_OURS, ERR_LAPACK, ERR_BETWEEN, NORMX

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    WRITE(6,*) 'Usage: ./test_block_lu N NB'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB

  WRITE(6,'(A,I0,A,I0)') ' N = ', N, ' , NB (taille de bloc) = ', NB

  ALLOCATE(A(N,N), A2(N,N), X(N,NRHS), B(N,NRHS), B2(N,NRHS))
  ALLOCATE(IPIV(N), IPIV2(N))
  ALLOCATE(RPART(N,N), IPART(N,N))

  ! ---------------------------------------------------------------------
  ! Matrice aleatoire a diagonale dominante (garantit l'existence d'une
  ! factorisation LU stable sans pivot exotique -- but ici : valider
  ! l'ORCHESTRATION par blocs, pas la robustesse numerique).
  ! ---------------------------------------------------------------------
  CALL SET_FIXED_SEED()   ! reproductible d'une execution a l'autre
  CALL RANDOM_NUMBER(RPART)
  CALL RANDOM_NUMBER(IPART)
  DO J = 1, N
    DO I = 1, N
      A(I,J) = CMPLX(RPART(I,J)-0.5_WP, IPART(I,J)-0.5_WP, KIND=WP)
    END DO
    A(J,J) = A(J,J) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)  ! dominance diagonale
  END DO
  A2 = A

  CALL RANDOM_NUMBER(RPART(:,1:NRHS))
  CALL RANDOM_NUMBER(IPART(:,1:NRHS))
  X(:,1) = CMPLX(RPART(:,1), IPART(:,1), KIND=WP)
  NORMX = SQRT(SUM(ABS(X(:,1))**2))

  CALL CGEMM('No transpose', 'No transpose', N, NRHS, N, CONE, A, N, &
              X, N, CZERO, B, N)
  B2 = B

  ! ---------------------------------------------------------------------
  ! Notre LU par blocs
  ! ---------------------------------------------------------------------
  CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO)
  IF (INFO .NE. 0) THEN
    WRITE(6,'(A,I0)') ' ERREUR BLOCK_LU_FACTOR, INFO = ', INFO
    STOP 1
  END IF
  CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, B, N)

  ! ---------------------------------------------------------------------
  ! Reference LAPACK directe (CGETRF/CGETRS non bloque "boite noire")
  ! ---------------------------------------------------------------------
  CALL CGETRF(N, N, A2, N, IPIV2, INFO)
  IF (INFO .NE. 0) THEN
    WRITE(6,'(A,I0)') ' ERREUR CGETRF (reference), INFO = ', INFO
    STOP 1
  END IF
  CALL CGETRS('No transpose', N, NRHS, A2, N, IPIV2, B2, N, INFO)

  ! ---------------------------------------------------------------------
  ! Comparaisons
  ! ---------------------------------------------------------------------
  ERR_OURS    = SQRT(SUM(ABS(B(:,1)  - X(:,1))**2)) / NORMX
  ERR_LAPACK  = SQRT(SUM(ABS(B2(:,1) - X(:,1))**2)) / NORMX
  ERR_BETWEEN = SQRT(SUM(ABS(B(:,1)  - B2(:,1))**2)) / NORMX

  WRITE(6,*) '--- Resultats ---'
  WRITE(6,'(A,ES12.4)') ' ||X_notre_LU   - X_exact||/||X_exact|| = ', ERR_OURS
  WRITE(6,'(A,ES12.4)') ' ||X_LAPACK     - X_exact||/||X_exact|| = ', ERR_LAPACK
  WRITE(6,'(A,ES12.4)') ' ||X_notre_LU   - X_LAPACK||/||X_exact|| = ', ERR_BETWEEN

  IF (ERR_OURS .LT. 1.0E-4_WP .AND. ERR_BETWEEN .LT. 1.0E-4_WP) THEN
    WRITE(6,*) 'TEST : OK'
  ELSE
    WRITE(6,*) 'TEST : ECHEC'
    STOP 1
  END IF

  DEALLOCATE(A, A2, X, B, B2, IPIV, IPIV2, RPART, IPART)

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

END PROGRAM TEST_BLOCK_LU
