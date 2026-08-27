!===============================================================================
! EXEMPLE CANONIQUE -- reference C (block_lu_mod.f90) : LU dense par blocs
! AVEC pivotage partiel complet (comme LAPACK CGETRF/CGETRS, mais
! reimplemente et decoupe en blocs "from scratch"), sequentiel/OpenMP
! (le parallelisme vient de la bibliotheque BLAS/LAPACK sous-jacente,
! MKL ou OpenBLAS selon la chaine de compilation -- voir INSTALL.md).
!
! C'est la reference la plus simple du projet : aucune compression, une
! seule matrice dense en memoire. Utile comme point de depart avant de
! passer a la version compressee (examples/example_blr_sequentiel.f90)
! ou distribuee (examples/example_mpi_dense.f90).
!
! Usage : ./example_dense [N] [NB]   (defauts : 2000 200)
!===============================================================================
PROGRAM EXAMPLE_DENSE_SEQUENTIEL
  USE BLOCK_LU_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, INFO, I, J, N_SWAPS
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A(:,:), B(:,:)
  INTEGER, ALLOCATABLE     :: IPIV(:)
  INTEGER, PARAMETER :: NRHS = 1

  N = 2000; NB = 200
  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  END IF

  WRITE(6,'(A,I0,A,I0)') ' N=', N, ' NB (taille de bloc) = ', NB

  ALLOCATE(A(N,N), B(N,NRHS), IPIV(N))

  ! Matrice arbitraire (a diagonale dominante ici, mais BLOCK_LU_FACTOR
  ! pivote reellement -- une matrice quelconque non singuliere convient).
  ! Dans un cas reel, A viendrait des donnees du probleme.
  DO J = 1, N
    DO I = 1, N
      A(I,J) = CMPLX(1.0_WP/(1.0_WP+REAL(ABS(I-J),WP)), 0.1_WP, KIND=WP)
    END DO
    A(J,J) = A(J,J) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
  END DO
  DO I = 1, N
    B(I,1) = CMPLX(REAL(I,WP)*0.01_WP, -0.5_WP, KIND=WP)
  END DO

  ! --- factorisation (A ecrase par L\U, IPIV = permutation des lignes) --
  CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO, N_SWAPS=N_SWAPS)
  IF (INFO .NE. 0) THEN
    WRITE(6,'(A,I0)') ' ERREUR BLOCK_LU_FACTOR, INFO=', INFO
    STOP 1
  END IF
  WRITE(6,'(A,I0)') ' Nombre de permutations de lignes (pivot) : ', N_SWAPS

  ! --- resolution, B ecrase par la solution ------------------------------
  CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, B, N)

  WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X :'
  DO I = 1, MIN(5, N)
    WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ') = ', REAL(B(I,1)), ' + i*', AIMAG(B(I,1))
  END DO

  DEALLOCATE(A, B, IPIV)

END PROGRAM EXAMPLE_DENSE_SEQUENTIEL
