!===============================================================================
! Test de MPI_MATVEC_2D (etape 6b) : matrice A et vecteur x connus par une
! formule deterministe (chaque processus ne calcule que SA part -- A comme
! en 6a, x seulement sur les processus MYROW=0), produit y=A*x distribue,
! rassemblement sur le rang 0 (GATHER_ROW_VECTOR), compare a une reference
! SEQUENTIELLE independante (A et x recalcules integralement sur le rang 0
! via les MEMES formules, CGEMV direct) -- meme methodologie que partout
! ailleurs dans ce projet. Tolerance non nulle cette fois (contrairement a
! 6a) : la sommation en virgule flottante n'est pas effectuee dans le meme
! ordre en distribue (reduction par blocs de colonnes) qu'en sequentiel
! (boucle simple), une difference de l'ordre de la precision machine est
! normale, pas un bug.
!
! Usage : mpirun -np <PR*PC> ./test_mpi_matvec PR PC N NB
!===============================================================================
PROGRAM TEST_MPI_MATVEC
  USE MPI_MATVEC_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, IERR, MLOC, NLOC, I, J, GI, GJ, MYID0
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), X_LOCAL(:), Y_LOCAL(:)
  COMPLEX(WP), ALLOCATABLE :: Y_GLOBAL(:), A_FULL(:,:), X_FULL(:), Y_REF(:)
  COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
  REAL(WP) :: ERR, NORMY

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 4) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_matvec PR PC N NB'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
  NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)

  ! --- A : chaque processus calcule SA part (comme en 6a), formule
  !     deterministe simple et bien definie ---
  ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
  DO J = 1, NLOC
    GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      A_LOCAL(I,J) = CMPLX(1.0_WP/(1.0_WP+ABS(REAL(GI-GJ,WP))), &
                            0.1_WP*REAL(MOD(GI+GJ,5),WP), KIND=WP)
    END DO
  END DO

  ! --- x : seuls les processus MYROW=0 calculent leur part (avant diffusion) ---
  ALLOCATE(X_LOCAL(MAX(NLOC,1)), Y_LOCAL(MAX(MLOC,1)))
  X_LOCAL = CZERO
  IF (GRID%MYROW .EQ. 0) THEN
    DO J = 1, NLOC
      GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
      X_LOCAL(J) = CMPLX(REAL(GJ,WP), -0.5_WP*REAL(GJ,WP), KIND=WP)
    END DO
  END IF

  CALL MPI_MATVEC_2D(GRID, MLOC, NLOC, A_LOCAL, X_LOCAL, Y_LOCAL)

  ALLOCATE(Y_GLOBAL(N))
  CALL GATHER_ROW_VECTOR(GRID, N, NB, Y_LOCAL, MAX(MLOC,1), 0, Y_GLOBAL)

  IF (MYID0 .EQ. 0) THEN
    ! --- Reference sequentielle independante (memes formules, recalculees
    !     integralement) ---
    ALLOCATE(A_FULL(N,N), X_FULL(N), Y_REF(N))
    DO GJ = 1, N
      X_FULL(GJ) = CMPLX(REAL(GJ,WP), -0.5_WP*REAL(GJ,WP), KIND=WP)
      DO GI = 1, N
        A_FULL(GI,GJ) = CMPLX(1.0_WP/(1.0_WP+ABS(REAL(GI-GJ,WP))), &
                               0.1_WP*REAL(MOD(GI+GJ,5),WP), KIND=WP)
      END DO
    END DO
    CALL CGEMV('No transpose', N, N, CONE, A_FULL, N, X_FULL, 1, CZERO, Y_REF, 1)

    NORMY = SQRT(SUM(ABS(Y_REF)**2))
    ERR = SQRT(SUM(ABS(Y_GLOBAL-Y_REF)**2)) / NORMY

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||y_distribue - y_reference||/||y_reference|| = ', ERR
    IF (ERR .LT. 1.0E-5_WP) THEN
      WRITE(6,*) '=== TEST MPI_MATVEC : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_MATVEC : ECHEC ==='
    END IF
    DEALLOCATE(A_FULL, X_FULL, Y_REF)
  END IF

  DEALLOCATE(A_LOCAL, X_LOCAL, Y_LOCAL, Y_GLOBAL)
  CALL MPI_FINALIZE(IERR)

END PROGRAM TEST_MPI_MATVEC
