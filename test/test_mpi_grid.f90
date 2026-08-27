!===============================================================================
! Test de l'infrastructure de distribution bloc-cyclique 2D (etape 6a) :
! chaque processus calcule INDEPENDAMMENT (sans aucune communication) sa
! propre part locale d'une matrice de test DETERMINISTE (A(i,j) =
! i*10000+j, en indices GLOBAUX 1-based -- valeur entiere exacte,
! representable exactement en simple precision jusqu'a N~1600, donc une
! comparaison d'EGALITE EXACTE est possible, pas juste "petite erreur").
! On rassemble ensuite sur le rang 0 (GATHER_2D) et on verifie que CHAQUE
! element reconstruit correspond EXACTEMENT a la formule -- toute erreur
! dans NUMROC/INDXL2G signalerait immediatement un decalage/trou/
! chevauchement.
!
! Usage : mpirun -np <PR*PC> ./test_mpi_grid <PR> <PC> <N> <NB>
!
! A tester avec plusieurs formes de grille pour un meme nombre de
! processus (ex. -np 6 : PR=2,PC=3 puis PR=3,PC=2 puis PR=1,PC=6...) et
! plusieurs N/NB (y compris N non multiple de NB, et NB > N/NPROCS pour
! forcer des blocs de tailles inegales entre processus).
!===============================================================================
PROGRAM TEST_MPI_GRID
  USE MPI_GRID_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, IERR, MLOC, NLOC, I, J, GI, GJ, MYID0
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_GLOBAL(:,:)
  LOGICAL :: OK
  INTEGER :: SUM_MLOC, SUM_NLOC, R, C

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 4) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_grid PR PC N NB'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (GRID%MYID .EQ. 0) WRITE(6,'(A,I0,A,I0,A,I0)') &
         ' ERREUR : PR*PC=', PR*PC, ' ne correspond pas au nombre de', &
         ' processus lances=', GRID%NPROCS
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
  NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)
  WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' [rang ', GRID%MYID, '] grille(', &
       GRID%MYROW, ',', GRID%MYCOL, ') MLOC=', MLOC, ' NLOC=', NLOC

  ! --- Chaque processus remplit SA PART, sans aucune communication ---
  ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
  DO J = 1, NLOC
    GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      A_LOCAL(I,J) = CMPLX(REAL(GI,WP)*10000.0_WP + REAL(GJ,WP), 0.0_WP, KIND=WP)
    END DO
  END DO

  ! --- Rassemblement sur le rang 0 et verification ---
  IF (GRID%MYID .EQ. 0) ALLOCATE(A_GLOBAL(N,N))
  CALL GATHER_2D(GRID, N, NB, A_LOCAL, MAX(MLOC,1), MAX(NLOC,1), 0, A_GLOBAL)

  IF (GRID%MYID .EQ. 0) THEN
    OK = .TRUE.
    DO GJ = 1, N
      DO GI = 1, N
        IF (A_GLOBAL(GI,GJ) .NE. CMPLX(REAL(GI,WP)*10000.0_WP+REAL(GJ,WP), &
                                        0.0_WP, KIND=WP)) THEN
          WRITE(6,'(A,I0,A,I0,A)') '  ECHEC a (', GI, ',', GJ, ') : valeur incorrecte'
          OK = .FALSE.
        END IF
      END DO
    END DO

    ! Sanity check supplementaire : NUMROC doit sommer exactement a N
    ! dans chaque dimension.
    SUM_MLOC = 0
    DO R = 0, GRID%PR-1
      SUM_MLOC = SUM_MLOC + NUMROC(N, NB, R, GRID%PR)
    END DO
    SUM_NLOC = 0
    DO C = 0, GRID%PC-1
      SUM_NLOC = SUM_NLOC + NUMROC(N, NB, C, GRID%PC)
    END DO
    IF (SUM_MLOC .NE. N .OR. SUM_NLOC .NE. N) THEN
      WRITE(6,'(A,I0,A,I0,A,I0)') '  ECHEC NUMROC : somme lignes=', SUM_MLOC, &
           ', somme colonnes=', SUM_NLOC, ', attendu N=', N
      OK = .FALSE.
    END IF

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB
    IF (OK) THEN
      WRITE(6,*) '=== TEST MPI_GRID : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_GRID : ECHEC ==='
    END IF
  END IF

  DEALLOCATE(A_LOCAL)
  IF (GRID%MYID .EQ. 0) DEALLOCATE(A_GLOBAL)
  CALL MPI_FINALIZE(IERR)

END PROGRAM TEST_MPI_GRID
