!===============================================================================
! EXEMPLE CANONIQUE -- reference MPI (mpi_lu_pivot_mod.f90 +
! mpi_lu_solve_2d_mod.f90) : LU dense DISTRIBUE (grille de processus
! PR x PC, distribution bloc-cyclique 2D via ScaLAPACK-like NUMROC/
! INDXL2G), avec pivotage partiel, resolution 2D complete (le second
! membre est lui aussi distribue par lignes ET par colonnes -- aucune
! replication residuelle sur un seul processus).
!
! But de ce fichier : montrer l'enchainement minimal (distribuer A et B,
! factoriser, resoudre, rassembler) sans toute l'instrumentation de
! mesure/validation croisee des drivers de test (voir
! test/test_mpi_lu_solve_2d.f90 pour la version qui verifie aussi la
! correction contre une reference CGESV independante sur le processus 0).
!
! Usage : mpirun -np <PR*PC> ./example_mpi_dense PR PC [N] [NB] [NRHS]
!         (defauts : N=2000 NB=200 NRHS=4)
!===============================================================================
PROGRAM EXAMPLE_MPI_DENSE
  USE MPI_LU_PIVOT_MOD
  USE MPI_LU_SOLVE_2D_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:)
  INTEGER :: MLOC, NLOC, NRHS_LOC, RSTART_MINE, I, J, GI, GJ, INFO, N_SWAPS

  N = 2000; NB = 200; NRHS = 4

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./example_mpi_dense PR PC [N] [NB] [NRHS]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS

  ! --- chaque processus ne possede que SA part locale de A et de B, selon
  !     la distribution bloc-cyclique 2D standard (NUMROC = nombre de
  !     lignes/colonnes locales, INDXL2G = indice global d'un indice local) --
  MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
  NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)
  NRHS_LOC = RHS_NUMROC(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ! Chaque processus calcule directement SA part locale via une formule
  ! deterministe sur les indices GLOBAUX -- pas de communication necessaire
  ! pour l'assemblage (dans un cas reel, cette formule viendrait des
  ! donnees du probleme).
  ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
  DO J = 1, NLOC
    GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      A_LOCAL(I,J) = CMPLX(1.0_WP/(1.0_WP+REAL(ABS(GI-GJ),WP)), 0.1_WP, KIND=WP)
      IF (GI .EQ. GJ) A_LOCAL(I,J) = A_LOCAL(I,J) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
    END DO
  END DO
  ! Empreinte THEORIQUE (8 octets/element, COMPLEX(WP) simple precision) --
  ! PAS mesuree via l'OS (voir CLUSTER.md pour la RSS reelle).
  CALL REPORT_MEMORY_MB(GRID, 'Matrice A (dense)', 8_8*INT(MLOC,8)*INT(NLOC,8))

  ALLOCATE(IPIV_GLOBAL(N))
  CALL MPI_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, INFO, N_SWAPS=N_SWAPS, &
                            SHOW_PROGRESS=.TRUE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_LU_PIVOT_FACTOR, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' Nombre de permutations de lignes (pivot) : ', N_SWAPS
  ! Facteurs L/U ecrasent A_LOCAL EN PLACE (meme stockage, meme taille --
  ! contrairement au BLR, il n'y a ici aucun phenomene de remplissage qui
  ! ferait grossir/retrecir l'empreinte entre la matrice et ses facteurs).
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU (dense)', 8_8*INT(MLOC,8)*INT(NLOC,8))

  ! --- second membre distribue PAR LIGNES (comme A) ET PAR COLONNES sur
  !     NRHS -- meme principe : formule deterministe, sans communication --
  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  DO J = 1, NRHS_LOC
    GJ = RSTART_MINE + J - 1
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      B_LOCAL(I,J) = CMPLX(REAL(GI,WP)*0.01_WP, -0.5_WP*REAL(GJ,WP), KIND=WP)
    END DO
  END DO
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B', 8_8*INT(MLOC,8)*INT(NRHS_LOC,8))

  CALL MPI_LU_SOLVE_DISTRIBUTED_2D(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                                    NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE.)

  ! --- rassemblement de la solution complete sur le processus 0 pour
  !     affichage (dans un cas reel, on garderait X distribue si la suite
  !     du calcul reste parallele) ---
  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
    DO I = 1, MIN(5, N)
      WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
    END DO
  END IF

  DEALLOCATE(A_LOCAL, IPIV_GLOBAL, B_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

END PROGRAM EXAMPLE_MPI_DENSE
