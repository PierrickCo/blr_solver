!===============================================================================
! Test de MPI_LU_SOLVE_DISTRIBUTED_2D : B distribue PAR LIGNES (comme A,
! NUMROC/PR) ET PAR COLONNES sur les seconds membres (RHS_NUMROC/PC) --
! AUCUNE replication residuelle. Meme matrice "anti-diagonale" que
! test_mpi_lu_solve.f90 (force un vrai pivotage, exerce SWAP_ROWS_RANGE).
!
! VALIDATION : rassemblement MANUEL (GATHER_2D existe pour un decoupage
! bloc-cyclique des DEUX cotes ; ici le cote "second membre" utilise
! RHS_LOCAL_RANGE, pas INDXG2L/NB -- rassemblement ecrit ici en reprenant
! le meme schema point-a-point que GATHER_2D, colonne mappee via
! RHS_LOCAL_RANGE au lieu de INDXL2G) puis comparaison a une reference
! INDEPENDANTE (CGESV).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_lu_solve_2d PR PC N NB NRHS
!===============================================================================
PROGRAM TEST_MPI_LU_SOLVE_2D
  USE MPI_LU_PIVOT_MOD
  USE MPI_LU_SOLVE_2D_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_REF(:,:), B_LOCAL(:,:), B_REF(:,:), X_GLOBAL(:,:)
  INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:), IPIV_REF(:)
  INTEGER :: MLOC, NLOC, NRHS_LOC, RSTART_MINE, I, J, GI, GJ, INFO, INFO_REF, R
  REAL(WP) :: ERR, NORMB, DF

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 5) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_lu_solve_2d PR PC N NB NRHS'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
  NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)
  NRHS_LOC = RHS_NUMROC(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ! --- meme matrice "anti-diagonale" que 6d : force un vrai pivotage ---
  ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
  DO J = 1, NLOC
    GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      IF (GJ .EQ. N+1-GI) THEN
        A_LOCAL(I,J) = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
      ELSE
        DF = 1.0_WP / (1.0_WP + ABS(REAL(N+1-GI-GJ,WP)))
        A_LOCAL(I,J) = CMPLX(DF, 0.3_WP*DF, KIND=WP)
      END IF
    END DO
  END DO

  ALLOCATE(IPIV_GLOBAL(N))
  CALL MPI_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, INFO)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_LU_PIVOT_FACTOR, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ! --- second membre distribue PAR LIGNES (comme A) ET PAR COLONNES sur
  !     NRHS (RSTART_MINE:RSTART_MINE+NRHS_LOC-1) -- chaque processus
  !     calcule SA part locale, formule deterministe sur l'indice RHS
  !     GLOBAL, sans communication ---
  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  DO R = 1, NRHS_LOC
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      GJ = RSTART_MINE + R - 1
      B_LOCAL(I,R) = CMPLX(REAL(GI*GJ,WP)*0.01_WP, -0.5_WP*REAL(GI,WP)+REAL(GJ,WP), KIND=WP)
    END DO
  END DO

  CALL MPI_LU_SOLVE_DISTRIBUTED_2D(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                                    NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))

  ! --- rassemblement manuel (meme schema que GATHER_2D, colonne mappee
  !     via RHS_LOCAL_RANGE au lieu de INDXL2G/NB) ---
  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)
  ! NOTE : GATHER_2D_RHS vient de mpi_lu_solve_2d_mod.f90 (USE'd plus haut)

  IF (MYID0 .EQ. 0) THEN
    ALLOCATE(A_REF(N,N), B_REF(N,NRHS), IPIV_REF(N))
    DO GJ = 1, N
      DO GI = 1, N
        IF (GJ .EQ. N+1-GI) THEN
          A_REF(GI,GJ) = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
        ELSE
          DF = 1.0_WP / (1.0_WP + ABS(REAL(N+1-GI-GJ,WP)))
          A_REF(GI,GJ) = CMPLX(DF, 0.3_WP*DF, KIND=WP)
        END IF
      END DO
    END DO
    DO GJ = 1, NRHS
      DO GI = 1, N
        B_REF(GI,GJ) = CMPLX(REAL(GI*GJ,WP)*0.01_WP, -0.5_WP*REAL(GI,WP)+REAL(GJ,WP), KIND=WP)
      END DO
    END DO
    CALL CGESV(N, NRHS, A_REF, N, IPIV_REF, B_REF, N, INFO_REF)

    NORMB = SQRT(SUM(ABS(B_REF)**2))
    ERR = SQRT(SUM(ABS(X_GLOBAL-B_REF)**2)) / NORMB

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_distribue - X_reference|| = ', ERR
    IF (ERR.LT.1.0E-4_WP .AND. INFO.EQ.0 .AND. INFO_REF.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_LU_SOLVE_2D : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_LU_SOLVE_2D : ECHEC ==='
    END IF
    DEALLOCATE(A_REF, B_REF, IPIV_REF)
  END IF

  DEALLOCATE(A_LOCAL, IPIV_GLOBAL, B_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

END PROGRAM TEST_MPI_LU_SOLVE_2D
