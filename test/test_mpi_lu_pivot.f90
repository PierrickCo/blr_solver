!===============================================================================
! Test de MPI_LU_PIVOT_FACTOR (etape 6d) : matrice construite pour forcer un
! VRAI pivotage -- BASE(R,C) est a diagonale dominante (diag=50, hors-diag
! decroissant en |R-C|, comme en 6c), mais on prend A(GI,GJ)=BASE(N+1-GI,GJ)
! (inversion des lignes). La diagonale NATURELLE de A (GI=GJ) tombe alors
! generiquement sur une entree FAIBLE de BASE (sauf au milieu) -- les
! entrees fortes de BASE se retrouvent sur l'ANTI-diagonale de A. Une
! factorisation SANS pivot serait donc numeriquement catastrophique ici ;
! un vrai pivotage doit "retrouver" ces entrees fortes. A reste bien
! conditionnee (c'est juste une permutation de lignes d'une matrice a
! diagonale dominante), donc la comparaison residuelle reste significative.
!
! VALIDATION (plus forte que 6a-6c) : au lieu de reconstruire L*U a la
! main, on reutilise BLOCK_LU_SOLVE (etape 1, deja valide independamment)
! sur le resultat rassemble (A_FACT + IPIV_GLOBAL) pour resoudre un systeme
! -- boucle une deuxieme fois par du code sequentiel deja teste plutot que
! par une nouvelle verification ad hoc. On compare la solution obtenue a
! une reference TOTALEMENT INDEPENDANTE (CGESV, sur une matrice/second
! membre recalcules a partir de zero).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_lu_pivot PR PC N NB
!===============================================================================
PROGRAM TEST_MPI_LU_PIVOT
  USE MPI_LU_PIVOT_MOD
  USE BLOCK_LU_MOD, ONLY: BLOCK_LU_SOLVE
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, IERR, MLOC, NLOC, I, J, GI, GJ, MYID0
  INTEGER :: INFO, NSWAPS, INFO_REF
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_FACT(:,:), A_REF(:,:)
  COMPLEX(WP), ALLOCATABLE :: B(:), B_REF(:)
  INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:), IPIV_REF(:)
  REAL(WP) :: ERR, NORMB, MINPIV, DF

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 4) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_lu_pivot PR PC N NB'
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

  ! --- A(GI,GJ) = BASE(N+1-GI,GJ) : diagonale NATURELLE de A faible,
  !     entrees fortes sur l'anti-diagonale -- force un vrai pivotage ---
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
  CALL MPI_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, &
                            INFO, NSWAPS, MINPIV)

  ALLOCATE(A_FACT(N,N))
  CALL GATHER_2D(GRID, N, NB, A_LOCAL, MAX(MLOC,1), MAX(NLOC,1), 0, A_FACT)

  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB
    WRITE(6,'(A,I0,A,ES12.4)') ' N_SWAPS=', NSWAPS, ' , MIN_PIVOT=', MINPIV
    IF (INFO .NE. 0) WRITE(6,*) 'ATTENTION : INFO=', INFO, ' (pivot proche de zero)'

    ! --- resolution via le code deja valide de l'etape 1 ---
    ALLOCATE(B(N))
    DO GI = 1, N
      B(GI) = CMPLX(REAL(GI,WP), -0.5_WP*REAL(GI,WP), KIND=WP)
    END DO
    CALL BLOCK_LU_SOLVE(N, NB, 1, A_FACT, N, IPIV_GLOBAL, B, N)

    ! --- reference totalement independante : A et B recalcules depuis
    !     zero, resolus par CGESV (appel LAPACK direct, un seul call) ---
    ALLOCATE(A_REF(N,N), B_REF(N), IPIV_REF(N))
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
    DO GI = 1, N
      B_REF(GI) = CMPLX(REAL(GI,WP), -0.5_WP*REAL(GI,WP), KIND=WP)
    END DO
    CALL CGESV(N, 1, A_REF, N, IPIV_REF, B_REF, N, INFO_REF)

    NORMB = SQRT(SUM(ABS(B_REF)**2))
    ERR = SQRT(SUM(ABS(B-B_REF)**2)) / NORMB

    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_distribue - X_reference||/||X_reference|| = ', ERR
    IF (ERR.LT.1.0E-4_WP .AND. INFO.EQ.0 .AND. INFO_REF.EQ.0 .AND. NSWAPS.GT.0) THEN
      WRITE(6,*) '=== TEST MPI_LU_PIVOT : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_LU_PIVOT : ECHEC ==='
    END IF
    DEALLOCATE(A_REF, B_REF, IPIV_REF, B)
  END IF

  DEALLOCATE(A_LOCAL, A_FACT, IPIV_GLOBAL)
  CALL MPI_FINALIZE(IERR)

END PROGRAM TEST_MPI_LU_PIVOT
