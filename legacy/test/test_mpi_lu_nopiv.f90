!===============================================================================
! Test de MPI_LU_NOPIV_FACTOR (etape 6c) : matrice A a diagonale dominante
! (indispensable, aucun pivot -- cf. commentaire du module), construite
! localement par une formule deterministe (meme principe qu'en 6a/6b :
! chaque processus ne calcule que SA part). Apres factorisation
! distribuee, on RASSEMBLE le resultat sur le rang 0 (GATHER_2D, deja
! valide en 6a), on y reconstruit L (triangulaire inferieure, diagonale
! unitaire implicite) et U (triangulaire superieure), on calcule L*U, et
! on compare au residu ||L*U - A_originale|| / ||A_originale|| avec
! A_originale RECALCULEE INDEPENDAMMENT (memes formules, en sequentiel
! complet) -- meme methodologie de validation que partout ailleurs dans
! ce projet.
!
! Usage : mpirun -np <PR*PC> ./test_mpi_lu_nopiv PR PC N NB
!===============================================================================
PROGRAM TEST_MPI_LU_NOPIV
  USE MPI_LU_NOPIV_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, IERR, MLOC, NLOC, I, J, GI, GJ, MYID0, INFO
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:)
  COMPLEX(WP), ALLOCATABLE :: A_FACT(:,:), A_REF(:,:), LMAT(:,:), UMAT(:,:), LU(:,:)
  COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
  REAL(WP) :: ERR, NORMA, DF

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 4) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_lu_nopiv PR PC N NB'
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

  ! --- A : diagonale fortement dominante (obligatoire, factorisation SANS
  !     pivot) ; partie hors-diagonale decroissante en |i-j| (reelle ET
  !     imaginaire) pour que la dominance diagonale reste vraie a tout N
  !     teste (cf. commentaire dans le module : une partie imaginaire non
  !     decroissante ferait grossir la somme de ligne comme O(N)) ---
  ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
  DO J = 1, NLOC
    GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      IF (GI .EQ. GJ) THEN
        A_LOCAL(I,J) = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
      ELSE
        DF = 1.0_WP / (1.0_WP + ABS(REAL(GI-GJ,WP)))
        A_LOCAL(I,J) = CMPLX(DF, 0.3_WP*DF, KIND=WP)
      END IF
    END DO
  END DO

  CALL MPI_LU_NOPIV_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, INFO)

  IF (MYID0 .EQ. 0) THEN
    IF (INFO .NE. 0) WRITE(6,*) 'ATTENTION : INFO=', INFO, ' (pivot proche de zero)'
  END IF

  ALLOCATE(A_FACT(N,N))
  CALL GATHER_2D(GRID, N, NB, A_LOCAL, MAX(MLOC,1), MAX(NLOC,1), 0, A_FACT)

  IF (MYID0 .EQ. 0) THEN
    ! --- Reference independante : A originale recalculee integralement ---
    ALLOCATE(A_REF(N,N))
    DO GJ = 1, N
      DO GI = 1, N
        IF (GI .EQ. GJ) THEN
          A_REF(GI,GJ) = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
        ELSE
          DF = 1.0_WP / (1.0_WP + ABS(REAL(GI-GJ,WP)))
          A_REF(GI,GJ) = CMPLX(DF, 0.3_WP*DF, KIND=WP)
        END IF
      END DO
    END DO

    ! --- Reconstruction de L (unitaire) et U depuis A_FACT (empaquetage
    !     style GETRF), puis produit L*U ---
    ALLOCATE(LMAT(N,N), UMAT(N,N), LU(N,N))
    LMAT = CZERO
    UMAT = CZERO
    DO J = 1, N
      LMAT(J,J) = CONE
      DO I = J+1, N
        LMAT(I,J) = A_FACT(I,J)
      END DO
      DO I = 1, J
        UMAT(I,J) = A_FACT(I,J)
      END DO
    END DO
    CALL CGEMM('No transpose', 'No transpose', N, N, N, CONE, LMAT, N, &
                UMAT, N, CZERO, LU, N)

    NORMA = SQRT(SUM(ABS(A_REF)**2))
    ERR = SQRT(SUM(ABS(LU-A_REF)**2)) / NORMA

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||L*U - A||/||A|| = ', ERR
    IF (ERR .LT. 1.0E-4_WP .AND. INFO.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_LU_NOPIV : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_LU_NOPIV : ECHEC ==='
    END IF
    DEALLOCATE(A_REF, LMAT, UMAT, LU)
  END IF

  DEALLOCATE(A_LOCAL, A_FACT)
  CALL MPI_FINALIZE(IERR)

END PROGRAM TEST_MPI_LU_NOPIV
