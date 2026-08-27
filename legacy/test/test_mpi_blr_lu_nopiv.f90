!===============================================================================
! Test de MPI_BLR_LU_NOPIV_FACTOR (case manquante de la grille dense/BLR x
! OpenMP/MPI-sans-pivot/MPI-avec-pivot). Matrice diagonale dominante (SANS
! pivot, comme 6c) ET compressible (meme recette rang-global-6 que
! test_mpi_blr_lu_pivot.f90/6e, formule deterministe -- chaque processus
! calcule sa part locale sans communication).
!
! VALIDATION : meme methode que 6c -- rassembler (GATHER_2D), reconstruire
! L*U, comparer a une reference recalculee independamment.
!
! Deux cas : TOL demande (compression active, residu de l'ordre de TOL) et
! TOL=0 (regression, doit redonner 6c pur -- precision machine).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_blr_lu_nopiv PR PC N NB TOL
!===============================================================================
PROGRAM TEST_MPI_BLR_LU_NOPIV
  USE MPI_BLR_LU_NOPIV_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, IERR, MYID0
  CHARACTER(LEN=32) :: ARG
  REAL(WP) :: TOL
  LOGICAL :: ALL_OK

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 5) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_lu_nopiv PR PC N NB TOL'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) TOL

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0 .EQ. 0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ALL_OK = .TRUE.
  CALL RUN_CASE('1) TOL demande (compression active)', TOL, &
                 .TRUE., MAX(TOL*100.0_WP, 1.0E-4_WP), ALL_OK)
  CALL RUN_CASE('2) Regression : TOL=0 (== etape 6c)', 0.0_WP, &
                 .FALSE., 1.0E-4_WP, ALL_OK)

  IF (MYID0 .EQ. 0) THEN
    IF (ALL_OK) THEN
      WRITE(6,*) '=== TOUS LES TESTS MPI_BLR_LU_NOPIV : OK ==='
    ELSE
      WRITE(6,*) '=== AU MOINS UN TEST MPI_BLR_LU_NOPIV A ECHOUE ==='
    END IF
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  SUBROUTINE RUN_CASE(LABEL, CASE_TOL, EXPECT_COMPRESSION, ERR_BOUND, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    REAL(WP), INTENT(IN)         :: CASE_TOL, ERR_BOUND
    LOGICAL, INTENT(IN)          :: EXPECT_COMPRESSION
    LOGICAL, INTENT(INOUT)       :: OK

    INTEGER, PARAMETER :: RGLOBAL = 6
    COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_FACT(:,:), A_REF(:,:), LMAT(:,:), UMAT(:,:), LU(:,:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: MLOC, NLOC, I, J, GI, GJ, R, INFO, N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8) :: RANK_SUM
    COMPLEX(WP) :: AVAL
    REAL(WP) :: ERR, NORMA, COMPRESS_FRAC

    MLOC = NUMROC(N, NB, GRID%MYROW, GRID%PR)
    NLOC = NUMROC(N, NB, GRID%MYCOL, GRID%PC)

    ALLOCATE(A_LOCAL(MAX(MLOC,1),MAX(NLOC,1)))
    DO J = 1, NLOC
      GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
      DO I = 1, MLOC
        GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
        END DO
        IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        A_LOCAL(I,J) = AVAL
      END DO
    END DO

    CALL MPI_BLR_LU_NOPIV_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, CASE_TOL, INFO, &
                                  N_COMPRESSED, N_TOTAL_PANELS, RANK_SUM)

    ALLOCATE(A_FACT(N,N))
    CALL GATHER_2D(GRID, N, NB, A_LOCAL, MAX(MLOC,1), MAX(NLOC,1), 0, A_FACT)

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'
      WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,ES10.3)') '  Grille ', PR, 'x', PC, &
           ' , N=', N, ' , NB=', NB, ' , TOL=', CASE_TOL
      IF (INFO .NE. 0) WRITE(6,*) '  ATTENTION : INFO=', INFO

      ALLOCATE(A_REF(N,N))
      DO GJ = 1, N
        DO GI = 1, N
          AVAL = (0.0_WP,0.0_WP)
          DO R = 1, RGLOBAL
            AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
          END DO
          IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
          A_REF(GI,GJ) = AVAL
        END DO
      END DO

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

      COMPRESS_FRAC = 0.0_WP
      IF (N_TOTAL_PANELS .GT. 0) COMPRESS_FRAC = 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(N_TOTAL_PANELS,WP)

      WRITE(6,'(A,ES12.4)') '  Erreur relative ||L*U - A||/||A|| = ', ERR
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Panneaux compresses : ', N_COMPRESSED, ' / ', &
           N_TOTAL_PANELS, ' (', COMPRESS_FRAC, ' %)'

      IF (ERR.GT.ERR_BOUND .OR. INFO.NE.0) THEN
        WRITE(6,*) '  -> ECHEC (residu trop grand)'
        OK = .FALSE.
      ELSE IF (EXPECT_COMPRESSION .AND. N_COMPRESSED.EQ.0) THEN
        WRITE(6,*) '  -> ECHEC (compression attendue mais aucun panneau compresse)'
        OK = .FALSE.
      ELSE IF (.NOT.EXPECT_COMPRESSION .AND. N_COMPRESSED.GT.0) THEN
        WRITE(6,*) '  -> ECHEC (compression inattendue)'
        OK = .FALSE.
      ELSE
        WRITE(6,*) '  -> OK'
      END IF
      DEALLOCATE(A_REF, LMAT, UMAT, LU)
    END IF

    DEALLOCATE(A_LOCAL, A_FACT)
  END SUBROUTINE RUN_CASE

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

END PROGRAM TEST_MPI_BLR_LU_NOPIV
