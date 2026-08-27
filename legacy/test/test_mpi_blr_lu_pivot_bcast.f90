!===============================================================================
! Test de MPI_BLR_LU_PIVOT_BCAST_FACTOR (proposition 3 : diffusion des
! panneaux DEJA COMPRESSES). Meme recette/validation que
! test_mpi_blr_lu_pivot.f90 (6e) -- meme formule deterministe "diagonale +
! correction basse-rang", meme reutilisation de BLOCK_LU_SOLVE (etape 1)
! pour resoudre, meme comparaison a CGESV independant.
!
! CE QUI EST NOUVEAU ICI : mesurer le volume de communication REEL envoye
! (BYTES_SENT) contre ce que 6e aurait envoye en diffusant dense
! (BYTES_DENSE_EQUIV) -- la promesse de cette proposition doit etre
! MESUREE, pas seulement supposee (meme discipline que l'etape 5 et la
! proposition 2).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_blr_lu_pivot_bcast PR PC N NB TOL
!===============================================================================
PROGRAM TEST_MPI_BLR_LU_PIVOT_BCAST
  USE MPI_BLR_LU_PIVOT_BCAST_MOD
  USE BLOCK_LU_MOD, ONLY: BLOCK_LU_SOLVE
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
    IF (MYID0 .EQ. 0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_lu_pivot_bcast PR PC N NB TOL'
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

  CALL RUN_CASE('1) TOL demande (compression + diffusion compressee)', TOL, &
                 .TRUE., MAX(TOL*100.0_WP, 1.0E-4_WP), ALL_OK)
  CALL RUN_CASE('2) Regression : TOL=0 (aucune compression)', 0.0_WP, &
                 .FALSE., 1.0E-4_WP, ALL_OK)

  IF (MYID0 .EQ. 0) THEN
    IF (ALL_OK) THEN
      WRITE(6,*) '=== TOUS LES TESTS MPI_BLR_LU_PIVOT_BCAST : OK ==='
    ELSE
      WRITE(6,*) '=== AU MOINS UN TEST MPI_BLR_LU_PIVOT_BCAST A ECHOUE ==='
    END IF
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  SUBROUTINE RUN_CASE(LABEL, CASE_TOL, EXPECT_SMALLER, ERR_BOUND, OK)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    REAL(WP), INTENT(IN)         :: CASE_TOL, ERR_BOUND
    LOGICAL, INTENT(IN)          :: EXPECT_SMALLER
    LOGICAL, INTENT(INOUT)       :: OK

    INTEGER, PARAMETER :: RGLOBAL = 6
    COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_FACT(:,:), A_REF(:,:)
    COMPLEX(WP), ALLOCATABLE :: B(:), B_REF(:)
    INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:), IPIV_REF(:)
    INTEGER :: MLOC, NLOC, I, J, GI, GJ, R, INFO, INFO_REF
    INTEGER :: N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8) :: RANK_SUM, BYTES_SENT, BYTES_DENSE_EQUIV
    COMPLEX(WP) :: AVAL
    REAL(WP) :: ERR, NORMB

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

    ALLOCATE(IPIV_GLOBAL(N))
    CALL MPI_BLR_LU_PIVOT_BCAST_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, CASE_TOL, &
                                        IPIV_GLOBAL, INFO, N_COMPRESSED, N_TOTAL_PANELS, &
                                        RANK_SUM, BYTES_SENT, BYTES_DENSE_EQUIV)

    ALLOCATE(A_FACT(N,N))
    CALL GATHER_2D(GRID, N, NB, A_LOCAL, MAX(MLOC,1), MAX(NLOC,1), 0, A_FACT)

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'
      WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,ES10.3)') '  Grille ', PR, 'x', PC, &
           ' , N=', N, ' , NB=', NB, ' , TOL=', CASE_TOL
      IF (INFO .NE. 0) WRITE(6,*) '  ATTENTION : INFO=', INFO

      ALLOCATE(B(N))
      DO GI = 1, N
        B(GI) = CMPLX(REAL(GI,WP), -0.5_WP*REAL(GI,WP), KIND=WP)
      END DO
      CALL BLOCK_LU_SOLVE(N, NB, 1, A_FACT, N, IPIV_GLOBAL, B, N)

      ALLOCATE(A_REF(N,N), B_REF(N), IPIV_REF(N))
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
      DO GI = 1, N
        B_REF(GI) = CMPLX(REAL(GI,WP), -0.5_WP*REAL(GI,WP), KIND=WP)
      END DO
      CALL CGESV(N, 1, A_REF, N, IPIV_REF, B_REF, N, INFO_REF)

      NORMB = SQRT(SUM(ABS(B_REF)**2))
      ERR = SQRT(SUM(ABS(B-B_REF)**2)) / NORMB

      WRITE(6,'(A,ES12.4)') '  Erreur relative ||X_distribue - X_reference|| = ', ERR
      WRITE(6,'(A,I0,A,I0)') '  Panneaux compresses : ', N_COMPRESSED, ' / ', N_TOTAL_PANELS
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Octets envoyes (compresse) = ', BYTES_SENT, &
           ' , equivalent dense (6e) = ', BYTES_DENSE_EQUIV, ' (', &
           100.0_WP*REAL(BYTES_SENT,WP)/REAL(MAX(BYTES_DENSE_EQUIV,1_8),WP), ' % du dense)'

      IF (ERR .GT. ERR_BOUND .OR. INFO.NE.0 .OR. INFO_REF.NE.0) THEN
        WRITE(6,*) '  -> ECHEC (residu trop grand)'
        OK = .FALSE.
      ELSE IF (EXPECT_SMALLER .AND. BYTES_SENT .GE. BYTES_DENSE_EQUIV) THEN
        WRITE(6,*) '  -> ECHEC (reduction du volume de communication attendue mais absente)'
        OK = .FALSE.
      ELSE
        WRITE(6,*) '  -> OK'
      END IF
      DEALLOCATE(A_REF, B_REF, IPIV_REF, B)
    END IF

    DEALLOCATE(A_LOCAL, A_FACT, IPIV_GLOBAL)
  END SUBROUTINE RUN_CASE

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

END PROGRAM TEST_MPI_BLR_LU_PIVOT_BCAST
