!===============================================================================
! Test de MPI_LU_SOLVE_DISTRIBUTED : B DISTRIBUE PAR LIGNES (comme A),
! replique seulement au sein d'une meme ligne de grille -- PAS replique
! integralement (evolution demandee par l'utilisateur : NRHS peut devenir
! tres grand). Meme matrice "anti-diagonale" que test_mpi_lu_pivot.f90
! (6d) -- force un vrai pivotage, exerce donc l'etape 0 (SWAP_ROWS_RANGE
! sur B, potentiellement inter-processus-ligne) pour de vrai.
!
! VALIDATION : rassemble B_LOCAL (GATHER_ROW_VECTOR, deja valide en 6b,
! appele NRHS fois -- une colonne a la fois) puis compare a une reference
! INDEPENDANTE (CGESV, matrice/second membre recalcules a partir de zero).
! Verifie AUSSI que la replication INTRA-LIGNE reste coherente (toutes les
! copies d'un meme MYROW identiques) -- diagnostic direct de la coherence
! du nouveau schema de diffusion ROW_COMM/COL_COMM.
!
! Usage : mpirun -np <PR*PC> ./test_mpi_lu_solve PR PC N NB NRHS
!===============================================================================
PROGRAM TEST_MPI_LU_SOLVE
  USE MPI_LU_PIVOT_MOD
  USE MPI_LU_SOLVE_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: A_LOCAL(:,:), A_REF(:,:), B_LOCAL(:,:), B_REF(:,:), X_GLOBAL(:,:)
  INTEGER, ALLOCATABLE :: IPIV_GLOBAL(:), IPIV_REF(:)
  COMPLEX(WP), PARAMETER :: CONE=(1.0_WP,0.0_WP), CZERO=(0.0_WP,0.0_WP)
  INTEGER :: MLOC, NLOC, I, J, GI, GJ, INFO, INFO_REF, R
  REAL(WP) :: ERR, ERR_REPL, NORMB, DF

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 5) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_lu_solve PR PC N NB NRHS'
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

  ! --- second membre DISTRIBUE PAR LIGNES (comme A) : chaque processus
  !     calcule SA part locale (formule deterministe -- meme donnee sur
  !     tous les MYCOL d'une meme ligne, automatiquement, sans
  !     communication, exactement comme A_LOCAL) ---
  ALLOCATE(B_LOCAL(MAX(MLOC,1),NRHS))
  DO R = 1, NRHS
    DO I = 1, MLOC
      GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
      B_LOCAL(I,R) = CMPLX(REAL(GI*R,WP)*0.01_WP, -0.5_WP*REAL(GI,WP)+REAL(R,WP), KIND=WP)
    END DO
  END DO

  CALL MPI_LU_SOLVE_DISTRIBUTED(GRID, N, NB, A_LOCAL, MLOC, NLOC, IPIV_GLOBAL, NRHS, B_LOCAL, MAX(MLOC,1))

  ! --- verifier la coherence de la replication INTRA-LIGNE : diffuser la
  !     copie du rang 0 DE CHAQUE LIGNE (root=0 dans ROW_COMM) et comparer
  !     -- plus simple et sans risque d'appariement qu'un anneau SENDRECV
  !     (une premiere version en anneau ne fonctionnait que pour PC=2,
  !     bloquait pour PC>=3 : chaque processus attendait un message d'un
  !     pair qui envoyait ailleurs -- bug de TEST, pas du module) ---
  ERR_REPL = 0.0_WP
  IF (GRID%PC .GT. 1) THEN
    BLOCK
      COMPLEX(WP), ALLOCATABLE :: B_ROW0(:,:)
      ALLOCATE(B_ROW0(MAX(MLOC,1),NRHS))
      B_ROW0 = B_LOCAL
      CALL MPI_BCAST(B_ROW0, MAX(MLOC,1)*NRHS, MPI_COMPLEX, 0, GRID%ROW_COMM, IERR)
      ERR_REPL = SQRT(SUM(ABS(B_LOCAL-B_ROW0)**2))
      DEALLOCATE(B_ROW0)
    END BLOCK
  END IF

  ! --- rassembler B_LOCAL (GATHER_ROW_VECTOR, deja valide en 6b, colonne
  !     par colonne) pour comparer a la reference ---
  ALLOCATE(X_GLOBAL(N,NRHS))
  DO R = 1, NRHS
    CALL GATHER_ROW_VECTOR(GRID, N, NB, B_LOCAL(:,R), MAX(MLOC,1), 0, X_GLOBAL(:,R))
  END DO

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
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = CMPLX(REAL(GI*R,WP)*0.01_WP, -0.5_WP*REAL(GI,WP)+REAL(R,WP), KIND=WP)
      END DO
    END DO
    CALL CGESV(N, NRHS, A_REF, N, IPIV_REF, B_REF, N, INFO_REF)

    NORMB = SQRT(SUM(ABS(B_REF)**2))
    ERR = SQRT(SUM(ABS(X_GLOBAL-B_REF)**2)) / NORMB

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_distribue - X_reference|| = ', ERR
    IF (ERR.LT.1.0E-4_WP .AND. INFO.EQ.0 .AND. INFO_REF.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_LU_SOLVE : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_LU_SOLVE : ECHEC ==='
    END IF
    DEALLOCATE(A_REF, B_REF, IPIV_REF)
  END IF

  IF (GRID%PC .GT. 1) THEN
    IF (ERR_REPL .GT. 1.0E-6_WP) THEN
      WRITE(6,'(A,I0,A,ES12.4)') ' [rang ',MYID0,'] ECHEC : divergence de la replication intra-ligne = ', ERR_REPL
    END IF
  END IF

  DEALLOCATE(A_LOCAL, IPIV_GLOBAL, B_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

END PROGRAM TEST_MPI_LU_SOLVE
