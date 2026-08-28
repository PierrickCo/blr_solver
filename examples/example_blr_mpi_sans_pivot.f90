!===============================================================================
! EXEMPLE CANONIQUE -- reference BLR+MPI, SANS pivot (mpi_blr_matrix_mod.f90,
! MPI_BLR_MATRIX_LU_FACTOR_NOPIV) : la piece la plus aboutie du projet --
! le stockage LOCAL de chaque processus reste compresse de bout en bout
! (jamais dense en dehors des blocs diagonaux), panneaux diffuses deja
! compresses. Reserve aux matrices a DIAGONALE DOMINANTE (pas de pivot ici
! -- voir example_blr_mpi_pivot.f90 pour la variante avec pivot confine).
!
! But de ce fichier : montrer l'enchainement minimal (assembler,
! factoriser, resoudre, rassembler) sans toute l'instrumentation de
! mesure/validation croisee des drivers de test (voir
! test/test_mpi_blr_matrix.f90 pour la version qui verifie aussi la
! correction contre une reference CGESV independante).
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_sans_pivot PR PC [N] [NB] [NRHS] [TOL]
!         (defauts : N=4000 NB=200 NRHS=4 TOL=1e-4)
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_SANS_PIVOT
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./example_blr_mpi_sans_pivot PR PC [N] [NB] [NRHS] [TOL]'
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
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL
  END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)') &
       ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL

  ! --- 1) assemblage BLR distribue : GET_BLOCK est appele par chaque
  !     processus UNIQUEMENT pour les blocs qu'il possede (distribution
  !     bloc-cyclique 2D sur la grille) -- jamais de matrice N x N dense,
  !     ni locale ni globale, meme transitoirement. ---
  ! SHOW_PROGRESS=.TRUE. : barre affichee par le processus 0 (approximative
  ! ici, l'assemblage n'est pas synchronise entre processus).
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE.)
  ! N_COMPRESSED/N_OFFDIAG sont deja des sommes GLOBALES (MPI_ALLREDUCE
  ! interne a MPI_ASSEMBLE_BLR_MATRIX) -- identiques sur tous les
  ! processus, un seul affichage (rang 0) suffit.
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Taux de compression (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF
  ! Empreinte THEORIQUE (blocs BLR_BLOCK_T reellement stockes -- dense pour
  ! les blocs diagonaux/repli, UFAC*VFAC pour les blocs basse-rang), PAS
  ! mesuree via l'OS.
  CALL REPORT_MEMORY_MB(GRID, 'Matrice A_BLR (compressee)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))

  ! --- 2) factorisation LU, SANS pivot -----------------------------------
  ! SHOW_PROGRESS=.TRUE. : ici K/NBLK est un pourcentage EXACT (boucle
  ! globalement synchronisee par les diffusions a chaque etape K).
  CALL MPI_BLR_MATRIX_LU_FACTOR_NOPIV(GRID, A_BLR, TOL, INFO, &
                                       N_COMPRESSED, N_OFFDIAG, RANK_SUM, T_FACTOR, &
                                       SHOW_PROGRESS=.TRUE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_NOPIV, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Taux de compression (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF
  ! A_BLR est factorisee EN PLACE (les facteurs L/U ecrasent DIAGBLK/OFFBLK) --
  ! cette valeur peut differer de celle de l'assemblage : le remplissage
  ! (fill-in) de la mise a jour de Schur peut faire grossir le rang de
  ! certains blocs (voire les faire tomber en repli dense).
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU (compresses)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))

  ! --- 3) second membre B, distribue par blocs-lignes (comme A_BLR, via
  !     NBLKR_LOC/BSTART/BEND) ET par colonnes sur NRHS -- formule
  !     deterministe sur l'indice global, sans communication. Dans un cas
  !     reel, B viendrait des donnees du probleme. ---
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_B(A_BLR, GRID, NRHS_LOC, RSTART_MINE, B_LOCAL, MLOC)
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B', 8_8*INT(MLOC,8)*INT(NRHS_LOC,8))

  ! --- 4) resolution 2D --------------------------------------------------
  ! SHOW_PROGRESS=.TRUE. : pourcentage exact, meme raison que la
  ! factorisation.
  CALL MPI_BLR_MATRIX_LU_SOLVE_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE.)

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

  DEALLOCATE(B_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Matrice "diagonale + basse-rang" a diagonale dominante -- SANS PIVOT,
  ! chaque processus calcule directement les blocs qu'il possede.
  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER, PARAMETER :: RGLOBAL = 6
    INTEGER :: II, JJ, RR
    COMPLEX(WP) :: AVAL
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        AVAL = (0.0_WP,0.0_WP)
        DO RR = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(II,RR)*LOWRANK_V(RR,JJ)
        END DO
        IF (II .EQ. JJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        BLK(II-ISTART+1, JJ-JSTART+1) = AVAL
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  SUBROUTINE FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: BL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, GRR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          GRR = RSTART + RR - 1
          BL(RS + (GII-AB%BSTART(II)), RR) = CMPLX(REAL(GII,WP)*0.01_WP, -0.3_WP*REAL(GRR,WP), KIND=WP)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_B

END PROGRAM EXAMPLE_BLR_MPI_SANS_PIVOT
