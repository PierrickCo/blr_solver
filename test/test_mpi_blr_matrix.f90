!===============================================================================
! Test du solveur LU MPI, BLR compresse des l'assemblage, SANS pivot :
! matrice "diagonale + basse-rang" deterministe (RGLOBAL=6, meme recette que
! partout ailleurs dans ce projet), stockage local en BLR_BLOCK_T (jamais
! dense en dehors des blocs diagonaux), resolution 2D complete.
!
! VALIDATION : rassemble le resultat (reutilise GATHER_2D_RHS, deja
! valide pour mpi_lu_solve_2d_mod.f90 -- meme convention de decoupage NRHS
! sur PC) puis compare a une reference INDEPENDANTE (CGESV, matrice/second
! membre recalcules a partir de zero).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_blr_matrix PR PC N NB NRHS TOL
!===============================================================================
PROGRAM TEST_MPI_BLR_MATRIX
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), XEX_LOCAL(:,:), X_GLOBAL(:,:)
  COMPLEX(WP), ALLOCATABLE :: B_REF(:,:)
  INTEGER, PARAMETER :: RGLOBAL = 6
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I, GI, GJ, R
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  DOUBLE PRECISION :: RSS_FACT_MB, RSS_SOLVE_MB, MEM_ANALYTIC_MB, TMP_MB
  DOUBLE PRECISION :: T_DIAG, T_BCAST_DKK, T_TRSM, T_BCAST_PANEL
  DOUBLE PRECISION :: T_ALLOC, T_MPICALL, T_SCHUR_PRODUCT, T_RECOMPRESS, T_SCHUR_DIAG
  REAL(WP) :: ERR

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 6) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_matrix PR PC N NB NRHS TOL'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL REPORT_THREAD_AFFINITY(MPI_COMM_WORLD)

  ! === 1) assemblage BLR distribue ===
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM)

  ! === 2) factorisation LU, sans pivot ===
  CALL MPI_BLR_MATRIX_LU_FACTOR_NOPIV(GRID, A_BLR, TOL, INFO, &
                                       N_COMPRESSED, N_OFFDIAG, RANK_SUM, T_FACTOR, &
                                       T_DIAG, T_BCAST_DKK, T_TRSM, T_BCAST_PANEL, &
                                       T_ALLOC, T_MPICALL, T_SCHUR_PRODUCT, T_RECOMPRESS, T_SCHUR_DIAG)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_NOPIV, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  ! RSS PAR PROCESSUS apres factorisation (max sur tous, meme convention
  ! que bench_mpi.f90) + memoire ANALYTIQUE (BLR_BLOCK_T reellement stocke).
  TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, RSS_FACT_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
  TMP_MB = REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR),8)/1024.0D0/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, MEM_ANALYTIC_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)

  ! === 3) second membre B = A*Xexact, distribue lignes ET colonnes ===
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(XEX_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)), B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_XEX(A_BLR, GRID, NRHS_LOC, RSTART_MINE, XEX_LOCAL, MLOC)

  ! B = A_BLR * Xexact -- utilise la matrice ORIGINALE, mais A_BLR est deja
  ! FACTORISEE a ce point : on reconstruit B directement via la formule
  ! (independante) plus bas pour la reference, et via une seconde matrice
  ! non factorisee ICI pour la coherence numerique locale -- plus simple :
  ! reassembler une copie fraiche juste pour construire B.
  BLOCK
    TYPE(MPI_BLR_MATRIX_T) :: A_BLR_FOR_B
    DOUBLE PRECISION :: T_ASM_DUMMY
    INTEGER :: NC_DUMMY, NO_DUMMY
    INTEGER(8) :: RS_DUMMY
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR_FOR_B, T_ASM_DUMMY, &
                                  NC_DUMMY, NO_DUMMY, RS_DUMMY)
    CALL APPLY_MPI_BLR_MATRIX(GRID, A_BLR_FOR_B, NRHS, NRHS_LOC, XEX_LOCAL, B_LOCAL, MLOC)
  END BLOCK

  ! === 4) resolution 2D ===
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_MATRIX_LU_SOLVE_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
  T1 = MPI_WTIME()
  T_SOLVE = T1 - T0
  TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, RSS_SOLVE_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)

  ! === validation : B a ete construit comme A_BLR*Xexact (APPLY_MPI_BLR_
  !     MATRIX ci-dessus, meme convention que bench_blr_pipeline.f90/"A")
  !     -- resoudre A*X=B doit donc redonner Xexact directement (PAS de
  !     CGESV separe necessaire ici : une premiere version comparait a
  !     tort au resultat de CGESV(A, Xexact_formule), ce qui n'a rien a
  !     voir avec ce qu'on cherche a verifier).
  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

  ERR = -1.0_WP
  IF (MYID0 .EQ. 0) THEN
    ALLOCATE(B_REF(N,NRHS))
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = XEX_FORMULA(GI,R)
      END DO
    END DO
    ERR = SQRT(SUM(ABS(X_GLOBAL-B_REF)**2)) / SQRT(SUM(ABS(B_REF)**2))

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)') ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, &
         ' , NRHS=', NRHS, ' , TOL=', TOL
    WRITE(6,'(A,F8.3,A,F8.3,A,F8.3,A)') '   Temps assemblage = ', T_ASM, ' s , factorisation = ', T_FACTOR, &
         ' s , resolution = ', T_SOLVE, ' s'
    WRITE(6,'(A)') '   --- Decomposition de la factorisation (max/processus) ---'
    WRITE(6,'(A,F8.3,A,F8.3,A,F8.3)') '     facteur diagonal = ', T_DIAG, ' s , bcast DKK = ', T_BCAST_DKK, &
         ' s , TRSM = ', T_TRSM
    WRITE(6,'(A,F8.3,A,F8.3,A,F8.3,A)') '     bcast panneaux L/U = ', T_BCAST_PANEL, ' s  (dont alloc/dealloc = ', &
         T_ALLOC, ' s , MPI_BCAST pur = ', T_MPICALL, ' s)'
    WRITE(6,'(A,F8.3,A,F8.3,A,F8.3)') '     Schur diag (dense) = ', T_SCHUR_DIAG, ' s , Schur hors-diag sans compression = ', &
         T_SCHUR_PRODUCT, ' s , recompression = ', T_RECOMPRESS
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Compression : ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', &
         100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F8.2)') '   Rang moyen = ', REAL(RANK_SUM,WP)/REAL(MAX(N_OFFDIAG,1),WP)
    WRITE(6,'(A,F10.1,A,F10.1,A)') '   Memoire ANALYTIQUE (BLR_BLOCK_T, max/processus) = ', MEM_ANALYTIC_MB, &
         ' Mo  ;  RSS apres factorisation = ', RSS_FACT_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)') '   RSS apres resolution (max/processus) = ', RSS_SOLVE_MB, ' Mo'
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_distribue - X_exact|| = ', ERR
    IF (ERR.LT.1.0E-4_WP .AND. INFO.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX : OK ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX : ECHEC ==='
    END IF
    DEALLOCATE(B_REF)
  END IF

  DEALLOCATE(B_LOCAL, XEX_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ, NN)
    INTEGER, INTENT(IN) :: GI, GJ, NN
    INTEGER :: RR
    LOWRANK_A = (0.0_WP,0.0_WP)
    DO RR = 1, RGLOBAL
      LOWRANK_A = LOWRANK_A + LOWRANK_U(GI,RR)*LOWRANK_V(RR,GJ)
    END DO
    IF (GI .EQ. GJ) LOWRANK_A = LOWRANK_A + CMPLX(REAL(NN,WP), 0.0_WP, KIND=WP)
  END FUNCTION LOWRANK_A

  ! Xexact deterministe (pas aleatoire -- chaque processus calcule sa
  ! propre part SANS communication, meme discipline que tous les tests
  ! MPI de ce projet).
  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        BLK(II-ISTART+1, JJ-JSTART+1) = LOWRANK_A(II, JJ, N)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  SUBROUTINE FILL_XEX(AB, GD, NRL, RSTART, XL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: XL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, GRR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          GRR = RSTART + RR - 1
          XL(RS + (GII-AB%BSTART(II)), RR) = XEX_FORMULA(GII, GRR)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_XEX

  ! Y = A_BLR * X (matrice COMPRESSEE, non factorisee), 2D complet --
  ! meme algorithme que la mise a jour de Schur/la resolution : diffuse
  ! chaque bloc-colonne le long de sa colonne de grille, applique
  ! localement via BLR_BLOCK_APPLY/CGEMM, accumule.
  SUBROUTINE APPLY_MPI_BLR_MATRIX(GD, AB, NR, NRL, XL, YL, ML)
    TYPE(GRID_T), INTENT(IN)           :: GD
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    INTEGER, INTENT(IN)                :: NR, NRL, ML
    COMPLEX(WP), INTENT(IN)            :: XL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), INTENT(OUT)           :: YL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), PARAMETER :: CONE=(1.0_WP,0.0_WP), CZERO=(0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: XJ_LOC(:,:), XJ_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: J, lj, li, I, NJ, MI, RS, RE, C, IERR2, PROWJ, PCOLJ, MLOC_ALL

    IF (NRL .GT. 0) YL = CZERO

    ALLOCATE(RSTARTS(0:GD%PC-1), RLOCS(0:GD%PC-1), RCOUNTS(0:GD%PC-1), RDISPLS(0:GD%PC-1))
    DO C = 0, GD%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NR, GD%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    DO J = 1, AB%NBLK
      NJ = AB%BEND(J)-AB%BSTART(J)+1
      PROWJ = MOD(J-1, GD%PR); PCOLJ = MOD(J-1, GD%PC)

      ! X(J) complet (NJ x NR), reconstruit sur la ligne PROWJ puis diffuse
      ALLOCATE(XJ_FULL(NJ,NR))
      IF (GD%MYROW .EQ. PROWJ) THEN
        RS = LOCAL_ROW_OFFSET(AB, GD, J)
        ALLOCATE(XJ_LOC(NJ,MAX(NRL,1)))
        IF (NRL .GT. 0) XJ_LOC = XL(RS:RS+NJ-1,1:NRL)
        DO C = 0, GD%PC-1
          RCOUNTS(C) = NJ*RLOCS(C); RDISPLS(C) = NJ*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XJ_LOC, NJ*NRL, MPI_COMPLEX, XJ_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                             GD%ROW_COMM, IERR2)
        DEALLOCATE(XJ_LOC)
      END IF
      CALL MPI_BCAST(XJ_FULL, NJ*NR, MPI_COMPLEX, PROWJ, GD%COL_COMM, IERR2)

      ! contribution de la colonne de bloc J a TOUS les blocs-lignes locaux
      ! -- MLOC_ALL ne depend que de MYROW, calcule par TOUS (necessaire
      ! pour que MPI_SCATTERV, appel COLLECTIF sur ROW_COMM, soit appele
      ! par TOUTE la ligne et pas seulement le proprietaire MYCOL==PCOLJ --
      ! meme bug/correction que MPI_BLR_MATRIX_LU_SOLVE_2D).
      MLOC_ALL = 0
      DO li = 1, AB%NBLKR_LOC
        I = (li-1)*GD%PR + GD%MYROW + 1
        MLOC_ALL = MLOC_ALL + (AB%BEND(I)-AB%BSTART(I)+1)
      END DO
      ALLOCATE(CONTRIB(MLOC_ALL,NR))
      IF (GD%MYCOL .EQ. PCOLJ) THEN
        lj = (J-1)/GD%PC + 1
        CONTRIB = CZERO
        RE = 0
        DO li = 1, AB%NBLKR_LOC
          I = (li-1)*GD%PR + GD%MYROW + 1
          MI = AB%BEND(I)-AB%BSTART(I)+1
          IF (I .EQ. J) THEN
            CALL CGEMM('No transpose','No transpose', MI, NR, NJ, CONE, AB%DIAGBLK(I)%D, MI, &
                        XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
          ELSE
            CALL BLR_BLOCK_APPLY(AB%OFFBLK(li,lj), NR, CONE, XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
          END IF
          RE = RE + MI
        END DO
      END IF
      ALLOCATE(CONTRIB_LOC(MLOC_ALL,MAX(NRL,1)))
      DO C = 0, GD%PC-1
        RCOUNTS(C) = MLOC_ALL*RLOCS(C); RDISPLS(C) = MLOC_ALL*(RSTARTS(C)-1)
      END DO
      CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, CONTRIB_LOC, MLOC_ALL*NRL, &
                         MPI_COMPLEX, PCOLJ, GD%ROW_COMM, IERR2)
      IF (NRL .GT. 0) YL(1:MLOC_ALL,1:NRL) = YL(1:MLOC_ALL,1:NRL) + CONTRIB_LOC(1:MLOC_ALL,1:NRL)
      DEALLOCATE(CONTRIB, CONTRIB_LOC)
      DEALLOCATE(XJ_FULL)
    END DO
    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE APPLY_MPI_BLR_MATRIX

END PROGRAM TEST_MPI_BLR_MATRIX
