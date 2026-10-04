!===============================================================================
! Test du solveur LU MPI, BLR compresse des l'assemblage, AVEC PIVOT PARTIEL
! CONFINE AU BLOC DIAGONAL (MPI_BLR_MATRIX_LU_FACTOR_PIVOT/_SOLVE_PIVOT_2D --
! etape 2 du plan valide avec l'utilisateur, apres le sans-pivot deja
! valide dans test_mpi_blr_matrix.f90).
!
! MATRICE : contrairement au test sans-pivot (diagonale dominante, recette
! LOWRANK_A -- ne force JAMAIS de pivot), on reutilise ICI la matrice
! "anti-diagonale forte / diagonale naturelle faible" DEJA utilisee par
! test_mpi_lu_pivot.f90 (6d, MPI dense avec pivot) : A(GI,GJ) = 50 si
! GJ=N+1-GI, sinon un noyau de Cauchy amorti 1/(1+|N+1-GI-GJ|) -- force un
! VRAI pivotage (diagonale naturelle petite face a l'anti-diagonale), et
! sert de test de non-regression pour le chemin DENSE de BLR_BLOCK_ROW_SWAP
! (cette matrice est chaotique, peu compressible -- la plupart des blocs
! hors-diagonale devraient tomber en repli dense, exercant precisement la
! branche non testee par le sans-pivot).
!
! VALIDATION : B = A_BLR*Xexact (formule deterministe), resoudre doit
! redonner Xexact -- meme methodologie que test_mpi_blr_matrix.f90. N_SWAPS
! (rapporte) doit etre > 0 pour confirmer que le pivot a reellement ete
! exerce (sinon le test ne prouverait rien sur la logique de swap).
!
! Usage : mpirun -np <PR*PC> ./test_mpi_blr_matrix_pivot PR PC N NB NRHS TOL
!===============================================================================
PROGRAM TEST_MPI_BLR_MATRIX_PIVOT
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
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I, GI, GJ, R
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  DOUBLE PRECISION :: RSS_FACT_MB, RSS_SOLVE_MB, MEM_ANALYTIC_MB, TMP_MB
  REAL(WP) :: ERR

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 6) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_matrix_pivot PR PC N NB NRHS TOL'
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

  ! === 2) factorisation LU, AVEC pivot confine au bloc diagonal ===
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                       N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, RSS_FACT_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
  TMP_MB = REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR),8)/1024.0D0/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, MEM_ANALYTIC_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)

  ! === VERIFICATION INDEPENDANTE : verifie directement P*A_orig == L*U (grille 1x1
  ! seulement, N raisonnable) -- isole si le bug est dans la FACTORISATION
  ! (propagation du pivot aux blocs hors-diagonale) ou dans la RESOLUTION. ===
  IF (PR.EQ.1 .AND. PC.EQ.1 .AND. N.LE.2000 .AND. MYID0.EQ.0) THEN
    CALL DEBUG_CHECK_FACTOR(A_BLR, N)
  END IF

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

  BLOCK
    TYPE(MPI_BLR_MATRIX_T) :: A_BLR_FOR_B
    DOUBLE PRECISION :: T_ASM_DUMMY
    INTEGER :: NC_DUMMY, NO_DUMMY
    INTEGER(8) :: RS_DUMMY
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR_FOR_B, T_ASM_DUMMY, &
                                  NC_DUMMY, NO_DUMMY, RS_DUMMY)
    CALL APPLY_MPI_BLR_MATRIX(GRID, A_BLR_FOR_B, NRHS, NRHS_LOC, XEX_LOCAL, B_LOCAL, MLOC)
  END BLOCK

  ! === VERIFICATION INDEPENDANTE : reference CGESV independante sur A dense vs B
  ! REEL (celui que la resolution va utiliser), AVANT que SOLVE_PIVOT_2D
  ! n'ecrase B_LOCAL en place -- isole si le bug est dans la resolution
  ! elle-meme (grille 1x1 uniquement). ===
  IF (PR.EQ.1 .AND. PC.EQ.1 .AND. N.LE.2000 .AND. MYID0.EQ.0) THEN
    CALL DEBUG_CHECK_SOLVE(N, NRHS, B_LOCAL, MLOC)
  END IF

  ! === 4) resolution 2D, AVEC pivot ===
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
  T1 = MPI_WTIME()
  T_SOLVE = T1 - T0
  TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, RSS_SOLVE_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)

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
    WRITE(6,'(A,I0)') '   Nombre de swaps de pivot effectues (total) = ', N_SWAPS
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Compression : ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', &
         100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F10.1,A,F10.1,A)') '   Memoire ANALYTIQUE (BLR_BLOCK_T, max/processus) = ', MEM_ANALYTIC_MB, &
         ' Mo  ;  RSS apres factorisation = ', RSS_FACT_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)') '   RSS apres resolution (max/processus) = ', RSS_SOLVE_MB, ' Mo'
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_distribue - X_exact|| = ', ERR
    IF (ERR.LT.1.0E-3_WP .AND. INFO.EQ.0 .AND. N_SWAPS.GT.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT : OK ==='
    ELSE IF (N_SWAPS.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT : SUSPECT (aucun pivot exerce, matrice pas assez agressive) ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT : ECHEC ==='
    END IF
    DEALLOCATE(B_REF)
  END IF

  DEALLOCATE(B_LOCAL, XEX_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! A(GI,GJ) = 50 sur une anti-diagonale CONFINEE AU BLOC DIAGONAL de GI
  ! (mirroir LOCAL a [BSTART,BEND] du bloc contenant GI, PAS a [1,N] global
  ! comme test_mpi_lu_pivot.f90/6d) -- diagonale NATURELLE (GI=GJ) faible
  ! face a cette anti-diagonale locale : force un VRAI pivotage QUE LE
  ! PIVOT CONFINE PEUT VOIR (contrairement a la version globale, dont
  ! l'element fort peut tomber HORS du bloc en cours de factorisation --
  ! injuste pour un pivot confine par construction, cf. memoire projet).
  ! Sinon, meme noyau de Cauchy amorti que 6d.
  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ, NN)
    INTEGER, INTENT(IN) :: GI, GJ, NN
    REAL(WP) :: DF
    INTEGER :: BSTART0, BEND0, GJ_STRONG
    BSTART0 = ((GI-1)/NB)*NB + 1
    BEND0 = MIN(BSTART0+NB-1, NN)
    GJ_STRONG = BEND0 - (GI-BSTART0)
    IF (GJ .EQ. GJ_STRONG) THEN
      LOWRANK_A = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
    ELSE
      DF = 1.0_WP / (1.0_WP + ABS(REAL(GJ_STRONG-GJ,WP)))
      LOWRANK_A = CMPLX(DF, 0.3_WP*DF, KIND=WP)
    END IF
  END FUNCTION LOWRANK_A

  ! VERIFICATION INDEPENDANTE : reconstruit P*A_orig (dense) et L*U (dense, a partir
  ! de A_BLR factorisee) et compare -- 1x1 uniquement, indexation OFFBLK
  ! directe (li=I, lj=J) donc pas d'ambiguite grille.
  SUBROUTINE DEBUG_CHECK_FACTOR(AB, NN)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    INTEGER, INTENT(IN) :: NN
    COMPLEX(WP), ALLOCATABLE :: A_ORIG(:,:), PA(:,:), LFULL(:,:), UFULL(:,:), LU(:,:)
    COMPLEX(WP), PARAMETER :: CONE=(1.0_WP,0.0_WP), CZERO=(0.0_WP,0.0_WP)
    INTEGER :: K, I, J, KW, NJ, MI, II, JJ
    REAL(WP) :: ERRDBG, NORMDBG

    ALLOCATE(A_ORIG(NN,NN), PA(NN,NN), LFULL(NN,NN), UFULL(NN,NN), LU(NN,NN))
    CALL GET_BLOCK(1, NN, 1, NN, A_ORIG)
    PA = A_ORIG

    DO K = 1, AB%NBLK
      KW = AB%BEND(K) - AB%BSTART(K) + 1
      DO I = 1, KW
        IF (AB%IPIV(I,K) .NE. I) THEN
          CALL CSWAP(NN, PA(AB%BSTART(K)+I-1,1), NN, PA(AB%BSTART(K)+AB%IPIV(I,K)-1,1), NN)
        END IF
      END DO
    END DO

    LFULL = CZERO; UFULL = CZERO
    DO K = 1, AB%NBLK
      KW = AB%BEND(K) - AB%BSTART(K) + 1
      DO II = 1, KW
        LFULL(AB%BSTART(K)+II-1, AB%BSTART(K)+II-1) = CONE
        DO JJ = 1, KW
          IF (JJ .LT. II) THEN
            LFULL(AB%BSTART(K)+II-1, AB%BSTART(K)+JJ-1) = AB%DIAGBLK(K)%D(II,JJ)
          ELSE
            UFULL(AB%BSTART(K)+II-1, AB%BSTART(K)+JJ-1) = AB%DIAGBLK(K)%D(II,JJ)
          END IF
        END DO
      END DO
      DO J = 1, AB%NBLK
        IF (J .EQ. K) CYCLE
        NJ = AB%BEND(J) - AB%BSTART(J) + 1
        IF (J .LT. K) THEN
          CALL EXPAND_BLOCK_LOCAL(AB%OFFBLK(K,J), LFULL(AB%BSTART(K):AB%BEND(K), AB%BSTART(J):AB%BEND(J)), KW, NJ)
        ELSE
          CALL EXPAND_BLOCK_LOCAL(AB%OFFBLK(K,J), UFULL(AB%BSTART(K):AB%BEND(K), AB%BSTART(J):AB%BEND(J)), KW, NJ)
        END IF
      END DO
    END DO

    CALL CGEMM('No transpose', 'No transpose', NN, NN, NN, CONE, LFULL, NN, UFULL, NN, CZERO, LU, NN)
    NORMDBG = SQRT(SUM(ABS(PA)**2))
    ERRDBG = SQRT(SUM(ABS(PA-LU)**2)) / MAX(NORMDBG, TINY(1.0_WP))
    WRITE(6,'(A,ES12.4)') ' [DEBUG] ||P*A - L*U|| / ||P*A|| = ', ERRDBG
    DEALLOCATE(A_ORIG, PA, LFULL, UFULL, LU)
  END SUBROUTINE DEBUG_CHECK_FACTOR

  ! VERIFICATION INDEPENDANTE : resout A_orig*X=B_reel (celui QUE la resolution va
  ! utiliser) via CGESV independant, sans passer par notre factorisation --
  ! ground truth pour isoler un bug dans MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D
  ! lui-meme (grille 1x1 uniquement, N raisonnable).
  SUBROUTINE DEBUG_CHECK_SOLVE(NN, NR, BLOC, MLC)
    INTEGER, INTENT(IN) :: NN, NR, MLC
    COMPLEX(WP), INTENT(IN) :: BLOC(MAX(MLC,1),MAX(NR,1))
    COMPLEX(WP), ALLOCATABLE :: A_ORIG(:,:), B_REF(:,:)
    INTEGER, ALLOCATABLE :: IPIV_REF(:)
    INTEGER :: INFO_REF
    REAL(WP) :: ERRDBG, NORMDBG

    ALLOCATE(A_ORIG(NN,NN), B_REF(NN,NR), IPIV_REF(NN))
    CALL GET_BLOCK(1, NN, 1, NN, A_ORIG)
    B_REF = BLOC(1:NN,1:NR)
    CALL CGESV(NN, NR, A_ORIG, NN, IPIV_REF, B_REF, NN, INFO_REF)
    NORMDBG = SQRT(SUM(ABS(B_REF)**2))
    ERRDBG = SQRT(SUM(ABS(B_REF-XEX_LOCAL(1:NN,1:NR))**2)) / MAX(NORMDBG, TINY(1.0_WP))
    WRITE(6,'(A,ES12.4)') ' [DEBUG] ||CGESV(A,B) - Xexact|| / ||CGESV(A,B)|| = ', ERRDBG
    DEALLOCATE(A_ORIG, B_REF, IPIV_REF)
  END SUBROUTINE DEBUG_CHECK_SOLVE

  SUBROUTINE EXPAND_BLOCK_LOCAL(BLK, OUT, M, N)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: BLK
    INTEGER, INTENT(IN) :: M, N
    COMPLEX(WP), INTENT(OUT) :: OUT(M,N)
    COMPLEX(WP), PARAMETER :: CONE=(1.0_WP,0.0_WP), CZERO=(0.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CGEMM('No transpose','No transpose', M, N, BLK%RANK, CONE, BLK%UFAC, M, BLK%VFAC, BLK%RANK, CZERO, OUT, M)
    ELSE
      OUT = BLK%DENSE
    END IF
  END SUBROUTINE EXPAND_BLOCK_LOCAL

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
      PROWJ = MOD(J-1, GD%PR); PCOLJ = BLR_COL_OWNER(AB%CMAP, J, GD)

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

      MLOC_ALL = 0
      DO li = 1, AB%NBLKR_LOC
        I = (li-1)*GD%PR + GD%MYROW + 1
        MLOC_ALL = MLOC_ALL + (AB%BEND(I)-AB%BSTART(I)+1)
      END DO
      ALLOCATE(CONTRIB(MLOC_ALL,NR))
      IF (GD%MYCOL .EQ. PCOLJ) THEN
        lj = BLR_COL_LOCAL(AB%CMAP, J, GD)
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

END PROGRAM TEST_MPI_BLR_MATRIX_PIVOT
