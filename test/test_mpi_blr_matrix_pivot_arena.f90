!===============================================================================
! Test du solveur LU MPI, BLR compresse des l'assemblage, AVEC PIVOT PARTIEL
! CONFINE AU BLOC DIAGONAL, VARIANTE ARENA (MPI_ASSEMBLE_BLR_MATRIX_ARENA /
! MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA / MPI_BLR_MATRIX_LU_SOLVE_PIVOT_ARENA_2D
! -- voir mpi_blr_matrix_mod.f90, type MPI_BLR_MATRIX_ARENA_T) -- MEME
! methodologie de validation que test_mpi_blr_matrix_pivot.f90 (la reference
! non-arena, JAMAIS modifiee), pour promouvoir l'arena au meme rang de
! "solveur de reference valide" -- voir README.md.
!
! POURQUOI CETTE VARIANTE EST UN SOLVEUR DE REFERENCE A PART ENTIERE, PAS UN
! SIMPLE EXEMPLE : elle resout le meme probleme (correctif de l'ecart
! RSS-vs-theorique cause par la fragmentation de tas, voir la memoire du
! projet) avec un compromis MEMOIRE/TEMPS mesure et documente -- a N=50000,
! NB=500 : RSS reelle inferieure a la reference classique (-11.5%, la
! fragmentation de l'existant devient le facteur dominant a cette echelle) et
! factorisation plus rapide (-30%), au prix d'un assemblage plus lent (+19%)
! et d'une empreinte THEORIQUE plus elevee (le compromis assume : capacite
! haute-marque + orphelins bas-rang->dense jamais recuperes, voir
! blr_arena_mod.f90).
!
! MATRICE ET VALIDATION : IDENTIQUES a test_mpi_blr_matrix_pivot.f90 --
! "anti-diagonale forte / diagonale naturelle faible" CONFINEE au bloc
! diagonal de GI (force un vrai pivotage QUE LE PIVOT CONFINE PEUT VOIR),
! B = A_BLR*Xexact (formule deterministe), resoudre doit redonner Xexact.
! N_SWAPS (rapporte) doit etre > 0 pour confirmer que le pivot a reellement
! ete exerce.
!
! Usage : mpirun -np <PR*PC> ./test_mpi_blr_matrix_pivot_arena PR PC N NB NRHS TOL
!===============================================================================
PROGRAM TEST_MPI_BLR_MATRIX_PIVOT_ARENA
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_ARENA_T), TARGET :: A_BLR
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
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./test_mpi_blr_matrix_pivot_arena PR PC N NB NRHS TOL'
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

  ! === 1) assemblage BLR distribue, arena ===
  CALL MPI_ASSEMBLE_BLR_MATRIX_ARENA(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                      N_COMPRESSED, N_OFFDIAG, RANK_SUM)

  ! === 2) factorisation LU, AVEC pivot confine au bloc diagonal, arena ===
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA(GRID, A_BLR, TOL, INFO, &
                                             N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
  CALL MPI_ALLREDUCE(TMP_MB, RSS_FACT_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
  TMP_MB = REAL(MPI_BLR_MATRIX_ARENA_MEMORY_BYTES(A_BLR),8)/1024.0D0/1024.0D0
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

  BLOCK
    TYPE(MPI_BLR_MATRIX_ARENA_T), TARGET :: A_BLR_FOR_B
    DOUBLE PRECISION :: T_ASM_DUMMY
    INTEGER :: NC_DUMMY, NO_DUMMY
    INTEGER(8) :: RS_DUMMY
    CALL MPI_ASSEMBLE_BLR_MATRIX_ARENA(GRID, N, NB, TOL, GET_BLOCK, A_BLR_FOR_B, T_ASM_DUMMY, &
                                        NC_DUMMY, NO_DUMMY, RS_DUMMY)
    CALL APPLY_MPI_BLR_MATRIX(GRID, A_BLR_FOR_B, NRHS, NRHS_LOC, XEX_LOCAL, B_LOCAL, MLOC)
  END BLOCK

  ! === 4) resolution 2D, AVEC pivot, arena ===
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_ARENA_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
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

    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)') ' [ARENA] Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, &
         ' , NRHS=', NRHS, ' , TOL=', TOL
    WRITE(6,'(A,F8.3,A,F8.3,A,F8.3,A)') '   Temps assemblage = ', T_ASM, ' s , factorisation = ', T_FACTOR, &
         ' s , resolution = ', T_SOLVE, ' s'
    WRITE(6,'(A,I0)') '   Nombre de swaps de pivot effectues (total) = ', N_SWAPS
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   Compression : ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', &
         100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F10.1,A,F10.1,A)') '   Memoire ANALYTIQUE (arena, max/processus) = ', MEM_ANALYTIC_MB, &
         ' Mo  ;  RSS apres factorisation = ', RSS_FACT_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)') '   RSS apres resolution (max/processus) = ', RSS_SOLVE_MB, ' Mo'
    WRITE(6,'(A,ES12.4)') ' Erreur relative ||X_distribue - X_exact|| = ', ERR
    IF (ERR.LT.1.0E-3_WP .AND. INFO.EQ.0 .AND. N_SWAPS.GT.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT_ARENA : OK ==='
    ELSE IF (N_SWAPS.EQ.0) THEN
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT_ARENA : SUSPECT (aucun pivot exerce, matrice pas assez agressive) ==='
    ELSE
      WRITE(6,*) '=== TEST MPI_BLR_MATRIX_PIVOT_ARENA : ECHEC ==='
    END IF
    DEALLOCATE(B_REF)
  END IF

  DEALLOCATE(B_LOCAL, XEX_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! A(GI,GJ) = 50 sur une anti-diagonale CONFINEE AU BLOC DIAGONAL de GI
  ! (mirroir LOCAL a [BSTART,BEND] du bloc contenant GI) -- diagonale
  ! NATURELLE (GI=GJ) faible face a cette anti-diagonale locale : force un
  ! VRAI pivotage QUE LE PIVOT CONFINE PEUT VOIR. Sinon, noyau de Cauchy
  ! amorti -- IDENTIQUE a test_mpi_blr_matrix_pivot.f90 (comparaison directe).
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
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)                 :: GD
    INTEGER, INTENT(IN)                      :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)                 :: XL(MAX(ML,1),MAX(NRL,1))
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
    TYPE(GRID_T), INTENT(IN)                 :: GD
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: AB
    INTEGER, INTENT(IN)                      :: NR, NRL, ML
    COMPLEX(WP), INTENT(IN)                  :: XL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), INTENT(OUT)                 :: YL(MAX(ML,1),MAX(NRL,1))
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
        RS = LOCAL_ROW_OFFSET_ARENA(AB, GD, J)
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
            CALL BLR_BLOCK_APPLY_ARENA(AB%OFFBLK(li,lj), NR, CONE, XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
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

END PROGRAM TEST_MPI_BLR_MATRIX_PIVOT_ARENA
