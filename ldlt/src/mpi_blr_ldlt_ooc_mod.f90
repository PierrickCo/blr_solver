!===============================================================================
! Variante OUT-OF-CORE du LDL^T BLR+MPI (mpi_blr_ldlt_mod.f90) -- meme
! algorithme numerique ; matrice, facteur et second membre vivent sur un
! fichier scratch PAR PROCESSUS, la RAM residente etant bornee par
! MAX_MB_PER_PROCESS via les deux caches LRU de blr_ooc_mod.f90 (CACHE_A
! pour les blocs hors-diagonale, CACHE_D pour les blocs diagonaux) -- meme
! infrastructure, memes conventions et memes mises en garde que
! MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC (voir le commentaire "DIMENSIONNER
! MAX_MB_PER_PROCESS" de MPI_BLR_MATRIX_OOC_T : un budget sous l'ensemble de
! travail reel fait thrasher le cache, il ne degrade pas gracieusement).
!
! Specificites LDL^T :
!   - seuls les blocs du triangle inferieur ont un descripteur rempli et une
!     place sur disque ;
!   - un bloc diagonal est stocke COMPACT (triangle inferieur, KW*(KW+1)/2
!     elements) sous la forme d'un bloc dense de taille NP x 1 : aucune
!     modification de blr_ooc_mod.f90 n'est necessaire, et disque comme
!     cache ne paient que la moitie d'un bloc plein ;
!   - un bloc du panneau K (colonne K) est evince des que L(I,K) est
!     calcule et diffuse (W) ; il sera relu une fois a l'etape I pour
!     l'echange de lignes de Bunch-Kaufman, puis a la resolution.
!
! Second membre : meme mecanisme que le LU OOC (OOC_RHS_T, paquets de
! NRHS_BATCH colonnes sur un fichier scratch separe, CONSUME_X_BATCH appele
! une fois par paquet). La remontee L^T (pleine largeur) garde, par paquet, la copie pleine
! largeur XFULL (MLOC x NRHS_BATCH) de la version classique.
!===============================================================================
MODULE MPI_BLR_LDLT_OOC_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE BLR_OOC_MOD
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE
  USE MPI_BLR_MATRIX_MOD, ONLY: BLOCK_NUMROC, BLR_COL_OWNER, BLR_COL_LOCAL, BLR_COL_GLOBAL, BLR_NBLKC_LOC, &
                                 STAGE_PROGRESS_CALLBACK_IFACE, STAGE_PROGRESS_SHOULD_FIRE, &
                                RHS_LOCAL_RANGE_MBM, OOC_RHS_T, OOC_RHS_FREE, FILL_RHS_IFACE, CONSUME_X_BATCH_IFACE
  USE MPI_BLR_LDLT_MOD, ONLY: SYM_PACK_LOWER, SYM_UNPACK_LOWER, BLR_BLOCK_COL_SWAP, &
                              BLR_TRSM_RIGHT_LOWER_TRANS_UNIT, BLR_BLOCK_DINV_RIGHT, BLR_BLOCK_TRANSPOSE, &
                              BLR_BLOCK_APPLY_T, APPLY_DINV_LEFT, COUNT_PIVOTS, &
                              SYM_BK_FACTOR, SYM_BK_LWORK, SYM_PACKED_SCHUR_UPDATE
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_BAR
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE

  TYPE :: MPI_BLR_SYM_OOC_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (BLR_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(OOC_BLOCK_DESC_T), ALLOCATABLE :: DIAGBLK(:)      ! bloc compact NP x 1
    TYPE(OOC_BLOCK_DESC_T), ALLOCATABLE :: OFFBLK(:,:)     ! (li,lj), rempli seulement si I>J
    INTEGER, ALLOCATABLE     :: IPIV(:,:)
    COMPLEX(WP), ALLOCATABLE :: DDIAG(:,:), DSUB(:,:)
    TYPE(OOC_STORE_T)      :: STORE_A
    TYPE(OOC_CACHE_T)      :: CACHE_A
    TYPE(OOC_DIAG_CACHE_T) :: CACHE_D
    REAL(WP) :: MAX_MB_PER_PROCESS = 0.0_WP
  END TYPE MPI_BLR_SYM_OOC_T

CONTAINS

  FUNCTION LOCAL_ROW_OFFSET_SYM_OOC(A, GRID, K) RESULT(RS)
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(IN) :: A
    TYPE(GRID_T), INTENT(IN)            :: GRID
    INTEGER, INTENT(IN)                 :: K
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (K-1)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + A%BEND(I) - A%BSTART(I) + 1
    END DO
  END FUNCTION LOCAL_ROW_OFFSET_SYM_OOC

  FUNCTION MPI_BLR_SYM_OOC_MLOC(A, GRID) RESULT(MLOC)
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(IN) :: A
    TYPE(GRID_T), INTENT(IN)            :: GRID
    INTEGER :: MLOC, li, I
    MLOC = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + A%BEND(I) - A%BSTART(I) + 1
    END DO
  END FUNCTION MPI_BLR_SYM_OOC_MLOC

  ! Octets REELLEMENT residents (les deux caches) -- a comparer au budget.
  FUNCTION MPI_BLR_SYM_OOC_RESIDENT_BYTES(A) RESULT(NBYTES)
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(IN) :: A
    INTEGER(8) :: NBYTES
    NBYTES = A%CACHE_A%USED_BYTES + A%CACHE_D%USED_BYTES
  END FUNCTION MPI_BLR_SYM_OOC_RESIDENT_BYTES

  ! Taille LOGIQUE des donnees sur disque (blocs diagonaux compacts + blocs
  ! I>J a leur taille courante) -- equivalent de MPI_BLR_SYM_MEMORY_BYTES.
  FUNCTION MPI_BLR_SYM_OOC_DATA_BYTES(A) RESULT(NBYTES)
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(IN) :: A
    INTEGER(8) :: NBYTES
    INTEGER :: K, li, lj
    NBYTES = 0_8
    DO K = 1, A%NBLK
      IF (A%DIAGBLK(K)%M .GT. 0) NBYTES = NBYTES + OOC_DESC_BYTES(A%DIAGBLK(K))
    END DO
    DO lj = 1, A%NBLKC_LOC
      DO li = 1, A%NBLKR_LOC
        IF (A%OFFBLK(li,lj)%M .GT. 0) NBYTES = NBYTES + OOC_DESC_BYTES(A%OFFBLK(li,lj))
      END DO
    END DO
  END FUNCTION MPI_BLR_SYM_OOC_DATA_BYTES

  SUBROUTINE MPI_BLR_SYM_OOC_FREE(A)
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(INOUT) :: A
    CALL OOC_CACHE_FREE(A%CACHE_A)
    CALL OOC_DIAG_CACHE_FREE(A%CACHE_D)
    CALL OOC_STORE_FREE(A%STORE_A)
  END SUBROUTINE MPI_BLR_SYM_OOC_FREE

  ! Diffuse le bloc diagonal K (compact) dans COMM depuis ROOT et le deplie
  ! dans DKK(KW,KW). Chez le proprietaire, lu via CACHE_D.
  SUBROUTINE BCAST_DIAG_LOWER_OOC(A, K, KW, OWNER, ROOT, COMM, DKK)
    INCLUDE 'mpif.h'
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(INOUT) :: A
    INTEGER, INTENT(IN)             :: K, KW, ROOT, COMM
    LOGICAL, INTENT(IN)             :: OWNER
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT) :: DKK(:,:)
    COMPLEX(WP), ALLOCATABLE :: APK(:)
    INTEGER :: NP, IERR, IDX_D
    NP = KW*(KW+1)/2
    ALLOCATE(APK(NP), DKK(KW,KW))
    IF (OWNER) THEN
      IDX_D = OOC_DIAG_CACHE_GET(A%CACHE_D, A%STORE_A, A%DIAGBLK, K)
      APK = A%CACHE_D%SLOTS(IDX_D)%BLK%DENSE(1:NP,1)
    END IF
    CALL MPI_BCAST(APK, NP, MPI_COMPLEX, ROOT, COMM, IERR)
    DKK = (0.0_WP,0.0_WP)
    CALL SYM_UNPACK_LOWER(KW, APK, DKK, KW)
    DEALLOCATE(APK)
  END SUBROUTINE BCAST_DIAG_LOWER_OOC

  !=======================================================================
  ! Assemblage OOC (triangle inferieur) : blocs diagonaux compacts ecrits
  ! directement sur STORE_A, blocs I>J inseres dans CACHE_A (evinces au fil
  ! de l'eau si le budget est atteint), tout est sur disque en sortie.
  !=======================================================================
  SUBROUTINE MPI_ASSEMBLE_BLR_SYM_OOC(GRID, N, NB, TOL, GET_BLOCK, A, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                      T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                      SHOW_PROGRESS, USE_RSVD, PROGRESS_CALLBACK, DIAG_SPREAD)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)              :: GRID
    INTEGER, INTENT(IN)                   :: N, NB
    REAL(WP), INTENT(IN)                  :: TOL
    PROCEDURE(GET_BLOCK_IFACE)            :: GET_BLOCK
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(OUT)  :: A
    CHARACTER(LEN=*), INTENT(IN)          :: SCRATCH_DIR
    REAL(WP), INTENT(IN)                  :: MAX_MB_PER_PROCESS
    DOUBLE PRECISION, INTENT(OUT)         :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)                  :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)               :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL         :: SHOW_PROGRESS, USE_RSVD
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir BLR_COL_OWNER) -- les blocs diagonaux visitent tous
    ! les processus ; .FALSE. : cyclique pur, comme le LU.
    LOGICAL, INTENT(IN), OPTIONAL            :: DIAG_SPREAD

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    TYPE(BLR_BLOCK_T) :: TMPBLK, DBLK
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, IDX, PBAR_STATE, NP
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, T_ASM_L
    LOGICAL :: PROGRESS, RSVD_FLAG

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD
    PBAR_STATE = -1
    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A%NBLK, A%BSTART, A%BEND)
    A%N = N; A%NB = NB
    A%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A%NBLKR_LOC = BLOCK_NUMROC(A%NBLK, GRID%PR, GRID%MYROW)
    A%NBLKC_LOC = BLR_NBLKC_LOC(A%CMAP, A%NBLK, GRID)
    A%MAX_MB_PER_PROCESS = MAX_MB_PER_PROCESS
    ALLOCATE(A%DIAGBLK(A%NBLK))
    ALLOCATE(A%OFFBLK(MAX(A%NBLKR_LOC,1), MAX(A%NBLKC_LOC,1)))
    NTOTAL = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .GE. J) NTOTAL = NTOTAL + 1
      END DO
    END DO
    CALL OOC_STORE_INIT(A%STORE_A, SCRATCH_DIR, 'ldlt_a', GRID%MYID)
    CALL OOC_CACHE_INIT(A%CACHE_A, MAX_MB_PER_PROCESS, NTOTAL+4, MAX(A%NBLKR_LOC,1), MAX(A%NBLKC_LOC,1))
    CALL OOC_DIAG_CACHE_INIT(A%CACHE_D, MAX_MB_PER_PROCESS, A%NBLKR_LOC+2, A%NBLK)

    NDONE = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A%BEND(I) - A%BSTART(I) + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LT. J) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        IF (I .EQ. J) THEN
          NP = MI*(MI+1)/2
          ALLOCATE(TMP(MI,MI))
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), TMP)
          DBLK%IS_LOWRANK = .FALSE.
          DBLK%M = NP; DBLK%N = 1; DBLK%RANK = 0
          ALLOCATE(DBLK%DENSE(NP,1))
          CALL SYM_PACK_LOWER(MI, TMP, MI, DBLK%DENSE)
          DEALLOCATE(TMP)
          CALL OOC_WRITE_BLOCK_FIXED(A%STORE_A, A%DIAGBLK(I), DBLK)
          DEALLOCATE(DBLK%DENSE)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), TMP)
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, TMPBLK)
          ELSE
            CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, TMPBLK)
          END IF
          DEALLOCATE(TMP)
          NO_LOCAL = NO_LOCAL + 1
          IF (TMPBLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + TMPBLK%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(MI,NJ)
          END IF
          A%OFFBLK(li,lj)%IS_LOWRANK = TMPBLK%IS_LOWRANK
          A%OFFBLK(li,lj)%M = TMPBLK%M; A%OFFBLK(li,lj)%N = TMPBLK%N; A%OFFBLK(li,lj)%RANK = TMPBLK%RANK
          IDX = OOC_CACHE_PUT_NEW(A%CACHE_A, A%STORE_A, A%OFFBLK, li, lj, TMPBLK)
        END IF
        NDONE = NDONE + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('1) Assemblage+compression LDLt ooc (rang 0)', NDONE, NTOTAL, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(NDONE, NTOTAL)) CALL PROGRESS_CALLBACK(NDONE, NTOTAL)
        END IF
      END DO
    END DO
    CALL OOC_CACHE_FLUSH_ALL(A%CACHE_A, A%STORE_A, A%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A%CACHE_D, A%STORE_A, A%DIAGBLK)
    T_ASM_L = MPI_WTIME() - T0
    CALL MPI_ALLREDUCE(T_ASM_L, T_ASSEMBLE_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_SYM_OOC

  !=======================================================================
  ! Second membre OOC : paquets de NRHS_BATCH colonnes ecrits sur un
  ! fichier scratch separe (meme format OOC_RHS_T que le LU OOC).
  !=======================================================================
  SUBROUTINE MPI_ASSEMBLE_RHS_SYM_OOC(GRID, A, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
    TYPE(GRID_T), INTENT(IN)              :: GRID
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(IN)   :: A
    INTEGER, INTENT(IN)                   :: NRHS, NRHS_BATCH
    PROCEDURE(FILL_RHS_IFACE)             :: FILL_RHS
    CHARACTER(LEN=*), INTENT(IN)          :: SCRATCH_DIR
    TYPE(OOC_RHS_T), INTENT(OUT)          :: RHS
    TYPE(BLR_BLOCK_T) :: BBLK
    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, I, RS, WIDTH, B, RSTART, REND, LOCAL_OFFSET, NRHS_LOC_B, GSTART, GEND

    RHS%NRHS = NRHS
    RHS%NRHS_BATCH = MAX(NRHS_BATCH,1)
    RHS%NBATCHES = (NRHS + RHS%NRHS_BATCH - 1) / RHS%NRHS_BATCH
    RHS%MLOC = MPI_BLR_SYM_OOC_MLOC(A, GRID)
    CALL OOC_STORE_INIT(RHS%STORE_B, SCRATCH_DIR, 'ldlt_b', GRID%MYID)
    ALLOCATE(RHS%BATCH(RHS%NBATCHES))
    DO B = 1, RHS%NBATCHES
      RSTART = (B-1)*RHS%NRHS_BATCH + 1
      REND = MIN(B*RHS%NRHS_BATCH, NRHS)
      WIDTH = REND - RSTART + 1
      CALL RHS_LOCAL_RANGE_MBM(WIDTH, GRID%PC, GRID%MYCOL, LOCAL_OFFSET, NRHS_LOC_B)
      GSTART = RSTART + LOCAL_OFFSET - 1
      GEND = GSTART + NRHS_LOC_B - 1
      ALLOCATE(TMP(MAX(RHS%MLOC,1),MAX(NRHS_LOC_B,1)))
      TMP = (0.0_WP,0.0_WP)
      RS = 1
      DO li = 1, A%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (NRHS_LOC_B .GT. 0) CALL FILL_RHS(A%BSTART(I), A%BEND(I), GSTART, GEND, &
                                             TMP(RS:RS+A%BEND(I)-A%BSTART(I),1:NRHS_LOC_B))
        RS = RS + (A%BEND(I)-A%BSTART(I)+1)
      END DO
      BBLK%IS_LOWRANK = .FALSE.
      BBLK%M = MAX(RHS%MLOC,1); BBLK%N = MAX(NRHS_LOC_B,1)
      CALL MOVE_ALLOC(TMP, BBLK%DENSE)
      CALL OOC_WRITE_BLOCK_RAW(RHS%STORE_B, RHS%BATCH(B), BBLK)
      DEALLOCATE(BBLK%DENSE)
    END DO
  END SUBROUTINE MPI_ASSEMBLE_RHS_SYM_OOC

  !=======================================================================
  ! Factorisation LDL^T OOC -- memes 6 sous-etapes que MPI_BLR_SYM_LDLT_FACTOR,
  ! chaque acces a un bloc passant par CACHE_A / CACHE_D. T_CACHE/T_COMPUTE/
  ! T_MPICALL (optionnels, MAX sur la grille) : meme decoupage que le LU OOC
  ! (T_CACHE regroupe lectures ET evictions de panneau).
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR_OOC(GRID, A, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                         N_SWAPS, N_2X2, T_FACTOR, SHOW_PROGRESS, SHOW_RSS_TRACE, &
                                         T_CACHE, T_COMPUTE, T_MPICALL, PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)               :: GRID
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(INOUT) :: A
    REAL(WP), INTENT(IN)                   :: TOL
    INTEGER, INTENT(OUT)                   :: INFO, N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)                :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL         :: N_SWAPS, N_2X2
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR, T_CACHE, T_COMPUTE, T_MPICALL
    LOGICAL, INTENT(IN), OPTIONAL          :: SHOW_PROGRESS, SHOW_RSS_TRACE
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), DBUF(:), E(:), WORK(:), DWORK(:,:)
    COMPLEX(WP) :: WQUERY(1)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: WRECV(:), LTRECV(:)
    TYPE(BLR_BLOCK_T) :: TMPB
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:), SWP(:), IPIV_Q(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP, LWORK, ROOT_K, IDX, IDX_D, NP
    INTEGER :: PROW_K, PCOL_K, PROW_J, IERR, NC_LOCAL, NO_LOCAL, li_K, lj_K, li_J, PBAR_STATE
    INTEGER :: NS_TOT, N2_TOT, NS_K, N2_K
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0, T_FACTOR_L, RSS_MAX_MB, CACHE_MAX_MB, TL(3), TG(3)
    LOGICAL :: PROGRESS, RSS_TRACE, OWNER

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    RSS_TRACE = .FALSE.
    IF (PRESENT(SHOW_RSS_TRACE)) RSS_TRACE = SHOW_RSS_TRACE .AND. (GRID%MYID.EQ.0)
    PBAR_STATE = -1
    RSS_MAX_MB = 0.0D0; CACHE_MAX_MB = 0.0D0
    TL = 0.0D0          ! 1 = cache (lectures + evictions), 2 = calcul, 3 = MPI
    T0 = MPI_WTIME()
    NBLK = A%NBLK
    INFO = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8; NS_TOT = 0; N2_TOT = 0
    ALLOCATE(WRECV(MAX(A%NBLKR_LOC,1)), LTRECV(MAX(A%NBLKC_LOC,1)))
    IF (ALLOCATED(A%IPIV))  DEALLOCATE(A%IPIV)
    IF (ALLOCATED(A%DDIAG)) DEALLOCATE(A%DDIAG)
    IF (ALLOCATED(A%DSUB))  DEALLOCATE(A%DSUB)
    ALLOCATE(A%IPIV(A%NB,NBLK), A%DDIAG(A%NB,NBLK), A%DSUB(A%NB,NBLK))
    A%IPIV = 0; A%DDIAG = (0.0_WP,0.0_WP); A%DSUB = (0.0_WP,0.0_WP)
    CALL BLR_WORKSPACE_ALLOC(WS, A%NB)
    ALLOCATE(E(A%NB), IPIV_Q(A%NB))
    LWORK = SYM_BK_LWORK(A%NB)
    ALLOCATE(WORK(LWORK), DWORK(A%NB,A%NB))
    DWORK = (0.0_WP,0.0_WP)

    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      NP = KW*(KW+1)/2
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A%CMAP, K, GRID)
      ROOT_K = PROW_K*GRID%PC + PCOL_K
      OWNER = (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K)
      li_K = (K-1)/GRID%PR + 1
      lj_K = BLR_COL_LOCAL(A%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW), SWP(KW), DBUF(2*KW))
      IPIV_K = 0; DBUF = (0.0_WP,0.0_WP)

      ! 1) Bunch-Kaufman borne sur le bloc diagonal (via CACHE_D)
      IF (OWNER) THEN
        TT0 = MPI_WTIME()
        IDX_D = OOC_DIAG_CACHE_GET(A%CACHE_D, A%STORE_A, A%DIAGBLK, K)
        TL(1) = TL(1) + (MPI_WTIME() - TT0)
        TT0 = MPI_WTIME()
        CALL SYM_UNPACK_LOWER(KW, A%CACHE_D%SLOTS(IDX_D)%BLK%DENSE, DWORK, A%NB)
        CALL SYM_BK_FACTOR(KW, DWORK, A%NB, E, IPIV_K, WORK, LWORK, IINFO)
        IF (IINFO .NE. 0) INFO_STEP = A%BSTART(K) - 1 + IINFO
        DO I = 1, KW
          DBUF(I) = DWORK(I,I)
        END DO
        DBUF(KW+1:2*KW) = E(1:KW)
        CALL SYM_PACK_LOWER(KW, DWORK, A%NB, A%CACHE_D%SLOTS(IDX_D)%BLK%DENSE)
        CALL OOC_DIAG_CACHE_MARK_DIRTY(A%CACHE_D, IDX_D)
        TL(2) = TL(2) + (MPI_WTIME() - TT0)
      END IF
      TT0 = MPI_WTIME()
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(WRECV, LTRECV, IPIV_K, SWP, DBUF, E, WORK, DWORK, IPIV_Q)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF
      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
      CALL MPI_BCAST(DBUF, 2*KW, MPI_COMPLEX, ROOT_K, GRID%COMM, IERR)
      TL(3) = TL(3) + (MPI_WTIME() - TT0)
      A%IPIV(1:KW,K)  = IPIV_K
      A%DDIAG(1:KW,K) = DBUF(1:KW)
      A%DSUB(1:KW,K)  = DBUF(KW+1:2*KW)
      SWP = ABS(IPIV_K)
      CALL COUNT_PIVOTS(KW, IPIV_K, DBUF(KW+1:2*KW), NS_K, N2_K)
      NS_TOT = NS_TOT + NS_K; N2_TOT = N2_TOT + N2_K

      ! 2) lignes des L(K,J<K) (via CACHE_A)
      IF (GRID%MYROW .EQ. PROW_K .AND. NS_K .GT. 0) THEN
        DO lj = 1, A%NBLKC_LOC
          J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
          IF (J .GE. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A%CACHE_A, A%STORE_A, A%OFFBLK, li_K, lj)
          TL(1) = TL(1) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BLR_BLOCK_ROW_SWAP(A%CACHE_A%SLOTS(IDX)%BLK, KW, SWP)
          CALL OOC_CACHE_MARK_DIRTY(A%CACHE_A, IDX)
          TL(2) = TL(2) + (MPI_WTIME() - TT0)
        END DO
      END IF

      ! 3+4) colonne PCOL_K : L_KK (compact) diffuse dans la colonne, puis par
      !      bloc I>K : echange de colonnes, W = A*P*L_KK^-T, diffusion de W
      !      le long de la ligne, L = W*D^-1 ecrit en cache puis EVINCE.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        TT0 = MPI_WTIME()
        CALL BCAST_DIAG_LOWER_OOC(A, K, KW, OWNER, PROW_K, GRID%COL_COMM, DKK)
        TL(3) = TL(3) + (MPI_WTIME() - TT0)
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A%CACHE_A, A%STORE_A, A%OFFBLK, li, lj_K)
          TL(1) = TL(1) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          IF (NS_K .GT. 0) CALL BLR_BLOCK_COL_SWAP(A%CACHE_A%SLOTS(IDX)%BLK, KW, SWP)
          CALL BLR_TRSM_RIGHT_LOWER_TRANS_UNIT(A%CACHE_A%SLOTS(IDX)%BLK, DKK, KW)
          TL(2) = TL(2) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, A%CACHE_A%SLOTS(IDX)%BLK)
          TL(3) = TL(3) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          WRECV(li) = A%CACHE_A%SLOTS(IDX)%BLK
          CALL BLR_BLOCK_DINV_RIGHT(A%CACHE_A%SLOTS(IDX)%BLK, KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          CALL OOC_CACHE_MARK_DIRTY(A%CACHE_A, IDX)
          NO_LOCAL = NO_LOCAL + 1
          IF (A%CACHE_A%SLOTS(IDX)%BLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A%CACHE_A%SLOTS(IDX)%BLK%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A%CACHE_A%SLOTS(IDX)%BLK%M, A%CACHE_A%SLOTS(IDX)%BLK%N)
          END IF
          TL(2) = TL(2) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL OOC_CACHE_EVICT_SLOT(A%CACHE_A, A%STORE_A, A%OFFBLK, IDX)
          TL(1) = TL(1) + (MPI_WTIME() - TT0)
        END DO
        DEALLOCATE(DKK)
      ELSE
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, WRECV(li))
          TL(3) = TL(3) + (MPI_WTIME() - TT0)
        END DO
      END IF

      ! 5) L(J,K)^T le long des colonnes depuis (ligne(J), MYCOL)
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        PROW_J = MOD(J-1, GRID%PR)
        IF (GRID%MYROW .EQ. PROW_J) THEN
          li_J = (J-1)/GRID%PR + 1
          TMPB = WRECV(li_J)
          CALL BLR_BLOCK_DINV_RIGHT(TMPB, KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          CALL BLR_BLOCK_TRANSPOSE(TMPB, LTRECV(lj))
        END IF
        TT0 = MPI_WTIME()
        CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_J, KW, NJ, LTRECV(lj))
        TL(3) = TL(3) + (MPI_WTIME() - TT0)
      END DO

      ! 6) Schur, triangle inferieur (via CACHE_D / CACHE_A)
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LT. J) CYCLE
          MI = A%BEND(I) - A%BSTART(I) + 1
          IF (I .EQ. J) THEN
            TT0 = MPI_WTIME()
            IDX_D = OOC_DIAG_CACHE_GET(A%CACHE_D, A%STORE_A, A%DIAGBLK, I)
            TL(1) = TL(1) + (MPI_WTIME() - TT0)
            TT0 = MPI_WTIME()
            CALL SYM_PACKED_SCHUR_UPDATE(MI, A%CACHE_D%SLOTS(IDX_D)%BLK%DENSE, KW, WRECV(li), LTRECV(lj), &
                                         DWORK, A%NB)
            CALL OOC_DIAG_CACHE_MARK_DIRTY(A%CACHE_D, IDX_D)
            TL(2) = TL(2) + (MPI_WTIME() - TT0)
          ELSE
            TT0 = MPI_WTIME()
            IDX = OOC_CACHE_GET(A%CACHE_A, A%STORE_A, A%OFFBLK, li, lj)
            TL(1) = TL(1) + (MPI_WTIME() - TT0)
            TT0 = MPI_WTIME()
            CALL SCHUR_UPDATE_BLR_BLOCK(A%CACHE_A%SLOTS(IDX)%BLK, MI, NJ, KW, WRECV(li), LTRECV(lj), TOL, WS)
            A%OFFBLK(li,lj)%IS_LOWRANK = A%CACHE_A%SLOTS(IDX)%BLK%IS_LOWRANK
            A%OFFBLK(li,lj)%RANK = A%CACHE_A%SLOTS(IDX)%BLK%RANK
            CALL OOC_CACHE_MARK_DIRTY(A%CACHE_A, IDX)
            TL(2) = TL(2) + (MPI_WTIME() - TT0)
          END IF
        END DO
      END DO

      DEALLOCATE(IPIV_K, SWP, DBUF)
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('2) Factorisation LDLt ooc (Bunch-Kaufman confine)', K, NBLK, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(K, NBLK)) CALL PROGRESS_CALLBACK(K, NBLK)
      END IF
      IF (RSS_TRACE) THEN
        RSS_MAX_MB = MAX(RSS_MAX_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
        CACHE_MAX_MB = MAX(CACHE_MAX_MB, REAL(A%CACHE_A%USED_BYTES+A%CACHE_D%USED_BYTES,8)/1048576.0D0)
      END IF
    END DO
    IF (RSS_TRACE) WRITE(6,'(A,F8.1,A,F8.1,A)') '     [RSS trace ooc] MAX observe pendant la factorisation : RSS = ', &
         RSS_MAX_MB, ' Mo , cache resident = ', CACHE_MAX_MB, ' Mo'

    TT0 = MPI_WTIME()
    CALL OOC_CACHE_FLUSH_ALL(A%CACHE_A, A%STORE_A, A%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A%CACHE_D, A%STORE_A, A%DIAGBLK)
    TL(1) = TL(1) + (MPI_WTIME() - TT0)
    DEALLOCATE(WRECV, LTRECV, E, WORK, DWORK, IPIV_Q)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) THEN
      T_FACTOR_L = MPI_WTIME() - T0
      CALL MPI_ALLREDUCE(T_FACTOR_L, T_FACTOR, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
    CALL MPI_ALLREDUCE(TL, TG, 3, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_CACHE))   T_CACHE = TG(1)
    IF (PRESENT(T_COMPUTE)) T_COMPUTE = TG(2)
    IF (PRESENT(T_MPICALL)) T_MPICALL = TG(3)
    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_TOT
    IF (PRESENT(N_2X2))   N_2X2 = N2_TOT
  END SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR_OOC

  !=======================================================================
  ! Resolution OOC, paquet par paquet : B relu depuis RHS%STORE_B, meme
  ! algorithme que MPI_BLR_SYM_LDLT_SOLVE_2D (blocs lus via CACHE_A/CACHE_D,
  ! cache persistant entre paquets), puis CONSUME_X_BATCH avec la solution
  ! locale du paquet.
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_OOC_2D(GRID, A, RHS, CONSUME_X_BATCH, SHOW_PROGRESS, PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)               :: GRID
    TYPE(MPI_BLR_SYM_OOC_T), INTENT(INOUT) :: A
    TYPE(OOC_RHS_T), INTENT(INOUT)         :: RHS
    PROCEDURE(CONSUME_X_BATCH_IFACE)       :: CONSUME_X_BATCH
    LOGICAL, INTENT(IN), OPTIONAL          :: SHOW_PROGRESS
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE :: XFULL(:,:), S(:,:), SSUM(:,:), SLOC(:,:), B_LOCAL(:,:)
    TYPE(BLR_BLOCK_T) :: DTMP
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, MK, RS, RE, NROWTRAIL, C, IERR, P, MLOC, LDB, IDX
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, PBAR_STATE
    INTEGER :: B, RSTART, REND, NRHS, NRHS_LOC, LOCAL_OFFSET, GSTART, GEND
    LOGICAL :: PROGRESS, OWNER

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NBLK = A%NBLK
    MLOC = MPI_BLR_SYM_OOC_MLOC(A, GRID)
    LDB = MAX(MLOC,1)
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))

    DO B = 1, RHS%NBATCHES
      RSTART = (B-1)*RHS%NRHS_BATCH + 1
      REND = MIN(B*RHS%NRHS_BATCH, RHS%NRHS)
      NRHS = REND - RSTART + 1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, LOCAL_OFFSET, NRHS_LOC)
      GSTART = RSTART + LOCAL_OFFSET - 1
      GEND = GSTART + NRHS_LOC - 1
      DO C = 0, GRID%PC-1
        CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
      END DO
      CALL OOC_READ_BLOCK_RAW(RHS%STORE_B, RHS%BATCH(B), DTMP)
      CALL MOVE_ALLOC(DTMP%DENSE, B_LOCAL)
      STEP = 0
      PBAR_STATE = -1
      NSTEPS = 2*NBLK

      ! 0) Q*B
      DO K = 1, NBLK
        KW = A%BEND(K) - A%BSTART(K) + 1
        IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
        RS = LOCAL_ROW_OFFSET_SYM_OOC(A, GRID, K)
        DO I = 1, KW
          P = ABS(A%IPIV(I,K))
          IF (P .NE. I) CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+P-1,1), LDB)
        END DO
      END DO

      ! 1) descente L
      DO K = 1, NBLK
        KW = A%BEND(K) - A%BSTART(K) + 1
        PROW_K = MOD(K-1, GRID%PR)
        PCOL_K = BLR_COL_OWNER(A%CMAP, K, GRID)
        IF (GRID%MYROW .EQ. PROW_K) THEN
          CALL BCAST_DIAG_LOWER_OOC(A, K, KW, GRID%MYCOL.EQ.PCOL_K, PCOL_K, GRID%ROW_COMM, DKK)
          RS = LOCAL_ROW_OFFSET_SYM_OOC(A, GRID, K)
          IF (NRHS_LOC .GT. 0) CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                                           CONE, DKK, KW, B_LOCAL(RS,1), LDB)
          DEALLOCATE(DKK)
          ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
          IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        END IF
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          ALLOCATE(XK_FULL(KW,NRHS))
        ELSE
          ALLOCATE(XK_FULL(1,1))
        END IF
        IF (GRID%MYROW .EQ. PROW_K) THEN
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
          END DO
          CALL MPI_GATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           PCOL_K, GRID%ROW_COMM, IERR)
          DEALLOCATE(XK_LOC)
        END IF
        IF (GRID%MYCOL .EQ. PCOL_K) CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
        lj = BLR_COL_LOCAL(A%CMAP, K, GRID)
        NROWTRAIL = 0
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .GT. K) NROWTRAIL = NROWTRAIL + A%BEND(I)-A%BSTART(I)+1
        END DO
        IF (NROWTRAIL .GT. 0) THEN
          IF (GRID%MYCOL .EQ. PCOL_K) THEN
            ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
            CONTRIB = CZERO
            RE = 0
            DO li = 1, A%NBLKR_LOC
              I = (li-1)*GRID%PR + GRID%MYROW + 1
              IF (I .LE. K) CYCLE
              MI = A%BEND(I)-A%BSTART(I)+1
              IDX = OOC_CACHE_GET(A%CACHE_A, A%STORE_A, A%OFFBLK, li, lj)
              CALL BLR_BLOCK_APPLY(A%CACHE_A%SLOTS(IDX)%BLK, NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
              RE = RE + MI
            END DO
          ELSE
            ALLOCATE(CONTRIB(1,1))
          END IF
          ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
          END DO
          CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, &
                            PCOL_K, GRID%ROW_COMM, IERR)
          IF (NRHS_LOC .GT. 0) THEN
            RE = 0; RS = 1
            DO li = 1, A%NBLKR_LOC
              I = (li-1)*GRID%PR + GRID%MYROW + 1
              MI = A%BEND(I)-A%BSTART(I)+1
              IF (I .GT. K) THEN
                B_LOCAL(RS:RS+MI-1,1:NRHS_LOC) = B_LOCAL(RS:RS+MI-1,1:NRHS_LOC) - CONTRIB_LOC(RE+1:RE+MI,1:NRHS_LOC)
                RE = RE + MI
              END IF
              RS = RS + MI
            END DO
          END IF
          DEALLOCATE(CONTRIB, CONTRIB_LOC)
        END IF
        DEALLOCATE(XK_FULL)
        STEP = STEP + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt ooc', STEP, NSTEPS, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
        END IF
      END DO

      ! 2) D
      DO K = 1, NBLK
        KW = A%BEND(K) - A%BSTART(K) + 1
        IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
        RS = LOCAL_ROW_OFFSET_SYM_OOC(A, GRID, K)
        CALL APPLY_DINV_LEFT(B_LOCAL(RS,1), LDB, KW, NRHS_LOC, A%IPIV(1:KW,K), A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
      END DO

      ! 3) L^T*Y = W, remontee "pleine largeur" : chaque processus recupere
      !    une fois ses lignes de W sur TOUTE la largeur NRHS (XFULL, un
      !    ALLGATHERV par ligne de grille). Etape K : la colonne PCOL_K calcule
      !    sum_(I>K) L(I,K)^T*Y_I et la reduit sur (PROW_K,PCOL_K), proprietaire
      !    du bloc diagonal, qui resout L_KK^T et diffuse Y_K sur sa ligne.
      !    2 communications par etape, aucune diffusion du bloc diagonal.
      ALLOCATE(XFULL(LDB,NRHS))
      XFULL = CZERO
      ALLOCATE(XK_LOC(LDB,MAX(NRHS_LOC,1)))
      IF (MLOC.GT.0 .AND. NRHS_LOC.GT.0) XK_LOC(1:MLOC,1:NRHS_LOC) = B_LOCAL(1:MLOC,1:NRHS_LOC)
      DO C = 0, GRID%PC-1
        RCOUNTS(C) = MLOC*RLOCS(C); RDISPLS(C) = MLOC*(RSTARTS(C)-1)
      END DO
      CALL MPI_ALLGATHERV(XK_LOC, MLOC*NRHS_LOC, MPI_COMPLEX, XFULL, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                          GRID%ROW_COMM, IERR)
      DEALLOCATE(XK_LOC)
      DO K = NBLK, 1, -1
        KW = A%BEND(K) - A%BSTART(K) + 1
        PROW_K = MOD(K-1, GRID%PR)
        PCOL_K = BLR_COL_OWNER(A%CMAP, K, GRID)
        OWNER = (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K)
        IF (GRID%MYCOL.EQ.PCOL_K .AND. K.LT.NBLK) THEN
          ALLOCATE(S(KW,NRHS))
          S = CZERO
          lj = BLR_COL_LOCAL(A%CMAP, K, GRID)
          RS = 1
          DO li = 1, A%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MK = A%BEND(I) - A%BSTART(I) + 1
            IF (I .GT. K) THEN
              IDX = OOC_CACHE_GET(A%CACHE_A, A%STORE_A, A%OFFBLK, li, lj)
              CALL BLR_BLOCK_APPLY_T(A%CACHE_A%SLOTS(IDX)%BLK, NRHS, CONE, XFULL(RS,1), LDB, CONE, S, KW)
            END IF
            RS = RS + MK
          END DO
          IF (GRID%MYROW .EQ. PROW_K) THEN
            ALLOCATE(SSUM(KW,NRHS))
          ELSE
            ALLOCATE(SSUM(1,1))
          END IF
          CALL MPI_REDUCE(S, SSUM, KW*NRHS, MPI_COMPLEX, MPI_SUM, PROW_K, GRID%COL_COMM, IERR)
          DEALLOCATE(S)
        END IF
        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(XK_FULL(KW,NRHS))
          RS = LOCAL_ROW_OFFSET_SYM_OOC(A, GRID, K)
          IF (OWNER) THEN
            XK_FULL = XFULL(RS:RS+KW-1,1:NRHS)
            IF (K .LT. NBLK) XK_FULL = XK_FULL - SSUM
            ALLOCATE(DKK(KW,KW))
            DKK = CZERO
            IDX = OOC_DIAG_CACHE_GET(A%CACHE_D, A%STORE_A, A%DIAGBLK, K)
            CALL SYM_UNPACK_LOWER(KW, A%CACHE_D%SLOTS(IDX)%BLK%DENSE, DKK, KW)
            CALL CTRSM('Left', 'Lower', 'Transpose', 'Unit', KW, NRHS, CONE, DKK, KW, XK_FULL, KW)
            DEALLOCATE(DKK)
          END IF
          CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
          XFULL(RS:RS+KW-1,1:NRHS) = XK_FULL
          DEALLOCATE(XK_FULL)
        END IF
        IF (ALLOCATED(SSUM)) DEALLOCATE(SSUM)
        STEP = STEP + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt ooc', STEP, NSTEPS, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
        END IF
      END DO
      ! part locale (colonnes RSTARTS(MYCOL)...) de la solution pleine largeur
      IF (MLOC.GT.0 .AND. NRHS_LOC.GT.0) &
        B_LOCAL(1:MLOC,1:NRHS_LOC) = XFULL(1:MLOC, RSTARTS(GRID%MYCOL):RSTARTS(GRID%MYCOL)+NRHS_LOC-1)
      DEALLOCATE(XFULL)

      ! 4) Q^T
      DO K = 1, NBLK
        KW = A%BEND(K) - A%BSTART(K) + 1
        IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
        RS = LOCAL_ROW_OFFSET_SYM_OOC(A, GRID, K)
        DO I = KW, 1, -1
          P = ABS(A%IPIV(I,K))
          IF (P .NE. I) CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+P-1,1), LDB)
        END DO
      END DO

      CALL CONSUME_X_BATCH(RSTART, REND, GSTART, GEND, B_LOCAL, LDB, NRHS_LOC)
      DEALLOCATE(B_LOCAL)
    END DO
    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
    ! Lecture seule : rien de DIRTY, le flush ne fait que liberer la RAM.
    CALL OOC_CACHE_FLUSH_ALL(A%CACHE_A, A%STORE_A, A%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A%CACHE_D, A%STORE_A, A%DIAGBLK)
  END SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_OOC_2D

END MODULE MPI_BLR_LDLT_OOC_MOD
