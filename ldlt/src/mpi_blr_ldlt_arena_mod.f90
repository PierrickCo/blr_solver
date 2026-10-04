!===============================================================================
! Variante ARENA du LDL^T BLR+MPI (mpi_blr_ldlt_mod.f90) -- meme algorithme
! numerique, mais les blocs hors-diagonale bas-rang sont carves dans un pool
! memoire dedie (BLR_ARENA_T, blr_arena_mod.f90) au lieu d'un ALLOCATE/
! DEALLOCATE individuel a chaque recompression de Schur : meme correctif de
! fragmentation de tas que MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA pour le LU.
!
! Ce qui change par rapport a la version classique :
!   - OFFBLK est un tableau de BLR_BLOCK_ARENA_T (vues pointeurs dans ARENA) ;
!   - mise a jour de Schur : SCHUR_UPDATE_BLR_BLOCK_ARENA ;
!   - operations de panneau (echange de colonnes, TRSM transpose, D^-1) :
!     faites sur une copie ORDINAIRE du bloc (deja necessaire pour la
!     diffusion), puis recopiees dans l'arena EN REUTILISANT le meme buffer
!     (COPY_REGULAR_TO_ARENA_REUSE -- ces operations ne changent pas le rang,
!     donc jamais de nouveau carve) ;
!   - echange de lignes des L(K,J<K) : BLR_BLOCK_ROW_SWAP_ARENA (en place).
! Blocs diagonaux, diffusions, resolution : identiques a la version
! classique (stockage compact DIAGP, memes primitives).
!
! Compromis (herite du LU arena, voir README principal) : une ancienne
! version de bloc remplacee par un rang plus grand n'est jamais rendue a
! l'arena -- l'empreinte analytique (ARENA_BYTES) peut depasser celle de la
! version classique ; le gain attendu est sur la RSS reelle a grande echelle.
!===============================================================================
MODULE MPI_BLR_LDLT_ARENA_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE BLR_ARENA_MOD
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE
  USE MPI_BLR_MATRIX_MOD, ONLY: BLOCK_NUMROC, STAGE_PROGRESS_CALLBACK_IFACE, STAGE_PROGRESS_SHOULD_FIRE, &
                                RHS_LOCAL_RANGE_MBM
  USE MPI_BLR_LDLT_MOD, ONLY: SYM_COL_OWNER, SYM_COL_LOCAL, SYM_COL_GLOBAL, SYM_NBLKC_LOC, &
                             PACKED_TRI_T, SYM_PACK_LOWER, SYM_UNPACK_LOWER, BLR_BLOCK_COL_SWAP, &
                              BLR_TRSM_RIGHT_LOWER_TRANS_UNIT, BLR_BLOCK_DINV_RIGHT, BLR_BLOCK_TRANSPOSE, &
                              APPLY_DINV_LEFT, COUNT_PIVOTS, &
                              SYM_BK_FACTOR, SYM_BK_LWORK, SYM_PACKED_SCHUR_UPDATE
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_BAR
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE

  TYPE :: MPI_BLR_SYM_ARENA_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (SYM_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(PACKED_TRI_T), ALLOCATABLE      :: DIAGP(:)
    TYPE(BLR_BLOCK_ARENA_T), ALLOCATABLE :: OFFBLK(:,:)
    INTEGER, ALLOCATABLE     :: IPIV(:,:)
    COMPLEX(WP), ALLOCATABLE :: DDIAG(:,:), DSUB(:,:)
    TYPE(BLR_ARENA_T) :: ARENA
  END TYPE MPI_BLR_SYM_ARENA_T

CONTAINS

  ! Recopie SRC (ordinaire) dans DST (arena) en REUTILISANT les buffers deja
  ! carves de DST quand leur capacite suffit (ARENA_CARVE_REUSE) -- DST est
  ! INTENT(INOUT), contrairement a COPY_REGULAR_TO_ARENA (assemblage).
  SUBROUTINE COPY_REGULAR_TO_ARENA_REUSE(SRC, DST, ARENA)
    TYPE(BLR_BLOCK_T), INTENT(IN)            :: SRC
    TYPE(BLR_BLOCK_ARENA_T), INTENT(INOUT)   :: DST
    TYPE(BLR_ARENA_T), INTENT(INOUT), TARGET :: ARENA
    COMPLEX(WP), POINTER :: FLAT(:)
    INTEGER :: M, N, R
    M = SRC%M; N = SRC%N; R = SRC%RANK
    DST%M = M; DST%N = N
    DST%IS_LOWRANK = SRC%IS_LOWRANK
    IF (SRC%IS_LOWRANK) THEN
      DST%RANK = R
      FLAT => ARENA_CARVE_REUSE(ARENA, DST%UFAC_FLAT, INT(M,8)*INT(R,8))
      DST%UFAC(1:M,1:R) => FLAT
      FLAT => ARENA_CARVE_REUSE(ARENA, DST%VFAC_FLAT, INT(R,8)*INT(N,8))
      DST%VFAC(1:R,1:N) => FLAT
      IF (R .GT. 0) THEN
        DST%UFAC = SRC%UFAC
        DST%VFAC = SRC%VFAC
      END IF
      NULLIFY(DST%DENSE)
    ELSE
      FLAT => ARENA_CARVE_REUSE(ARENA, DST%DENSE_FLAT, INT(M,8)*INT(N,8))
      DST%DENSE(1:M,1:N) => FLAT
      IF (M.GT.0 .AND. N.GT.0) DST%DENSE = SRC%DENSE
      NULLIFY(DST%UFAC, DST%VFAC)
    END IF
  END SUBROUTINE COPY_REGULAR_TO_ARENA_REUSE

  ! Y(BLK%N x NRHS) <- ALPHA*BLK^T*X + BETA*Y, bloc arena (transposee simple).
  SUBROUTINE BLR_BLOCK_APPLY_T_ARENA(BLK, NRHS, ALPHA, X, LDX, BETA, Y, LDY)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(IN) :: BLK
    INTEGER, INTENT(IN)           :: NRHS, LDX, LDY
    COMPLEX(WP), INTENT(IN)       :: ALPHA, BETA
    COMPLEX(WP), INTENT(IN)       :: X(LDX,NRHS)
    COMPLEX(WP), INTENT(INOUT)    :: Y(LDY,NRHS)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: T(:,:)
    IF (BLK%M.EQ.0 .OR. BLK%N.EQ.0 .OR. NRHS.EQ.0) RETURN
    IF (BLK%IS_LOWRANK) THEN
      IF (BLK%RANK .EQ. 0) RETURN
      ALLOCATE(T(BLK%RANK,NRHS))
      CALL CGEMM('Transpose', 'No transpose', BLK%RANK, NRHS, BLK%M, CONE, BLK%UFAC, BLK%M, X, LDX, CZERO, T, BLK%RANK)
      CALL CGEMM('Transpose', 'No transpose', BLK%N, NRHS, BLK%RANK, ALPHA, BLK%VFAC, BLK%RANK, T, BLK%RANK, BETA, Y, LDY)
      DEALLOCATE(T)
    ELSE
      CALL CGEMM('Transpose', 'No transpose', BLK%N, NRHS, BLK%M, ALPHA, BLK%DENSE, BLK%M, X, LDX, BETA, Y, LDY)
    END IF
  END SUBROUTINE BLR_BLOCK_APPLY_T_ARENA

  SUBROUTINE BCAST_DIAG_LOWER_ARENA(A, K, KW, OWNER, ROOT, COMM, DKK)
    INCLUDE 'mpif.h'
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(IN) :: A
    INTEGER, INTENT(IN)             :: K, KW, ROOT, COMM
    LOGICAL, INTENT(IN)             :: OWNER
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT) :: DKK(:,:)
    COMPLEX(WP), ALLOCATABLE :: APK(:)
    INTEGER :: NP, IERR
    NP = KW*(KW+1)/2
    ALLOCATE(APK(NP), DKK(KW,KW))
    IF (OWNER) APK = A%DIAGP(K)%AP
    CALL MPI_BCAST(APK, NP, MPI_COMPLEX, ROOT, COMM, IERR)
    DKK = (0.0_WP,0.0_WP)
    CALL SYM_UNPACK_LOWER(KW, APK, DKK, KW)
    DEALLOCATE(APK)
  END SUBROUTINE BCAST_DIAG_LOWER_ARENA

  FUNCTION LOCAL_ROW_OFFSET_SYM_ARENA(A, GRID, K) RESULT(RS)
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(IN) :: A
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: K
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (K-1)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + A%BEND(I) - A%BSTART(I) + 1
    END DO
  END FUNCTION LOCAL_ROW_OFFSET_SYM_ARENA

  FUNCTION MPI_BLR_SYM_ARENA_MLOC(A, GRID) RESULT(MLOC)
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(IN) :: A
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER :: MLOC, li, I
    MLOC = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + A%BEND(I) - A%BSTART(I) + 1
    END DO
  END FUNCTION MPI_BLR_SYM_ARENA_MLOC

  ! Empreinte analytique LOCALE : diagonale compacte + ARENA_BYTES (inclut
  ! l'espace des anciennes versions de blocs jamais recupere -- chiffre
  ! honnete, comme MPI_BLR_MATRIX_ARENA_MEMORY_BYTES).
  FUNCTION MPI_BLR_SYM_ARENA_MEMORY_BYTES(A) RESULT(NBYTES)
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(IN) :: A
    INTEGER(8) :: NBYTES
    INTEGER :: K
    NBYTES = 0_8
    DO K = 1, A%NBLK
      IF (ALLOCATED(A%DIAGP(K)%AP)) NBYTES = NBYTES + 8_8*INT(SIZE(A%DIAGP(K)%AP),8)
    END DO
    NBYTES = NBYTES + ARENA_BYTES(A%ARENA)
  END FUNCTION MPI_BLR_SYM_ARENA_MEMORY_BYTES

  SUBROUTINE MPI_BLR_SYM_ARENA_FREE(A)
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(INOUT) :: A
    IF (ALLOCATED(A%OFFBLK)) DEALLOCATE(A%OFFBLK)      ! vues seulement : la memoire est dans ARENA
    CALL ARENA_FREE(A%ARENA)
    IF (ALLOCATED(A%DIAGP))  DEALLOCATE(A%DIAGP)
    IF (ALLOCATED(A%BSTART)) DEALLOCATE(A%BSTART)
    IF (ALLOCATED(A%BEND))   DEALLOCATE(A%BEND)
    IF (ALLOCATED(A%IPIV))   DEALLOCATE(A%IPIV)
    IF (ALLOCATED(A%DDIAG))  DEALLOCATE(A%DDIAG)
    IF (ALLOCATED(A%DSUB))   DEALLOCATE(A%DSUB)
    A%N = 0; A%NB = 0; A%NBLK = 0; A%NBLKR_LOC = 0; A%NBLKC_LOC = 0
  END SUBROUTINE MPI_BLR_SYM_ARENA_FREE

  !=======================================================================
  ! Assemblage (triangle inferieur), blocs hors-diagonale carves dans ARENA.
  !=======================================================================
  SUBROUTINE MPI_ASSEMBLE_BLR_SYM_ARENA(GRID, N, NB, TOL, GET_BLOCK, A, T_ASSEMBLE_COMPRESS, &
                                        N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS, USE_RSVD, &
                                        PROGRESS_CALLBACK, DIAG_SPREAD)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)                 :: GRID
    INTEGER, INTENT(IN)                      :: N, NB
    REAL(WP), INTENT(IN)                     :: TOL
    PROCEDURE(GET_BLOCK_IFACE)               :: GET_BLOCK
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(OUT), TARGET :: A
    DOUBLE PRECISION, INTENT(OUT)            :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)                     :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)                  :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL            :: SHOW_PROGRESS, USE_RSVD
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir SYM_COL_OWNER) -- les blocs diagonaux visitent tous
    ! les processus ; .FALSE. : cyclique pur, comme le LU.
    LOGICAL, INTENT(IN), OPTIONAL            :: DIAG_SPREAD

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    TYPE(BLR_BLOCK_T) :: TMPBLK
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, PBAR_STATE
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, T_ASM_L
    LOGICAL :: PROGRESS, RSVD_FLAG

    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    PBAR_STATE = -1
    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A%NBLK, A%BSTART, A%BEND)
    A%N = N; A%NB = NB
    A%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A%NBLKR_LOC = BLOCK_NUMROC(A%NBLK, GRID%PR, GRID%MYROW)
    A%NBLKC_LOC = SYM_NBLKC_LOC(A%CMAP, A%NBLK, GRID)
    ALLOCATE(A%DIAGP(A%NBLK))
    ALLOCATE(A%OFFBLK(MAX(A%NBLKR_LOC,1), MAX(A%NBLKC_LOC,1)))
    NTOTAL = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = SYM_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .GE. J) NTOTAL = NTOTAL + 1
      END DO
    END DO
    NDONE = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A%BEND(I) - A%BSTART(I) + 1
      DO lj = 1, A%NBLKC_LOC
        J = SYM_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LT. J) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(TMP(MI,MI), A%DIAGP(I)%AP(MI*(MI+1)/2))
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), TMP)
          CALL SYM_PACK_LOWER(MI, TMP, MI, A%DIAGP(I)%AP)
          DEALLOCATE(TMP)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), TMP)
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, TMPBLK)
            CALL COPY_REGULAR_TO_ARENA(TMPBLK, A%OFFBLK(li,lj), A%ARENA)
          ELSE
            CALL COMPRESS_BLOCK_ARENA(MI, NJ, TMP, MI, TOL, A%OFFBLK(li,lj), A%ARENA)
          END IF
          NO_LOCAL = NO_LOCAL + 1
          IF (A%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(MI,NJ)
          END IF
          DEALLOCATE(TMP)
        END IF
        NDONE = NDONE + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('1) Assemblage+compression arena (rang 0)', NDONE, NTOTAL, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(NDONE, NTOTAL)) CALL PROGRESS_CALLBACK(NDONE, NTOTAL)
        END IF
      END DO
    END DO
    T_ASM_L = MPI_WTIME() - T0
    CALL MPI_ALLREDUCE(T_ASM_L, T_ASSEMBLE_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_SYM_ARENA

  !=======================================================================
  ! Factorisation LDL^T arena -- memes 6 sous-etapes que
  ! MPI_BLR_SYM_LDLT_FACTOR (voir son commentaire).
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR_ARENA(GRID, A, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                           N_SWAPS, N_2X2, T_FACTOR, SHOW_PROGRESS, SHOW_RSS_TRACE, &
                                           PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)                 :: GRID
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(INOUT), TARGET :: A
    REAL(WP), INTENT(IN)                     :: TOL
    INTEGER, INTENT(OUT)                     :: INFO, N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)                  :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL           :: N_SWAPS, N_2X2
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL  :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL            :: SHOW_PROGRESS, SHOW_RSS_TRACE
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), DBUF(:), E(:), WORK(:), DWORK(:,:)
    COMPLEX(WP) :: WQUERY(1)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: WRECV(:), LTRECV(:)
    TYPE(BLR_BLOCK_T) :: TMPB
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:), SWP(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP, LWORK, ROOT_K
    INTEGER :: PROW_K, PCOL_K, PROW_J, IERR, NC_LOCAL, NO_LOCAL, li_K, lj_K, li_J, PBAR_STATE
    INTEGER :: NS_TOT, N2_TOT, NS_K, N2_K
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, T_FACTOR_L, RSS_MAX_MB
    LOGICAL :: PROGRESS, RSS_TRACE, OWNER

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    RSS_TRACE = .FALSE.
    IF (PRESENT(SHOW_RSS_TRACE)) RSS_TRACE = SHOW_RSS_TRACE .AND. (GRID%MYID.EQ.0)
    PBAR_STATE = -1
    RSS_MAX_MB = 0.0D0
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
    ALLOCATE(E(A%NB))
    LWORK = SYM_BK_LWORK(A%NB)
    ALLOCATE(WORK(LWORK), DWORK(A%NB,A%NB))
    DWORK = (0.0_WP,0.0_WP)

    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = SYM_COL_OWNER(A%CMAP, K, GRID)
      ROOT_K = PROW_K*GRID%PC + PCOL_K
      OWNER = (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K)
      li_K = (K-1)/GRID%PR + 1
      lj_K = SYM_COL_LOCAL(A%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW), SWP(KW), DBUF(2*KW))
      IPIV_K = 0; DBUF = (0.0_WP,0.0_WP)

      ! 1) Bunch-Kaufman borne sur le bloc diagonal (deplie le temps du calcul)
      IF (OWNER) THEN
        CALL SYM_UNPACK_LOWER(KW, A%DIAGP(K)%AP, DWORK, A%NB)
        CALL SYM_BK_FACTOR(KW, DWORK, A%NB, E, IPIV_K, WORK, LWORK, IINFO)
        IF (IINFO .NE. 0) INFO_STEP = A%BSTART(K) - 1 + IINFO
        DO I = 1, KW
          DBUF(I) = DWORK(I,I)
        END DO
        DBUF(KW+1:2*KW) = E(1:KW)
        CALL SYM_PACK_LOWER(KW, DWORK, A%NB, A%DIAGP(K)%AP)
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(WRECV, LTRECV, IPIV_K, SWP, DBUF, E, WORK, DWORK)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF
      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
      CALL MPI_BCAST(DBUF, 2*KW, MPI_COMPLEX, ROOT_K, GRID%COMM, IERR)
      A%IPIV(1:KW,K)  = IPIV_K
      A%DDIAG(1:KW,K) = DBUF(1:KW)
      A%DSUB(1:KW,K)  = DBUF(KW+1:2*KW)
      SWP = ABS(IPIV_K)
      CALL COUNT_PIVOTS(KW, IPIV_K, DBUF(KW+1:2*KW), NS_K, N2_K)
      NS_TOT = NS_TOT + NS_K; N2_TOT = N2_TOT + N2_K

      ! 2) lignes des L(K,J<K) : echange en place dans l'arena (UFAC seul)
      IF (GRID%MYROW .EQ. PROW_K .AND. NS_K .GT. 0) THEN
        DO lj = 1, A%NBLKC_LOC
          J = SYM_COL_GLOBAL(A%CMAP, lj, GRID)
          IF (J .GE. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP_ARENA(A%OFFBLK(li_K,lj), KW, SWP)
        END DO
      END IF

      ! 3+4) colonne PCOL_K : copie ordinaire du bloc, echange de colonnes,
      !      W = A*P*L_KK^-T, diffusion de W le long de la ligne, puis
      !      L = W*D^-1 recopie dans l'arena (meme buffer, rang inchange).
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        CALL BCAST_DIAG_LOWER_ARENA(A, K, KW, OWNER, PROW_K, GRID%COL_COMM, DKK)
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL COPY_ARENA_TO_REGULAR(A%OFFBLK(li,lj_K), WRECV(li))
          IF (NS_K .GT. 0) CALL BLR_BLOCK_COL_SWAP(WRECV(li), KW, SWP)
          CALL BLR_TRSM_RIGHT_LOWER_TRANS_UNIT(WRECV(li), DKK, KW)
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, WRECV(li))
          TMPB = WRECV(li)
          CALL BLR_BLOCK_DINV_RIGHT(TMPB, KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          CALL COPY_REGULAR_TO_ARENA_REUSE(TMPB, A%OFFBLK(li,lj_K), A%ARENA)
          NO_LOCAL = NO_LOCAL + 1
          IF (TMPB%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + TMPB%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(TMPB%M, TMPB%N)
          END IF
        END DO
        DEALLOCATE(DKK)
      ELSE
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, WRECV(li))
        END DO
      END IF

      ! 5) L(J,K)^T le long des colonnes depuis (ligne(J), MYCOL)
      DO lj = 1, A%NBLKC_LOC
        J = SYM_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        PROW_J = MOD(J-1, GRID%PR)
        IF (GRID%MYROW .EQ. PROW_J) THEN
          li_J = (J-1)/GRID%PR + 1
          TMPB = WRECV(li_J)
          CALL BLR_BLOCK_DINV_RIGHT(TMPB, KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          CALL BLR_BLOCK_TRANSPOSE(TMPB, LTRECV(lj))
        END IF
        CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_J, KW, NJ, LTRECV(lj))
      END DO

      ! 6) Schur, triangle inferieur (blocs hors-diagonale carves dans l'arena)
      DO lj = 1, A%NBLKC_LOC
        J = SYM_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LT. J) CYCLE
          MI = A%BEND(I) - A%BSTART(I) + 1
          IF (I .EQ. J) THEN
            CALL SYM_PACKED_SCHUR_UPDATE(MI, A%DIAGP(I)%AP, KW, WRECV(li), LTRECV(lj), DWORK, A%NB)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK_ARENA(A%OFFBLK(li,lj), MI, NJ, KW, WRECV(li), LTRECV(lj), TOL, &
                                              A%ARENA, WS)
          END IF
        END DO
      END DO

      DEALLOCATE(IPIV_K, SWP, DBUF)
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('2) Factorisation LDLt arena (Bunch-Kaufman confine)', K, NBLK, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(K, NBLK)) CALL PROGRESS_CALLBACK(K, NBLK)
      END IF
      IF (RSS_TRACE) RSS_MAX_MB = MAX(RSS_MAX_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
    END DO
    IF (RSS_TRACE) WRITE(6,'(A,F8.1,A)') '     [RSS trace] MAX observe pendant la factorisation = ', RSS_MAX_MB, ' Mo'

    DEALLOCATE(WRECV, LTRECV, E, WORK, DWORK)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) THEN
      T_FACTOR_L = MPI_WTIME() - T0
      CALL MPI_ALLREDUCE(T_FACTOR_L, T_FACTOR, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_TOT
    IF (PRESENT(N_2X2))   N_2X2 = N2_TOT
  END SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR_ARENA

  !=======================================================================
  ! Resolution -- meme algorithme que MPI_BLR_SYM_LDLT_SOLVE_2D, blocs arena.
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_ARENA_2D(GRID, A, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS, PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)               :: GRID
    TYPE(MPI_BLR_SYM_ARENA_T), INTENT(IN)  :: A
    INTEGER, INTENT(IN)                    :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)             :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL          :: SHOW_PROGRESS
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE :: XFULL(:,:), S(:,:), SSUM(:,:), SLOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, MK, RS, RE, NROWTRAIL, C, IERR, P, MLOC
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, PBAR_STATE
    LOGICAL :: PROGRESS, OWNER

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    STEP = 0
    PBAR_STATE = -1
    NBLK = A%NBLK
    NSTEPS = 2*NBLK
    MLOC = MPI_BLR_SYM_ARENA_MLOC(A, GRID)
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    ! 0) Q*B
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_SYM_ARENA(A, GRID, K)
      DO I = 1, KW
        P = ABS(A%IPIV(I,K))
        IF (P .NE. I) CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+P-1,1), LDB)
      END DO
    END DO

    ! 1) descente L
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = SYM_COL_OWNER(A%CMAP, K, GRID)
      IF (GRID%MYROW .EQ. PROW_K) THEN
        CALL BCAST_DIAG_LOWER_ARENA(A, K, KW, GRID%MYCOL.EQ.PCOL_K, PCOL_K, GRID%ROW_COMM, DKK)
        RS = LOCAL_ROW_OFFSET_SYM_ARENA(A, GRID, K)
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
      lj = SYM_COL_LOCAL(A%CMAP, K, GRID)
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
            CALL BLR_BLOCK_APPLY_ARENA(A%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
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
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt arena', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO

    ! 2) D
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_SYM_ARENA(A, GRID, K)
      CALL APPLY_DINV_LEFT(B_LOCAL(RS,1), LDB, KW, NRHS_LOC, A%IPIV(1:KW,K), A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
    END DO

    ! 3) L^T*Y = W, remontee "pleine largeur" : chaque processus recupere
    !    une fois ses lignes de W sur TOUTE la largeur NRHS (XFULL, un
    !    ALLGATHERV par ligne de grille). Etape K : la colonne PCOL_K calcule
    !    sum_(I>K) L(I,K)^T*Y_I et la reduit sur (PROW_K,PCOL_K), proprietaire
    !    du bloc diagonal, qui resout L_KK^T et diffuse Y_K sur sa ligne.
    !    2 communications par etape, aucune diffusion du bloc diagonal.
    ALLOCATE(XFULL(MAX(MLOC,1),NRHS))
    XFULL = CZERO
    ALLOCATE(XK_LOC(MAX(MLOC,1),MAX(NRHS_LOC,1)))
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
      PCOL_K = SYM_COL_OWNER(A%CMAP, K, GRID)
      OWNER = (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K)
      IF (GRID%MYCOL.EQ.PCOL_K .AND. K.LT.NBLK) THEN
        ALLOCATE(S(KW,NRHS))
        S = CZERO
        lj = SYM_COL_LOCAL(A%CMAP, K, GRID)
        RS = 1
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          MK = A%BEND(I) - A%BSTART(I) + 1
          IF (I .GT. K) THEN
            CALL BLR_BLOCK_APPLY_T_ARENA(A%OFFBLK(li,lj), NRHS, CONE, XFULL(RS,1), MAX(MLOC,1), CONE, S, KW)
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
        RS = LOCAL_ROW_OFFSET_SYM_ARENA(A, GRID, K)
        IF (OWNER) THEN
          XK_FULL = XFULL(RS:RS+KW-1,1:NRHS)
          IF (K .LT. NBLK) XK_FULL = XK_FULL - SSUM
          ALLOCATE(DKK(KW,KW))
          DKK = CZERO
          CALL SYM_UNPACK_LOWER(KW, A%DIAGP(K)%AP, DKK, KW)
          CALL CTRSM('Left', 'Lower', 'Transpose', 'Unit', KW, NRHS, CONE, DKK, KW, XK_FULL, KW)
          DEALLOCATE(DKK)
        END IF
        CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        XFULL(RS:RS+KW-1,1:NRHS) = XK_FULL
        DEALLOCATE(XK_FULL)
      END IF
      IF (ALLOCATED(SSUM)) DEALLOCATE(SSUM)
      STEP = STEP + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt arena', STEP, NSTEPS, PBAR_STATE)
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
      RS = LOCAL_ROW_OFFSET_SYM_ARENA(A, GRID, K)
      DO I = KW, 1, -1
        P = ABS(A%IPIV(I,K))
        IF (P .NE. I) CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+P-1,1), LDB)
      END DO
    END DO
    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_ARENA_2D

END MODULE MPI_BLR_LDLT_ARENA_MOD
