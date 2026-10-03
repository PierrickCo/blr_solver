!===============================================================================
! Solveur LDL^T MPI+BLR pour matrices COMPLEXES SYMETRIQUES (A^T = A, PAS
! hermitiennes : transposee simple partout, jamais de conjugaison) --
! calque sur le LU BLR+MPI avec pivot confine (mpi_blr_matrix_mod.f90,
! MPI_BLR_MATRIX_LU_FACTOR_PIVOT / MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D).
!
! STOCKAGE : seulement le triangle INFERIEUR par blocs -- blocs diagonaux
! (denses, DIAGBLK) + blocs hors-diagonale I>J (OFFBLK, dense ou bas-rang).
! Les blocs I<J ne sont jamais assembles ni alloues : ~2x moins de memoire
! que le LU pour la meme matrice. Meme distribution bloc-cyclique 2D que
! MPI_BLR_MATRIX_T (BLOCK_NUMROC, bloc-ligne local li <-> global
! I=(li-1)*PR+MYROW+1).
!
! PIVOT : Bunch-Kaufman borne (CSYTRF_RK, pivots 1x1 et 2x2) CONFINE au
! bloc diagonal -- meme philosophie que le pivot partiel confine du LU.
! Le bloc diagonal K est factorise EN ENTIER, une fois toutes ses mises a
! jour de Schur recues (right-looking au niveau bloc) : aucune mise a jour
! differee ne traverse la compression (cause de l'abandon de la premiere
! tentative LDLt+BLR, voir legacy/src/blr_ldlt_mod.f90).
!   CSYTRF_RK donne P^T*A_KK*P = L_KK*D_K*L_KK^T, ou P^T = sequence
!   d'echanges (k,|IPIV(k)|) appliques dans l'ordre k=1..KW (convention
!   VERIFIEE numeriquement, ldlt/test/test_ldlt_kernels.f90), D_K diagonal
!   par blocs 1x1/2x2 (IPIV(k)<0 : bloc 2x2 en (k,k+1), terme hors-
!   diagonal dans E(k)).
!   La permutation SYMETRIQUE du bloc K se propage sans decompression :
!     - lignes des L(K,J), J<K, deja calcules  -> echange de lignes de UFAC
!     - colonnes des A(I,K), I>K, pas encore traites -> echange de colonnes
!       de VFAC
!   (les blocs A(K,J>K) de la ligne K ne sont pas stockes : ce sont les
!   transposes des A(J,K), deja couverts par le second point).
!   Factorisation globale obtenue : Q*A*Q^T = L*D*L^T, Q = diag(P_K^T).
!
! ETAPE K (factorisation) :
!   1) proprietaire de (K,K) : CSYTRF_RK ; IPIV_K + D_K diffuses a TOUTE la
!      grille (petits, reutilises a la resolution).
!   2) ligne PROW_K : echange de lignes des L(K,J<K).
!   3) colonne PCOL_K : echange de colonnes puis W(I,K) = A(I,K)*P_K*
!      L_KK^-T (CTRSM sur VFAC seulement si bas-rang).
!   4) W diffuse le long des lignes de grille (deja compresse) ; le
!      proprietaire stocke L(I,K) = W(I,K)*D_K^-1.
!   5) L(J,K)^T diffuse le long des COLONNES de grille depuis le processus
!      (ligne(J), colonne(J)) -- qui a recu W(J,K) a l'etape 4 (schema
!      classique du Cholesky ScaLAPACK).
!   6) Schur, triangle inferieur seulement : A(I,J) -= W(I,K)*L(J,K)^T,
!      I>=J>K (bloc diagonal I=J mis a jour en entier, dense).
!
! RESOLUTION (A*X=B, B distribue 2D comme le LU) :
!   Q*B -> L*Z=Q*B (descente, identique au LU) -> D*W=Z (local, aucune
!   communication) -> L^T*Y=W (remontee "left-looking" : une reduction le
!   long d'une colonne de grille + une diffusion par etape) -> X=Q^T*Y.
!   La remontee garde une copie PLEINE LARGEUR (NRHS colonnes) des
!   inconnues deja calculees de SES blocs-lignes (XFULL, MLOC x NRHS) :
!   cout memoire xPC sur le second membre pendant cette phase -- borne en
!   pratique par le traitement par lots de colonnes (--nrhs-batch).
!===============================================================================
MODULE MPI_BLR_LDLT_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE, DENSE_PANEL_T
  USE MPI_BLR_MATRIX_MOD, ONLY: BLOCK_NUMROC, STAGE_PROGRESS_CALLBACK_IFACE, STAGE_PROGRESS_SHOULD_FIRE, &
                                RHS_LOCAL_RANGE_MBM, RHS_NUMROC_MBM
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_BAR
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE

  TYPE :: MPI_BLR_SYM_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)           ! decoupage GLOBAL, partage
    TYPE(DENSE_PANEL_T), ALLOCATABLE :: DIAGBLK(:)       ! bloc GLOBAL K -- alloue seulement chez le proprietaire
    TYPE(BLR_BLOCK_T), ALLOCATABLE   :: OFFBLK(:,:)      ! (li,lj) LOCAUX -- rempli seulement si I>J
    ! Remplis par MPI_BLR_SYM_LDLT_FACTOR, diffuses a TOUTE la grille :
    INTEGER, ALLOCATABLE     :: IPIV(:,:)                ! (NB,NBLK) IPIV brut de CSYTRF_RK, bloc K
    COMPLEX(WP), ALLOCATABLE :: DDIAG(:,:), DSUB(:,:)    ! (NB,NBLK) diagonale / sous-diagonale (2x2) de D_K
  END TYPE MPI_BLR_SYM_T

CONTAINS

  !=======================================================================
  ! Primitives de bloc (dense ou bas-rang) propres au LDL^T
  !=======================================================================

  ! BLK <- BLK*P, colonnes echangees (i,SWP(i)) dans l'ordre i=1..KW
  ! (BLK%N = KW). Bas-rang : seules les colonnes de VFAC bougent, O(rang).
  SUBROUTINE BLR_BLOCK_COL_SWAP(BLK, KW, SWP)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN) :: KW, SWP(KW)
    INTEGER :: I
    IF (BLK%M .EQ. 0 .OR. BLK%N .EQ. 0) RETURN
    IF (BLK%IS_LOWRANK) THEN
      IF (BLK%RANK .EQ. 0) RETURN
      DO I = 1, KW
        IF (SWP(I) .NE. I) CALL CSWAP(BLK%RANK, BLK%VFAC(1,I), 1, BLK%VFAC(1,SWP(I)), 1)
      END DO
    ELSE
      DO I = 1, KW
        IF (SWP(I) .NE. I) CALL CSWAP(BLK%M, BLK%DENSE(1,I), 1, BLK%DENSE(1,SWP(I)), 1)
      END DO
    END IF
  END SUBROUTINE BLR_BLOCK_COL_SWAP

  ! BLK <- BLK * L^-T, L (KWxKW, DKK) triangulaire inferieure unitaire.
  ! Bas-rang : (U*V)*L^-T = U*(V*L^-T), CTRSM sur VFAC (RANGxKW) seulement.
  SUBROUTINE BLR_TRSM_RIGHT_LOWER_TRANS_UNIT(BLK, DKK, KW)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN)              :: KW
    COMPLEX(WP), INTENT(IN)          :: DKK(KW,KW)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%M .EQ. 0) RETURN
    IF (BLK%IS_LOWRANK) THEN
      IF (BLK%RANK .GT. 0) CALL CTRSM('Right', 'Lower', 'Transpose', 'Unit', BLK%RANK, KW, CONE, DKK, KW, BLK%VFAC, BLK%RANK)
    ELSE
      CALL CTRSM('Right', 'Lower', 'Transpose', 'Unit', BLK%M, KW, CONE, DKK, KW, BLK%DENSE, BLK%M)
    END IF
  END SUBROUTINE BLR_TRSM_RIGHT_LOWER_TRANS_UNIT

  ! X(1:NR,1:KW) <- X * D^-1, D bloc-diagonale 1x1/2x2 (format CSYTRF_RK).
  SUBROUTINE APPLY_DINV_RIGHT(X, LDX, NR, KW, IPIV, DD, DS)
    INTEGER, INTENT(IN)        :: LDX, NR, KW, IPIV(KW)
    COMPLEX(WP), INTENT(INOUT) :: X(LDX,KW)
    COMPLEX(WP), INTENT(IN)    :: DD(KW), DS(KW)
    COMPLEX(WP) :: A11, A21, A22, DET, X1, X2
    INTEGER :: K, R
    IF (NR .LE. 0) RETURN
    K = 1
    DO WHILE (K .LE. KW)
      IF (IPIV(K) .LT. 0 .AND. K .LT. KW) THEN
        A11 = DD(K); A21 = DS(K); A22 = DD(K+1)
        DET = A11*A22 - A21*A21
        DO R = 1, NR
          X1 = X(R,K); X2 = X(R,K+1)
          X(R,K)   = (X1*A22 - X2*A21) / DET
          X(R,K+1) = (X2*A11 - X1*A21) / DET
        END DO
        K = K + 2
      ELSE
        X(1:NR,K) = X(1:NR,K) / DD(K)
        K = K + 1
      END IF
    END DO
  END SUBROUTINE APPLY_DINV_RIGHT

  ! Y(1:KW,1:NC) <- D^-1 * Y (D symetrique : meme inverse 2x2 qu'a droite).
  SUBROUTINE APPLY_DINV_LEFT(Y, LDY, KW, NC, IPIV, DD, DS)
    INTEGER, INTENT(IN)        :: LDY, KW, NC, IPIV(KW)
    COMPLEX(WP), INTENT(INOUT) :: Y(LDY,*)
    COMPLEX(WP), INTENT(IN)    :: DD(KW), DS(KW)
    COMPLEX(WP) :: A11, A21, A22, DET, Y1, Y2
    INTEGER :: K, C
    IF (NC .LE. 0) RETURN
    K = 1
    DO WHILE (K .LE. KW)
      IF (IPIV(K) .LT. 0 .AND. K .LT. KW) THEN
        A11 = DD(K); A21 = DS(K); A22 = DD(K+1)
        DET = A11*A22 - A21*A21
        DO C = 1, NC
          Y1 = Y(K,C); Y2 = Y(K+1,C)
          Y(K,C)   = (A22*Y1 - A21*Y2) / DET
          Y(K+1,C) = (A11*Y2 - A21*Y1) / DET
        END DO
        K = K + 2
      ELSE
        Y(K,1:NC) = Y(K,1:NC) / DD(K)
        K = K + 1
      END IF
    END DO
  END SUBROUTINE APPLY_DINV_LEFT

  ! BLK <- BLK * D^-1 (bas-rang : sur VFAC seulement).
  SUBROUTINE BLR_BLOCK_DINV_RIGHT(BLK, KW, IPIV, DD, DS)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN)              :: KW, IPIV(KW)
    COMPLEX(WP), INTENT(IN)          :: DD(KW), DS(KW)
    IF (BLK%M .EQ. 0) RETURN
    IF (BLK%IS_LOWRANK) THEN
      IF (BLK%RANK .GT. 0) CALL APPLY_DINV_RIGHT(BLK%VFAC, BLK%RANK, BLK%RANK, KW, IPIV, DD, DS)
    ELSE
      CALL APPLY_DINV_RIGHT(BLK%DENSE, BLK%M, BLK%M, KW, IPIV, DD, DS)
    END IF
  END SUBROUTINE BLR_BLOCK_DINV_RIGHT

  ! BOUT <- BIN^T (transposee SIMPLE). Bas-rang : (U*V)^T = V^T*U^T.
  SUBROUTINE BLR_BLOCK_TRANSPOSE(BIN, BOUT)
    TYPE(BLR_BLOCK_T), INTENT(IN)  :: BIN
    TYPE(BLR_BLOCK_T), INTENT(OUT) :: BOUT
    BOUT%M = BIN%N; BOUT%N = BIN%M
    BOUT%IS_LOWRANK = BIN%IS_LOWRANK
    BOUT%RANK = BIN%RANK
    IF (BIN%IS_LOWRANK) THEN
      ALLOCATE(BOUT%UFAC(BIN%N,BIN%RANK), BOUT%VFAC(BIN%RANK,BIN%M))
      IF (BIN%RANK .GT. 0) THEN
        BOUT%UFAC = TRANSPOSE(BIN%VFAC)
        BOUT%VFAC = TRANSPOSE(BIN%UFAC)
      END IF
    ELSE
      ALLOCATE(BOUT%DENSE(BIN%N,BIN%M))
      BOUT%DENSE = TRANSPOSE(BIN%DENSE)
    END IF
  END SUBROUTINE BLR_BLOCK_TRANSPOSE

  ! Y(BLK%N x NRHS) <- ALPHA*BLK^T*X + BETA*Y, X de taille BLK%M x NRHS
  ! (transposee SIMPLE). Bas-rang : V^T*(U^T*X), jamais d'expansion.
  SUBROUTINE BLR_BLOCK_APPLY_T(BLK, NRHS, ALPHA, X, LDX, BETA, Y, LDY)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: BLK
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
  END SUBROUTINE BLR_BLOCK_APPLY_T

  ! Nombre d'echanges effectifs et de blocs 2x2 dans un IPIV CSYTRF_RK.
  SUBROUTINE COUNT_PIVOTS(KW, IPIV, NSWAP, N2X2)
    INTEGER, INTENT(IN)  :: KW, IPIV(KW)
    INTEGER, INTENT(OUT) :: NSWAP, N2X2
    INTEGER :: K
    NSWAP = 0; N2X2 = 0
    DO K = 1, KW
      IF (ABS(IPIV(K)) .NE. K) NSWAP = NSWAP + 1
    END DO
    K = 1
    DO WHILE (K .LE. KW)
      IF (IPIV(K) .LT. 0 .AND. K .LT. KW) THEN
        N2X2 = N2X2 + 1
        K = K + 2
      ELSE
        K = K + 1
      END IF
    END DO
  END SUBROUTINE COUNT_PIVOTS

  ! Decalage (1-base) de la premiere ligne du bloc-ligne GLOBAL K dans le
  ! stockage LOCAL des lignes (B_LOCAL) -- K doit appartenir a MYROW.
  FUNCTION LOCAL_ROW_OFFSET_SYM(A, GRID, K) RESULT(RS)
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: K
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (K-1)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + A%BEND(I) - A%BSTART(I) + 1
    END DO
  END FUNCTION LOCAL_ROW_OFFSET_SYM

  ! Nombre de lignes locales (somme des blocs-lignes possedes).
  FUNCTION MPI_BLR_SYM_MLOC(A, GRID) RESULT(MLOC)
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER :: MLOC, li, I
    MLOC = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + A%BEND(I) - A%BSTART(I) + 1
    END DO
  END FUNCTION MPI_BLR_SYM_MLOC

  !=======================================================================
  ! Assemblage : triangle inferieur seulement (I>=J), compression immediate
  ! de chaque bloc hors-diagonale (CGESVD ou SVD randomisee adaptative).
  ! GET_BLOCK n'est JAMAIS appele pour I<J -- l'appelant n'a a fournir que
  ! le triangle inferieur (la symetrie est supposee, pas verifiee).
  ! Memes options/diagnostics que MPI_ASSEMBLE_BLR_MATRIX.
  !=======================================================================
  SUBROUTINE MPI_ASSEMBLE_BLR_SYM(GRID, N, NB, TOL, GET_BLOCK, A, T_ASSEMBLE_COMPRESS, &
                                  N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS, &
                                  T_GET_BLOCK, T_COMPRESS, USE_RSVD, RSVD_OVERSAMPLE, RSVD_NPOWER, &
                                  PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)          :: GRID
    INTEGER, INTENT(IN)               :: N, NB
    REAL(WP), INTENT(IN)              :: TOL
    PROCEDURE(GET_BLOCK_IFACE)        :: GET_BLOCK
    TYPE(MPI_BLR_SYM_T), INTENT(OUT)  :: A
    DOUBLE PRECISION, INTENT(OUT)     :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)              :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)           :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL     :: SHOW_PROGRESS
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_GET_BLOCK, T_COMPRESS
    LOGICAL, INTENT(IN), OPTIONAL     :: USE_RSVD
    INTEGER, INTENT(IN), OPTIONAL     :: RSVD_OVERSAMPLE, RSVD_NPOWER
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, PBAR_STATE
    INTEGER :: RSVD_OVS, RSVD_NPOW
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0, T_GETBLK_L, T_COMPRESS_L, T_ASM_L
    LOGICAL :: PROGRESS, RSVD_FLAG

    RSVD_FLAG = .FALSE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD
    RSVD_OVS = 10; IF (PRESENT(RSVD_OVERSAMPLE)) RSVD_OVS = RSVD_OVERSAMPLE
    RSVD_NPOW = 1; IF (PRESENT(RSVD_NPOWER)) RSVD_NPOW = RSVD_NPOWER
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    PBAR_STATE = -1
    T_GETBLK_L = 0.0D0; T_COMPRESS_L = 0.0D0

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A%NBLK, A%BSTART, A%BEND)
    A%N = N
    A%NB = NB
    A%NBLKR_LOC = BLOCK_NUMROC(A%NBLK, GRID%PR, GRID%MYROW)
    A%NBLKC_LOC = BLOCK_NUMROC(A%NBLK, GRID%PC, GRID%MYCOL)
    ALLOCATE(A%DIAGBLK(A%NBLK))
    ALLOCATE(A%OFFBLK(MAX(A%NBLKR_LOC,1), MAX(A%NBLKC_LOC,1)))

    ! Barre de progression : part LOCALE du rang 0 (triangle inferieur).
    NTOTAL = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        IF (I .GE. J) NTOTAL = NTOTAL + 1
      END DO
    END DO
    NDONE = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A%BEND(I) - A%BSTART(I) + 1
      DO lj = 1, A%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        IF (I .LT. J) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A%DIAGBLK(I)%D(MI,MI))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), A%DIAGBLK(I)%D)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), TMP)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, A%OFFBLK(li,lj), RSVD_OVS, RSVD_NPOW)
          ELSE
            CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, A%OFFBLK(li,lj))
          END IF
          T_COMPRESS_L = T_COMPRESS_L + (MPI_WTIME() - TT0)
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
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('1) Assemblage+compression (rang 0)', NDONE, NTOTAL, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(NDONE, NTOTAL)) CALL PROGRESS_CALLBACK(NDONE, NTOTAL)
        END IF
      END DO
    END DO
    T_ASM_L = MPI_WTIME() - T0
    CALL MPI_ALLREDUCE(T_ASM_L, T_ASSEMBLE_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_GET_BLOCK)) CALL MPI_ALLREDUCE(T_GETBLK_L, T_GET_BLOCK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_COMPRESS))  CALL MPI_ALLREDUCE(T_COMPRESS_L, T_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_SYM

  ! Empreinte memoire LOCALE (octets) : blocs diagonaux + blocs I>J.
  FUNCTION MPI_BLR_SYM_MEMORY_BYTES(A) RESULT(NBYTES)
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    INTEGER(8) :: NBYTES
    INTEGER :: li, lj, K
    NBYTES = 0_8
    DO K = 1, A%NBLK
      IF (ALLOCATED(A%DIAGBLK(K)%D)) NBYTES = NBYTES + 8_8*INT(SIZE(A%DIAGBLK(K)%D),8)
    END DO
    DO lj = 1, A%NBLKC_LOC
      DO li = 1, A%NBLKR_LOC
        NBYTES = NBYTES + BLR_BLOCK_BYTES(A%OFFBLK(li,lj))
      END DO
    END DO
  END FUNCTION MPI_BLR_SYM_MEMORY_BYTES

  SUBROUTINE MPI_BLR_SYM_FREE(A)
    TYPE(MPI_BLR_SYM_T), INTENT(INOUT) :: A
    IF (ALLOCATED(A%DIAGBLK)) DEALLOCATE(A%DIAGBLK)
    IF (ALLOCATED(A%OFFBLK))  DEALLOCATE(A%OFFBLK)
    IF (ALLOCATED(A%BSTART))  DEALLOCATE(A%BSTART)
    IF (ALLOCATED(A%BEND))    DEALLOCATE(A%BEND)
    IF (ALLOCATED(A%IPIV))    DEALLOCATE(A%IPIV)
    IF (ALLOCATED(A%DDIAG))   DEALLOCATE(A%DDIAG)
    IF (ALLOCATED(A%DSUB))    DEALLOCATE(A%DSUB)
    A%N = 0; A%NB = 0; A%NBLK = 0; A%NBLKR_LOC = 0; A%NBLKC_LOC = 0
  END SUBROUTINE MPI_BLR_SYM_FREE

  !=======================================================================
  ! Factorisation LDL^T en place (voir le commentaire de tete pour les 6
  ! sous-etapes de chaque etape K). N_COMPRESSED/N_OFFDIAG/RANK_SUM
  ! caracterisent le facteur L final (blocs I>J, sommes sur la grille).
  ! N_SWAPS/N_2X2 (optionnels) : echanges effectifs / pivots 2x2 au total.
  ! Memes options que MPI_BLR_MATRIX_LU_FACTOR_PIVOT.
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR(GRID, A, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                     N_SWAPS, N_2X2, T_FACTOR, SHOW_PROGRESS, SHOW_RSS_TRACE, &
                                     USE_MOVE_ALLOC, PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)           :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(INOUT) :: A
    REAL(WP), INTENT(IN)               :: TOL
    INTEGER, INTENT(OUT)               :: INFO
    INTEGER, INTENT(OUT)               :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)            :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL     :: N_SWAPS, N_2X2
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS, SHOW_RSS_TRACE, USE_MOVE_ALLOC
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), DBUF(:), E(:), WORK(:)
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

    ! Espace de travail CSYTRF_RK dimensionne une fois pour NB (le plus
    ! grand bloc diagonal possible).
    ALLOCATE(E(A%NB))
    CALL CSYTRF_RK('L', A%NB, WQUERY, A%NB, E, A%IPIV(:,1), WQUERY, -1, IINFO)
    LWORK = MAX(1, INT(REAL(WQUERY(1))))
    ALLOCATE(WORK(LWORK))

    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = MOD(K-1, GRID%PC)
      ROOT_K = PROW_K*GRID%PC + PCOL_K
      OWNER = (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K)
      li_K = (K-1)/GRID%PR + 1
      lj_K = (K-1)/GRID%PC + 1
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW), SWP(KW), DBUF(2*KW))
      IPIV_K = 0; DBUF = (0.0_WP,0.0_WP)

      ! 1) Bunch-Kaufman borne sur le bloc diagonal, chez le proprietaire.
      IF (OWNER) THEN
        CALL CSYTRF_RK('L', KW, A%DIAGBLK(K)%D, KW, E, IPIV_K, WORK, LWORK, IINFO)
        IF (IINFO .NE. 0) INFO_STEP = A%BSTART(K) - 1 + IINFO
        DO I = 1, KW
          DBUF(I) = A%DIAGBLK(K)%D(I,I)
        END DO
        DBUF(KW+1:2*KW) = E(1:KW)
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(WRECV, LTRECV, IPIV_K, SWP, DBUF, E, WORK)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF
      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
      CALL MPI_BCAST(DBUF, 2*KW, MPI_COMPLEX, ROOT_K, GRID%COMM, IERR)
      A%IPIV(1:KW,K)  = IPIV_K
      A%DDIAG(1:KW,K) = DBUF(1:KW)
      A%DSUB(1:KW,K)  = DBUF(KW+1:2*KW)
      SWP = ABS(IPIV_K)
      CALL COUNT_PIVOTS(KW, IPIV_K, NS_K, N2_K)
      NS_TOT = NS_TOT + NS_K; N2_TOT = N2_TOT + N2_K

      ! 2) permutation des LIGNES des L(K,J<K) deja finalises (ligne PROW_K).
      IF (GRID%MYROW .EQ. PROW_K .AND. NS_K .GT. 0) THEN
        DO lj = 1, A%NBLKC_LOC
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .GE. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP(A%OFFBLK(li_K,lj), KW, SWP)
        END DO
      END IF

      ! 3) L_KK diffuse dans la colonne PCOL_K ; permutation des COLONNES
      !    des A(I>K,K) puis W = A(I,K)*P_K*L_KK^-T.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (OWNER) DKK = A%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          IF (NS_K .GT. 0) CALL BLR_BLOCK_COL_SWAP(A%OFFBLK(li,lj_K), KW, SWP)
          CALL BLR_TRSM_RIGHT_LOWER_TRANS_UNIT(A%OFFBLK(li,lj_K), DKK, KW)
        END DO
        DEALLOCATE(DKK)
      END IF

      ! 4) W(I,K) diffuse le long des lignes de grille (deja compresse) ;
      !    le proprietaire garde ensuite L(I,K) = W(I,K)*D_K^-1.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, A%OFFBLK(li,lj_K))
          WRECV(li) = A%OFFBLK(li,lj_K)
          CALL BLR_BLOCK_DINV_RIGHT(A%OFFBLK(li,lj_K), KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          NO_LOCAL = NO_LOCAL + 1
          IF (A%OFFBLK(li,lj_K)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A%OFFBLK(li,lj_K)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A%OFFBLK(li,lj_K)%M, A%OFFBLK(li,lj_K)%N)
          END IF
        END DO
      ELSE
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, WRECV(li))
        END DO
      END IF

      ! 5) L(J,K)^T diffuse le long des colonnes de grille, racine =
      !    processus (ligne(J), MYCOL), qui possede WRECV(J) depuis l'etape 4.
      !    Tous les membres d'un COL_COMM (meme MYCOL) parcourent la meme
      !    liste lj : condition uniforme sur le communicateur.
      DO lj = 1, A%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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

      ! 6) Schur, triangle inferieur : A(I,J) -= W(I,K)*L(J,K)^T, I>=J>K.
      DO lj = 1, A%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LT. J) CYCLE
          MI = A%BEND(I) - A%BSTART(I) + 1
          IF (I .EQ. J) THEN
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, WRECV(li), LTRECV(lj), A%DIAGBLK(I)%D, MI)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK(A%OFFBLK(li,lj), MI, NJ, KW, WRECV(li), LTRECV(lj), TOL, WS, &
                                        USE_MOVE_ALLOC=USE_MOVE_ALLOC)
          END IF
        END DO
      END DO

      DEALLOCATE(IPIV_K, SWP, DBUF)
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('2) Factorisation LDLt (Bunch-Kaufman confine)', K, NBLK, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(K, NBLK)) CALL PROGRESS_CALLBACK(K, NBLK)
      END IF
      IF (RSS_TRACE) RSS_MAX_MB = MAX(RSS_MAX_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
    END DO
    IF (RSS_TRACE) WRITE(6,'(A,F8.1,A)') '     [RSS trace] MAX observe pendant la factorisation = ', RSS_MAX_MB, ' Mo'

    DEALLOCATE(WRECV, LTRECV, E, WORK)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) THEN
      T_FACTOR_L = MPI_WTIME() - T0
      CALL MPI_ALLREDUCE(T_FACTOR_L, T_FACTOR, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    ! NS_TOT/N2_TOT deja identiques partout (IPIV_K diffuse a toute la grille).
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_TOT
    IF (PRESENT(N_2X2))   N_2X2 = N2_TOT
  END SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR

  !=======================================================================
  ! Resolution A*X=B apres MPI_BLR_SYM_LDLT_FACTOR. B_LOCAL (MLOC lignes
  ! locales x NRHS_LOC colonnes locales, decoupage RHS_LOCAL_RANGE_MBM) est
  ! ecrase par X. Generique sur la largeur NRHS (appelable par lot).
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, A, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS, PROGRESS_CALLBACK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)         :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN)  :: A
    INTEGER, INTENT(IN)              :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)       :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL    :: SHOW_PROGRESS
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE :: XFULL(:,:), S(:,:), SSUM(:,:), SLOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, MK, RS, RE, NROWTRAIL, C, IERR, P, MLOC
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, PBAR_STATE
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    STEP = 0
    PBAR_STATE = -1
    NBLK = A%NBLK
    NSTEPS = 2*NBLK
    MLOC = MPI_BLR_SYM_MLOC(A, GRID)
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    ! === 0) B <- Q*B : echanges de chaque bloc, dans l'ordre, sur le B
    !     ORIGINAL (avant toute soustraction croisee -- meme piege que le LU).
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
      DO I = 1, KW
        P = ABS(A%IPIV(I,K))
        IF (P .NE. I) CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+P-1,1), LDB)
      END DO
    END DO

    ! === 1) L*Z = Q*B (descente) -- meme algorithme que le LU pivote ===
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = MOD(K-1, GRID%PC)
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
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

      lj = (K-1)/GRID%PC + 1
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
            CALL BLR_BLOCK_APPLY(A%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
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
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO

    ! === 2) D*W = Z : local a chaque ligne de grille, aucune communication ===
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
      CALL APPLY_DINV_LEFT(B_LOCAL(RS,1), LDB, KW, NRHS_LOC, A%IPIV(1:KW,K), A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
    END DO

    ! === 3) L^T*Y = W (remontee "left-looking") ===
    !     Y_I = L_II^-T * (W_I - sum_{K>I} L(K,I)^T * Y_K). Les L(K,I) sont
    !     dans la colonne de grille PCOL_I : chaque processus de cette
    !     colonne somme ses termes locaux (Y_K pleine largeur depuis XFULL),
    !     MPI_REDUCE vers la ligne PROW_I, puis SCATTERV par colonnes de RHS.
    ALLOCATE(XFULL(MAX(MLOC,1),NRHS))
    XFULL = CZERO
    DO K = NBLK, 1, -1
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = MOD(K-1, GRID%PC)
      IF (K .LT. NBLK) THEN
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          ALLOCATE(S(KW,NRHS))
          S = CZERO
          lj = (K-1)/GRID%PC + 1
          RS = 1
          DO li = 1, A%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MK = A%BEND(I) - A%BSTART(I) + 1
            IF (I .GT. K) CALL BLR_BLOCK_APPLY_T(A%OFFBLK(li,lj), NRHS, CONE, XFULL(RS,1), MAX(MLOC,1), CONE, S, KW)
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
          IF (.NOT. ALLOCATED(SSUM)) ALLOCATE(SSUM(1,1))
          ALLOCATE(SLOC(KW,MAX(NRHS_LOC,1)))
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
          END DO
          CALL MPI_SCATTERV(SSUM, RCOUNTS, RDISPLS, MPI_COMPLEX, SLOC, KW*NRHS_LOC, MPI_COMPLEX, &
                            PCOL_K, GRID%ROW_COMM, IERR)
          IF (NRHS_LOC .GT. 0) THEN
            RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
            B_LOCAL(RS:RS+KW-1,1:NRHS_LOC) = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC) - SLOC(:,1:NRHS_LOC)
          END IF
          DEALLOCATE(SLOC)
        END IF
        IF (ALLOCATED(SSUM)) DEALLOCATE(SSUM)
      END IF

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
        IF (NRHS_LOC .GT. 0) CALL CTRSM('Left', 'Lower', 'Transpose', 'Unit', KW, NRHS_LOC, &
                                         CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        DEALLOCATE(DKK)
        ! Y_K pleine largeur pour toute la ligne de grille (termes K>I des
        ! etapes suivantes).
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)), XK_FULL(KW,NRHS))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                            GRID%ROW_COMM, IERR)
        XFULL(RS:RS+KW-1,1:NRHS) = XK_FULL
        DEALLOCATE(XK_LOC, XK_FULL)
      END IF
      STEP = STEP + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO
    DEALLOCATE(XFULL)

    ! === 4) X = Q^T*Y : echanges de chaque bloc dans l'ordre INVERSE ===
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
      DO I = KW, 1, -1
        P = ABS(A%IPIV(I,K))
        IF (P .NE. I) CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+P-1,1), LDB)
      END DO
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_2D

END MODULE MPI_BLR_LDLT_MOD
