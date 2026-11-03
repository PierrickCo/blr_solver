!===============================================================================
! Solveur LDL^T MPI+BLR pour matrices COMPLEXES SYMETRIQUES (A^T = A, PAS
! hermitiennes : transposee simple partout, jamais de conjugaison) --
! calque sur le LU BLR+MPI avec pivot confine (mpi_blr_matrix_mod.f90,
! MPI_BLR_MATRIX_LU_FACTOR_PIVOT / MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D).
!
! STOCKAGE : seulement le triangle INFERIEUR par blocs -- blocs diagonaux
! en STOCKAGE COMPACT (triangle inferieur seul, KW*(KW+1)/2 elements,
! format LAPACK 'L', DIAGP) + blocs hors-diagonale I>J (OFFBLK, dense ou
! bas-rang). Un bloc diagonal n'est deplie en tableau plein que le temps
! d'une mise a jour de Schur ou de CSYTRF_RK (un seul espace de travail
! NB x NB), et il est DIFFUSE compact (moitie du volume d'un bloc plein).
! Les blocs I<J ne sont jamais assembles ni alloues : ~2x moins de memoire
! que le LU pour la meme matrice. Meme distribution bloc-cyclique 2D que
! MPI_BLR_MATRIX_T pour les blocs-lignes (BLOCK_NUMROC, local li <->
! global I=(li-1)*PR+MYROW+1) ; blocs-colonnes par defaut en cyclique
! DECALE pour repartir les blocs diagonaux sur tous les processus
! (BLR_COL_OWNER / BLR_COL_LOCAL / BLR_COL_GLOBAL, champ CMAP).
!
! PIVOT : Bunch-Kaufman borne (CSYTRF_RK, pivots 1x1 et 2x2) CONFINE au
! bloc diagonal -- meme philosophie que le pivot partiel confine du LU.
! Le bloc diagonal K est factorise EN ENTIER, une fois toutes ses mises a
! jour de Schur recues (right-looking au niveau bloc) : aucune mise a jour
! differee ne traverse la compression (cause de l'abandon de la premiere
! tentative LDLt+BLR, voir legacy/src/blr_ldlt_mod.f90).
!   CSYTRF puis CSYCONVF (SYM_BK_FACTOR, meme resultat que CSYTRF_RK, plus
!   rapide) donnent P^T*A_KK*P = L_KK*D_K*L_KK^T, ou P^T = sequence
!   d'echanges (k,|IPIV(k)|) appliques dans l'ordre k=1..KW (convention
!   VERIFIEE numeriquement, ldlt/test/test_ldlt_kernels.f90), D_K diagonal
!   par blocs 1x1/2x2 (bloc 2x2 en (k,k+1) <=> terme hors-diagonal E(k)
!   non nul -- critere commun aux formats _RK et CSYCONVF).
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
!   communication) -> L^T*Y=W (remontee "pleine largeur" : une reduction le
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
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE
  USE MPI_BLR_MATRIX_MOD, ONLY: BLOCK_NUMROC, BLR_COL_OWNER, BLR_COL_LOCAL, BLR_COL_GLOBAL, BLR_NBLKC_LOC, &
                                 STAGE_PROGRESS_CALLBACK_IFACE, STAGE_PROGRESS_SHOULD_FIRE, &
                                RHS_LOCAL_RANGE_MBM, RHS_NUMROC_MBM
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_BAR
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB, BLR_TRIM_HEAP
  IMPLICIT NONE

  ! Triangle inferieur (diagonale comprise) d'un bloc diagonal KW x KW,
  ! stocke compact colonne par colonne (CTRTTP/CTPTTR, UPLO='L').
  TYPE :: PACKED_TRI_T
    COMPLEX(WP), ALLOCATABLE :: AP(:)
  END TYPE PACKED_TRI_T

  TYPE :: MPI_BLR_SYM_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (BLR_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)           ! decoupage GLOBAL, partage
    TYPE(PACKED_TRI_T), ALLOCATABLE  :: DIAGP(:)         ! bloc GLOBAL K, compact -- alloue seulement chez le proprietaire
    TYPE(BLR_BLOCK_T), ALLOCATABLE   :: OFFBLK(:,:)      ! (li,lj) LOCAUX -- rempli seulement si I>J
    ! Remplis par MPI_BLR_SYM_LDLT_FACTOR, diffuses a TOUTE la grille :
    INTEGER, ALLOCATABLE     :: IPIV(:,:)                ! (NB,NBLK) IPIV brut de CSYTRF_RK, bloc K
    COMPLEX(WP), ALLOCATABLE :: DDIAG(:,:), DSUB(:,:)    ! (NB,NBLK) diagonale / sous-diagonale (2x2) de D_K
  END TYPE MPI_BLR_SYM_T

CONTAINS

  ! Distribution des blocs-colonnes : BLR_COL_OWNER/LOCAL/GLOBAL et
  ! BLR_NBLKC_LOC (mpi_blr_matrix_mod.f90), partagees avec le LU.

  !=======================================================================
  ! Primitives de bloc (dense ou bas-rang) propres au LDL^T
  !=======================================================================

  ! Triangle inferieur de FULL(1:KW,1:KW) -> stockage compact AP.
  SUBROUTINE SYM_PACK_LOWER(KW, FULL, LDF, AP)
    INTEGER, INTENT(IN)      :: KW, LDF
    COMPLEX(WP), INTENT(IN)  :: FULL(LDF,*)
    COMPLEX(WP), INTENT(OUT) :: AP(*)
    INTEGER :: INFO
    IF (KW .GT. 0) CALL CTRTTP('L', KW, FULL, LDF, AP, INFO)
  END SUBROUTINE SYM_PACK_LOWER

  ! Stockage compact AP -> triangle inferieur de FULL(1:KW,1:KW) (le
  ! triangle superieur n'est pas touche : jamais lu par CSYTRF_RK('L'),
  ! CTRSM('Lower'), ni repacke).
  SUBROUTINE SYM_UNPACK_LOWER(KW, AP, FULL, LDF)
    INTEGER, INTENT(IN)        :: KW, LDF
    COMPLEX(WP), INTENT(IN)    :: AP(*)
    COMPLEX(WP), INTENT(INOUT) :: FULL(LDF,*)
    INTEGER :: INFO
    IF (KW .GT. 0) CALL CTPTTR('L', KW, AP, FULL, LDF, INFO)
  END SUBROUTINE SYM_UNPACK_LOWER

  ! Diffuse le bloc diagonal K COMPACT dans COMM (racine ROOT, OWNER =
  ! proprietaire) et le deplie dans DKK(KW,KW) chez tous les membres.
  SUBROUTINE BCAST_DIAG_LOWER(A, K, KW, OWNER, ROOT, COMM, DKK)
    INCLUDE 'mpif.h'
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
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
  END SUBROUTINE BCAST_DIAG_LOWER

  !=======================================================================
  ! Factorisation de Bunch-Kaufman d'un bloc diagonal : CSYTRF (classique)
  ! puis CSYCONVF('C'), qui convertit le resultat au format de CSYTRF_RK
  ! (L unitaire explicite, D diagonale + E sous-diagonale, IPIV au format
  ! _RK) -- toute la propagation de pivot reste inchangee. Mesure
  ! (test/bench_diag_kernels.f90) : 5 a 15 % plus rapide que CSYTRF_RK a
  ! 1 thread, 10 a 30 % a 2 threads.
  !=======================================================================
  SUBROUTINE SYM_BK_FACTOR(KW, A, LDA, E, IPIV, WORK, LWORK, INFO)
    INTEGER, INTENT(IN)        :: KW, LDA, LWORK
    COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
    COMPLEX(WP), INTENT(OUT)   :: E(*)
    INTEGER, INTENT(OUT)       :: IPIV(*), INFO
    COMPLEX(WP), INTENT(INOUT) :: WORK(*)
    INTEGER :: INFO2
    CALL CSYTRF('L', KW, A, LDA, IPIV, WORK, LWORK, INFO)
    CALL CSYCONVF('L', 'C', KW, A, LDA, E, IPIV, INFO2)
  END SUBROUTINE SYM_BK_FACTOR

  ! Taille d'espace de travail pour SYM_BK_FACTOR sur un bloc NB x NB.
  INTEGER FUNCTION SYM_BK_LWORK(NB)
    INTEGER, INTENT(IN) :: NB
    COMPLEX(WP) :: WQ(1), AQ(1)
    INTEGER :: IPQ(1), INFO
    CALL CSYTRF('L', NB, AQ, MAX(NB,1), IPQ, WQ, -1, INFO)
    SYM_BK_LWORK = MAX(1, INT(REAL(WQ(1))))
  END FUNCTION SYM_BK_LWORK

  ! AP(triangle inferieur compact M x M) -= X(M x R) * Y(R x M), triangle
  ! inferieur seul : une CGEMV par colonne, directement dans le stockage
  ! compact (aucun depliage).
  SUBROUTINE PACKED_LOWER_UPDATE(M, R, X, LDX, Y, LDY, AP)
    INTEGER, INTENT(IN)        :: M, R, LDX, LDY
    COMPLEX(WP), INTENT(IN)    :: X(LDX,*), Y(LDY,*)
    COMPLEX(WP), INTENT(INOUT) :: AP(*)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    INTEGER :: J, OFF, LEN
    OFF = 1
    DO J = 1, M
      LEN = M - J + 1
      CALL CGEMV('N', LEN, R, -CONE, X(J,1), LDX, Y(1,J), 1, CONE, AP(OFF), 1)
      OFF = OFF + LEN
    END DO
  END SUBROUTINE PACKED_LOWER_UPDATE

  !=======================================================================
  ! Mise a jour de Schur d'un bloc diagonal COMPACT : AP -= LBLK * UBLK
  ! (LBLK = W(J,K), M x KW ; UBLK = L(J,K)^T, KW x M), triangle inferieur
  ! seulement. Si l'un des deux facteurs est bas-rang, le produit s'ecrit
  ! X(M x R) * Y(R x M) avec R petit et la mise a jour se fait directement
  ! dans le stockage compact (PACKED_LOWER_UPDATE) ; sinon (deux blocs
  ! denses, ou R > RMAX_PACK) on deplie le bloc dans DWORK. Mesure a
  ! M=800, R=8 : 0.21 ms contre 0.77 ms pour depliage + CGEMM + repliage,
  ! resultat identique.
  !=======================================================================
  SUBROUTINE SYM_PACKED_SCHUR_UPDATE(M, AP, KW, LBLK, UBLK, DWORK, LDW)
    INTEGER, INTENT(IN)           :: M, KW, LDW
    COMPLEX(WP), INTENT(INOUT)    :: AP(*)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: LBLK, UBLK
    COMPLEX(WP), INTENT(INOUT)    :: DWORK(LDW,*)
    INTEGER, PARAMETER :: RMAX_PACK = 48
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: X(:,:), Y(:,:), T1(:,:)
    INTEGER :: RL, RU, R
    IF (M .EQ. 0 .OR. KW .EQ. 0) RETURN
    RL = KW; RU = KW
    IF (LBLK%IS_LOWRANK) RL = LBLK%RANK
    IF (UBLK%IS_LOWRANK) RU = UBLK%RANK
    IF (RL .EQ. 0 .OR. RU .EQ. 0) RETURN
    R = MIN(RL, RU)
    IF ((.NOT.LBLK%IS_LOWRANK .AND. .NOT.UBLK%IS_LOWRANK) .OR. R .GT. RMAX_PACK) THEN
      CALL SYM_UNPACK_LOWER(M, AP, DWORK, LDW)
      CALL SCHUR_UPDATE_BLOCK(M, M, KW, LBLK, UBLK, DWORK, LDW)
      CALL SYM_PACK_LOWER(M, DWORK, LDW, AP)
      RETURN
    END IF
    IF (LBLK%IS_LOWRANK .AND. UBLK%IS_LOWRANK) THEN
      ALLOCATE(T1(RL,RU), Y(RL,M))
      CALL CGEMM('N', 'N', RL, RU, KW, CONE, LBLK%VFAC, RL, UBLK%UFAC, KW, CZERO, T1, RL)
      CALL CGEMM('N', 'N', RL, M, RU, CONE, T1, RL, UBLK%VFAC, RU, CZERO, Y, RL)
      CALL PACKED_LOWER_UPDATE(M, RL, LBLK%UFAC, M, Y, RL, AP)
    ELSE IF (LBLK%IS_LOWRANK) THEN
      ALLOCATE(Y(RL,M))
      CALL CGEMM('N', 'N', RL, M, KW, CONE, LBLK%VFAC, RL, UBLK%DENSE, KW, CZERO, Y, RL)
      CALL PACKED_LOWER_UPDATE(M, RL, LBLK%UFAC, M, Y, RL, AP)
    ELSE
      ALLOCATE(X(M,RU))
      CALL CGEMM('N', 'N', M, RU, KW, CONE, LBLK%DENSE, M, UBLK%UFAC, KW, CZERO, X, M)
      CALL PACKED_LOWER_UPDATE(M, RU, X, M, UBLK%VFAC, RU, AP)
    END IF
  END SUBROUTINE SYM_PACKED_SCHUR_UPDATE

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
    BLK%ORTHO = .FALSE.
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
      IF (DS(K) .NE. (0.0_WP,0.0_WP) .AND. K .LT. KW) THEN      ! bloc 2x2 (voir COUNT_PIVOTS)
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
      IF (DS(K) .NE. (0.0_WP,0.0_WP) .AND. K .LT. KW) THEN      ! bloc 2x2 (voir COUNT_PIVOTS)
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
    BLK%ORTHO = .FALSE.
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

  ! Nombre d'echanges effectifs et de blocs 2x2. Un bloc 2x2 en (k,k+1) est
  ! repere par son terme sous-diagonal DS(k) = E(k) non nul : c'est le seul
  ! critere valable a la fois pour la sortie de CSYTRF_RK (IPIV(k) et
  ! IPIV(k+1) negatifs) et pour celle de CSYTRF + CSYCONVF (IPIV(k) = k,
  ! IPIV(k+1) < 0) -- verifie numeriquement. Les echanges, eux, s'appliquent
  ! de la meme facon dans les deux cas : (k,|IPIV(k)|) pour k = 1..KW.
  SUBROUTINE COUNT_PIVOTS(KW, IPIV, DS, NSWAP, N2X2)
    INTEGER, INTENT(IN)     :: KW, IPIV(KW)
    COMPLEX(WP), INTENT(IN) :: DS(KW)
    INTEGER, INTENT(OUT)    :: NSWAP, N2X2
    INTEGER :: K
    NSWAP = 0; N2X2 = 0
    DO K = 1, KW
      IF (ABS(IPIV(K)) .NE. K) NSWAP = NSWAP + 1
    END DO
    K = 1
    DO WHILE (K .LE. KW)
      IF (DS(K) .NE. (0.0_WP,0.0_WP) .AND. K .LT. KW) THEN
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
                                  PROGRESS_CALLBACK, DIAG_SPREAD)
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
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir BLR_COL_OWNER) -- les blocs diagonaux visitent tous
    ! les processus ; .FALSE. : cyclique pur, comme le LU.
    LOGICAL, INTENT(IN), OPTIONAL            :: DIAG_SPREAD

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, PBAR_STATE
    INTEGER :: RSVD_OVS, RSVD_NPOW
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0, T_GETBLK_L, T_COMPRESS_L, T_ASM_L
    LOGICAL :: PROGRESS, RSVD_FLAG

    RSVD_FLAG = .TRUE.     ! defaut : SVD randomisee adaptative (comme le LU)
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
    A%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A%NBLKR_LOC = BLOCK_NUMROC(A%NBLK, GRID%PR, GRID%MYROW)
    A%NBLKC_LOC = BLR_NBLKC_LOC(A%CMAP, A%NBLK, GRID)
    ALLOCATE(A%DIAGP(A%NBLK))
    ALLOCATE(A%OFFBLK(MAX(A%NBLKR_LOC,1), MAX(A%NBLKC_LOC,1)))

    ! Barre de progression : part LOCALE du rang 0 (triangle inferieur).
    NTOTAL = 0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .GE. J) NTOTAL = NTOTAL + 1
      END DO
    END DO
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
          ALLOCATE(TMP(MI,MI), A%DIAGP(I)%AP(MI*(MI+1)/2))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A%BSTART(I), A%BEND(I), A%BSTART(J), A%BEND(J), TMP)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
          CALL SYM_PACK_LOWER(MI, TMP, MI, A%DIAGP(I)%AP)
          DEALLOCATE(TMP)
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
      IF (ALLOCATED(A%DIAGP(K)%AP)) NBYTES = NBYTES + 8_8*INT(SIZE(A%DIAGP(K)%AP),8)
    END DO
    DO lj = 1, A%NBLKC_LOC
      DO li = 1, A%NBLKR_LOC
        NBYTES = NBYTES + BLR_BLOCK_BYTES(A%OFFBLK(li,lj))
      END DO
    END DO
  END FUNCTION MPI_BLR_SYM_MEMORY_BYTES

  SUBROUTINE MPI_BLR_SYM_FREE(A)
    TYPE(MPI_BLR_SYM_T), INTENT(INOUT) :: A
    IF (ALLOCATED(A%DIAGP))   DEALLOCATE(A%DIAGP)
    IF (ALLOCATED(A%OFFBLK))  DEALLOCATE(A%OFFBLK)
    IF (ALLOCATED(A%BSTART))  DEALLOCATE(A%BSTART)
    IF (ALLOCATED(A%BEND))    DEALLOCATE(A%BEND)
    IF (ALLOCATED(A%IPIV))    DEALLOCATE(A%IPIV)
    IF (ALLOCATED(A%DDIAG))   DEALLOCATE(A%DDIAG)
    IF (ALLOCATED(A%DSUB))    DEALLOCATE(A%DSUB)
    A%N = 0; A%NB = 0; A%NBLK = 0; A%NBLKR_LOC = 0; A%NBLKC_LOC = 0
  END SUBROUTINE MPI_BLR_SYM_FREE

  !=======================================================================
  ! Sauvegarde / reprise (checkpoint) de A : un fichier par processus
  ! (PREFIX_rankNNNNNN.dat), STREAM non formate, aucune communication sauf
  ! le verdict OK reduit sur la grille (un rang qui echoue arrete TOUS les
  ! rangs, jamais un sous-ensemble qui continuerait seul vers un appel
  ! collectif). WITH_FACTOR : ecrit/relit aussi IPIV/DDIAG/DSUB (apres
  ! factorisation). Meme conventions que les checkpoints de l'exemple LU.
  ! En-tete verifie a la relecture : N, NB, PR, PC, position dans la grille,
  ! rang, et la nature (matrice / facteur) du fichier ; le numero magique
  ! (7311 + CMAP) porte la distribution des blocs-colonnes, reprise telle quelle.
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_SAVE(GRID, A, PREFIX, WITH_FACTOR, OK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    CHARACTER(LEN=*), INTENT(IN)    :: PREFIX
    LOGICAL, INTENT(IN)             :: WITH_FACTOR
    LOGICAL, INTENT(OUT)            :: OK
    INTEGER :: U, li, lj, I, J, ND, IOS, IERR, KIND_TAG
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16)  :: RSTR
    CHARACTER(LEN=512) :: FNAME
    WRITE(RSTR,'(I6.6)') GRID%MYID
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    OK_LOCAL = (IOS .EQ. 0)
    IF (.NOT. OK_LOCAL) WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible de creer ', TRIM(FNAME), &
         ' (rang ', GRID%MYID, ') -- le repertoire de checkpoint doit deja exister.'
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U, STATUS='DELETE')
      RETURN
    END IF
    KIND_TAG = MERGE(2, 1, WITH_FACTOR)       ! 1 = matrice assemblee, 2 = facteur LDLt
    WRITE(U) 7311 + A%CMAP, KIND_TAG, A%N, A%NB, A%NBLK, GRID%PR, GRID%PC, GRID%MYROW, GRID%MYCOL, GRID%MYID, &
             A%NBLKR_LOC, A%NBLKC_LOC
    ND = 0
    DO I = 1, A%NBLK
      IF (ALLOCATED(A%DIAGP(I)%AP)) ND = ND + 1
    END DO
    WRITE(U) ND
    DO I = 1, A%NBLK
      IF (ALLOCATED(A%DIAGP(I)%AP)) THEN
        WRITE(U) I, SIZE(A%DIAGP(I)%AP)
        WRITE(U) A%DIAGP(I)%AP
      END IF
    END DO
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LE. J) CYCLE                    ! triangle inferieur strict seulement
        WRITE(U) I, J, A%OFFBLK(li,lj)%IS_LOWRANK, A%OFFBLK(li,lj)%M, A%OFFBLK(li,lj)%N, A%OFFBLK(li,lj)%RANK
        IF (A%OFFBLK(li,lj)%IS_LOWRANK) THEN
          IF (A%OFFBLK(li,lj)%RANK .GT. 0) THEN
            WRITE(U) A%OFFBLK(li,lj)%UFAC
            WRITE(U) A%OFFBLK(li,lj)%VFAC
          END IF
        ELSE
          WRITE(U) A%OFFBLK(li,lj)%DENSE
        END IF
      END DO
    END DO
    IF (WITH_FACTOR) WRITE(U) A%IPIV, A%DDIAG, A%DSUB
    CLOSE(U)
  END SUBROUTINE MPI_BLR_SYM_SAVE

  SUBROUTINE MPI_BLR_SYM_LOAD(GRID, A, N, NB, PREFIX, WITH_FACTOR, OK)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)         :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(OUT) :: A
    INTEGER, INTENT(IN)              :: N, NB
    CHARACTER(LEN=*), INTENT(IN)     :: PREFIX
    LOGICAL, INTENT(IN)              :: WITH_FACTOR
    LOGICAL, INTENT(OUT)             :: OK
    INTEGER :: U, li, lj, I, J, ND, K, NP, IOS, IERR
    INTEGER :: F_MAGIC, F_KIND, F_N, F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
    INTEGER :: R_I, R_J, R_M, R_N, R_RANK
    LOGICAL :: R_LOWRANK, OK_LOCAL
    CHARACTER(LEN=16)  :: RSTR
    CHARACTER(LEN=512) :: FNAME
    OK_LOCAL = .TRUE.
    WRITE(RSTR,'(I6.6)') GRID%MYID
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), ' (rang ', GRID%MYID, &
           ') -- sauvegarde absente ?'
      OK_LOCAL = .FALSE.
    ELSE
      READ(U, IOSTAT=IOS) F_MAGIC, F_KIND, F_N, F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
      IF (IOS.NE.0 .OR. (F_MAGIC.NE.7311 .AND. F_MAGIC.NE.7312) .OR. F_KIND.NE.MERGE(2,1,WITH_FACTOR) .OR. F_N.NE.N .OR. F_NB.NE.NB .OR. &
          F_PR.NE.GRID%PR .OR. F_PC.NE.GRID%PC .OR. F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. &
          F_MYID.NE.GRID%MYID) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), ' incompatible avec le run courant (rang ', &
             GRID%MYID, ') -- N/NB/grille ou nature (matrice/facteur LDLt) differente.'
        OK_LOCAL = .FALSE.
        CLOSE(U)
      END IF
    END IF
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U)
      RETURN
    END IF

    A%N = N; A%NB = NB
    A%CMAP = F_MAGIC - 7311                ! distribution des colonnes de la sauvegarde
    CALL BLOCK_PARTITION(N, NB, A%NBLK, A%BSTART, A%BEND)
    A%NBLKR_LOC = BLOCK_NUMROC(A%NBLK, GRID%PR, GRID%MYROW)
    A%NBLKC_LOC = BLR_NBLKC_LOC(A%CMAP, A%NBLK, GRID)
    ALLOCATE(A%DIAGP(A%NBLK))
    ALLOCATE(A%OFFBLK(MAX(A%NBLKR_LOC,1), MAX(A%NBLKC_LOC,1)))
    READ(U) ND
    DO K = 1, ND
      READ(U) I, NP
      ALLOCATE(A%DIAGP(I)%AP(NP))
      READ(U) A%DIAGP(I)%AP
    END DO
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LE. J) CYCLE
        READ(U) R_I, R_J, R_LOWRANK, R_M, R_N, R_RANK
        A%OFFBLK(li,lj)%IS_LOWRANK = R_LOWRANK
        A%OFFBLK(li,lj)%M = R_M; A%OFFBLK(li,lj)%N = R_N; A%OFFBLK(li,lj)%RANK = R_RANK
        IF (R_LOWRANK) THEN
          ALLOCATE(A%OFFBLK(li,lj)%UFAC(R_M,R_RANK), A%OFFBLK(li,lj)%VFAC(R_RANK,R_N))
          IF (R_RANK .GT. 0) THEN
            READ(U) A%OFFBLK(li,lj)%UFAC
            READ(U) A%OFFBLK(li,lj)%VFAC
          END IF
        ELSE
          ALLOCATE(A%OFFBLK(li,lj)%DENSE(R_M,R_N))
          READ(U) A%OFFBLK(li,lj)%DENSE
        END IF
      END DO
    END DO
    IF (WITH_FACTOR) THEN
      ALLOCATE(A%IPIV(NB,A%NBLK), A%DDIAG(NB,A%NBLK), A%DSUB(NB,A%NBLK))
      READ(U) A%IPIV, A%DDIAG, A%DSUB
    END IF
    CLOSE(U)
  END SUBROUTINE MPI_BLR_SYM_LOAD

  ! Statistiques de compression du triangle inferieur strict (apres une
  ! reprise, qui saute l'appel qui les produit normalement).
  SUBROUTINE MPI_BLR_SYM_COMPRESSION_STATS(GRID, A, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    INTEGER :: li, lj, I, J, NC, NO, IERR
    INTEGER(8) :: RS
    NC = 0; NO = 0; RS = 0_8
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LE. J) CYCLE
        NO = NO + 1
        IF (A%OFFBLK(li,lj)%IS_LOWRANK) THEN
          NC = NC + 1
          RS = RS + A%OFFBLK(li,lj)%RANK
        ELSE
          RS = RS + MIN(A%OFFBLK(li,lj)%M, A%OFFBLK(li,lj)%N)
        END IF
      END DO
    END DO
    CALL MPI_ALLREDUCE(NC, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_BLR_SYM_COMPRESSION_STATS

  !=======================================================================
  ! Factorisation LDL^T en place (voir le commentaire de tete pour les 6
  ! sous-etapes de chaque etape K). N_COMPRESSED/N_OFFDIAG/RANK_SUM
  ! caracterisent le facteur L final (blocs I>J, sommes sur la grille).
  ! N_SWAPS/N_2X2 (optionnels) : echanges effectifs / pivots 2x2 au total.
  ! Memes options que MPI_BLR_MATRIX_LU_FACTOR_PIVOT.
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR(GRID, A, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                     N_SWAPS, N_2X2, T_FACTOR, SHOW_PROGRESS, SHOW_RSS_TRACE, &
                                     USE_MOVE_ALLOC, PROGRESS_CALLBACK, T_BREAKDOWN, LOOKAHEAD, T_BREAKDOWN_MIN)
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
    ! T_BREAKDOWN (optionnel, MAX sur la grille, s) : (1) factorisation des
    ! blocs diagonaux, (2) echanges de lignes/colonnes du pivot, (3) operations
    ! de panneau (TRSM, D^-1, transposees), (4) communications MPI (attentes
    ! comprises), (5) Schur des blocs diagonaux, (6) Schur hors-diagonale,
    ! (7) boucle complete. Chaque categorie est un MAX independant : leur
    ! somme peut differer de (7).
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_BREAKDOWN(7)
    ! T_BREAKDOWN_MIN (optionnel) : memes categories, MIN sur la grille
    ! (avec T_BREAKDOWN : desequilibre entre processus).
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_BREAKDOWN_MIN(7)
    ! LOOKAHEAD (defaut .TRUE.) : anticipation de profondeur 1. A l'etape K,
    ! des que la colonne K+1 a recu sa mise a jour de Schur, le proprietaire
    ! du bloc (K+1,K+1) le factorise et lance des diffusions NON BLOQUANTES
    ! (MPI_IBCAST) de IPIV/D a toute la grille et du bloc compact a la
    ! colonne de grille ; le reste du Schur de l'etape K recouvre ces
    ! diffusions, et l'etape K+1 commence par un simple MPI_WAITALL. Resultat
    ! numerique identique (memes operations, seul l'ordre change).
    LOGICAL, INTENT(IN), OPTIONAL      :: LOOKAHEAD

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), DBUF(:), E(:), WORK(:), DWORK(:,:)
    ! tampons des diffusions anticipees (en vol entre deux etapes)
    INTEGER, ALLOCATABLE, ASYNCHRONOUS :: IBUF_N(:)
    COMPLEX(WP), ALLOCATABLE, ASYNCHRONOUS :: DBUF_N(:), APK_N(:)
    INTEGER :: REQ(3), NREQ
    LOGICAL :: LA, PENDING
    DOUBLE PRECISION :: TB(7), TT0, TLOOP0
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

    ! Espace de travail de SYM_BK_FACTOR dimensionne une fois pour NB (le
    ! plus grand bloc diagonal possible).
    ALLOCATE(E(A%NB))
    LWORK = SYM_BK_LWORK(A%NB)
    ALLOCATE(WORK(LWORK))
    ! Espace de travail unique ou un bloc diagonal compact est deplie le
    ! temps de CSYTRF_RK ou d'une mise a jour de Schur.
    ALLOCATE(DWORK(A%NB,A%NB))
    DWORK = (0.0_WP,0.0_WP)
    LA = .TRUE.
    IF (PRESENT(LOOKAHEAD)) LA = LOOKAHEAD
    PENDING = .FALSE.
    NREQ = 0
    IF (LA) THEN
      ALLOCATE(IBUF_N(0:A%NB), DBUF_N(2*A%NB), APK_N(A%NB*(A%NB+1)/2))
    END IF
    TB = 0.0D0
    TLOOP0 = MPI_WTIME()

    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A%CMAP, K, GRID)
      ROOT_K = PROW_K*GRID%PC + PCOL_K
      OWNER = (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K)
      li_K = (K-1)/GRID%PR + 1
      lj_K = BLR_COL_LOCAL(A%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW), SWP(KW), DBUF(2*KW))
      IPIV_K = 0; DBUF = (0.0_WP,0.0_WP)

      ! 1) Bunch-Kaufman borne sur le bloc diagonal, chez le proprietaire --
      !    ou, en anticipation, deja fait a l'etape K-1 : on attend alors
      !    seulement la fin des diffusions non bloquantes.
      IF (PENDING) THEN
        TT0 = MPI_WTIME()
        CALL MPI_WAITALL(NREQ, REQ, MPI_STATUSES_IGNORE, IERR)
        TB(4) = TB(4) + (MPI_WTIME() - TT0)
        PENDING = .FALSE.
        INFO_STEP = IBUF_N(0)
        IPIV_K = IBUF_N(1:KW)
        DBUF = DBUF_N(1:2*KW)
      ELSE
        IF (OWNER) THEN
          TT0 = MPI_WTIME()
          CALL SYM_DIAG_FACTOR(K, KW, IPIV_K, DBUF, INFO_STEP)
          TB(1) = TB(1) + (MPI_WTIME() - TT0)
        END IF
        TT0 = MPI_WTIME()
        CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
        IF (INFO_STEP .EQ. 0) THEN
          CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, ROOT_K, GRID%COMM, IERR)
          CALL MPI_BCAST(DBUF, 2*KW, MPI_COMPLEX, ROOT_K, GRID%COMM, IERR)
        END IF
        TB(4) = TB(4) + (MPI_WTIME() - TT0)
      END IF
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(WRECV, LTRECV, IPIV_K, SWP, DBUF, E, WORK, DWORK)
        IF (ALLOCATED(IBUF_N)) DEALLOCATE(IBUF_N, DBUF_N, APK_N)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF
      A%IPIV(1:KW,K)  = IPIV_K
      A%DDIAG(1:KW,K) = DBUF(1:KW)
      A%DSUB(1:KW,K)  = DBUF(KW+1:2*KW)
      SWP = ABS(IPIV_K)
      CALL COUNT_PIVOTS(KW, IPIV_K, DBUF(KW+1:2*KW), NS_K, N2_K)
      NS_TOT = NS_TOT + NS_K; N2_TOT = N2_TOT + N2_K

      ! 2) permutation des LIGNES des L(K,J<K) deja finalises (ligne PROW_K).
      TT0 = MPI_WTIME()
      IF (GRID%MYROW .EQ. PROW_K .AND. NS_K .GT. 0) THEN
        DO lj = 1, A%NBLKC_LOC
          J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
          IF (J .GE. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP(A%OFFBLK(li_K,lj), KW, SWP)
        END DO
      END IF
      TB(2) = TB(2) + (MPI_WTIME() - TT0)

      ! 3) L_KK diffuse (compact) dans la colonne PCOL_K ; permutation des
      !    COLONNES des A(I>K,K) puis W = A(I,K)*P_K*L_KK^-T. Deja fait a
      !    l'etape K-1 en anticipation (POST_NEXT_DIAG).
      IF (GRID%MYCOL .EQ. PCOL_K .AND. .NOT. (LA .AND. K .GT. 1)) THEN
        TT0 = MPI_WTIME()
        CALL BCAST_DIAG_LOWER(A, K, KW, OWNER, PROW_K, GRID%COL_COMM, DKK)
        TB(4) = TB(4) + (MPI_WTIME() - TT0)
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          IF (NS_K .GT. 0) CALL BLR_BLOCK_COL_SWAP(A%OFFBLK(li,lj_K), KW, SWP)
          TB(2) = TB(2) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BLR_TRSM_RIGHT_LOWER_TRANS_UNIT(A%OFFBLK(li,lj_K), DKK, KW)
          TB(3) = TB(3) + (MPI_WTIME() - TT0)
        END DO
        DEALLOCATE(DKK)
      END IF

      ! 4) W(I,K) diffuse le long des lignes de grille (deja compresse) ;
      !    le proprietaire garde ensuite L(I,K) = W(I,K)*D_K^-1.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, A%OFFBLK(li,lj_K))
          TB(4) = TB(4) + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          WRECV(li) = A%OFFBLK(li,lj_K)
          CALL BLR_BLOCK_DINV_RIGHT(A%OFFBLK(li,lj_K), KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          TB(3) = TB(3) + (MPI_WTIME() - TT0)
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
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A%BEND(I)-A%BSTART(I)+1, KW, WRECV(li))
          TB(4) = TB(4) + (MPI_WTIME() - TT0)
        END DO
      END IF

      ! 5) L(J,K)^T diffuse le long des colonnes de grille, racine =
      !    processus (ligne(J), MYCOL), qui possede WRECV(J) depuis l'etape 4.
      !    Tous les membres d'un COL_COMM (meme MYCOL) parcourent la meme
      !    liste lj : condition uniforme sur le communicateur.
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        PROW_J = MOD(J-1, GRID%PR)
        IF (GRID%MYROW .EQ. PROW_J) THEN
          TT0 = MPI_WTIME()
          li_J = (J-1)/GRID%PR + 1
          TMPB = WRECV(li_J)
          CALL BLR_BLOCK_DINV_RIGHT(TMPB, KW, IPIV_K, A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
          CALL BLR_BLOCK_TRANSPOSE(TMPB, LTRECV(lj))
          TB(3) = TB(3) + (MPI_WTIME() - TT0)
        END IF
        TT0 = MPI_WTIME()
        CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_J, KW, NJ, LTRECV(lj))
        TB(4) = TB(4) + (MPI_WTIME() - TT0)
      END DO

      ! 6) Schur, triangle inferieur : A(I,J) -= W(I,K)*L(J,K)^T, I>=J>K.
      !    En anticipation, la colonne K+1 (premier lj traite par la colonne
      !    de grille PCOL_(K+1)) est mise a jour d'abord, puis le bloc
      !    (K+1,K+1) est factorise et diffuse sans bloquer ; les autres
      !    colonnes de grille postent leur reception avant leur Schur.
      IF (LA .AND. K .LT. NBLK) THEN
        IF (GRID%MYCOL .NE. BLR_COL_OWNER(A%CMAP, K+1, GRID)) CALL POST_NEXT_DIAG(K+1)
      END IF
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A%BEND(J) - A%BSTART(J) + 1
        DO li = 1, A%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LT. J) CYCLE
          MI = A%BEND(I) - A%BSTART(I) + 1
          TT0 = MPI_WTIME()
          IF (I .EQ. J) THEN
            CALL SYM_PACKED_SCHUR_UPDATE(MI, A%DIAGP(I)%AP, KW, WRECV(li), LTRECV(lj), DWORK, A%NB)
            TB(5) = TB(5) + (MPI_WTIME() - TT0)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK(A%OFFBLK(li,lj), MI, NJ, KW, WRECV(li), LTRECV(lj), TOL, WS, &
                                        USE_MOVE_ALLOC=USE_MOVE_ALLOC)
            TB(6) = TB(6) + (MPI_WTIME() - TT0)
          END IF
        END DO
        IF (LA .AND. J .EQ. K+1) CALL POST_NEXT_DIAG(K+1)
      END DO

      DEALLOCATE(IPIV_K, SWP, DBUF)
      CALL BLR_TRIM_HEAP()        ! fin d'etape : memoire liberee rendue au systeme
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('2) Factorisation LDLt (Bunch-Kaufman confine)', K, NBLK, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(K, NBLK)) CALL PROGRESS_CALLBACK(K, NBLK)
      END IF
      IF (RSS_TRACE) RSS_MAX_MB = MAX(RSS_MAX_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
    END DO
    IF (RSS_TRACE) WRITE(6,'(A,F8.1,A)') '     [RSS trace] MAX observe pendant la factorisation = ', RSS_MAX_MB, ' Mo'
    TB(7) = MPI_WTIME() - TLOOP0
    IF (PRESENT(T_BREAKDOWN)) CALL MPI_ALLREDUCE(TB, T_BREAKDOWN, 7, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_BREAKDOWN_MIN)) CALL MPI_ALLREDUCE(TB, T_BREAKDOWN_MIN, 7, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, &
                                                     IERR)

    DEALLOCATE(WRECV, LTRECV, E, WORK, DWORK)
    IF (ALLOCATED(IBUF_N)) DEALLOCATE(IBUF_N, DBUF_N, APK_N)
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

  CONTAINS

    ! Bunch-Kaufman sur le bloc diagonal compact KD (proprietaire seulement) :
    ! IP = IPIV, DB = (D, E), INFO_D = indice global du premier pivot nul.
    SUBROUTINE SYM_DIAG_FACTOR(KD, KWD, IP, DB, INFO_D)
      INTEGER, INTENT(IN)        :: KD, KWD
      INTEGER, INTENT(OUT)       :: IP(KWD), INFO_D
      COMPLEX(WP), INTENT(OUT)   :: DB(2*KWD)
      INTEGER :: II, IINFO_D
      INFO_D = 0
      CALL SYM_UNPACK_LOWER(KWD, A%DIAGP(KD)%AP, DWORK, A%NB)
      CALL SYM_BK_FACTOR(KWD, DWORK, A%NB, E, IP, WORK, LWORK, IINFO_D)
      IF (IINFO_D .NE. 0) INFO_D = A%BSTART(KD) - 1 + IINFO_D
      DO II = 1, KWD
        DB(II) = DWORK(II,II)
      END DO
      CALL SYM_PACK_LOWER(KWD, DWORK, A%NB, A%DIAGP(KD)%AP)
      DB(KWD+1:2*KWD) = E(1:KWD)
    END SUBROUTINE SYM_DIAG_FACTOR

    ! Anticipation : factorise (proprietaire) le bloc KN deja a jour et
    ! lance les diffusions non bloquantes de (INFO, IPIV), (D, E) sur la
    ! grille et du bloc compact sur la colonne de grille PCOL_KN. La colonne
    ! de grille PCOL_KN attend ensuite le bloc et traite tout de suite son
    ! panneau (sous-etape 3 de l'etape KN : echange de colonnes + TRSM), qui
    ! sort ainsi du chemin critique. Appelee une fois par etape par CHAQUE
    ! processus, dans le meme ordre relatif aux autres collectives de chaque
    ! communicateur. N'utilise que des variables locales (appelee dans les
    ! boucles de l'appelant).
    SUBROUTINE POST_NEXT_DIAG(KN)
      INTEGER, INTENT(IN) :: KN
      INTEGER :: KWN, PRN, PCN, NPN, IERRN, LIN, IN, LJN, NSN, N2N
      INTEGER, ALLOCATABLE :: SWPN(:)
      COMPLEX(WP), ALLOCATABLE :: DKN(:,:)
      DOUBLE PRECISION :: TF0
      KWN = A%BEND(KN) - A%BSTART(KN) + 1
      PRN = MOD(KN-1, GRID%PR)
      PCN = BLR_COL_OWNER(A%CMAP, KN, GRID)
      NPN = KWN*(KWN+1)/2
      IF (GRID%MYROW.EQ.PRN .AND. GRID%MYCOL.EQ.PCN) THEN
        TF0 = MPI_WTIME()
        CALL SYM_DIAG_FACTOR(KN, KWN, IBUF_N(1:KWN), DBUF_N(1:2*KWN), IBUF_N(0))
        APK_N(1:NPN) = A%DIAGP(KN)%AP(1:NPN)
        TB(1) = TB(1) + (MPI_WTIME() - TF0)
      END IF
      CALL MPI_IBCAST(IBUF_N, KWN+1, MPI_INTEGER, PRN*GRID%PC + PCN, GRID%COMM, REQ(1), IERRN)
      CALL MPI_IBCAST(DBUF_N, 2*KWN, MPI_COMPLEX, PRN*GRID%PC + PCN, GRID%COMM, REQ(2), IERRN)
      NREQ = 2
      IF (GRID%MYCOL .EQ. PCN) THEN
        CALL MPI_IBCAST(APK_N, NPN, MPI_COMPLEX, PRN, GRID%COL_COMM, REQ(3), IERRN)
        NREQ = 3
        TF0 = MPI_WTIME()
        CALL MPI_WAITALL(NREQ, REQ, MPI_STATUSES_IGNORE, IERRN)
        TB(4) = TB(4) + (MPI_WTIME() - TF0)
        IF (IBUF_N(0) .EQ. 0) THEN
          ALLOCATE(SWPN(KWN), DKN(KWN,KWN))
          SWPN = ABS(IBUF_N(1:KWN))
          CALL COUNT_PIVOTS(KWN, IBUF_N(1:KWN), DBUF_N(KWN+1:2*KWN), NSN, N2N)
          DKN = (0.0_WP,0.0_WP)
          CALL SYM_UNPACK_LOWER(KWN, APK_N, DKN, KWN)
          LJN = BLR_COL_LOCAL(A%CMAP, KN, GRID)
          DO LIN = 1, A%NBLKR_LOC
            IN = (LIN-1)*GRID%PR + GRID%MYROW + 1
            IF (IN .LE. KN) CYCLE
            TF0 = MPI_WTIME()
            IF (NSN .GT. 0) CALL BLR_BLOCK_COL_SWAP(A%OFFBLK(LIN,LJN), KWN, SWPN)
            TB(2) = TB(2) + (MPI_WTIME() - TF0)
            TF0 = MPI_WTIME()
            CALL BLR_TRSM_RIGHT_LOWER_TRANS_UNIT(A%OFFBLK(LIN,LJN), DKN, KWN)
            TB(3) = TB(3) + (MPI_WTIME() - TF0)
          END DO
          DEALLOCATE(SWPN, DKN)
        END IF
      END IF
      PENDING = .TRUE.
    END SUBROUTINE POST_NEXT_DIAG

  END SUBROUTINE MPI_BLR_SYM_LDLT_FACTOR

  !=======================================================================
  ! Resolution A*X=B apres MPI_BLR_SYM_LDLT_FACTOR. B_LOCAL (MLOC lignes
  ! locales x NRHS_LOC colonnes locales, decoupage RHS_LOCAL_RANGE_MBM) est
  ! ecrase par X. Generique sur la largeur NRHS (appelable par lot).
  !=======================================================================
  SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, A, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS, PROGRESS_CALLBACK, &
                                       T_PHASES, DESCENT_MODE)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)         :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN)  :: A
    INTEGER, INTENT(IN)              :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)       :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL    :: SHOW_PROGRESS
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK
    ! T_PHASES (optionnel) : temps MAX sur la grille de (1) Q*B + descente L,
    ! (2) D, (3) remontee L^T + Q^T.
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_PHASES(3)
    ! DESCENT_MODE (optionnel, defaut 0 = automatique) : algorithme de la
    ! descente L ET de la remontee L^T (memes deux schemas, symetriques ;
    ! en mode 2 la remontee reduit S = sum L(I,K)^T*Y_I sur MES seconds
    ! membres le long de la colonne de grille, sans copie pleine largeur).
    !   1 = "colonne" : X_K rassemble pleine largeur sur la colonne de grille
    !       PCOL_K, qui calcule seule les contributions L(I,K)*X_K pour TOUS
    !       les seconds membres puis les disperse le long des lignes.
    !       Communications proportionnelles a NRHS ; une seule colonne de
    !       grille calcule a chaque etape.
    !   2 = "distribuee" : X_K (tranche locale de seconds membres) diffuse le
    !       long de la colonne de grille, colonne L(.,K) de mes lignes diffusee
    !       COMPRESSEE le long de la ligne de grille (un paquet : en-tetes +
    !       donnees), chaque processus applique L a SES seconds membres --
    !       toutes les colonnes de grille calculent ; volume de L independant
    !       de NRHS.
    !   0 = choix automatique (voir DESCENT_AUTO).
    INTEGER, INTENT(IN), OPTIONAL :: DESCENT_MODE

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE :: LBUF(:), TR(:,:)
    INTEGER, ALLOCATABLE :: HDR(:,:)
    INTEGER :: DMODE, NT, IT, OFF, R, RSI
    INTEGER(8) :: NTOT8
    DOUBLE PRECISION :: TPH(3), TPH0
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
    MLOC = MPI_BLR_SYM_MLOC(A, GRID)
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO
    DMODE = 0
    IF (PRESENT(DESCENT_MODE)) DMODE = DESCENT_MODE
    IF (DMODE .NE. 1 .AND. DMODE .NE. 2) DMODE = DESCENT_AUTO(GRID, A, NRHS)

    TPH0 = MPI_WTIME()
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
      PCOL_K = BLR_COL_OWNER(A%CMAP, K, GRID)
      IF (GRID%MYROW .EQ. PROW_K) THEN
        CALL BCAST_DIAG_LOWER(A, K, KW, GRID%MYCOL.EQ.PCOL_K, PCOL_K, GRID%ROW_COMM, DKK)
        RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
        IF (NRHS_LOC .GT. 0) CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                                         CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
      END IF
      IF (DMODE .EQ. 2) THEN
        CALL DESCENT_STEP_DISTRIBUTED()
        CYCLE
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

    TPH(1) = MPI_WTIME() - TPH0; TPH0 = MPI_WTIME()
    ! === 2) D*W = Z : local a chaque ligne de grille, aucune communication ===
    DO K = 1, NBLK
      KW = A%BEND(K) - A%BSTART(K) + 1
      IF (GRID%MYROW.NE.MOD(K-1,GRID%PR) .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
      CALL APPLY_DINV_LEFT(B_LOCAL(RS,1), LDB, KW, NRHS_LOC, A%IPIV(1:KW,K), A%DDIAG(1:KW,K), A%DSUB(1:KW,K))
    END DO

    TPH(2) = MPI_WTIME() - TPH0; TPH0 = MPI_WTIME()
    IF (DMODE .EQ. 2) THEN
      ! 3) L^T*Y = W, remontee distribuee : voir UPWARD_STEP_DISTRIBUTED.
      DO K = NBLK, 1, -1
        CALL UPWARD_STEP_DISTRIBUTED()
      END DO
    ELSE
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
            CALL BLR_BLOCK_APPLY_T(A%OFFBLK(li,lj), NRHS, CONE, XFULL(RS,1), MAX(MLOC,1), CONE, S, KW)
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
        RS = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
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
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO
    ! part locale (colonnes RSTARTS(MYCOL)...) de la solution pleine largeur
    IF (MLOC.GT.0 .AND. NRHS_LOC.GT.0) &
      B_LOCAL(1:MLOC,1:NRHS_LOC) = XFULL(1:MLOC, RSTARTS(GRID%MYCOL):RSTARTS(GRID%MYCOL)+NRHS_LOC-1)
    DEALLOCATE(XFULL)
    END IF

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
    TPH(3) = MPI_WTIME() - TPH0
    IF (PRESENT(T_PHASES)) CALL MPI_ALLREDUCE(TPH, T_PHASES, 3, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)

  CONTAINS

    ! Etape K de la descente "distribuee" (DESCENT_MODE=2), apres la
    ! resolution L_KK sur la ligne PROW_K (XK_LOC = tranche locale de X_K).
    SUBROUTINE DESCENT_STEP_DISTRIBUTED()
      INTEGER :: LID, ID, MID, RSD, IERRD
      ! X_K, tranche locale de seconds membres, le long de la colonne de grille
      IF (.NOT. ALLOCATED(XK_LOC)) ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
      IF (NRHS_LOC .GT. 0) CALL MPI_BCAST(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERRD)
      CALL BCAST_LCOL()
      IF (NT .GT. 0) THEN
        ! B(I) -= L(I,K) * X_K sur mes seconds membres
        IF (NRHS_LOC .GT. 0) THEN
          OFF = 1; IT = 0; RSD = 1
          DO LID = 1, A%NBLKR_LOC
            ID = (LID-1)*GRID%PR + GRID%MYROW + 1
            MID = A%BEND(ID) - A%BSTART(ID) + 1
            IF (ID .GT. K) THEN
              IT = IT + 1
              IF (HDR(1,IT) .EQ. 1) THEN
                R = HDR(2,IT)
                IF (R .GT. 0) THEN
                  ALLOCATE(TR(R,NRHS_LOC))
                  CALL CGEMM('N', 'N', R, NRHS_LOC, KW, CONE, LBUF(OFF+MID*R), R, XK_LOC, KW, CZERO, TR, R)
                  CALL CGEMM('N', 'N', MID, NRHS_LOC, R, -CONE, LBUF(OFF), MID, TR, R, CONE, B_LOCAL(RSD,1), LDB)
                  DEALLOCATE(TR)
                END IF
                OFF = OFF + (MID+KW)*R
              ELSE
                CALL CGEMM('N', 'N', MID, NRHS_LOC, KW, -CONE, LBUF(OFF), MID, XK_LOC, KW, CONE, B_LOCAL(RSD,1), LDB)
                OFF = OFF + MID*KW
              END IF
            END IF
            RSD = RSD + MID
          END DO
        END IF
        DEALLOCATE(HDR, LBUF)
      END IF
      DEALLOCATE(XK_LOC)
      STEP = STEP + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END SUBROUTINE DESCENT_STEP_DISTRIBUTED

    ! Diffuse le long de la ligne de grille, depuis PCOL_K, la colonne
    ! L(.,K) restreinte a MES blocs-lignes I>K : en-tetes HDR(2,NT)
    ! (bas-rang ?, rang) puis donnees empilees dans LBUF (UFAC puis VFAC,
    ! ou DENSE), deux diffusions. Meme liste de lignes pour tout ROW_COMM
    ! (meme MYROW) : condition uniforme. NT = 0 : rien n'est alloue.
    SUBROUTINE BCAST_LCOL()
      INTEGER :: LID, ID, MID, LJD, IERRD
      NT = 0
      DO LID = 1, A%NBLKR_LOC
        ID = (LID-1)*GRID%PR + GRID%MYROW + 1
        IF (ID .GT. K) NT = NT + 1
      END DO
      IF (NT .EQ. 0) RETURN
      LJD = BLR_COL_LOCAL(A%CMAP, K, GRID)
      ALLOCATE(HDR(2,NT))
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        IT = 0
        DO LID = 1, A%NBLKR_LOC
          ID = (LID-1)*GRID%PR + GRID%MYROW + 1
          IF (ID .LE. K) CYCLE
          IT = IT + 1
          HDR(1,IT) = MERGE(1, 0, A%OFFBLK(LID,LJD)%IS_LOWRANK)
          HDR(2,IT) = A%OFFBLK(LID,LJD)%RANK
        END DO
      END IF
      CALL MPI_BCAST(HDR, 2*NT, MPI_INTEGER, PCOL_K, GRID%ROW_COMM, IERRD)
      NTOT8 = 0_8; IT = 0
      DO LID = 1, A%NBLKR_LOC
        ID = (LID-1)*GRID%PR + GRID%MYROW + 1
        IF (ID .LE. K) CYCLE
        IT = IT + 1
        MID = A%BEND(ID) - A%BSTART(ID) + 1
        IF (HDR(1,IT) .EQ. 1) THEN
          NTOT8 = NTOT8 + INT(MID+KW,8)*INT(HDR(2,IT),8)
        ELSE
          NTOT8 = NTOT8 + INT(MID,8)*INT(KW,8)
        END IF
      END DO
      ALLOCATE(LBUF(MAX(NTOT8,1_8)))
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        OFF = 1
        DO LID = 1, A%NBLKR_LOC
          ID = (LID-1)*GRID%PR + GRID%MYROW + 1
          IF (ID .LE. K) CYCLE
          MID = A%BEND(ID) - A%BSTART(ID) + 1
          ASSOCIATE (BK => A%OFFBLK(LID,LJD))
            IF (BK%IS_LOWRANK) THEN
              R = BK%RANK
              IF (R .GT. 0) THEN
                LBUF(OFF:OFF+MID*R-1) = RESHAPE(BK%UFAC, [MID*R])
                LBUF(OFF+MID*R:OFF+(MID+KW)*R-1) = RESHAPE(BK%VFAC, [R*KW])
              END IF
              OFF = OFF + (MID+KW)*R
            ELSE
              LBUF(OFF:OFF+MID*KW-1) = RESHAPE(BK%DENSE, [MID*KW])
              OFF = OFF + MID*KW
            END IF
          END ASSOCIATE
        END DO
      END IF
      IF (NTOT8 .GT. 0) CALL MPI_BCAST(LBUF, INT(NTOT8), MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERRD)
    END SUBROUTINE BCAST_LCOL

    ! Etape K de la remontee "distribuee" (DESCENT_MODE=2) : chaque processus
    ! calcule sur SES seconds membres S = sum_(I>K, mes lignes) L(I,K)^T*Y_I
    ! (colonne L(.,K) recue par BCAST_LCOL), S est reduit le long de la
    ! colonne de grille sur la ligne PROW_K, qui fait Y_K = L_KK^-T*(W_K - S)
    ! en place dans B_LOCAL. Aucune copie pleine largeur (pas de XFULL), et
    ! toutes les colonnes de grille calculent.
    SUBROUTINE UPWARD_STEP_DISTRIBUTED()
      INTEGER :: LID, ID, MID, RSD, IERRD
      COMPLEX(WP), ALLOCATABLE :: SU(:,:), SR(:,:)
      KW = A%BEND(K) - A%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A%CMAP, K, GRID)
      CALL BCAST_LCOL()
      ALLOCATE(SU(KW,MAX(NRHS_LOC,1)))
      SU = CZERO
      IF (NT .GT. 0) THEN
        IF (NRHS_LOC .GT. 0) THEN
          OFF = 1; IT = 0; RSD = 1
          DO LID = 1, A%NBLKR_LOC
            ID = (LID-1)*GRID%PR + GRID%MYROW + 1
            MID = A%BEND(ID) - A%BSTART(ID) + 1
            IF (ID .GT. K) THEN
              IT = IT + 1
              IF (HDR(1,IT) .EQ. 1) THEN
                R = HDR(2,IT)
                IF (R .GT. 0) THEN
                  ALLOCATE(TR(R,NRHS_LOC))
                  ! L^T*Y = V^T*(U^T*Y)
                  CALL CGEMM('T', 'N', R, NRHS_LOC, MID, CONE, LBUF(OFF), MID, B_LOCAL(RSD,1), LDB, CZERO, TR, R)
                  CALL CGEMM('T', 'N', KW, NRHS_LOC, R, CONE, LBUF(OFF+MID*R), R, TR, R, CONE, SU, KW)
                  DEALLOCATE(TR)
                END IF
                OFF = OFF + (MID+KW)*R
              ELSE
                CALL CGEMM('T', 'N', KW, NRHS_LOC, MID, CONE, LBUF(OFF), MID, B_LOCAL(RSD,1), LDB, CONE, SU, KW)
                OFF = OFF + MID*KW
              END IF
            END IF
            RSD = RSD + MID
          END DO
        END IF
        DEALLOCATE(HDR, LBUF)
      END IF
      ! Reduction sur la ligne PROW_K (tous les membres de la colonne de
      ! grille, meme ceux sans ligne I>K : condition uniforme K < NBLK).
      IF (K .LT. NBLK .AND. NRHS_LOC .GT. 0) THEN
        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(SR(KW,NRHS_LOC))
        ELSE
          ALLOCATE(SR(1,1))
        END IF
        CALL MPI_REDUCE(SU, SR, KW*NRHS_LOC, MPI_COMPLEX, MPI_SUM, PROW_K, GRID%COL_COMM, IERRD)
      END IF
      IF (GRID%MYROW .EQ. PROW_K) THEN
        CALL BCAST_DIAG_LOWER(A, K, KW, GRID%MYCOL.EQ.PCOL_K, PCOL_K, GRID%ROW_COMM, DKK)
        IF (NRHS_LOC .GT. 0) THEN
          RSD = LOCAL_ROW_OFFSET_SYM(A, GRID, K)
          IF (K .LT. NBLK) B_LOCAL(RSD:RSD+KW-1,1:NRHS_LOC) = B_LOCAL(RSD:RSD+KW-1,1:NRHS_LOC) - SR
          CALL CTRSM('Left', 'Lower', 'Transpose', 'Unit', KW, NRHS_LOC, CONE, DKK, KW, B_LOCAL(RSD,1), LDB)
        END IF
        DEALLOCATE(DKK)
      END IF
      IF (ALLOCATED(SR)) DEALLOCATE(SR)
      DEALLOCATE(SU)
      STEP = STEP + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution LDLt', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END SUBROUTINE UPWARD_STEP_DISTRIBUTED

  END SUBROUTINE MPI_BLR_SYM_LDLT_SOLVE_2D

  ! Largeur effective moyenne de L : stockage total des blocs I>J (dense
  ! M*N, bas-rang R*(M+N)) divise par la somme des M, i.e. nombre moyen de
  ! colonnes equivalentes par ligne de L. C'est le volume de L diffuse par
  ! ligne de matrice et par etape en mode 2, a comparer a NRHS (mode 1).
  REAL(WP) FUNCTION SYM_L_EFFECTIVE_WIDTH(GRID, A)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    INTEGER :: li, lj, I, J, IERR
    DOUBLE PRECISION :: LOC(2), GLOB(2)
    LOC = 0.0D0
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LE. J) CYCLE
        ASSOCIATE (BK => A%OFFBLK(li,lj))
          IF (BK%IS_LOWRANK) THEN
            LOC(1) = LOC(1) + DBLE(BK%RANK)*DBLE(BK%M + BK%N)
          ELSE
            LOC(1) = LOC(1) + DBLE(BK%M)*DBLE(BK%N)
          END IF
          LOC(2) = LOC(2) + DBLE(BK%M)
        END ASSOCIATE
      END DO
    END DO
    CALL MPI_ALLREDUCE(LOC, GLOB, 2, MPI_DOUBLE_PRECISION, MPI_SUM, GRID%COMM, IERR)
    SYM_L_EFFECTIVE_WIDTH = REAL(GLOB(1)/MAX(GLOB(2),1.0D0), WP)
  END FUNCTION SYM_L_EFFECTIVE_WIDTH

  ! Choix automatique de la descente/remontee (DESCENT_MODE=0). Par etape
  ! et par ligne de matrice, le mode 1 deplace et calcule de l'ordre de NRHS
  ! colonnes sur une seule colonne de grille, le mode 2 diffuse la largeur
  ! effective W de L (independante de NRHS) puis calcule sur toutes les
  ! colonnes de grille. Mesure (N=16000, NB=400, TOL=1e-4, W=4.5, grilles
  ! 2x2 et 2x4) : mode 1 meilleur a NRHS=1 (x1.5 a x1.8), egalite vers
  ! NRHS=8, mode 2 x1.7 a NRHS=32 et x3 a NRHS=512-1024 -> mode 2 des que
  ! NRHS > 1.5*W. Decision identique sur toute la grille (W reduit).
  INTEGER FUNCTION DESCENT_AUTO(GRID, A, NRHS)
    TYPE(GRID_T), INTENT(IN)        :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    INTEGER, INTENT(IN)             :: NRHS
    REAL(WP) :: W
    W = SYM_L_EFFECTIVE_WIDTH(GRID, A)
    IF (REAL(NRHS,WP) .GT. 1.5_WP*W) THEN
      DESCENT_AUTO = 2
    ELSE
      DESCENT_AUTO = 1
    END IF
  END FUNCTION DESCENT_AUTO

END MODULE MPI_BLR_LDLT_MOD
