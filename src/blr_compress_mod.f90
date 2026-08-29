!===============================================================================
! ETAPE 2 (partie 1/2) : primitives de compression basse-rang d'un bloc
! dense, complexe simple precision.
!
! METHODE : SVD tronquee (CGESVD), pas QR avec pivotage de colonnes (ce que
! font MUMPS/STRUMPACK en production, plus rapide) -- choix DELIBERE pour
! cette premiere version : la SVD donne directement le rang numerique (a
! une tolerance donnee) sans bookkeeping de permutation de colonnes a
! gerer, ce qui reduit le risque de bug pendant qu'on met l'algorithme au
! point (cf. le bug de pivot trouve a l'etape 1). On pourra basculer vers
! une QR tronquee plus tard une fois la logique BLR validee.
!
! REPRESENTATION : un bloc M x N est stocke soit DENSE (tel quel), soit
! BASSE-RANG comme A ~= UFAC * VFAC avec UFAC (MxR) et VFAC (RxN) -- le
! facteur singulier Sigma est absorbe dans VFAC (VFAC(i,:) = S(i)*Vh(i,:)),
! donc UFAC*VFAC reconstruit directement l'approximation tronquee, sans
! manipulation supplementaire.
!
! SEUIL DE RANG : on ne compresse que si le rang numerique trouve (a la
! tolerance TOL, relative a la plus grande valeur singuliere) est
! significativement plus petit que MIN(M,N) -- sinon la compression
! n'apporte rien (voire coute plus cher a exploiter qu'un GEMM dense direct)
! et on garde le bloc dense.
!===============================================================================
! GENERALISATION (demande explicite : une librairie d'algebre lineaire sur
! blocs full/low-rank qui tire au mieux parti de MKL/LAPACK) : ce module,
! deja le point de convergence naturel de toutes les primitives bloc du
! projet, en devient le catalogue COMPLET et le point d'entree UNIQUE --
! chaque routine choisit systematiquement le chemin de calcul le moins cher
! selon que ses operandes sont denses ou basse-rang, et s'appuie sur les
! memes noyaux BLAS3/LAPACK threades par MKL que le reste du projet
! (CGEMM, CTRSM, CGESVD, CGEQRF/CUNGQR) -- aucune reimplementation "maison"
! de ces operations, uniquement de la logique de DISPATCH et de
! recombinaison des facteurs low-rank autour d'elles. Les routines
! partagees entre plusieurs solveurs (BLOCK_PARTITION, taille memoire d'un
! bloc, LU sans pivot, TRSM applique a un facteur low-rank seul) ont ete
! consolidees ICI (au lieu d'etre dupliquees dans blr_lu_mod.f90,
! blr_lu_storage_mod.f90, blr_matrix_mod.f90 comme precedemment) : tout
! module qui fait deja "USE BLR_COMPRESS_MOD" les recoit automatiquement,
! aucun changement de nom de fichier/module necessaire ailleurs.
!
! PISTES MKL NON ADOPTEES ICI (documentees, pas implementees -- eviter la
! sur-ingenierie tant que le besoin n'est pas mesure) : CGEMM3M (variante
! MKL de CGEMM, parfois plus rapide sur gros GEMM complexes via l'algorithme
! 3M) pourrait remplacer CGEMM en UN SEUL point (ici) si un profilage futur
! le justifie ; le "batching" de petits GEMM (CGEMM_BATCH) pourrait accelerer
! les tres nombreuses petites operations de rang faible pendant une
! factorisation BLR, mais demande une refonte des boucles appelantes pour
! grouper les appels, hors scope de cette generalisation.
MODULE BLR_COMPRESS_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE OMP_LIB, ONLY: OMP_GET_WTIME
  IMPLICIT NONE

  ! Fraction de MIN(M,N) en-dessous de laquelle la compression est jugee
  ! "rentable". Expose en PARAMETER modifiable si besoin de retuner plus
  ! tard sans casser l'interface des routines.
  REAL(WP), PARAMETER :: COMPRESS_WORTHWHILE_RATIO = 0.8_WP

  TYPE :: BLR_BLOCK_T
    LOGICAL :: IS_LOWRANK = .FALSE.
    INTEGER :: M = 0, N = 0, RANK = 0
    COMPLEX(WP), ALLOCATABLE :: DENSE(:,:)         ! si IS_LOWRANK = .FALSE.
    COMPLEX(WP), ALLOCATABLE :: UFAC(:,:), VFAC(:,:)  ! si IS_LOWRANK = .TRUE.
  END TYPE BLR_BLOCK_T

  !-----------------------------------------------------------------------
  ! REUTILISATION DES TABLEAUX DE TRAVAIL TRANSITOIRES (demande explicite :
  ! reduire l'allocation/desallocation repetee des scratch buffers de
  ! COMPRESS_BLOCK/FAST_RECOMPRESS -- SVD, QR -- au lieu d'un ALLOCATE/
  ! DEALLOCATE a CHAQUE appel). Toutes les matrices 2D sont dimensionnees
  ! UNE FOIS a (MAXDIM,MAXDIM) (MAXDIM = plus grande taille de bloc/rang
  ! total attendue, typiquement NB) et REUTILISEES telles quelles via leur
  ! parametre LDA/LDU/LDVT (convention LAPACK standard : LDA >= M, on peut
  ! passer un buffer plus grand que la sous-partie logiquement utilisee)
  ! -- AUCUNE copie de section n'est faite (on passe le tableau ENTIER, pas
  ! une tranche WS%X(1:M,1:N), qui forcerait une copie temporaire pour un
  ! tableau non contigu). Le WORK complexe (LWORK) est dimensionne une fois
  ! au pire cas (MAXDIM,MAXDIM) : LAPACK accepte un LWORK plus grand que le
  ! minimum requis pour un probleme plus petit (documente, sans erreur).
  !
  ! THREAD-SAFETY : UNE instance de ce type ne doit JAMAIS etre partagee
  ! entre taches concurrentes (threads OpenMP) -- chaque tache (ou
  ! processus MPI, naturellement isole en memoire) doit avoir SA PROPRE
  ! instance. Voir BLR_WORKSPACE_ALLOC.
  !-----------------------------------------------------------------------
  TYPE :: BLR_WORKSPACE_T
    INTEGER :: MAXDIM = 0
    ! -- COMPRESS_BLOCK (SVD tronquee) --
    COMPLEX(WP), ALLOCATABLE :: BCOPY(:,:), UFULL(:,:), VTFULL(:,:), CWORK1(:)
    REAL(WP), ALLOCATABLE    :: S1(:), RWORK1(:)
    ! -- FAST_RECOMPRESS (QR + petite SVD) --
    COMPLEX(WP), ALLOCATABLE :: UCAT(:,:), VCATT(:,:), TAUU(:), TAUV(:)
    COMPLEX(WP), ALLOCATABLE :: RU(:,:), RV(:,:), WMAT(:,:), UW(:,:), VWH(:,:), SVWH(:,:)
    COMPLEX(WP), ALLOCATABLE :: CWORK2(:)
    REAL(WP), ALLOCATABLE    :: SW(:), RWORK2(:)
  END TYPE BLR_WORKSPACE_T

CONTAINS

  !-----------------------------------------------------------------------
  ! Alloue un workspace reutilisable pour des blocs de taille au plus
  ! MAXDIM x MAXDIM (typiquement NB, la taille de bloc du solveur
  ! appelant). Le LWORK complexe est determine UNE FOIS via les requetes
  ! standard (LWORK=-1) de CGESVD/CGEQRF/CUNGQR au pire cas (MAXDIM,MAXDIM),
  ! puis reutilise pour tous les appels ulterieurs quelle que soit leur
  ! taille reelle (toujours <= MAXDIM).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_WORKSPACE_ALLOC(WS, MAXDIM)
    TYPE(BLR_WORKSPACE_T), INTENT(OUT) :: WS
    INTEGER, INTENT(IN) :: MAXDIM
    INTEGER :: LWORK, INFO
    COMPLEX(WP) :: WQ(1)

    WS%MAXDIM = MAX(MAXDIM, 1)

    ALLOCATE(WS%BCOPY(WS%MAXDIM,WS%MAXDIM), WS%UFULL(WS%MAXDIM,WS%MAXDIM), WS%VTFULL(WS%MAXDIM,WS%MAXDIM))
    ALLOCATE(WS%S1(WS%MAXDIM), WS%RWORK1(5*WS%MAXDIM))

    ALLOCATE(WS%UCAT(WS%MAXDIM,WS%MAXDIM), WS%VCATT(WS%MAXDIM,WS%MAXDIM))
    ALLOCATE(WS%TAUU(WS%MAXDIM), WS%TAUV(WS%MAXDIM))
    ALLOCATE(WS%RU(WS%MAXDIM,WS%MAXDIM), WS%RV(WS%MAXDIM,WS%MAXDIM), WS%WMAT(WS%MAXDIM,WS%MAXDIM))
    ALLOCATE(WS%UW(WS%MAXDIM,WS%MAXDIM), WS%VWH(WS%MAXDIM,WS%MAXDIM), WS%SVWH(WS%MAXDIM,WS%MAXDIM))
    ALLOCATE(WS%SW(WS%MAXDIM), WS%RWORK2(5*WS%MAXDIM))

    LWORK = 1
    CALL CGESVD('S', 'S', WS%MAXDIM, WS%MAXDIM, WS%BCOPY, WS%MAXDIM, WS%S1, &
                 WS%UFULL, WS%MAXDIM, WS%VTFULL, WS%MAXDIM, WQ, -1, WS%RWORK1, INFO)
    LWORK = MAX(LWORK, INT(REAL(WQ(1))))
    CALL CGEQRF(WS%MAXDIM, WS%MAXDIM, WS%UCAT, WS%MAXDIM, WS%TAUU, WQ, -1, INFO)
    LWORK = MAX(LWORK, INT(REAL(WQ(1))))
    CALL CUNGQR(WS%MAXDIM, WS%MAXDIM, WS%MAXDIM, WS%UCAT, WS%MAXDIM, WS%TAUU, WQ, -1, INFO)
    LWORK = MAX(LWORK, INT(REAL(WQ(1))))
    CALL CGESVD('S', 'S', WS%MAXDIM, WS%MAXDIM, WS%WMAT, WS%MAXDIM, WS%SW, &
                 WS%UW, WS%MAXDIM, WS%VWH, WS%MAXDIM, WQ, -1, WS%RWORK2, INFO)
    LWORK = MAX(LWORK, INT(REAL(WQ(1))))

    ALLOCATE(WS%CWORK1(LWORK), WS%CWORK2(LWORK))
  END SUBROUTINE BLR_WORKSPACE_ALLOC

  SUBROUTINE BLR_WORKSPACE_FREE(WS)
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT) :: WS
    IF (ALLOCATED(WS%BCOPY))  DEALLOCATE(WS%BCOPY)
    IF (ALLOCATED(WS%UFULL))  DEALLOCATE(WS%UFULL)
    IF (ALLOCATED(WS%VTFULL)) DEALLOCATE(WS%VTFULL)
    IF (ALLOCATED(WS%CWORK1)) DEALLOCATE(WS%CWORK1)
    IF (ALLOCATED(WS%S1))     DEALLOCATE(WS%S1)
    IF (ALLOCATED(WS%RWORK1)) DEALLOCATE(WS%RWORK1)
    IF (ALLOCATED(WS%UCAT))   DEALLOCATE(WS%UCAT)
    IF (ALLOCATED(WS%VCATT))  DEALLOCATE(WS%VCATT)
    IF (ALLOCATED(WS%TAUU))   DEALLOCATE(WS%TAUU)
    IF (ALLOCATED(WS%TAUV))   DEALLOCATE(WS%TAUV)
    IF (ALLOCATED(WS%RU))     DEALLOCATE(WS%RU)
    IF (ALLOCATED(WS%RV))     DEALLOCATE(WS%RV)
    IF (ALLOCATED(WS%WMAT))   DEALLOCATE(WS%WMAT)
    IF (ALLOCATED(WS%UW))     DEALLOCATE(WS%UW)
    IF (ALLOCATED(WS%VWH))    DEALLOCATE(WS%VWH)
    IF (ALLOCATED(WS%SVWH))   DEALLOCATE(WS%SVWH)
    IF (ALLOCATED(WS%CWORK2)) DEALLOCATE(WS%CWORK2)
    IF (ALLOCATED(WS%SW))     DEALLOCATE(WS%SW)
    IF (ALLOCATED(WS%RWORK2)) DEALLOCATE(WS%RWORK2)
    WS%MAXDIM = 0
  END SUBROUTINE BLR_WORKSPACE_FREE

  !-----------------------------------------------------------------------
  ! Compresse le bloc dense B (M x N) en BLK, a la tolerance relative TOL
  ! (TOL=0 force le mode dense ; typiquement TOL~1e-4 a 1e-6 en simple
  ! precision). B n'est PAS modifie (copie interne avant SVD, qui est
  ! destructive).
  !-----------------------------------------------------------------------
  ! WS (optionnel) : workspace pre-alloue (BLR_WORKSPACE_ALLOC) reutilise
  ! pour BCOPY/UFULL/VTFULL/CWORK1/S1/RWORK1 -- AUCUNE allocation dans cette
  ! routine si fourni (M,N doivent etre <= WS%MAXDIM). Absent : comportement
  ! d'origine (allocation/desallocation locale a chaque appel).
  SUBROUTINE COMPRESS_BLOCK(M, N, B, LDB, TOL, BLK, WS)
    INTEGER, INTENT(IN)          :: M, N, LDB
    COMPLEX(WP), INTENT(IN)      :: B(LDB,N)
    REAL(WP), INTENT(IN)         :: TOL
    TYPE(BLR_BLOCK_T), INTENT(OUT) :: BLK
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS

    COMPLEX(WP), ALLOCATABLE, TARGET :: BCOPY_LOC(:,:), UFULL_LOC(:,:), VTFULL_LOC(:,:), CWORK_LOC(:)
    REAL(WP), ALLOCATABLE, TARGET :: S_LOC(:), RWORK_LOC(:)
    COMPLEX(WP), POINTER :: BCOPY(:,:), UFULL(:,:), VTFULL(:,:), CWORK(:)
    REAL(WP), POINTER :: S(:), RWORK(:)
    INTEGER :: MN, LDW, LDV, LWORK, INFO, R, I
    REAL(WP) :: THRESH

    BLK%M = M
    BLK%N = N
    MN = MIN(M,N)

    IF (TOL .LE. 0.0_WP .OR. MN .EQ. 0) THEN
      BLK%IS_LOWRANK = .FALSE.
      ALLOCATE(BLK%DENSE(M,N))
      IF (M.GT.0 .AND. N.GT.0) BLK%DENSE = B(1:M,1:N)
      RETURN
    END IF

    IF (PRESENT(WS)) THEN
      LDW = WS%MAXDIM
      LDV = WS%MAXDIM
      BCOPY  => WS%BCOPY
      UFULL  => WS%UFULL
      VTFULL => WS%VTFULL
      S      => WS%S1
      RWORK  => WS%RWORK1
      CWORK  => WS%CWORK1
      LWORK  = SIZE(CWORK)
    ELSE
      LDW = M
      LDV = MN
      ALLOCATE(BCOPY_LOC(M,N), UFULL_LOC(M,MN), VTFULL_LOC(MN,N), S_LOC(MN), RWORK_LOC(5*MN))
      BCOPY => BCOPY_LOC; UFULL => UFULL_LOC; VTFULL => VTFULL_LOC; S => S_LOC; RWORK => RWORK_LOC
      ALLOCATE(CWORK_LOC(1))
      LWORK = -1
      CALL CGESVD('S', 'S', M, N, BCOPY, LDW, S, UFULL, LDW, VTFULL, LDV, CWORK_LOC, LWORK, RWORK, INFO)
      LWORK = INT(REAL(CWORK_LOC(1)))
      DEALLOCATE(CWORK_LOC)
      ALLOCATE(CWORK_LOC(MAX(LWORK,1)))
      CWORK => CWORK_LOC
    END IF

    BCOPY(1:M,1:N) = B(1:M,1:N)
    CALL CGESVD('S', 'S', M, N, BCOPY, LDW, S, UFULL, LDW, VTFULL, LDV, &
                 CWORK, LWORK, RWORK, INFO)

    IF (INFO .NE. 0) THEN
      ! SVD non convergee (rarissime) : repli securise sur le mode dense
      BLK%IS_LOWRANK = .FALSE.
      ALLOCATE(BLK%DENSE(M,N))
      BLK%DENSE = B(1:M,1:N)
      IF (.NOT. PRESENT(WS)) DEALLOCATE(BCOPY_LOC, UFULL_LOC, VTFULL_LOC, S_LOC, RWORK_LOC, CWORK_LOC)
      RETURN
    END IF

    THRESH = TOL * MAX(S(1), TINY(1.0_WP))
    R = 0
    DO I = 1, MN
      IF (S(I) .GT. THRESH) R = R + 1
    END DO
    R = MAX(R, 1)

    IF (REAL(R,WP) .LT. COMPRESS_WORTHWHILE_RATIO*REAL(MN,WP)) THEN
      BLK%IS_LOWRANK = .TRUE.
      BLK%RANK = R
      ALLOCATE(BLK%UFAC(M,R), BLK%VFAC(R,N))
      BLK%UFAC = UFULL(1:M,1:R)
      DO I = 1, R
        BLK%VFAC(I,1:N) = S(I) * VTFULL(I,1:N)
      END DO
    ELSE
      BLK%IS_LOWRANK = .FALSE.
      ALLOCATE(BLK%DENSE(M,N))
      BLK%DENSE = B(1:M,1:N)
    END IF

    IF (.NOT. PRESENT(WS)) DEALLOCATE(BCOPY_LOC, UFULL_LOC, VTFULL_LOC, S_LOC, RWORK_LOC, CWORK_LOC)
  END SUBROUTINE COMPRESS_BLOCK

  !-----------------------------------------------------------------------
  ! Reconstruit la forme dense d'un bloc (compresse ou non) -- utilise pour
  ! les tests et pour EXPAND explicite si besoin.
  !-----------------------------------------------------------------------
  SUBROUTINE EXPAND_BLOCK(BLK, OUT, LDOUT)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: BLK
    INTEGER, INTENT(IN)           :: LDOUT
    COMPLEX(WP), INTENT(OUT)      :: OUT(LDOUT,BLK%N)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)

    IF (BLK%IS_LOWRANK) THEN
      CALL CGEMM('No transpose', 'No transpose', BLK%M, BLK%N, BLK%RANK, &
                  CONE, BLK%UFAC, BLK%M, BLK%VFAC, BLK%RANK, CZERO, OUT, LDOUT)
    ELSE
      OUT(1:BLK%M,1:BLK%N) = BLK%DENSE
    END IF
  END SUBROUTINE EXPAND_BLOCK

  !-----------------------------------------------------------------------
  ! Mise a jour de Schur d'un bloc de la sous-matrice restante :
  !     A22 <- A22 - L*U
  ! avec L (MI x KW, decrit par LBLK) et U (KW x NJ, decrit par UBLK),
  ! chacun dense OU basse-rang -- 4 cas geres explicitement, chacun choisi
  ! pour minimiser le nombre d'operations quand un ou deux facteurs sont
  ! basse-rang (c'est le coeur du gain de performance BLR).
  !-----------------------------------------------------------------------
  SUBROUTINE SCHUR_UPDATE_BLOCK(MI, NJ, KW, LBLK, UBLK, A22, LDA)
    INTEGER, INTENT(IN)             :: MI, NJ, KW, LDA
    TYPE(BLR_BLOCK_T), INTENT(IN)   :: LBLK, UBLK
    COMPLEX(WP), INTENT(INOUT)      :: A22(LDA,NJ)

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: T(:,:), T1(:,:), T2(:,:)

    IF (MI.EQ.0 .OR. NJ.EQ.0 .OR. KW.EQ.0) RETURN

    IF (.NOT.LBLK%IS_LOWRANK .AND. .NOT.UBLK%IS_LOWRANK) THEN
      ! Dense x Dense : GEMM classique (identique a l'etape 1)
      CALL CGEMM('No transpose', 'No transpose', MI, NJ, KW, -CONE, &
                  LBLK%DENSE, MI, UBLK%DENSE, KW, CONE, A22, LDA)

    ELSE IF (LBLK%IS_LOWRANK .AND. .NOT.UBLK%IS_LOWRANK) THEN
      ! (UL*VL) x U_dense : T = VL*U_dense (RxNJ), A22 -= UL*T
      ALLOCATE(T(LBLK%RANK,NJ))
      CALL CGEMM('No transpose', 'No transpose', LBLK%RANK, NJ, KW, &
                  CONE, LBLK%VFAC, LBLK%RANK, UBLK%DENSE, KW, CZERO, T, LBLK%RANK)
      CALL CGEMM('No transpose', 'No transpose', MI, NJ, LBLK%RANK, &
                  -CONE, LBLK%UFAC, MI, T, LBLK%RANK, CONE, A22, LDA)
      DEALLOCATE(T)

    ELSE IF (.NOT.LBLK%IS_LOWRANK .AND. UBLK%IS_LOWRANK) THEN
      ! L_dense x (UU*VU) : T = L_dense*UU (MIxR), A22 -= T*VU
      ALLOCATE(T(MI,UBLK%RANK))
      CALL CGEMM('No transpose', 'No transpose', MI, UBLK%RANK, KW, &
                  CONE, LBLK%DENSE, MI, UBLK%UFAC, KW, CZERO, T, MI)
      CALL CGEMM('No transpose', 'No transpose', MI, NJ, UBLK%RANK, &
                  -CONE, T, MI, UBLK%VFAC, UBLK%RANK, CONE, A22, LDA)
      DEALLOCATE(T)

    ELSE
      ! (UL*VL) x (UU*VU) : T1 = VL*UU (R1xR2, PETIT), T2 = T1*VU (R1xNJ),
      ! A22 -= UL*T2 -- jamais de matrice KW-large formee explicitement.
      ALLOCATE(T1(LBLK%RANK,UBLK%RANK), T2(LBLK%RANK,NJ))
      CALL CGEMM('No transpose', 'No transpose', LBLK%RANK, UBLK%RANK, KW, &
                  CONE, LBLK%VFAC, LBLK%RANK, UBLK%UFAC, KW, CZERO, T1, LBLK%RANK)
      CALL CGEMM('No transpose', 'No transpose', LBLK%RANK, NJ, UBLK%RANK, &
                  CONE, T1, LBLK%RANK, UBLK%VFAC, UBLK%RANK, CZERO, T2, LBLK%RANK)
      CALL CGEMM('No transpose', 'No transpose', MI, NJ, LBLK%RANK, &
                  -CONE, LBLK%UFAC, MI, T2, LBLK%RANK, CONE, A22, LDA)
      DEALLOCATE(T1, T2)
    END IF
  END SUBROUTINE SCHUR_UPDATE_BLOCK

  !-----------------------------------------------------------------------
  ! OPTIMISATION (suite au constat des ~835s de factorisation a N=10000
  ! dans blr_matrix_mod.f90, causees par ~NBLK^3/3 SVD MIxNJ COMPLETES --
  ! une par mise a jour de Schur hors-diagonale, via l'ancien schema
  ! EXPAND -> SCHUR_UPDATE_BLOCK -> COMPRESS_BLOCK systematique) :
  !
  ! PRINCIPE : le produit L*U (L=LBLK MIxKW, U=UBLK KWxNJ) admet DEJA une
  ! representation basse-rang bon marche, SANS jamais former le MIxNJ
  ! complet, des que L OU U est basse-rang (COMPUTE_LU_PRODUCT -- seul le
  ! cas dense x dense force un produit plein, cout MIxKWxNJ inevitable).
  ! Si en plus le bloc A cible est LUI-MEME basse-rang, alors
  ! A - L*U = concatenation de deux termes basse-rang (rang RA+RP), qu'on
  ! recompresse a bas cout via QR + SVD du "coeur" (RA+RP)x(RA+RP) SEUL
  ! (FAST_RECOMPRESS) -- jamais de SVD MIxNJ.
  !
  ! DEUX CAS geres par SCHUR_UPDATE_BLR_BLOCK, comme demande :
  !   1) A OU L*U non basse-rang (ou rang recompresse trop grand par
  !      rapport a MIN(MI,NJ), voir FAST_RECOMPRESS/RMAX) -> bascule en
  !      dense (expand + soustraction, LR-aware quand possible cote
  !      SUBTRACT_PRODUCT_DENSE).
  !   2) A ET L*U basse-rang -> reste en low-rank (concatenation des
  !      facteurs + compression rapide).
  !-----------------------------------------------------------------------
  SUBROUTINE COMPUTE_LU_PRODUCT(MI, NJ, KW, LBLK, UBLK, PBLK)
    INTEGER, INTENT(IN)             :: MI, NJ, KW
    TYPE(BLR_BLOCK_T), INTENT(IN)   :: LBLK, UBLK
    TYPE(BLR_BLOCK_T), INTENT(OUT)  :: PBLK
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: T(:,:)

    PBLK%M = MI
    PBLK%N = NJ

    IF (.NOT.LBLK%IS_LOWRANK .AND. .NOT.UBLK%IS_LOWRANK) THEN
      ! dense x dense : produit force dense (GEMM plein, cout inevitable)
      PBLK%IS_LOWRANK = .FALSE.
      ALLOCATE(PBLK%DENSE(MI,NJ))
      CALL CGEMM('No transpose', 'No transpose', MI, NJ, KW, CONE, &
                  LBLK%DENSE, MI, UBLK%DENSE, KW, CZERO, PBLK%DENSE, MI)

    ELSE IF (LBLK%IS_LOWRANK .AND. .NOT.UBLK%IS_LOWRANK) THEN
      ! (UL*VL) x U_dense : produit basse-rang de rang RL, T=VL*U_dense
      PBLK%IS_LOWRANK = .TRUE.
      PBLK%RANK = LBLK%RANK
      ALLOCATE(PBLK%UFAC(MI,LBLK%RANK), PBLK%VFAC(LBLK%RANK,NJ))
      PBLK%UFAC = LBLK%UFAC
      CALL CGEMM('No transpose', 'No transpose', LBLK%RANK, NJ, KW, CONE, &
                  LBLK%VFAC, LBLK%RANK, UBLK%DENSE, KW, CZERO, PBLK%VFAC, LBLK%RANK)

    ELSE IF (.NOT.LBLK%IS_LOWRANK .AND. UBLK%IS_LOWRANK) THEN
      ! L_dense x (UU*VU) : produit basse-rang de rang RU, T=L_dense*UU
      PBLK%IS_LOWRANK = .TRUE.
      PBLK%RANK = UBLK%RANK
      ALLOCATE(PBLK%UFAC(MI,UBLK%RANK), PBLK%VFAC(UBLK%RANK,NJ))
      CALL CGEMM('No transpose', 'No transpose', MI, UBLK%RANK, KW, CONE, &
                  LBLK%DENSE, MI, UBLK%UFAC, KW, CZERO, PBLK%UFAC, MI)
      PBLK%VFAC = UBLK%VFAC

    ELSE
      ! (UL*VL) x (UU*VU) : produit basse-rang, rang = min(RL,RU) choisi
      ! pour minimiser le cout (jamais de matrice KW-large formee)
      PBLK%IS_LOWRANK = .TRUE.
      ALLOCATE(T(LBLK%RANK,UBLK%RANK))
      CALL CGEMM('No transpose', 'No transpose', LBLK%RANK, UBLK%RANK, KW, CONE, &
                  LBLK%VFAC, LBLK%RANK, UBLK%UFAC, KW, CZERO, T, LBLK%RANK)
      IF (LBLK%RANK .LE. UBLK%RANK) THEN
        PBLK%RANK = LBLK%RANK
        ALLOCATE(PBLK%UFAC(MI,LBLK%RANK), PBLK%VFAC(LBLK%RANK,NJ))
        PBLK%UFAC = LBLK%UFAC
        CALL CGEMM('No transpose', 'No transpose', LBLK%RANK, NJ, UBLK%RANK, CONE, &
                    T, LBLK%RANK, UBLK%VFAC, UBLK%RANK, CZERO, PBLK%VFAC, LBLK%RANK)
      ELSE
        PBLK%RANK = UBLK%RANK
        ALLOCATE(PBLK%UFAC(MI,UBLK%RANK), PBLK%VFAC(UBLK%RANK,NJ))
        PBLK%VFAC = UBLK%VFAC
        CALL CGEMM('No transpose', 'No transpose', MI, UBLK%RANK, LBLK%RANK, CONE, &
                    LBLK%UFAC, MI, T, LBLK%RANK, CZERO, PBLK%UFAC, MI)
      END IF
      DEALLOCATE(T)
    END IF
  END SUBROUTINE COMPUTE_LU_PRODUCT

  !-----------------------------------------------------------------------
  ! Applique a BLK la MEME sequence de permutations de lignes que celle
  ! produite par un CGETRF sur les KW premieres lignes de son bloc-ligne
  ! (convention LAPACK/CLASWP : pour I=1..KW, si IPIV(I).NE.I, echanger la
  ! ligne I et la ligne IPIV(I) -- PAS un vecteur de permutation a appliquer
  ! d'un coup, une SEQUENCE ordonnee de swaps).
  !
  ! POURQUOI CA RESTE BON MARCHE EN BLR (le point cle du pivot confine au
  ! bloc diagonal, cf. discussion avec l'utilisateur) : permuter des LIGNES
  ! d'un bloc bas-rang A~=UFAC*VFAC ne touche QUE UFAC (VFAC n'a aucune
  ! notion de ligne) -- cout O(rang) par swap, JAMAIS de decompression.
  ! Cas dense : swap direct sur DENSE, cout O(N) par swap (standard).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_BLOCK_ROW_SWAP(BLK, KW, IPIV)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN) :: KW, IPIV(KW)
    INTEGER :: I
    IF (BLK%M .EQ. 0 .OR. BLK%N .EQ. 0) RETURN
    IF (BLK%IS_LOWRANK) THEN
      DO I = 1, KW
        IF (IPIV(I) .NE. I) CALL CSWAP(BLK%RANK, BLK%UFAC(I,1), KW, BLK%UFAC(IPIV(I),1), KW)
      END DO
    ELSE
      DO I = 1, KW
        IF (IPIV(I) .NE. I) CALL CSWAP(BLK%N, BLK%DENSE(I,1), KW, BLK%DENSE(IPIV(I),1), KW)
      END DO
    END IF
  END SUBROUTINE BLR_BLOCK_ROW_SWAP

  !-----------------------------------------------------------------------
  ! A(MIxNJ) -= PBLK (dense ou basse-rang) -- jamais besoin d'expanser
  ! PBLK si low-rank (CGEMM direct avec ses facteurs).
  !-----------------------------------------------------------------------
  SUBROUTINE SUBTRACT_PRODUCT_DENSE(MI, NJ, PBLK, A)
    INTEGER, INTENT(IN)           :: MI, NJ
    TYPE(BLR_BLOCK_T), INTENT(IN) :: PBLK
    COMPLEX(WP), INTENT(INOUT)    :: A(MI,NJ)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)

    IF (PBLK%IS_LOWRANK) THEN
      CALL CGEMM('No transpose', 'No transpose', MI, NJ, PBLK%RANK, -CONE, &
                  PBLK%UFAC, MI, PBLK%VFAC, PBLK%RANK, CONE, A, MI)
    ELSE
      A = A - PBLK%DENSE
    END IF
  END SUBROUTINE SUBTRACT_PRODUCT_DENSE

  !-----------------------------------------------------------------------
  ! Recompresse C = UA*VA - PU*PV (MIxNJ, rang RA+RP en entree) SANS
  ! jamais former le MIxNJ complet : QR de [UA|-PU] (MIx(RA+RP)) et QR de
  ! la transposee de [VA;PV] (NJx(RA+RP)), puis SVD COMPLETE du "coeur"
  ! (RA+RP)x(RA+RP) SEULEMENT -- typiquement quelques dizaines de lignes/
  ! colonnes, donc bien moins cher qu'une SVD MIxNJ complete.
  !
  ! Note sur la transposee : Vcat = [VA;PV] n'est jamais formee comme
  ! telle -- on construit directement VCATT = transpose(Vcat) (NJxRTOT),
  ! on en prend la QR (VCATT = QV*RV), d'ou Vcat = RV^T * QV^T par simple
  ! transposition de l'egalite -- transposee ORDINAIRE (CGEMM 'Transpose'),
  ! PAS conjuguee, puisqu'aucune conjugaison n'intervient dans cette
  ! derivation purement algebrique.
  !
  ! OK=.FALSE. si le rang tronque a TOL depasse RMAX (l'appelant doit
  ! alors passer le bloc en dense) OU si RA+RP >= MIN(MI,NJ) (recompression
  ! sans objet -- sortie immediate avant tout calcul).
  !-----------------------------------------------------------------------
  ! WS (optionnel) : workspace pre-alloue reutilise pour TOUS les tableaux
  ! de travail (QR + petite SVD) -- AUCUNE allocation dans cette routine si
  ! fourni (RA+RP doivent etre <= WS%MAXDIM, deja garanti par le test RMAX
  ! ci-dessous des lors que RMAX <= WS%MAXDIM). Absent : comportement
  ! d'origine (allocation/desallocation locale a chaque appel).
  SUBROUTINE FAST_RECOMPRESS(MI, NJ, RA, UA, VA, RP, PU, PV, TOL, RMAX, NEWBLK, OK, WS)
    INTEGER, INTENT(IN)            :: MI, NJ, RA, RP, RMAX
    COMPLEX(WP), INTENT(IN)        :: UA(MI,RA), VA(RA,NJ), PU(MI,RP), PV(RP,NJ)
    REAL(WP), INTENT(IN)           :: TOL
    TYPE(BLR_BLOCK_T), INTENT(OUT) :: NEWBLK
    LOGICAL, INTENT(OUT)           :: OK
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS

    INTEGER :: RTOT, R, I, LWORK, INFO, LDU, LDVQ, LDSQ, LDSV
    COMPLEX(WP), ALLOCATABLE, TARGET :: UCAT_LOC(:,:), VCATT_LOC(:,:), TAUU_LOC(:), TAUV_LOC(:)
    COMPLEX(WP), ALLOCATABLE, TARGET :: RU_LOC(:,:), RV_LOC(:,:), W_LOC(:,:), UW_LOC(:,:), VWH_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE, TARGET :: SVWH_LOC(:,:), CWORK_LOC(:)
    REAL(WP), ALLOCATABLE, TARGET    :: SW_LOC(:), RWORK_LOC(:)
    COMPLEX(WP), POINTER :: UCAT(:,:), VCATT(:,:), TAUU(:), TAUV(:)
    COMPLEX(WP), POINTER :: RU(:,:), RV(:,:), WMAT(:,:), UW(:,:), VWH(:,:), SVWH(:,:), CWORK(:)
    REAL(WP), POINTER :: SW(:), RWORK(:)
    REAL(WP) :: THRESH
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)

    RTOT = RA + RP
    IF (RTOT .GE. MIN(MI,NJ) .OR. RTOT .GT. RMAX) THEN
      OK = .FALSE.
      RETURN
    END IF

    IF (PRESENT(WS)) THEN
      LDU = WS%MAXDIM; LDVQ = WS%MAXDIM; LDSQ = WS%MAXDIM
      UCAT => WS%UCAT; VCATT => WS%VCATT; TAUU => WS%TAUU; TAUV => WS%TAUV
      RU => WS%RU; RV => WS%RV; WMAT => WS%WMAT; UW => WS%UW; VWH => WS%VWH
      SW => WS%SW; RWORK => WS%RWORK2; CWORK => WS%CWORK2
      LWORK = SIZE(CWORK)
    ELSE
      LDU = MI; LDVQ = NJ; LDSQ = RTOT
      ALLOCATE(UCAT_LOC(MI,RTOT), VCATT_LOC(NJ,RTOT), TAUU_LOC(RTOT), TAUV_LOC(RTOT))
      ALLOCATE(RU_LOC(RTOT,RTOT), RV_LOC(RTOT,RTOT), W_LOC(RTOT,RTOT), UW_LOC(RTOT,RTOT), VWH_LOC(RTOT,RTOT))
      ALLOCATE(SW_LOC(RTOT), RWORK_LOC(5*RTOT))
      UCAT => UCAT_LOC; VCATT => VCATT_LOC; TAUU => TAUU_LOC; TAUV => TAUV_LOC
      RU => RU_LOC; RV => RV_LOC; WMAT => W_LOC; UW => UW_LOC; VWH => VWH_LOC
      SW => SW_LOC; RWORK => RWORK_LOC
      ALLOCATE(CWORK_LOC(1))
      LWORK = -1
      CALL CGEQRF(MI, RTOT, UCAT, LDU, TAUU, CWORK_LOC, LWORK, INFO)
      LWORK = INT(REAL(CWORK_LOC(1)))
      CALL CUNGQR(MI, RTOT, RTOT, UCAT, LDU, TAUU, CWORK_LOC, -1, INFO)
      LWORK = MAX(LWORK, INT(REAL(CWORK_LOC(1))))
      CALL CGEQRF(NJ, RTOT, VCATT, LDVQ, TAUV, CWORK_LOC, -1, INFO)
      LWORK = MAX(LWORK, INT(REAL(CWORK_LOC(1))))
      CALL CUNGQR(NJ, RTOT, RTOT, VCATT, LDVQ, TAUV, CWORK_LOC, -1, INFO)
      LWORK = MAX(LWORK, INT(REAL(CWORK_LOC(1))))
      CALL CGESVD('S', 'S', RTOT, RTOT, WMAT, LDSQ, SW, UW, LDSQ, VWH, LDSQ, CWORK_LOC, -1, RWORK, INFO)
      LWORK = MAX(LWORK, INT(REAL(CWORK_LOC(1))))
      DEALLOCATE(CWORK_LOC)
      ALLOCATE(CWORK_LOC(MAX(LWORK,1)))
      CWORK => CWORK_LOC
    END IF

    UCAT(1:MI,1:RA)      = UA
    UCAT(1:MI,RA+1:RTOT) = -PU

    DO I = 1, RA
      VCATT(1:NJ,I) = VA(I,1:NJ)
    END DO
    DO I = 1, RP
      VCATT(1:NJ,RA+I) = PV(I,1:NJ)
    END DO

    ! QR de UCAT -> QU (ecrase UCAT), RU
    CALL CGEQRF(MI, RTOT, UCAT, LDU, TAUU, CWORK, LWORK, INFO)
    RU(1:RTOT,1:RTOT) = CZERO
    DO I = 1, RTOT
      RU(1:I,I) = UCAT(1:I,I)
    END DO
    CALL CUNGQR(MI, RTOT, RTOT, UCAT, LDU, TAUU, CWORK, LWORK, INFO)
    ! UCAT contient maintenant QU (MI x RTOT)

    ! QR de VCATT -> QV (ecrase VCATT), RV
    CALL CGEQRF(NJ, RTOT, VCATT, LDVQ, TAUV, CWORK, LWORK, INFO)
    RV(1:RTOT,1:RTOT) = CZERO
    DO I = 1, RTOT
      RV(1:I,I) = VCATT(1:I,I)
    END DO
    CALL CUNGQR(NJ, RTOT, RTOT, VCATT, LDVQ, TAUV, CWORK, LWORK, INFO)
    ! VCATT contient maintenant QV (NJ x RTOT)

    ! W = RU * RV^T (transposee ORDINAIRE -- voir note ci-dessus)
    CALL CGEMM('No transpose', 'Transpose', RTOT, RTOT, RTOT, CONE, &
                RU, LDSQ, RV, LDSQ, CZERO, WMAT, LDSQ)

    CALL CGESVD('S', 'S', RTOT, RTOT, WMAT, LDSQ, SW, UW, LDSQ, VWH, LDSQ, CWORK, LWORK, RWORK, INFO)

    IF (INFO .NE. 0) THEN
      OK = .FALSE.
      IF (.NOT. PRESENT(WS)) DEALLOCATE(UCAT_LOC, VCATT_LOC, TAUU_LOC, TAUV_LOC, RU_LOC, RV_LOC, &
                                         W_LOC, UW_LOC, VWH_LOC, SW_LOC, RWORK_LOC, CWORK_LOC)
      RETURN
    END IF

    THRESH = TOL * MAX(SW(1), TINY(1.0_WP))
    R = 0
    DO I = 1, RTOT
      IF (SW(I) .GT. THRESH) R = R + 1
    END DO
    R = MAX(R, 1)

    IF (R .GT. RMAX) THEN
      OK = .FALSE.
    ELSE
      OK = .TRUE.
      NEWBLK%IS_LOWRANK = .TRUE.
      NEWBLK%M = MI; NEWBLK%N = NJ; NEWBLK%RANK = R
      ALLOCATE(NEWBLK%UFAC(MI,R), NEWBLK%VFAC(R,NJ))
      CALL CGEMM('No transpose', 'No transpose', MI, R, RTOT, CONE, &
                  UCAT, LDU, UW, LDSQ, CZERO, NEWBLK%UFAC, MI)

      IF (PRESENT(WS)) THEN
        LDSV = WS%MAXDIM
        SVWH => WS%SVWH
      ELSE
        LDSV = R
        ALLOCATE(SVWH_LOC(R,RTOT))
        SVWH => SVWH_LOC
      END IF
      DO I = 1, R
        SVWH(I,1:RTOT) = SW(I) * VWH(I,1:RTOT)
      END DO
      CALL CGEMM('No transpose', 'Transpose', R, NJ, RTOT, CONE, &
                  SVWH, LDSV, VCATT, LDVQ, CZERO, NEWBLK%VFAC, R)
      IF (.NOT. PRESENT(WS)) DEALLOCATE(SVWH_LOC)
    END IF

    IF (.NOT. PRESENT(WS)) DEALLOCATE(UCAT_LOC, VCATT_LOC, TAUU_LOC, TAUV_LOC, RU_LOC, RV_LOC, &
                                       W_LOC, UW_LOC, VWH_LOC, SW_LOC, RWORK_LOC, CWORK_LOC)
  END SUBROUTINE FAST_RECOMPRESS

  !-----------------------------------------------------------------------
  ! Point d'entree : ABLK <- ABLK - L*U, ABLK restant COMPRESSE si possible
  ! (voir le commentaire de tete de bloc ci-dessus pour les 2 cas geres).
  !
  ! T_PRODUCT/T_RECOMP (optionnels, INTENT(INOUT), CUMULES par l'appelant --
  ! initialiser a 0 avant la premiere boucle) : instrumentation demandee
  ! pour distinguer, a l'interieur d'UNE mise a jour de Schur, le temps
  ! passe a CALCULER le produit L*U ("Schur sans compression" -- inclut
  ! aussi le repli dense EXPAND_BLOCK/SUBTRACT_PRODUCT_DENSE quand il a
  ! lieu) du temps passe a le RECOMPRESSER (FAST_RECOMPRESS seul). Non
  ! fournis (defaut) : zero surcout, comportement inchange -- sans danger
  ! pour les appelants existants (blr_matrix_mod.f90 sous !$OMP PARALLEL DO,
  ! qui ne les passe pas).
  !-----------------------------------------------------------------------
  ! USE_MOVE_ALLOC (optionnel, defaut .FALSE. = comportement INCHANGE) :
  ! DIAGNOSTIC/OPTIMISATION EN COURS D'EVALUATION (voir la memoire du
  ! projet -- ecart RSS-vs-theorique de MPI_BLR_MATRIX_LU_FACTOR_PIVOT).
  ! Par defaut, "ABLK = NEWBLK" (assignation Fortran de type derive a
  ! composantes ALLOCATABLE) declenche une REALLOCATION IMPLICITE de
  ! ABLK%UFAC/VFAC a la forme de NEWBLK PUIS une COPIE -- meme si NEWBLK
  ! vient d'etre alloue a la bonne taille, c'est un DEUXIEME cycle
  ! allocation/copie en plus du premier (dans FAST_RECOMPRESS). Avec
  ! USE_MOVE_ALLOC=.TRUE., on desalloue explicitement l'ancien contenu de
  ! ABLK puis on TRANSFERE la propriete des tableaux de NEWBLK/DENSE_TMP
  ! vers ABLK via MOVE_ALLOC (aucune copie, aucune reallocation
  ! supplementaire) -- resultat numerique IDENTIQUE, seul le nombre de
  ! cycles allocate/deallocate change.
  SUBROUTINE SCHUR_UPDATE_BLR_BLOCK(ABLK, MI, NJ, KW, LBLK, UBLK, TOL, WS, T_PRODUCT, T_RECOMP, USE_MOVE_ALLOC)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: ABLK
    INTEGER, INTENT(IN)              :: MI, NJ, KW
    TYPE(BLR_BLOCK_T), INTENT(IN)    :: LBLK, UBLK
    REAL(WP), INTENT(IN)             :: TOL
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: T_PRODUCT, T_RECOMP
    LOGICAL, INTENT(IN), OPTIONAL     :: USE_MOVE_ALLOC

    TYPE(BLR_BLOCK_T) :: PBLK, NEWBLK, DENSEBLK
    COMPLEX(WP), ALLOCATABLE :: DENSE_TMP(:,:)
    INTEGER :: RMAX
    LOGICAL :: OK, MOVEALLOC_ON
    DOUBLE PRECISION :: TT0

    IF (MI.EQ.0 .OR. NJ.EQ.0 .OR. KW.EQ.0) RETURN

    MOVEALLOC_ON = .FALSE.
    IF (PRESENT(USE_MOVE_ALLOC)) MOVEALLOC_ON = USE_MOVE_ALLOC

    IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
    CALL COMPUTE_LU_PRODUCT(MI, NJ, KW, LBLK, UBLK, PBLK)
    IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)

    IF (ABLK%IS_LOWRANK .AND. PBLK%IS_LOWRANK) THEN
      ! cas 2 : A et L*U basse-rang -> tenter de rester en low-rank
      RMAX = MAX(1, INT(COMPRESS_WORTHWHILE_RATIO*REAL(MIN(MI,NJ),WP)) - 1)
      IF (PRESENT(T_RECOMP)) TT0 = OMP_GET_WTIME()
      IF (PRESENT(WS)) THEN
        CALL FAST_RECOMPRESS(MI, NJ, ABLK%RANK, ABLK%UFAC, ABLK%VFAC, &
                              PBLK%RANK, PBLK%UFAC, PBLK%VFAC, TOL, RMAX, NEWBLK, OK, WS)
      ELSE
        CALL FAST_RECOMPRESS(MI, NJ, ABLK%RANK, ABLK%UFAC, ABLK%VFAC, &
                              PBLK%RANK, PBLK%UFAC, PBLK%VFAC, TOL, RMAX, NEWBLK, OK)
      END IF
      IF (PRESENT(T_RECOMP)) T_RECOMP = T_RECOMP + (OMP_GET_WTIME() - TT0)
      IF (OK) THEN
        IF (MOVEALLOC_ON) THEN
          ! transfert de propriete (NEWBLK%UFAC/VFAC deja alloues a la
          ! bonne taille dans FAST_RECOMPRESS) -- AUCUNE copie, AUCUNE
          ! reallocation supplementaire, contrairement a "ABLK = NEWBLK".
          IF (ALLOCATED(ABLK%UFAC))  DEALLOCATE(ABLK%UFAC)
          IF (ALLOCATED(ABLK%VFAC))  DEALLOCATE(ABLK%VFAC)
          IF (ALLOCATED(ABLK%DENSE)) DEALLOCATE(ABLK%DENSE)
          CALL MOVE_ALLOC(NEWBLK%UFAC, ABLK%UFAC)
          CALL MOVE_ALLOC(NEWBLK%VFAC, ABLK%VFAC)
          ABLK%IS_LOWRANK = .TRUE.
          ABLK%M = MI; ABLK%N = NJ; ABLK%RANK = NEWBLK%RANK
        ELSE
          ABLK = NEWBLK
        END IF
        RETURN
      END IF
      ! rang recompresse trop grand -> tombe dans le repli dense ci-dessous
    END IF

    ! cas 1 (ou repli du cas 2) : A ou L*U non basse-rang -> dense
    IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
    ALLOCATE(DENSE_TMP(MI,NJ))
    CALL EXPAND_BLOCK(ABLK, DENSE_TMP, MI)
    CALL SUBTRACT_PRODUCT_DENSE(MI, NJ, PBLK, DENSE_TMP)
    IF (MOVEALLOC_ON) THEN
      ! meme principe : DENSE_TMP est deja a la bonne taille, transfert
      ! direct vers ABLK%DENSE au lieu de passer par un DENSEBLK
      ! intermediaire puis "ABLK = DENSEBLK".
      IF (ALLOCATED(ABLK%UFAC))  DEALLOCATE(ABLK%UFAC)
      IF (ALLOCATED(ABLK%VFAC))  DEALLOCATE(ABLK%VFAC)
      IF (ALLOCATED(ABLK%DENSE)) DEALLOCATE(ABLK%DENSE)
      CALL MOVE_ALLOC(DENSE_TMP, ABLK%DENSE)
      ABLK%IS_LOWRANK = .FALSE.
      ABLK%M = MI; ABLK%N = NJ
    ELSE
      DENSEBLK%IS_LOWRANK = .FALSE.
      DENSEBLK%M = MI; DENSEBLK%N = NJ
      CALL MOVE_ALLOC(DENSE_TMP, DENSEBLK%DENSE)
      ABLK = DENSEBLK
    END IF
    IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)
  END SUBROUTINE SCHUR_UPDATE_BLR_BLOCK

  !-----------------------------------------------------------------------
  ! PROPOSITION 2 (post-etape 6) : Y <- ALPHA*BLK*X + BETA*Y, ou BLK est
  ! un bloc M x N (dense OU basse-rang) et X est N x NRHS. Utilise pour la
  ! substitution par blocs (avant/arriere) d'une factorisation dont les
  ! facteurs L21/U12 sont GARDES compresses en memoire (pas seulement
  ! compresses le temps de la mise a jour de Schur, comme a l'etape 2/6e)
  ! -- voir blr_lu_storage_mod.f90. Cas basse-rang : Y <- ALPHA*UFAC*
  ! (VFAC*X) + BETA*Y (deux CGEMM avec une matrice intermediaire RANG x
  ! NRHS, PETITE -- jamais besoin de reconstruire le bloc dense complet,
  ! c'est exactement le gain memoire ET calcul vise par cette proposition).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_BLOCK_APPLY(BLK, NRHS, ALPHA, X, LDX, BETA, Y, LDY)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: BLK
    INTEGER, INTENT(IN)           :: NRHS, LDX, LDY
    COMPLEX(WP), INTENT(IN)       :: ALPHA, BETA
    COMPLEX(WP), INTENT(IN)       :: X(LDX,NRHS)
    COMPLEX(WP), INTENT(INOUT)    :: Y(LDY,NRHS)

    COMPLEX(WP), PARAMETER :: CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: T(:,:)

    IF (BLK%M.EQ.0 .OR. BLK%N.EQ.0) RETURN

    IF (BLK%IS_LOWRANK) THEN
      ALLOCATE(T(BLK%RANK,NRHS))
      CALL CGEMM('No transpose', 'No transpose', BLK%RANK, NRHS, BLK%N, &
                  (1.0_WP,0.0_WP), BLK%VFAC, BLK%RANK, X, LDX, CZERO, T, BLK%RANK)
      CALL CGEMM('No transpose', 'No transpose', BLK%M, NRHS, BLK%RANK, &
                  ALPHA, BLK%UFAC, BLK%M, T, BLK%RANK, BETA, Y, LDY)
      DEALLOCATE(T)
    ELSE
      CALL CGEMM('No transpose', 'No transpose', BLK%M, NRHS, BLK%N, &
                  ALPHA, BLK%DENSE, BLK%M, X, LDX, BETA, Y, LDY)
    END IF
  END SUBROUTINE BLR_BLOCK_APPLY

  !-----------------------------------------------------------------------
  ! Decoupage de [1:N] en blocs de taille reguliere NB (dernier bloc
  ! eventuellement plus petit). BSTART/BEND sont alloues ici. Consolide
  ! depuis blr_lu_mod.f90 ET blr_matrix_mod.f90 (code identique duplique
  ! dans les deux avant cette generalisation).
  !-----------------------------------------------------------------------
  SUBROUTINE BLOCK_PARTITION(N, NB, NBLK, BSTART, BEND)
    INTEGER, INTENT(IN)               :: N, NB
    INTEGER, INTENT(OUT)              :: NBLK
    INTEGER, ALLOCATABLE, INTENT(OUT) :: BSTART(:), BEND(:)
    INTEGER :: MAXBLK, K

    MAXBLK = (N + NB - 1) / NB
    ALLOCATE(BSTART(MAXBLK), BEND(MAXBLK))
    NBLK = 0
    K = 1
    DO WHILE (K .LE. N)
      NBLK = NBLK + 1
      BSTART(NBLK) = K
      BEND(NBLK)   = MIN(K+NB-1, N)
      K = BEND(NBLK) + 1
    END DO
  END SUBROUTINE BLOCK_PARTITION

  !-----------------------------------------------------------------------
  ! Empreinte memoire (octets) d'un bloc dense ou basse-rang. Consolide
  ! depuis blr_lu_storage_mod.f90 ET blr_matrix_mod.f90 (BLOCK_BYTES,
  ! code identique duplique dans les deux -- renomme BLR_BLOCK_BYTES ici
  ! pour eviter toute ambiguite si un futur fichier importe plusieurs
  ! modules definissant chacun un BLOCK_BYTES local).
  !-----------------------------------------------------------------------
  PURE FUNCTION BLR_BLOCK_BYTES(BLK) RESULT(NBYTES)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: BLK
    INTEGER(8) :: NBYTES
    IF (BLK%IS_LOWRANK) THEN
      NBYTES = 8_8 * (INT(BLK%M,8)*INT(BLK%RANK,8) + INT(BLK%RANK,8)*INT(BLK%N,8))
    ELSE
      NBYTES = 8_8 * INT(BLK%M,8) * INT(BLK%N,8)
    END IF
  END FUNCTION BLR_BLOCK_BYTES

  !-----------------------------------------------------------------------
  ! LU sans pivot, non blocuee (CSCAL+CGERU, LAPACK n'a pas de "GETRF sans
  ! pivot"). Consolide depuis blr_matrix_mod.f90 -- nommee AVEC prefixe
  ! BLR_ (et pas juste LU_NOPIV_UNBLOCKED) car mpi_lu_nopiv_mod.f90 exporte
  ! DEJA une routine de ce nom, et mpi_blr_lu_nopiv_mod.f90 importe les
  ! deux modules simultanement -- un nom identique aurait cree une
  ! ambiguite de resolution a la compilation.
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_LU_NOPIV_UNBLOCKED(KB, A, LDA, INFO)
    INTEGER, INTENT(IN)          :: KB, LDA
    COMPLEX(WP), INTENT(INOUT)   :: A(LDA,KB)
    INTEGER, INTENT(OUT)         :: INFO
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    REAL(WP), PARAMETER :: PIVTOL = 1.0E-6_WP
    INTEGER :: J
    INFO = 0
    DO J = 1, KB
      IF (ABS(A(J,J)) .LT. PIVTOL) THEN
        INFO = J
        RETURN
      END IF
      IF (J .LT. KB) THEN
        CALL CSCAL(KB-J, CONE/A(J,J), A(J+1,J), 1)
        CALL CGERU(KB-J, KB-J, -CONE, A(J+1,J), 1, A(J,J+1), LDA, A(J+1,J+1), LDA)
      END IF
    END DO
  END SUBROUTINE BLR_LU_NOPIV_UNBLOCKED

  !-----------------------------------------------------------------------
  ! LU sans pivot, BLOCKEE (right-looking, meme structure que BLOCK_LU_
  ! FACTOR/block_lu_mod.f90, mais SANS CGETRF/pivot -- ecrite pour remplacer
  ! BLR_LU_NOPIV_UNBLOCKED sur le bloc diagonal de mpi_blr_matrix_mod.f90,
  ! qui domine 40-50% du temps de factorisation a NB>=750 (voir memoire
  ! projet, decomposition fine 2026-08-27) : la version non blockee ne fait
  ! que du BLAS niveau 1/2 (CSCAL/CGERU), intrinsequement sequentiel et
  ! sans reutilisation de cache -- ici, seul le PANNEAU (INNER_NB colonnes,
  ! petit) reste non blocke (CALL a BLR_LU_NOPIV_UNBLOCKED ci-dessus), le
  ! GROS du travail (mise a jour de la sous-matrice restante) passe par
  ! CTRSM/CGEMM (BLAS niveau 3 -- threadable par MKL, bien meilleure
  ! intensite arithmetique).
  !
  ! DELIBEREMENT PAS CGETRF (LAPACK) : CGETRF fait du pivot partiel PAR
  ! CONSTRUCTION, ce qui exigerait de propager les permutations de lignes
  ! aux blocs L21/U12 de la meme ligne/colonne pour rester correct --
  ! exactement le travail de la Phase 2 ("pivot qui ne decompresse
  ! jamais") explicitement reportee par l'utilisateur. Cette version reste
  ! SANS PIVOT DE BOUT EN BOUT, seulement blockee -- meme resultat
  ! numerique que BLR_LU_NOPIV_UNBLOCKED (mêmes operations, juste
  ! reordonnees/groupees), pas une nouvelle methode.
  !
  ! INFO : indice (1-based, GLOBAL au bloc KBxKB) du premier pivot sous
  ! PIVTOL rencontre, 0 si aucun -- meme convention que la version non
  ! blockee.
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_LU_NOPIV_BLOCKED(KB, A, LDA, INFO)
    INTEGER, INTENT(IN)          :: KB, LDA
    COMPLEX(WP), INTENT(INOUT)   :: A(LDA,KB)
    INTEGER, INTENT(OUT)         :: INFO
    INTEGER, PARAMETER :: INNER_NB = 64
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    INTEGER :: J, JB, IINFO, NREST

    INFO = 0
    DO J = 1, KB, INNER_NB
      JB = MIN(INNER_NB, KB-J+1)

      ! Panneau diagonal JBxJB : trop petit pour valoir un BLAS niveau 3
      ! dedie -- reste non blocke (reutilise la routine existante).
      CALL BLR_LU_NOPIV_UNBLOCKED(JB, A(J,J), LDA, IINFO)
      IF (IINFO .NE. 0) THEN
        INFO = J - 1 + IINFO
        RETURN
      END IF

      NREST = KB - J - JB + 1
      IF (NREST .GT. 0) THEN
        ! U12 <- L11^-1 * A12 (TRSM, triangulaire inferieure UNITAIRE)
        CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', JB, NREST, CONE, &
                    A(J,J), LDA, A(J,J+JB), LDA)
        ! L21 <- A21 * U11^-1 (TRSM, triangulaire superieure NON unitaire)
        CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', NREST, JB, CONE, &
                    A(J,J), LDA, A(J+JB,J), LDA)
        ! A22 <- A22 - L21*U12 (GEMM, BLAS niveau 3 -- le gros du travail)
        CALL CGEMM('No transpose', 'No transpose', NREST, NREST, JB, -CONE, &
                    A(J+JB,J), LDA, A(J,J+JB), LDA, CONE, A(J+JB,J+JB), LDA)
      END IF
    END DO
  END SUBROUTINE BLR_LU_NOPIV_BLOCKED

  !-----------------------------------------------------------------------
  ! BLK <- L11^-1 * BLK, ou L11 (KWxKW, dense DKK) est unitaire triangulaire
  ! inferieure, et BLK est un panneau KWxN dense OU basse-rang. Cas
  ! basse-rang : resoudre sur le SEUL facteur UFAC (KWxRANG) donne
  ! directement le bon resultat (L11^-1*(UFAC*VFAC) = (L11^-1*UFAC)*VFAC) --
  ! un vrai gain de calcul (CTRSM sur KWxRANG au lieu de KWxN), pas une
  ! simplification. Consolide depuis blr_matrix_mod.f90 (TRANSFORM_U12).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_TRSM_LEFT_LOWER_UNIT(DKK, KW, BLK)
    INTEGER, INTENT(IN)             :: KW
    COMPLEX(WP), INTENT(IN)         :: DKK(KW,KW)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, BLK%RANK, CONE, DKK, KW, BLK%UFAC, KW)
    ELSE
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, BLK%N, CONE, DKK, KW, BLK%DENSE, KW)
    END IF
  END SUBROUTINE BLR_TRSM_LEFT_LOWER_UNIT

  !-----------------------------------------------------------------------
  ! BLK <- BLK * U11^-1, ou U11 (KWxKW, dense DKK) est triangulaire
  ! superieure non-unitaire, et BLK est un panneau MxKW dense OU
  ! basse-rang. Cas basse-rang : resoudre sur le SEUL facteur VFAC
  ! (RANGxKW). Consolide depuis blr_matrix_mod.f90 (TRANSFORM_L21).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_TRSM_RIGHT_UPPER_NONUNIT(BLK, DKK, KW)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN)             :: KW
    COMPLEX(WP), INTENT(IN)         :: DKK(KW,KW)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', BLK%RANK, KW, CONE, DKK, KW, BLK%VFAC, BLK%RANK)
    ELSE
      CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', BLK%M, KW, CONE, DKK, KW, BLK%DENSE, BLK%M)
    END IF
  END SUBROUTINE BLR_TRSM_RIGHT_UPPER_NONUNIT

  !-----------------------------------------------------------------------
  ! Norme de Frobenius d'un bloc, SANS jamais l'expanser en dense si
  ! basse-rang : ||UV||_F^2 = trace((U^H U)(V V^H)), calcule sur deux
  ! matrices RANGxRANG (PETITES) plutot que sur le bloc MxN complet --
  ! meme philosophie que le reste de cette librairie (ne jamais payer un
  ! cout MxN quand un cout RANGxRANG suffit). NOUVELLE primitive (n'existait
  ! sous aucune forme avant cette generalisation).
  !-----------------------------------------------------------------------
  FUNCTION BLR_BLOCK_NORM(BLK) RESULT(NRM)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: BLK
    REAL(WP) :: NRM
    COMPLEX(WP), ALLOCATABLE :: G(:,:), H(:,:)
    COMPLEX(WP) :: TR
    INTEGER :: I, J

    IF (.NOT. BLK%IS_LOWRANK) THEN
      NRM = SQRT(SUM(ABS(BLK%DENSE)**2))
      RETURN
    END IF

    ALLOCATE(G(BLK%RANK,BLK%RANK), H(BLK%RANK,BLK%RANK))
    CALL CGEMM('Conjugate transpose', 'No transpose', BLK%RANK, BLK%RANK, BLK%M, &
                (1.0_WP,0.0_WP), BLK%UFAC, BLK%M, BLK%UFAC, BLK%M, (0.0_WP,0.0_WP), G, BLK%RANK)
    CALL CGEMM('No transpose', 'Conjugate transpose', BLK%RANK, BLK%RANK, BLK%N, &
                (1.0_WP,0.0_WP), BLK%VFAC, BLK%RANK, BLK%VFAC, BLK%RANK, (0.0_WP,0.0_WP), H, BLK%RANK)

    TR = (0.0_WP,0.0_WP)
    DO J = 1, BLK%RANK
      DO I = 1, BLK%RANK
        TR = TR + G(I,J)*H(J,I)
      END DO
    END DO
    NRM = SQRT(MAX(REAL(TR,WP), 0.0_WP))
    DEALLOCATE(G, H)
  END FUNCTION BLR_BLOCK_NORM

END MODULE BLR_COMPRESS_MOD
