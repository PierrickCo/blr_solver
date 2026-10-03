!===============================================================================
! Bibliotheque d'algebre lineaire par blocs full/low-rank : primitives de
! compression basse-rang d'un bloc dense (complexe simple precision) et
! catalogue complet des operations BLR/LU qui s'appuient dessus.
!
! METHODE DE COMPRESSION : SVD tronquee (CGESVD), pas QR avec pivotage de
! colonnes (ce que font MUMPS/STRUMPACK en production, plus rapide) --
! choix delibere : la SVD donne directement le rang numerique (a une
! tolerance donnee) sans bookkeeping de permutation de colonnes a gerer.
!
! REPRESENTATION : un bloc M x N est stocke soit dense (tel quel), soit
! basse-rang comme A ~= UFAC * VFAC avec UFAC (MxR) et VFAC (RxN) -- le
! facteur singulier Sigma est absorbe dans VFAC (VFAC(i,:) = S(i)*Vh(i,:)),
! donc UFAC*VFAC reconstruit directement l'approximation tronquee, sans
! manipulation supplementaire.
!
! SEUIL DE RANG : la compression n'a lieu que si le rang numerique trouve
! (a la tolerance TOL, relative a la plus grande valeur singuliere) est
! significativement plus petit que MIN(M,N) -- sinon elle n'apporte rien
! (voire coute plus cher a exploiter qu'un GEMM dense direct) et le bloc
! reste dense.
!
! ROLE DE CE MODULE : point d'entree unique pour l'algebre bloc full/low-
! rank du solveur -- chaque routine choisit systematiquement le chemin de
! calcul le moins cher selon que ses operandes sont denses ou basse-rang,
! et s'appuie sur les memes noyaux BLAS3/LAPACK threades par MKL (CGEMM,
! CTRSM, CGESVD, CGEQRF/CUNGQR) -- aucune reimplementation "maison" de ces
! operations, uniquement de la logique de dispatch et de recombinaison des
! facteurs low-rank autour d'elles. Les routines partagees entre plusieurs
! solveurs (BLOCK_PARTITION, taille memoire d'un bloc, LU sans pivot, TRSM
! applique a un facteur low-rank seul) sont consolidees ici : tout module
! qui fait "USE BLR_COMPRESS_MOD" les recoit automatiquement.
!
! PISTES MKL NON ADOPTEES ICI (documentees, pas implementees -- eviter la
! sur-ingenierie tant que le besoin n'est pas mesure) : CGEMM3M (variante
! MKL de CGEMM, parfois plus rapide sur gros GEMM complexes via l'algorithme
! 3M) pourrait remplacer CGEMM en un seul point (ici) si un profilage futur
! le justifie ; le "batching" de petits GEMM (CGEMM_BATCH) pourrait accelerer
! les tres nombreuses petites operations de rang faible pendant une
! factorisation BLR, mais demande une refonte des boucles appelantes pour
! grouper les appels, hors scope de cette generalisation.
MODULE BLR_COMPRESS_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE OMP_LIB, ONLY: OMP_GET_WTIME
  USE BLR_ARENA_MOD
  IMPLICIT NONE

  ! Fraction de MIN(M,N) en-dessous de laquelle la compression est jugee
  ! "rentable". Expose en PARAMETER modifiable si besoin de retuner plus
  ! tard sans casser l'interface des routines.
  !
  ! 0.5 = seuil d'EQUILIBRE MEMOIRE bas-rang/dense (stockage bas-rang =
  ! R*(M+N), dense = M*N ; egaux pour R = M*N/(M+N), soit M/2 pour un
  ! bloc carre -- voir blr_compress_mod.f90/SCHUR_UPDATE_BLR_BLOCK_ARENA,
  ! meme raisonnement deja applique au plafond de croissance de l'arena).
  ! Un plafond a 0.8 (choix tune calcul : un GEMM bas-rang a rang R coute
  ! O(M*N*R), toujours moins cher qu'un GEMM dense O(M*N*NB) tant que
  ! R<NB, donc accepter un rang plus eleve que le seuil memoire reste, en
  ! theorie, un gain FLOPs) a ete mesure (bench/sweep_compressibility.sh)
  ! moins bon en pratique : la zone de compressibilite marginale (rang
  ! final proche de l'ancien RMAX, ex. ~121 sur 200) est la plus lente de
  ! tout le balayage (16.3s, plus de 2x le cas totalement dense ~7-8s) --
  ! l'algorithme tente la QR+SVD complete de FAST_RECOMPRESS pour un
  ! resultat qui, une fois accepte, coute de toute facon plus de memoire
  ! que le repli dense (R>M/2). 0.5 evite ce "pire des deux mondes" : le
  ! test rapide RTOT.GT.RMAX (avant la QR+SVD couteuse) rejette plus tot
  ! les blocs qui n'auraient de toute facon pas ete memoire-rentables.
  ! Variable de module (PAS PARAMETER) : valeur par defaut utilisee par
  ! TOUS les appelants existants (tous les benchmarks/tests, sauf
  ! test/bench_blr_mpi_pivot_ratio.f90 qui la surcharge explicitement pour
  ! son propre balayage), modifiable a l'execution via
  ! SET_COMPRESS_WORTHWHILE_RATIO ci-dessous.
  !
  ! 1/3 (mesure via bench/sweep_compress_ratio.sh au point marginal
  ! EPS=1e-6) : plus prudent que 0.5 (l'equilibre memoire bas-rang/dense
  ! strict, R=MIN(M,N)/2) -- accepte un rang encore plus bas avant de
  ! basculer en dense, ce qui reste dans le plateau "rapide" du temps de
  ! factorisation (~3.7-3.9s mesures, RATIO 0.1-0.5) au lieu de s'approcher
  ! de la zone ou le temps commence deja a grimper (~4.6s a 0.6, ~12.7-20s
  ! a 0.8-0.9) -- marge de securite deliberee sur le temps, au prix de
  ! quelques opportunites de compression memoire-rentables (rang entre M/3
  ! et M/2) rejetees.
  REAL(WP) :: COMPRESS_WORTHWHILE_RATIO = 1.0_WP/3.0_WP

  TYPE :: BLR_BLOCK_T
    LOGICAL :: IS_LOWRANK = .FALSE.
    INTEGER :: M = 0, N = 0, RANK = 0
    COMPLEX(WP), ALLOCATABLE :: DENSE(:,:)         ! si IS_LOWRANK = .FALSE.
    COMPLEX(WP), ALLOCATABLE :: UFAC(:,:), VFAC(:,:)  ! si IS_LOWRANK = .TRUE.
  END TYPE BLR_BLOCK_T

  ! Variante ARENA de BLR_BLOCK_T (composantes POINTEUR, vues carvees dans
  ! un BLR_ARENA_T partage au lieu d'ALLOCATABLE individuel) -- voir la
  ! section "VARIANTES ARENA" en fin de fichier pour le detail et le
  ! compromis assume (blr_arena_mod.f90).
  TYPE :: BLR_BLOCK_ARENA_T
    LOGICAL :: IS_LOWRANK = .FALSE.
    INTEGER :: M = 0, N = 0, RANK = 0
    COMPLEX(WP), POINTER :: DENSE(:,:) => NULL()
    COMPLEX(WP), POINTER :: UFAC(:,:) => NULL(), VFAC(:,:) => NULL()
    ! Handles PERSISTANTS (1D) vers le buffer arena REELLEMENT carve pour
    ! DENSE/UFAC/VFAC -- SIZE(xxx_FLAT) est la CAPACITE actuelle (peut etre
    ! PLUS GRANDE que la forme logique ci-dessus, ex. UFAC(1:M,1:R) est une
    ! vue sur un PREFIXE de UFAC_FLAT si R a deja ete plus grand par le
    ! passe). Permet a ARENA_CARVE_REUSE (blr_arena_mod.f90) de reutiliser
    ! le meme buffer d'une recompression a l'autre au lieu d'en carver un
    ! nouveau a chaque fois -- reduit fortement le gaspillage par rapport a
    ! un carve systematique (voir SCHUR_UPDATE_BLR_BLOCK_ARENA plus bas).
    COMPLEX(WP), POINTER :: DENSE_FLAT(:) => NULL()
    COMPLEX(WP), POINTER :: UFAC_FLAT(:) => NULL(), VFAC_FLAT(:) => NULL()
  END TYPE BLR_BLOCK_ARENA_T

  !-----------------------------------------------------------------------
  ! Variante "CAPACITE FIXE" : les blocs bas-rang stockent en permanence
  ! le rang maximal que COMPRESS_WORTHWHILE_RATIO autorise ; les echanges
  ! MPI et les operations associees ne portent, eux, que sur le rang reel.
  ! Approche differente de l'arena pour le meme probleme (fragmentation de
  ! tas due aux ALLOCATE/DEALLOCATE de taille variable a chaque
  ! recompression).
  !
  ! PRINCIPE : UFAC/VFAC sont alloues UNE SEULE FOIS par bloc, a la
  ! capacite MAXIMALE RMAX = MAX(1, INT(COMPRESS_WORTHWHILE_RATIO*
  ! MIN(M,N))-1) plutot qu'a la taille EXACTE du rang courant -- puisque
  ! FAST_RECOMPRESS garantit RANK<=RMAX en permanence (sinon bascule en
  ! dense, voir SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP), cette capacite reste
  ! TOUJOURS suffisante pour TOUTE recompression future du meme bloc :
  ! JAMAIS de reallocation apres la premiere, contrairement a
  ! BLR_BLOCK_T (taille exacte, realloue a CHAQUE etape K) ou meme a
  ! l'arena (carve-ou-reutilise, mais peut encore grandir par doublement).
  ! Row-major/column-major : seul le SECOND indice de UFAC (le rang) est
  ! "gaspille" au-dela du rang reel -- le PREMIER indice (M, jamais
  ! paddes) rend tous les calculs sur UFAC directement reutilisables tels
  ! quels (le stride BLAS reste M). VFAC, a l'inverse, a son PREMIER
  ! indice (le rang) paddé a RMAX -- toute operation BLAS qui le
  ! reference doit passer la VUE tranchee VFAC(1:RANK,:) explicitement
  ! (Fortran gere alors automatiquement la copie temporaire correcte,
  ! aucune gestion manuelle de leading dimension necessaire).
  !
  ! MPI/calcul restent bornes au rang REEL (jamais RMAX) : la diffusion
  ! (COPY_FIXEDCAP_TO_REGULAR + BCAST_BLR_BLOCK, INCHANGEE, meme pattern
  ! que l'arena) copie SEULEMENT les RANK premieres colonnes/lignes vers
  ! un BLR_BLOCK_T ordinaire EXACTEMENT dimensionne avant l'envoi -- le
  ! volume MPI et le cout des TRSM/GEMM ne dependent QUE du rang reel,
  ! jamais de la capacite allouee.
  !-----------------------------------------------------------------------
  TYPE :: BLR_BLOCK_FIXEDCAP_T
    LOGICAL :: IS_LOWRANK = .FALSE.
    INTEGER :: M = 0, N = 0, RANK = 0, RMAX = 0
    COMPLEX(WP), ALLOCATABLE :: DENSE(:,:)
    COMPLEX(WP), ALLOCATABLE :: UFAC(:,:)   ! (M,RMAX) une fois alloue, jamais realloue ensuite
    COMPLEX(WP), ALLOCATABLE :: VFAC(:,:)   ! (RMAX,N) idem
  END TYPE BLR_BLOCK_FIXEDCAP_T

  !-----------------------------------------------------------------------
  ! REUTILISATION DES TABLEAUX DE TRAVAIL TRANSITOIRES : reduit
  ! l'allocation/desallocation repetee des scratch buffers de
  ! COMPRESS_BLOCK/FAST_RECOMPRESS (SVD, QR) au lieu d'un ALLOCATE/
  ! DEALLOCATE a chaque appel. Toutes les matrices 2D sont dimensionnees
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
  ! Surcharge COMPRESS_WORTHWHILE_RATIO (voir sa declaration ci-dessus) --
  ! a appeler AVANT assemblage/factorisation pour un run donne (ex. un
  ! benchmark qui balaie ce ratio, voir bench/sweep_compress_ratio.sh). Ne
  ! touche a AUCUN autre etat -- un run qui n'appelle jamais cette
  ! routine garde la valeur par defaut, comportement IDENTIQUE a avant
  ! que COMPRESS_WORTHWHILE_RATIO devienne modifiable.
  !-----------------------------------------------------------------------
  SUBROUTINE SET_COMPRESS_WORTHWHILE_RATIO(R)
    REAL(WP), INTENT(IN) :: R
    COMPRESS_WORTHWHILE_RATIO = R
  END SUBROUTINE SET_COMPRESS_WORTHWHILE_RATIO

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
  ! EXPERIMENTAL -- orthonormalise en place les COLONNES de A (ROWS x COLS,
  ! ROWS>=COLS) via CGEQRF+CUNGQR. Appelee a CHAQUE demi-etape de RSVD_
  ! RANGE_SVD (pas seulement a la fin) -- voir son commentaire pour
  ! pourquoi c'est necessaire en simple precision des que NPOWER>1.
  !-----------------------------------------------------------------------
  SUBROUTINE ORTHONORMALIZE_COLS(A, ROWS, COLS)
    INTEGER, INTENT(IN)          :: ROWS, COLS
    COMPLEX(WP), INTENT(INOUT)   :: A(ROWS,COLS)
    COMPLEX(WP), ALLOCATABLE :: TAU(:), CWORK(:)
    INTEGER :: INFO, LWORK

    ALLOCATE(TAU(COLS))
    ALLOCATE(CWORK(1)); LWORK = -1
    CALL CGEQRF(ROWS, COLS, A, ROWS, TAU, CWORK, LWORK, INFO)
    LWORK = INT(REAL(CWORK(1))); DEALLOCATE(CWORK); ALLOCATE(CWORK(MAX(LWORK,1)))
    CALL CGEQRF(ROWS, COLS, A, ROWS, TAU, CWORK, LWORK, INFO)
    DEALLOCATE(CWORK); ALLOCATE(CWORK(1)); LWORK = -1
    CALL CUNGQR(ROWS, COLS, COLS, A, ROWS, TAU, CWORK, LWORK, INFO)
    LWORK = INT(REAL(CWORK(1))); DEALLOCATE(CWORK); ALLOCATE(CWORK(MAX(LWORK,1)))
    CALL CUNGQR(ROWS, COLS, COLS, A, ROWS, TAU, CWORK, LWORK, INFO)
    DEALLOCATE(CWORK, TAU)
  END SUBROUTINE ORTHONORMALIZE_COLS

  !-----------------------------------------------------------------------
  ! EXPERIMENTAL -- coeur PARTAGE de la SVD randomisee (schema Halko-
  ! Martinsson-Tropp : projection aleatoire + iterations de puissance +
  ! QR + petite SVD), factorise ICI pour etre appele soit UNE FOIS a taille
  ! de sous-espace L FIXE (COMPRESS_BLOCK_RSVD, rang cible connu), soit EN
  ! BOUCLE a L CROISSANT (COMPRESS_BLOCK_RSVD_ADAPTIVE, rang inconnu) --
  ! une seule implementation des GEMM/QR/SVD, aucune duplication entre les
  ! deux routines publiques.
  !
  ! ALGORITHME (L colonnes de sous-espace) :
  !   1) Omega (N x L) matrice ALEATOIRE (parties reelle/imaginaire
  !      uniformes independantes, PAS besoin de gaussien strict pour un
  !      bon sous-espace d'echantillonnage en pratique).
  !   2) Y = B*Omega (M x L), puis NPOWER iterations de puissance
  !      Y <- B*(B^H*Y) (ameliore la qualite quand les valeurs
  !      singulieres decroissent lentement -- coute 2 GEMM(M,N,L)
  !      supplementaires par iteration).
  !   3) Q = orthonormalisation de Y (CGEQRF+CUNGQR) -- base du sous-
  !      espace image approxime de B.
  !   4) B_proj = Q^H*B (L x N, PETIT) -- SVD COMPLETE de B_proj (CGESVD)
  !      coute O(L*N*L) au lieu de O(M*N*min(M,N)) pour B directement :
  !      c'est le remplacement du cout domine, PROPORTIONNEL a L/min(M,N).
  !   5) UFULL = Q*Uhat (M x L) -- S(1:L) et VTHAT(1:L,:) approximent les L
  !      premieres valeurs/vecteurs singuliers de B (de plus en plus mal
  !      approximes en approchant l'indice L, d'ou la marge d'echantillon-
  !      nage OVERSAMPLE consommee par les appelants).
  !-----------------------------------------------------------------------
  SUBROUTINE RSVD_RANGE_SVD(M, N, B, LDB, L, NPOWER, S, UFULL, VTHAT)
    INTEGER, INTENT(IN)         :: M, N, LDB, L, NPOWER
    COMPLEX(WP), INTENT(IN)     :: B(LDB,N)
    REAL(WP), INTENT(OUT)       :: S(L)
    COMPLEX(WP), INTENT(OUT)    :: UFULL(M,L), VTHAT(L,N)

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: OMEGA(:,:), Y(:,:), Z(:,:), BPROJ(:,:)
    COMPLEX(WP), ALLOCATABLE :: UHAT(:,:), CWORK(:)
    REAL(WP), ALLOCATABLE    :: RWORK(:), RE_PART(:,:), IM_PART(:,:)
    INTEGER :: IP, INFO, LWORK

    ! 1) Omega (N x L) aleatoire, parties reelle/imaginaire uniformes sur
    ! [-1,1) independantes.
    ALLOCATE(OMEGA(N,L), RE_PART(N,L), IM_PART(N,L))
    CALL RANDOM_NUMBER(RE_PART)
    CALL RANDOM_NUMBER(IM_PART)
    OMEGA = CMPLX(2.0_WP*RE_PART-1.0_WP, 2.0_WP*IM_PART-1.0_WP, KIND=WP)
    DEALLOCATE(RE_PART, IM_PART)

    ! 2) Y = B*Omega, orthonormalise immediatement (voir ORTHONORMALIZE_
    ! COLS), puis NPOWER iterations de puissance -- chaque demi-iteration
    ! (Z=B^H*Y puis Y=B*Z) est re-orthonormalisee aussitot, jamais
    ! seulement a la toute fin. Cette re-orthonormalisation intermediaire
    ! est necessaire : sans elle, la plage dynamique de Y/Z explose au fil
    ! des iterations (le rapport des valeurs singulieres extremes est
    ! porte a la puissance 2*IP+1), ce qui degrade catastrophiquement la
    ! precision en simple precision (WP) des qu'on depasse NPOWER=1
    ! (erreur relative x30 pire entre NPOWER=1 et NPOWER=3 sur un cas reel,
    ! voir test/bench_rsvd_compress.f90) plutot que de l'ameliorer comme
    ! prevu par la theorie (qui suppose une arithmetique a precision
    ! infinie ou re-orthonormalisee a chaque etape).
    ALLOCATE(Y(M,L))
    CALL CGEMM('No transpose', 'No transpose', M, L, N, CONE, B, LDB, OMEGA, N, CZERO, Y, M)
    DEALLOCATE(OMEGA)
    CALL ORTHONORMALIZE_COLS(Y, M, L)
    DO IP = 1, NPOWER
      ALLOCATE(Z(N,L))
      CALL CGEMM('Conjugate transpose', 'No transpose', N, L, M, CONE, B, LDB, Y, M, CZERO, Z, N)
      CALL ORTHONORMALIZE_COLS(Z, N, L)
      CALL CGEMM('No transpose', 'No transpose', M, L, N, CONE, B, LDB, Z, N, CZERO, Y, M)
      CALL ORTHONORMALIZE_COLS(Y, M, L)
      DEALLOCATE(Z)
    END DO
    ! Y contient maintenant Q (M x L), colonnes orthonormees (par
    ! construction : derniere operation de la boucle, ou de l'etape
    ! initiale si NPOWER=0).

    ! 4) B_proj = Q^H*B (L x N, PETIT) -- SVD COMPLETE de ce petit systeme.
    ALLOCATE(BPROJ(L,N))
    CALL CGEMM('Conjugate transpose', 'No transpose', L, N, M, CONE, Y, M, B, LDB, CZERO, BPROJ, L)

    ALLOCATE(UHAT(L,L), RWORK(5*L))
    ALLOCATE(CWORK(1)); LWORK = -1
    CALL CGESVD('S', 'S', L, N, BPROJ, L, S, UHAT, L, VTHAT, L, CWORK, LWORK, RWORK, INFO)
    LWORK = INT(REAL(CWORK(1))); DEALLOCATE(CWORK); ALLOCATE(CWORK(MAX(LWORK,1)))
    CALL CGESVD('S', 'S', L, N, BPROJ, L, S, UHAT, L, VTHAT, L, CWORK, LWORK, RWORK, INFO)
    DEALLOCATE(CWORK, BPROJ, RWORK)

    ! 5) UFULL = Q*Uhat (M x L).
    CALL CGEMM('No transpose', 'No transpose', M, L, L, CONE, Y, M, UHAT, L, CZERO, UFULL, M)
    DEALLOCATE(Y, UHAT)
  END SUBROUTINE RSVD_RANGE_SVD

  !-----------------------------------------------------------------------
  ! EXPERIMENTAL -- SVD randomisee a RANG CIBLE FIXE, PREMIERE PASSE
  ! demandee explicitement pour mesurer l'ordre de grandeur du gain face a
  ! COMPRESS_BLOCK (CGESVD complet, mesure : ~90-95% du temps d'assemblage
  ! a N=20000, voir bench/bench_blr_mpi_pivot_ooc_compare.f90) -- RANK_
  ! TARGET choisi par l'APPELANT (typiquement en copiant le rang deja
  ! trouve par COMPRESS_BLOCK pour une comparaison a rang EGAL). Voir
  ! COMPRESS_BLOCK_RSVD_ADAPTIVE ci-dessous pour la version SANS rang
  ! connu a l'avance (le cas d'usage reel). Routine INDEPENDANTE de
  ! COMPRESS_BLOCK, n'affecte AUCUN code existant.
  !
  ! GARANTIE : PROBABILISTE, pas deterministe (contrairement a CGESVD) --
  ! l'erreur d'approximation depend du tirage aleatoire d'Omega ET de la
  ! decroissance reelle des valeurs singulieres de B. Ne PAS utiliser tel
  ! quel en production sans validation de la marge d'erreur (voir
  ! test/bench_rsvd_compress.f90 pour la mesure temps ET precision).
  !-----------------------------------------------------------------------
  SUBROUTINE COMPRESS_BLOCK_RSVD(M, N, B, LDB, RANK_TARGET, OVERSAMPLE, NPOWER, BLK)
    INTEGER, INTENT(IN)            :: M, N, LDB, RANK_TARGET, OVERSAMPLE, NPOWER
    COMPLEX(WP), INTENT(IN)        :: B(LDB,N)
    TYPE(BLR_BLOCK_T), INTENT(OUT) :: BLK

    COMPLEX(WP), ALLOCATABLE :: UFULL(:,:), VTHAT(:,:)
    REAL(WP), ALLOCATABLE    :: S(:)
    INTEGER :: L, R, I

    L = MIN(RANK_TARGET + OVERSAMPLE, MIN(M,N))
    R = MIN(RANK_TARGET, L)
    BLK%M = M; BLK%N = N

    ALLOCATE(UFULL(M,L), VTHAT(L,N), S(L))
    CALL RSVD_RANGE_SVD(M, N, B, LDB, L, NPOWER, S, UFULL, VTHAT)

    BLK%IS_LOWRANK = .TRUE.
    BLK%RANK = R
    ALLOCATE(BLK%UFAC(M,R), BLK%VFAC(R,N))
    BLK%UFAC = UFULL(1:M,1:R)
    DO I = 1, R
      BLK%VFAC(I,1:N) = S(I) * VTHAT(I,1:N)
    END DO
    DEALLOCATE(UFULL, VTHAT, S)
  END SUBROUTINE COMPRESS_BLOCK_RSVD

  !-----------------------------------------------------------------------
  ! EXPERIMENTAL -- SVD randomisee ADAPTATIVE : rang PAS connu a l'avance,
  ! decouvert en faisant CROITRE la taille du sous-espace d'echantillon
  ! (doublement) jusqu'a ce qu'un CRITERE D'ARRET base sur la marge de
  ! sur-echantillonnage restante indique une confiance suffisante -- voir
  ! COMPRESS_BLOCK_RSVD ci-dessus pour la version a rang deja connu.
  !
  ! CRITERE D'ARRET (heuristique PRATIQUE, PAS le critere probabiliste
  ! rigoureux de Halko-Martinsson-Tropp Algorithme 4.2, qui necessiterait
  ! un Gram-Schmidt incremental colonne par colonne avec estimateur
  ! d'erreur par vecteurs de sondage supplementaires -- plus complexe a
  ! implementer correctement, laisse comme etape ulterieure si celui-ci
  ! s'avere insuffisant en pratique) : a sous-espace de taille L, on compte
  ! R = nombre de valeurs singulieres approximees S(1:L) au-dessus du
  ! seuil TOL*S(1) (MEME regle que COMPRESS_BLOCK). Si R <= L-OVERSAMPLE,
  ! il reste une marge de sur-echantillonnage NON consommee -- signe que
  ! le sous-espace est assez grand pour avoir capture toutes les
  ! directions significatives (les L-R dernieres colonnes, dans la zone de
  ! marge, ont ete jugees NON significatives). Si R depasse ce seuil, la
  ! marge est "epuisee" -- signe qu'il pourrait exister des directions
  ! significatives AU-DELA de L qu'on n'a pas encore vues : on DOUBLE L et
  ! on recommence. Arret aussi si L atteint LCAP (= seuil COMPRESS_
  ! WORTHWHILE_RATIO*min(M,N), le MEME plafond que COMPRESS_BLOCK utilise
  ! pour son repli dense -- au-dela, la compression ne vaudrait de toute
  ! facon pas le coup) : reste alors DENSE, comme COMPRESS_BLOCK.
  !
  ! GARANTIE : PROBABILISTE ET HEURISTIQUE (le critere d'arret lui-meme
  ! n'a pas de borne de confiance formelle, contrairement au vrai
  ! algorithme adaptatif HMT). MESURE (test/bench_rsvd_compress.f90, cas
  ! limite volontairement difficile -- rang juste en-dessous du seuil de
  ! COMPRESS_WORTHWHILE_RATIO) : le rang trouve correspond EXACTEMENT au
  ! rang de reference (COMPRESS_BLOCK) et l'erreur de reconstruction est
  ! QUASI IDENTIQUE (ratio 1.00), y compris avec un OVERSAMPLE reduit a 2
  ! (2 doublements de sous-espace necessaires dans ce cas). Ce resultat
  ! depend d'un point CRITIQUE (voir RSVD_RANGE_SVD/ORTHONORMALIZE_COLS) :
  ! chaque demi-etape d'iteration de puissance doit etre RE-orthonormalisee
  ! immediatement -- un bogue initial (orthonormalisation UNIQUEMENT en
  ! fin de calcul) donnait une erreur x30 a x2000 PIRE des que NPOWER>1 en
  ! simple precision, corrige depuis. Reste malgre tout une garantie
  ! HEURISTIQUE, pas une preuve -- valider par mesure sur le cas d'usage
  ! reel avant tout usage en production.
  !-----------------------------------------------------------------------
  SUBROUTINE COMPRESS_BLOCK_RSVD_ADAPTIVE(M, N, B, LDB, TOL, BLK, OVERSAMPLE, NPOWER, N_GROWTHS)
    INTEGER, INTENT(IN)            :: M, N, LDB
    COMPLEX(WP), INTENT(IN)        :: B(LDB,N)
    REAL(WP), INTENT(IN)           :: TOL
    TYPE(BLR_BLOCK_T), INTENT(OUT) :: BLK
    INTEGER, INTENT(IN), OPTIONAL  :: OVERSAMPLE, NPOWER
    ! N_GROWTHS (optionnel, INTENT(OUT)) : nombre de doublements de L avant
    ! arret -- diagnostic pratique pour juger si le critere d'arret
    ! converge vite (0-1 doublement, cas courant) ou tatonne (plusieurs
    ! doublements, signe que le sous-espace initial etait mal dimensionne).
    INTEGER, INTENT(OUT), OPTIONAL :: N_GROWTHS

    COMPLEX(WP), ALLOCATABLE :: UFULL(:,:), VTHAT(:,:)
    REAL(WP), ALLOCATABLE    :: S(:)
    INTEGER :: OVS, NPOW, MN, L, LCAP, LNEW, R, I, IGROW
    REAL(WP) :: THRESH
    LOGICAL :: MARGIN_OK

    OVS = 10; IF (PRESENT(OVERSAMPLE)) OVS = OVERSAMPLE
    NPOW = 1; IF (PRESENT(NPOWER)) NPOW = NPOWER
    MN = MIN(M,N)
    BLK%M = M; BLK%N = N

    ! Meme plafond que COMPRESS_BLOCK (blr_arena_mod-like RMAX) : au-dela,
    ! la compression n'est de toute facon pas rentable, inutile de faire
    ! grossir le sous-espace plus loin.
    LCAP = MAX(1, INT(COMPRESS_WORTHWHILE_RATIO*REAL(MN,WP)) - 1)
    LCAP = MAX(LCAP, OVS+1)   ! garde-fou : LCAP doit rester > OVS pour que le critere d'arret ait un sens

    ! ZONE MORTE CARACTERISEE (mesure sur bloc a spectre PRESCRIT, rang
    ! exact connu par construction, voir test/bench_rsvd_compress.f90
    ! RANK_WANT) : pour tout rang REEL dans (LCAP-OVERSAMPLE, seuil reel de
    ! COMPRESS_BLOCK] (ex. M=N=400, OVERSAMPLE=10 -> rangs 123 a 133), la
    ! marge R<=L-OVERSAMPLE ne peut JAMAIS etre satisfaite avant que L
    ! n'atteigne LCAP (L ne grossit jamais au-dela) -- l'adaptatif abandonne
    ! et replie SYSTEMATIQUEMENT en dense (0/20 tirages, PAS une question de
    ! chance), alors que COMPRESS_BLOCK compresserait correctement ces
    ! blocs. Repli dense = reconstruction EXACTE, donc pas de corruption,
    ! seulement une occasion de gain RSVD manquee sur cette bande etroite
    ! (~OVERSAMPLE valeurs de rang, pres du seuil COMPRESS_WORTHWHILE_
    ! RATIO). Comportement accepte tel quel (pas de correctif demande) --
    ! reduire OVERSAMPLE reduirait cette bande mais aussi la marge de
    ! securite statistique du critere d'arret ailleurs.

    ! MARGIN_OK distingue explicitement "la marge R<=L-OVERSAMPLE a ete
    ! satisfaite" (rang digne de confiance) de "le plafond LCAP/MN a ete
    ! atteint sans jamais satisfaire la marge" (la boucle a simplement
    ! abandonne, R n'a alors aucune valeur probante) : comme R est de
    ! toute facon toujours <= LCAP < RATIO*MN par construction, un test
    ! final qui se contenterait de R < RATIO*MN serait toujours vrai des
    ! qu'on atteint le plafond, meme pour un bloc franchement plein rang,
    ! et inventerait un rang = LCAP au lieu du repli dense correct.
    ! MARGIN_OK trace donc explicitement laquelle des deux conditions
    ! d'arret a declenche la sortie de boucle, pour que la decision finale
    ! s'appuie dessus plutot que sur une relecture de R seul.
    L = MIN(2*OVS, MN, LCAP)
    IGROW = 0
    MARGIN_OK = .FALSE.
    DO
      IF (ALLOCATED(UFULL)) DEALLOCATE(UFULL, VTHAT, S)
      ALLOCATE(UFULL(M,L), VTHAT(L,N), S(L))
      CALL RSVD_RANGE_SVD(M, N, B, LDB, L, NPOW, S, UFULL, VTHAT)

      THRESH = TOL * MAX(S(1), TINY(1.0_WP))
      R = 0
      DO I = 1, L
        IF (S(I) .GT. THRESH) R = R + 1
      END DO

      IF (R .LE. L-OVS) THEN
        MARGIN_OK = .TRUE.
        EXIT
      END IF
      IF (L .GE. LCAP .OR. L .GE. MN) EXIT   ! plafond atteint, marge JAMAIS satisfaite -- abandon
      LNEW = MIN(2*L, MN, LCAP)
      IF (LNEW .LE. L) EXIT   ! plus de marge de croissance possible
      L = LNEW
      IGROW = IGROW + 1
    END DO
    IF (PRESENT(N_GROWTHS)) N_GROWTHS = IGROW

    IF (MARGIN_OK) THEN
      ! Rang DIGNE DE CONFIANCE (marge de sur-echantillonnage non
      ! consommee) -- R est de toute facon <= L <= LCAP < RATIO*MN a ce
      ! stade, la compression est automatiquement rentable.
      BLK%IS_LOWRANK = .TRUE.
      BLK%RANK = R
      ALLOCATE(BLK%UFAC(M,R), BLK%VFAC(R,N))
      BLK%UFAC = UFULL(1:M,1:R)
      DO I = 1, R
        BLK%VFAC(I,1:N) = S(I) * VTHAT(I,1:N)
      END DO
    ELSE
      ! Repli dense, EXACTEMENT comme COMPRESS_BLOCK -- pas de SVD
      ! supplementaire : soit la marge n'a jamais ete satisfaite avant le
      ! plafond (bloc pas assez compressible, R n'est PAS fiable), soit le
      ! bloc est trop petit pour meme atteindre LCAP.
      BLK%IS_LOWRANK = .FALSE.
      ALLOCATE(BLK%DENSE(M,N))
      BLK%DENSE = B(1:M,1:N)
    END IF
    DEALLOCATE(UFULL, VTHAT, S)
  END SUBROUTINE COMPRESS_BLOCK_RSVD_ADAPTIVE

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
  ! POURQUOI CA RESTE BON MARCHE EN BLR (le point cle qui rend le pivot
  ! confine au bloc diagonal compatible avec la compression) : permuter
  ! des lignes d'un bloc bas-rang A~=UFAC*VFAC ne touche que UFAC (VFAC
  ! n'a aucune notion de ligne) -- cout O(rang) par swap, jamais de
  ! decompression. Cas dense : swap direct sur DENSE, cout O(N) par swap
  ! (standard).
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
  ! LU sans pivot, blockee (right-looking, meme structure que BLOCK_LU_
  ! FACTOR/block_lu_mod.f90, mais sans CGETRF/pivot) -- remplace
  ! BLR_LU_NOPIV_UNBLOCKED sur le bloc diagonal de mpi_blr_matrix_mod.f90,
  ! dont le cout domine 40-50% du temps de factorisation a NB>=750 : la
  ! version non blockee ne fait que du BLAS niveau 1/2 (CSCAL/CGERU),
  ! intrinsequement sequentiel et sans reutilisation de cache -- ici, seul
  ! le panneau (INNER_NB colonnes, petit) reste non blocke (CALL a
  ! BLR_LU_NOPIV_UNBLOCKED ci-dessus), le gros du travail (mise a jour de
  ! la sous-matrice restante) passe par CTRSM/CGEMM (BLAS niveau 3 --
  ! threadable par MKL, bien meilleure intensite arithmetique).
  !
  ! DELIBEREMENT PAS CGETRF (LAPACK) : CGETRF fait du pivot partiel par
  ! construction, ce qui exigerait de propager les permutations de lignes
  ! aux blocs L21/U12 de la meme ligne/colonne pour rester correct --
  ! un pivot qui ne decompresse jamais un bloc bas-rang, hors de la portee
  ! de cette routine. Cette version reste sans pivot de bout en bout,
  ! seulement blockee -- meme resultat numerique que BLR_LU_NOPIV_
  ! UNBLOCKED (memes operations, juste reordonnees/groupees), pas une
  ! nouvelle methode.
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

  !===========================================================================
  ! VARIANTES "ARENA" -- voir blr_arena_mod.f90 pour le principe general et
  ! le compromis assume (pas de compaction, memoire des anciennes versions
  ! de blocs jamais recuperee). Type PARALLELE a BLR_BLOCK_T (composantes
  ! POINTEUR, vues carvees dans un BLR_ARENA_T partage, au lieu de
  ! ALLOCATABLE) + copies des fonctions qui touchent le stockage d'un bloc
  ! -- Fortran exige une correspondance de type EXACTE entre argument reel
  ! et factice, meme quand le corps est identique a la version non-arena.
  !
  ! Utilisees UNIQUEMENT par MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA
  ! (mpi_blr_matrix_mod.f90, nouveau) -- AUCUNE fonction/type existant(e)
  ! ci-dessus n'est modifie(e), les solveurs actuels (A/B/MPI/C/BLR+MPI,
  ! avec et sans pivot) restent inchanges et continuent d'utiliser
  ! BLR_BLOCK_T tel quel.
  !
  ! Les panneaux TRANSITOIRES (diffuses via BCAST_BLR_BLOCK, duree de vie
  ! d'UNE etape K -- LRECV/URECV cote appelant, PBLK ici) restent des
  ! BLR_BLOCK_T ORDINAIRES (ALLOCATABLE) : seul le stockage PERSISTANT
  ! (A_BLR%OFFBLK, qui doit survivre a toutes les etapes futures) a besoin
  ! d'etre arena-compatible -- voir COPY_ARENA_TO_REGULAR pour le point de
  ! passage entre les deux mondes.
  !===========================================================================

  !-----------------------------------------------------------------------
  ! Copie une vue arena (POINTEUR) vers un BLR_BLOCK_T ALLOCATABLE
  ! independant -- point de passage vers le monde "transitoire" ordinaire
  ! (ex. avant de diffuser un panneau via l'existant BCAST_BLR_BLOCK, qui
  ! attend un BLR_BLOCK_T). Vraie copie de donnees (pas un partage de
  ! pointeur) -- MEME COUT que ce que le code non-arena fait deja
  ! aujourd'hui ("LRECV(li) = A_BLR%OFFBLK(...)", une copie identique),
  ! aucun surcout supplementaire introduit.
  !-----------------------------------------------------------------------
  SUBROUTINE COPY_ARENA_TO_REGULAR(SRC, DST)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(IN) :: SRC
    TYPE(BLR_BLOCK_T), INTENT(OUT)      :: DST
    DST%IS_LOWRANK = SRC%IS_LOWRANK
    DST%M = SRC%M; DST%N = SRC%N; DST%RANK = SRC%RANK
    IF (SRC%IS_LOWRANK) THEN
      ALLOCATE(DST%UFAC(SRC%M,SRC%RANK), DST%VFAC(SRC%RANK,SRC%N))
      DST%UFAC = SRC%UFAC
      DST%VFAC = SRC%VFAC
    ELSE
      ALLOCATE(DST%DENSE(SRC%M,SRC%N))
      DST%DENSE = SRC%DENSE
    END IF
  END SUBROUTINE COPY_ARENA_TO_REGULAR

  !-----------------------------------------------------------------------
  ! Equivalent COMPRESS_BLOCK -- ecrit dans une vue carvee de ARENA au
  ! lieu d'ALLOUER BLK%DENSE/UFAC/VFAC directement. Utilise a l'assemblage
  ! initial (MPI_ASSEMBLE_BLR_MATRIX_ARENA). Corps IDENTIQUE a
  ! COMPRESS_BLOCK a l'exception des 3 endroits qui ecrivent le resultat
  ! final (marques ci-dessous).
  !-----------------------------------------------------------------------
  SUBROUTINE COMPRESS_BLOCK_ARENA(M, N, B, LDB, TOL, BLK, ARENA, WS)
    INTEGER, INTENT(IN)          :: M, N, LDB
    COMPLEX(WP), INTENT(IN)      :: B(LDB,N)
    REAL(WP), INTENT(IN)         :: TOL
    TYPE(BLR_BLOCK_ARENA_T), INTENT(OUT) :: BLK
    TYPE(BLR_ARENA_T), INTENT(INOUT), TARGET :: ARENA
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS

    COMPLEX(WP), ALLOCATABLE, TARGET :: BCOPY_LOC(:,:), UFULL_LOC(:,:), VTFULL_LOC(:,:), CWORK_LOC(:)
    REAL(WP), ALLOCATABLE, TARGET :: S_LOC(:), RWORK_LOC(:)
    COMPLEX(WP), POINTER :: BCOPY(:,:), UFULL(:,:), VTFULL(:,:), CWORK(:)
    REAL(WP), POINTER :: S(:), RWORK(:)
    COMPLEX(WP), POINTER :: FLAT(:)
    INTEGER :: MN, LDW, LDV, LWORK, INFO, R, I
    REAL(WP) :: THRESH

    BLK%M = M
    BLK%N = N
    MN = MIN(M,N)

    IF (TOL .LE. 0.0_WP .OR. MN .EQ. 0) THEN
      BLK%IS_LOWRANK = .FALSE.
      FLAT => ARENA_CARVE_REUSE(ARENA, BLK%DENSE_FLAT, INT(M,8)*INT(N,8))   ! <-- carve/reutilise
      BLK%DENSE(1:M,1:N) => FLAT
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
      FLAT => ARENA_CARVE_REUSE(ARENA, BLK%DENSE_FLAT, INT(M,8)*INT(N,8))   ! <-- carve/reutilise
      BLK%DENSE(1:M,1:N) => FLAT
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
      FLAT => ARENA_CARVE_REUSE(ARENA, BLK%UFAC_FLAT, INT(M,8)*INT(R,8))   ! <-- carve/reutilise UFAC
      BLK%UFAC(1:M,1:R) => FLAT
      FLAT => ARENA_CARVE_REUSE(ARENA, BLK%VFAC_FLAT, INT(R,8)*INT(N,8))   ! <-- carve/reutilise VFAC
      BLK%VFAC(1:R,1:N) => FLAT
      BLK%UFAC = UFULL(1:M,1:R)
      DO I = 1, R
        BLK%VFAC(I,1:N) = S(I) * VTFULL(I,1:N)
      END DO
    ELSE
      BLK%IS_LOWRANK = .FALSE.
      FLAT => ARENA_CARVE_REUSE(ARENA, BLK%DENSE_FLAT, INT(M,8)*INT(N,8))   ! <-- carve/reutilise
      BLK%DENSE(1:M,1:N) => FLAT
      BLK%DENSE = B(1:M,1:N)
    END IF

    IF (.NOT. PRESENT(WS)) DEALLOCATE(BCOPY_LOC, UFULL_LOC, VTFULL_LOC, S_LOC, RWORK_LOC, CWORK_LOC)
  END SUBROUTINE COMPRESS_BLOCK_ARENA

  !-----------------------------------------------------------------------
  ! Copie un BLR_BLOCK_T ordinaire (ex. issu de COMPRESS_BLOCK_RSVD_
  ! ADAPTIVE) dans un bloc ARENA, en carvant DENSE ou UFAC/VFAC dans ARENA
  ! (memes appels ARENA_CARVE_REUSE que COMPRESS_BLOCK_ARENA).
  !-----------------------------------------------------------------------
  SUBROUTINE COPY_REGULAR_TO_ARENA(SRC, DST, ARENA)
    TYPE(BLR_BLOCK_T), INTENT(IN)            :: SRC
    TYPE(BLR_BLOCK_ARENA_T), INTENT(OUT)     :: DST
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
    ELSE
      FLAT => ARENA_CARVE_REUSE(ARENA, DST%DENSE_FLAT, INT(M,8)*INT(N,8))
      DST%DENSE(1:M,1:N) => FLAT
      IF (M.GT.0 .AND. N.GT.0) DST%DENSE = SRC%DENSE
    END IF
  END SUBROUTINE COPY_REGULAR_TO_ARENA

  !-----------------------------------------------------------------------
  ! Equivalents EXPAND_BLOCK / BLR_BLOCK_ROW_SWAP / BLR_TRSM_*/
  ! BLR_BLOCK_APPLY, corps STRICTEMENT IDENTIQUE aux versions non-arena --
  ! seule la declaration de type de BLK differe (Fortran l'exige).
  !-----------------------------------------------------------------------
  SUBROUTINE EXPAND_BLOCK_ARENA(BLK, OUT, LDOUT)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(IN) :: BLK
    INTEGER, INTENT(IN)           :: LDOUT
    COMPLEX(WP), INTENT(OUT)      :: OUT(LDOUT,BLK%N)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)

    IF (BLK%IS_LOWRANK) THEN
      CALL CGEMM('No transpose', 'No transpose', BLK%M, BLK%N, BLK%RANK, &
                  CONE, BLK%UFAC, BLK%M, BLK%VFAC, BLK%RANK, CZERO, OUT, LDOUT)
    ELSE
      OUT(1:BLK%M,1:BLK%N) = BLK%DENSE
    END IF
  END SUBROUTINE EXPAND_BLOCK_ARENA

  SUBROUTINE BLR_BLOCK_ROW_SWAP_ARENA(BLK, KW, IPIV)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(INOUT) :: BLK
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
  END SUBROUTINE BLR_BLOCK_ROW_SWAP_ARENA

  SUBROUTINE BLR_TRSM_LEFT_LOWER_UNIT_ARENA(DKK, KW, BLK)
    INTEGER, INTENT(IN)             :: KW
    COMPLEX(WP), INTENT(IN)         :: DKK(KW,KW)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(INOUT) :: BLK
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, BLK%RANK, CONE, DKK, KW, BLK%UFAC, KW)
    ELSE
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, BLK%N, CONE, DKK, KW, BLK%DENSE, KW)
    END IF
  END SUBROUTINE BLR_TRSM_LEFT_LOWER_UNIT_ARENA

  SUBROUTINE BLR_TRSM_RIGHT_UPPER_NONUNIT_ARENA(BLK, DKK, KW)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN)             :: KW
    COMPLEX(WP), INTENT(IN)         :: DKK(KW,KW)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', BLK%RANK, KW, CONE, DKK, KW, BLK%VFAC, BLK%RANK)
    ELSE
      CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', BLK%M, KW, CONE, DKK, KW, BLK%DENSE, BLK%M)
    END IF
  END SUBROUTINE BLR_TRSM_RIGHT_UPPER_NONUNIT_ARENA

  SUBROUTINE BLR_BLOCK_APPLY_ARENA(BLK, NRHS, ALPHA, X, LDX, BETA, Y, LDY)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(IN) :: BLK
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
  END SUBROUTINE BLR_BLOCK_APPLY_ARENA

  !-----------------------------------------------------------------------
  ! Equivalent FAST_RECOMPRESS -- corps IDENTIQUE (meme algorithme QR+SVD
  ! du "coeur", memes buffers de travail via WS) a l'exception du calcul
  ! du resultat final NEWBLK%UFAC/VFAC, carve dans ARENA au lieu d'etre
  ! ALLOUE directement (seul endroit modifie, marque ci-dessous). Les
  ! arguments d'entree UA/VA/PU/PV restent de simples tableaux (pas de
  ! type BLR_BLOCK_T/BLR_BLOCK_ARENA_T ici) -- deja reutilisable tel quel
  ! par l'appelant, qu'il lise depuis un bloc ordinaire ou arena.
  !-----------------------------------------------------------------------
  SUBROUTINE FAST_RECOMPRESS_ARENA(MI, NJ, RA, UA, VA, RP, PU, PV, TOL, RMAX, NEWBLK, OK, ARENA, &
                                    UFAC_FLAT, VFAC_FLAT, WS)
    INTEGER, INTENT(IN)            :: MI, NJ, RA, RP, RMAX
    COMPLEX(WP), INTENT(IN)        :: UA(MI,RA), VA(RA,NJ), PU(MI,RP), PV(RP,NJ)
    REAL(WP), INTENT(IN)           :: TOL
    TYPE(BLR_BLOCK_ARENA_T), INTENT(OUT) :: NEWBLK
    LOGICAL, INTENT(OUT)           :: OK
    TYPE(BLR_ARENA_T), INTENT(INOUT), TARGET :: ARENA
    ! Handles PERSISTANTS du bloc ABLK appelant (voir BLR_BLOCK_ARENA_T) --
    ! passes ici pour permettre a ARENA_CARVE_REUSE de reutiliser le buffer
    ! DEJA carve pour ABLK%UFAC/VFAC lors d'un appel PRECEDENT, au lieu
    ! d'en carver un nouveau a chaque recompression (voir blr_arena_mod.f90
    ! /ARENA_CARVE_REUSE). SANS RISQUE d'aliasing : UA/VA (== l'ANCIEN
    ! contenu de ABLK%UFAC/VFAC) sont entierement recopies dans UCAT/VCATT
    ! (buffers de travail SEPARES) des le debut de cette routine, AVANT
    ! que le resultat final ne soit ecrit dans ce meme buffer reutilise --
    ! l'ancien contenu n'est plus necessaire au moment de l'ecriture.
    COMPLEX(WP), POINTER, INTENT(INOUT) :: UFAC_FLAT(:), VFAC_FLAT(:)
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS

    INTEGER :: RTOT, R, I, LWORK, INFO, LDU, LDVQ, LDSQ, LDSV
    COMPLEX(WP), ALLOCATABLE, TARGET :: UCAT_LOC(:,:), VCATT_LOC(:,:), TAUU_LOC(:), TAUV_LOC(:)
    COMPLEX(WP), ALLOCATABLE, TARGET :: RU_LOC(:,:), RV_LOC(:,:), W_LOC(:,:), UW_LOC(:,:), VWH_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE, TARGET :: SVWH_LOC(:,:), CWORK_LOC(:)
    REAL(WP), ALLOCATABLE, TARGET    :: SW_LOC(:), RWORK_LOC(:)
    COMPLEX(WP), POINTER :: UCAT(:,:), VCATT(:,:), TAUU(:), TAUV(:)
    COMPLEX(WP), POINTER :: RU(:,:), RV(:,:), WMAT(:,:), UW(:,:), VWH(:,:), SVWH(:,:), CWORK(:)
    REAL(WP), POINTER :: SW(:), RWORK(:)
    COMPLEX(WP), POINTER :: FLAT(:)
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

    CALL CGEQRF(MI, RTOT, UCAT, LDU, TAUU, CWORK, LWORK, INFO)
    RU(1:RTOT,1:RTOT) = CZERO
    DO I = 1, RTOT
      RU(1:I,I) = UCAT(1:I,I)
    END DO
    CALL CUNGQR(MI, RTOT, RTOT, UCAT, LDU, TAUU, CWORK, LWORK, INFO)

    CALL CGEQRF(NJ, RTOT, VCATT, LDVQ, TAUV, CWORK, LWORK, INFO)
    RV(1:RTOT,1:RTOT) = CZERO
    DO I = 1, RTOT
      RV(1:I,I) = VCATT(1:I,I)
    END DO
    CALL CUNGQR(NJ, RTOT, RTOT, VCATT, LDVQ, TAUV, CWORK, LWORK, INFO)

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
      ! RCAP = MIN(MI,NJ)/2, PAS RMAX : plafond de croissance amortie (voir
      ! ARENA_CARVE_REUSE/blr_arena_mod.f90) choisi au point d'EQUILIBRE
      ! MEMOIRE bas-rang/dense (stockage bas-rang = R*(M+N), dense = M*N ;
      ! egaux pour R = M*N/(M+N), soit M/2 pour un bloc carre) -- PAS
      ! RMAX (~0.8*MIN(MI,NJ), plafond ALGORITHMIQUE de
      ! COMPRESS_WORTHWHILE_RATIO, tune pour le cout CALCUL, pas memoire).
      ! Un rang R.LE.RCAP est de toute facon TOUJOURS moins cher que le
      ! bloc dense final -- garantit qu'un bloc qui finit par basculer en
      ! dense (voir SCHUR_UPDATE_BLR_BLOCK_ARENA) n'aura jamais engage,
      ! via ce plafond de croissance, plus que ce que coute le bloc dense
      ! lui-meme. Au-dela de RCAP (rang REELLEMENT necessaire plus grand,
      ! toujours possible jusqu'a RMAX), ARENA_CARVE_REUSE carve alors
      ! EXACTEMENT le besoin (le garde-fou CEILING=MAX(CAP_MAX,NELEMS) de
      ! ARENA_CARVE_REUSE l'assure), sans jamais refuser la place.
      FLAT => ARENA_CARVE_REUSE(ARENA, UFAC_FLAT, INT(MI,8)*INT(R,8), INT(MI,8)*INT(MAX(1,MIN(MI,NJ)/2),8))
      NEWBLK%UFAC(1:MI,1:R) => FLAT
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
      FLAT => ARENA_CARVE_REUSE(ARENA, VFAC_FLAT, INT(R,8)*INT(NJ,8), INT(MAX(1,MIN(MI,NJ)/2),8)*INT(NJ,8))
      NEWBLK%VFAC(1:R,1:NJ) => FLAT
      CALL CGEMM('No transpose', 'Transpose', R, NJ, RTOT, CONE, &
                  SVWH, LDSV, VCATT, LDVQ, CZERO, NEWBLK%VFAC, R)
      IF (.NOT. PRESENT(WS)) DEALLOCATE(SVWH_LOC)
    END IF

    IF (.NOT. PRESENT(WS)) DEALLOCATE(UCAT_LOC, VCATT_LOC, TAUU_LOC, TAUV_LOC, RU_LOC, RV_LOC, &
                                       W_LOC, UW_LOC, VWH_LOC, SW_LOC, RWORK_LOC, CWORK_LOC)
  END SUBROUTINE FAST_RECOMPRESS_ARENA

  !-----------------------------------------------------------------------
  ! Equivalent SCHUR_UPDATE_BLR_BLOCK -- ABLK est ici BLR_BLOCK_ARENA_T
  ! (stockage PERSISTANT). LBLK/UBLK (panneaux DIFFUSES, transitoires,
  ! duree de vie d'une etape K) restent des BLR_BLOCK_T ORDINAIRES, tout
  ! comme PBLK (produit L*U, purement local a cet appel) -- voir le
  ! commentaire de tete de la section arena. Cas succes bas-rang : la
  ! "reaffectation" de ABLK est une SIMPLE COPIE DE POINTEUR (zero copie
  ! de donnees, zero (re)allocation) -- c'est le gain recherche par
  ! rapport a "ABLK = NEWBLK" (reallocation+copie implicite, voir
  ! SCHUR_UPDATE_BLR_BLOCK/USE_MOVE_ALLOC). L'ancienne vue de ABLK devient
  ! simplement orpheline dans l'arena (jamais recuperee -- compromis
  ! assume, voir blr_arena_mod.f90).
  !-----------------------------------------------------------------------
  SUBROUTINE SCHUR_UPDATE_BLR_BLOCK_ARENA(ABLK, MI, NJ, KW, LBLK, UBLK, TOL, ARENA, WS, T_PRODUCT, T_RECOMP)
    TYPE(BLR_BLOCK_ARENA_T), INTENT(INOUT) :: ABLK
    INTEGER, INTENT(IN)              :: MI, NJ, KW
    TYPE(BLR_BLOCK_T), INTENT(IN)    :: LBLK, UBLK
    REAL(WP), INTENT(IN)             :: TOL
    TYPE(BLR_ARENA_T), INTENT(INOUT), TARGET :: ARENA
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: T_PRODUCT, T_RECOMP

    TYPE(BLR_BLOCK_T) :: PBLK
    TYPE(BLR_BLOCK_ARENA_T) :: NEWBLK
    COMPLEX(WP), POINTER :: FLAT(:), DENSE_VIEW(:,:)
    INTEGER :: RMAX
    LOGICAL :: OK
    DOUBLE PRECISION :: TT0

    IF (MI.EQ.0 .OR. NJ.EQ.0 .OR. KW.EQ.0) RETURN

    IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
    CALL COMPUTE_LU_PRODUCT(MI, NJ, KW, LBLK, UBLK, PBLK)
    IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)

    IF (ABLK%IS_LOWRANK .AND. PBLK%IS_LOWRANK) THEN
      RMAX = MAX(1, INT(COMPRESS_WORTHWHILE_RATIO*REAL(MIN(MI,NJ),WP)) - 1)
      IF (PRESENT(T_RECOMP)) TT0 = OMP_GET_WTIME()
      IF (PRESENT(WS)) THEN
        CALL FAST_RECOMPRESS_ARENA(MI, NJ, ABLK%RANK, ABLK%UFAC, ABLK%VFAC, &
                              PBLK%RANK, PBLK%UFAC, PBLK%VFAC, TOL, RMAX, NEWBLK, OK, ARENA, &
                              ABLK%UFAC_FLAT, ABLK%VFAC_FLAT, WS)
      ELSE
        CALL FAST_RECOMPRESS_ARENA(MI, NJ, ABLK%RANK, ABLK%UFAC, ABLK%VFAC, &
                              PBLK%RANK, PBLK%UFAC, PBLK%VFAC, TOL, RMAX, NEWBLK, OK, ARENA, &
                              ABLK%UFAC_FLAT, ABLK%VFAC_FLAT)
      END IF
      IF (PRESENT(T_RECOMP)) T_RECOMP = T_RECOMP + (OMP_GET_WTIME() - TT0)
      IF (OK) THEN
        ABLK%IS_LOWRANK = .TRUE.
        ABLK%M = MI; ABLK%N = NJ; ABLK%RANK = NEWBLK%RANK
        ABLK%UFAC => NEWBLK%UFAC   ! copie de POINTEUR seulement
        ABLK%VFAC => NEWBLK%VFAC
        NULLIFY(ABLK%DENSE)
        RETURN
      END IF
      ! rang recompresse trop grand -> tombe dans le repli dense ci-dessous
    END IF

    ! cas 1 (ou repli du cas 2) : A ou L*U non basse-rang -> dense.
    !
    ! ABLK DEJA dense (le cas DE LOIN le plus frequent en pratique une fois
    ! qu'un bloc est tombe en repli -- ce meme algorithme, arena ou non,
    ! n'y revient JAMAIS ensuite : la branche bas-rang ci-dessus exige
    ! ABLK%IS_LOWRANK, donc un bloc dense le reste pour toujours) : mise a
    ! jour EN PLACE via SUBTRACT_PRODUCT_DENSE (meme principe que
    ! SCHUR_UPDATE_BLOCK pour les blocs diagonaux, CGEMM BETA=1) --
    ! ABLK%DENSE_FLAT (deja de la BONNE taille MI*NJ, fixe pour ce slot de
    ! bloc) est reutilise TEL QUEL, ZERO carve, ZERO copie intermediaire.
    ! Sans cette branche, EXPAND_BLOCK_ARENA recopierait ABLK%DENSE dans
    ! un nouveau buffer a chaque etape K ou ce bloc est touche, alors que
    ! son contenu est deja exactement la ou il doit etre.
    IF (.NOT. ABLK%IS_LOWRANK) THEN
      IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
      CALL SUBTRACT_PRODUCT_DENSE(MI, NJ, PBLK, ABLK%DENSE)
      ABLK%M = MI; ABLK%N = NJ
      IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)
      RETURN
    END IF

    ! ABLK actuellement bas-rang -> PREMIERE bascule vers dense pour ce
    ! bloc (n'arrivera plus jamais ensuite, voir ci-dessus). IMPORTANT
    ! (ordre) : EXPAND_BLOCK_ARENA doit lire l'ANCIEN contenu de ABLK
    ! (UFAC/VFAC) AVANT qu'on ne redirige ses pointeurs -- DENSE_VIEW est
    ! carve/reutilise depuis ABLK%DENSE_FLAT, un buffer TOUJOURS DISTINCT
    ! de UFAC/VFAC (jamais utilise tant que le bloc etait bas-rang), donc
    ! sans aucun risque d'aliasing avec la lecture qui suit.
    IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
    FLAT => ARENA_CARVE_REUSE(ARENA, ABLK%DENSE_FLAT, INT(MI,8)*INT(NJ,8))
    DENSE_VIEW(1:MI,1:NJ) => FLAT
    CALL EXPAND_BLOCK_ARENA(ABLK, DENSE_VIEW, MI)
    CALL SUBTRACT_PRODUCT_DENSE(MI, NJ, PBLK, DENSE_VIEW)
    ABLK%IS_LOWRANK = .FALSE.
    ABLK%M = MI; ABLK%N = NJ
    NULLIFY(ABLK%UFAC, ABLK%VFAC)
    ABLK%DENSE => DENSE_VIEW
    IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)
  END SUBROUTINE SCHUR_UPDATE_BLR_BLOCK_ARENA

  !=========================================================================
  ! VARIANTES FIXEDCAP (BLR_BLOCK_FIXEDCAP_T) -- voir le type ci-dessus
  ! pour la justification complete. Contrairement a l'arena (composantes
  ! POINTER, necessitant l'attribut TARGET propage), UFAC/VFAC/DENSE sont
  ! ici de simples ALLOCATABLE -- aucune contrainte TARGET, affectation de
  ! tableau ordinaire.
  !=========================================================================

  !-----------------------------------------------------------------------
  ! Copie un bloc COMPRESS_BLOCK ordinaire (taille exacte) vers un bloc
  ! FIXEDCAP (capacite RMAX) -- utilise a l'assemblage (COMPRESS_BLOCK
  ! INCHANGEE, puis cette copie carve la capacite definitive du bloc pour
  ! TOUTE la suite de la factorisation). RMAX recalcule ICI (pas encore
  ! connu a l'assemblage) via COMPRESS_WORTHWHILE_RATIO -- MAX(...,
  ! SRC%RANK) en garde-fou : le seuil d'acceptation de COMPRESS_BLOCK
  ! (R < RATIO*MN, inegalite stricte) et la formule RMAX (INT(RATIO*MN)-1)
  ! peuvent differer de 1 par arrondi ; la capacite allouee doit rester AU
  ! MOINS aussi grande que le rang REELLEMENT accepte.
  !-----------------------------------------------------------------------
  SUBROUTINE COPY_REGULAR_TO_FIXEDCAP(SRC, DST)
    TYPE(BLR_BLOCK_T), INTENT(IN) :: SRC
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(OUT) :: DST
    DST%IS_LOWRANK = SRC%IS_LOWRANK
    DST%M = SRC%M; DST%N = SRC%N; DST%RANK = SRC%RANK
    IF (SRC%IS_LOWRANK) THEN
      DST%RMAX = MAX(1, INT(COMPRESS_WORTHWHILE_RATIO*REAL(MIN(SRC%M,SRC%N),WP)) - 1)
      DST%RMAX = MAX(DST%RMAX, SRC%RANK)
      ALLOCATE(DST%UFAC(SRC%M,DST%RMAX), DST%VFAC(DST%RMAX,SRC%N))
      DST%UFAC(1:SRC%M,1:SRC%RANK) = SRC%UFAC
      DST%VFAC(1:SRC%RANK,1:SRC%N) = SRC%VFAC
    ELSE
      ALLOCATE(DST%DENSE(SRC%M,SRC%N))
      DST%DENSE = SRC%DENSE
    END IF
  END SUBROUTINE COPY_REGULAR_TO_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Copie profonde d'un bloc FIXEDCAP vers un BLR_BLOCK_T ordinaire,
  ! EXACTEMENT dimensionne au rang REEL (ignore la capacite RMAX) --
  ! utilise avant diffusion MPI (BCAST_BLR_BLOCK, INCHANGEE, ne sait pas
  ! lire un BLR_BLOCK_FIXEDCAP_T) : le volume envoye ne depend QUE du rang
  ! reel, jamais de la capacite allouee.
  !-----------------------------------------------------------------------
  SUBROUTINE COPY_FIXEDCAP_TO_REGULAR(SRC, DST)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(IN) :: SRC
    TYPE(BLR_BLOCK_T), INTENT(OUT)         :: DST
    DST%IS_LOWRANK = SRC%IS_LOWRANK
    DST%M = SRC%M; DST%N = SRC%N; DST%RANK = SRC%RANK
    IF (SRC%IS_LOWRANK) THEN
      ALLOCATE(DST%UFAC(SRC%M,SRC%RANK), DST%VFAC(SRC%RANK,SRC%N))
      DST%UFAC = SRC%UFAC(1:SRC%M,1:SRC%RANK)
      DST%VFAC = SRC%VFAC(1:SRC%RANK,1:SRC%N)
    ELSE
      ALLOCATE(DST%DENSE(SRC%M,SRC%N))
      DST%DENSE = SRC%DENSE
    END IF
  END SUBROUTINE COPY_FIXEDCAP_TO_REGULAR

  !-----------------------------------------------------------------------
  ! Equivalent EXPAND_BLOCK -- BLK%UFAC/VFAC sont declares (M,RMAX)/
  ! (RMAX,N) mais seules les RANK premieres colonnes/lignes sont valides :
  ! la vue tranchee VFAC(1:RANK,:) force Fortran a passer une copie
  ! temporaire CORRECTEMENT dimensionnee/stridee a CGEMM (le premier
  ! indice de VFAC, paddé a RMAX, ne peut pas etre utilise tel quel comme
  ! leading dimension) ; UFAC(1:M,1:RANK) est deja contigu (seul le
  ! SECOND indice est paddé), la vue est gratuite (aucune copie).
  !-----------------------------------------------------------------------
  SUBROUTINE EXPAND_BLOCK_FIXEDCAP(BLK, OUT, LDOUT)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(IN) :: BLK
    INTEGER, INTENT(IN)           :: LDOUT
    COMPLEX(WP), INTENT(OUT)      :: OUT(LDOUT,BLK%N)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)

    IF (BLK%IS_LOWRANK) THEN
      CALL CGEMM('No transpose', 'No transpose', BLK%M, BLK%N, BLK%RANK, &
                  CONE, BLK%UFAC(1:BLK%M,1:BLK%RANK), BLK%M, BLK%VFAC(1:BLK%RANK,1:BLK%N), BLK%RANK, &
                  CZERO, OUT, LDOUT)
    ELSE
      OUT(1:BLK%M,1:BLK%N) = BLK%DENSE
    END IF
  END SUBROUTINE EXPAND_BLOCK_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Equivalent BLR_BLOCK_ROW_SWAP -- corps STRICTEMENT identique cote
  ! bas-rang : le swap porte sur les LIGNES de UFAC (premier indice,
  ! JAMAIS paddé), le stride BLAS (KW = M) est inchange que UFAC soit
  ! declare (M,RANK) ou (M,RMAX). Cote dense, DENSE n'est de toute facon
  ! jamais paddé (taille fixe M*N, voir le type).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_BLOCK_ROW_SWAP_FIXEDCAP(BLK, KW, IPIV)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(INOUT) :: BLK
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
  END SUBROUTINE BLR_BLOCK_ROW_SWAP_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Equivalent BLR_TRSM_LEFT_LOWER_UNIT -- opere sur UFAC (premier indice
  ! M jamais paddé) : passer BLK%UFAC tel quel avec N=BLK%RANK suffit,
  ! CTRSM ne touche que les BLK%RANK premieres colonnes quel que soit
  ! RMAX -- AUCUN changement necessaire par rapport a la version
  ! BLR_BLOCK_T (seule la declaration de type de BLK differe).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_TRSM_LEFT_LOWER_UNIT_FIXEDCAP(DKK, KW, BLK)
    INTEGER, INTENT(IN)             :: KW
    COMPLEX(WP), INTENT(IN)         :: DKK(KW,KW)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(INOUT) :: BLK
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, BLK%RANK, CONE, DKK, KW, BLK%UFAC, KW)
    ELSE
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, BLK%N, CONE, DKK, KW, BLK%DENSE, KW)
    END IF
  END SUBROUTINE BLR_TRSM_LEFT_LOWER_UNIT_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Equivalent BLR_TRSM_RIGHT_UPPER_NONUNIT -- ICI, contrairement a UFAC,
  ! VFAC a son PREMIER indice (le rang) paddé a RMAX : passer BLK%VFAC
  ! TEL QUEL avec LDB=BLK%RANK serait FAUX (le stride physique reel est
  ! RMAX, pas RANK -- corruption memoire). La vue tranchee
  ! BLK%VFAC(1:BLK%RANK,:) force Fortran a passer une copie temporaire
  ! CORRECTEMENT dimensionnee (LDB=BLK%RANK devient alors exact pour cette
  ! copie), et a la recopier dans BLK%VFAC(1:BLK%RANK,:) en sortie (BLK
  ! est INTENT(INOUT), CTRSM modifie ce facteur en place).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_TRSM_RIGHT_UPPER_NONUNIT_FIXEDCAP(BLK, DKK, KW)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(INOUT) :: BLK
    INTEGER, INTENT(IN)             :: KW
    COMPLEX(WP), INTENT(IN)         :: DKK(KW,KW)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    IF (BLK%IS_LOWRANK) THEN
      CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', BLK%RANK, KW, CONE, DKK, KW, &
                  BLK%VFAC(1:BLK%RANK,1:KW), BLK%RANK)
    ELSE
      CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', BLK%M, KW, CONE, DKK, KW, BLK%DENSE, BLK%M)
    END IF
  END SUBROUTINE BLR_TRSM_RIGHT_UPPER_NONUNIT_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Equivalent BLR_BLOCK_APPLY -- lecture SEULE de UFAC/VFAC (INTENT(IN)),
  ! meme raisonnement que EXPAND_BLOCK_FIXEDCAP pour la vue tranchee de
  ! VFAC.
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_BLOCK_APPLY_FIXEDCAP(BLK, NRHS, ALPHA, X, LDX, BETA, Y, LDY)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(IN) :: BLK
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
                  (1.0_WP,0.0_WP), BLK%VFAC(1:BLK%RANK,1:BLK%N), BLK%RANK, X, LDX, CZERO, T, BLK%RANK)
      CALL CGEMM('No transpose', 'No transpose', BLK%M, NRHS, BLK%RANK, &
                  ALPHA, BLK%UFAC(1:BLK%M,1:BLK%RANK), BLK%M, T, BLK%RANK, BETA, Y, LDY)
      DEALLOCATE(T)
    ELSE
      CALL CGEMM('No transpose', 'No transpose', BLK%M, NRHS, BLK%N, &
                  ALPHA, BLK%DENSE, BLK%M, X, LDX, BETA, Y, LDY)
    END IF
  END SUBROUTINE BLR_BLOCK_APPLY_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Equivalent SCHUR_UPDATE_BLR_BLOCK -- ABLK est ici BLR_BLOCK_FIXEDCAP_T
  ! (stockage PERSISTANT, capacite RMAX fixee a la premiere compression).
  ! LBLK/UBLK (panneaux diffuses, transitoires) et PBLK/NEWBLK (locaux a
  ! cet appel) restent des BLR_BLOCK_T ORDINAIRES -- FAST_RECOMPRESS
  ! (INCHANGEE) est reutilisee TELLE QUELLE en lui passant des VUES
  ! tranchees de ABLK%UFAC/VFAC comme UA/VA (dummy arguments a forme
  ! EXPLICITE UA(MI,RA)/VA(RA,NJ) : Fortran gere automatiquement la copie
  ! temporaire correcte pour VA, aucune gestion manuelle de leading
  ! dimension necessaire) -- elle produit un NEWBLK ORDINAIRE (taille
  ! exacte, transitoire), qu'on RECOPIE ENSUITE dans ABLK%UFAC/VFAC EN
  ! PLACE (ECRASEMENT, PAS de nouvelle ALLOCATE) : NEWBLK%RANK<=ABLK%RMAX
  ! est GARANTI par le test RTOT.GT.RMAX de FAST_RECOMPRESS -- c'est la
  ! PROPRIETE CLE de cette variante, AUCUNE reallocation de la memoire
  ! persistante apres la toute premiere compression du bloc.
  !
  ! Repli dense : une fois dense, TOUJOURS dense ensuite (meme garde
  ! ABLK%IS_LOWRANK que la version non-FIXEDCAP) -- mise a jour directe EN
  ! PLACE via SUBTRACT_PRODUCT_DENSE des le 2e passage (ZERO allocation,
  ! comme SCHUR_UPDATE_BLOCK pour les blocs diagonaux) ; la PREMIERE
  ! bascule alloue ABLK%DENSE(MI,NJ) (taille fixe, jamais realloue
  ! ensuite) et LIBERE UFAC/VFAC (rendus a glibc -- un chunk de taille
  ! COMMUNE a tous les blocs du meme calibre, immediatement reutilisable
  ! par la PROCHAINE premiere compression d'un autre bloc, contrairement a
  ! l'arena qui ne recupere jamais cet espace).
  !-----------------------------------------------------------------------
  SUBROUTINE SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP(ABLK, MI, NJ, KW, LBLK, UBLK, TOL, WS, T_PRODUCT, T_RECOMP)
    TYPE(BLR_BLOCK_FIXEDCAP_T), INTENT(INOUT) :: ABLK
    INTEGER, INTENT(IN)              :: MI, NJ, KW
    TYPE(BLR_BLOCK_T), INTENT(IN)    :: LBLK, UBLK
    REAL(WP), INTENT(IN)             :: TOL
    TYPE(BLR_WORKSPACE_T), INTENT(INOUT), OPTIONAL, TARGET :: WS
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: T_PRODUCT, T_RECOMP

    TYPE(BLR_BLOCK_T) :: PBLK, NEWBLK
    LOGICAL :: OK
    DOUBLE PRECISION :: TT0

    IF (MI.EQ.0 .OR. NJ.EQ.0 .OR. KW.EQ.0) RETURN

    IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
    CALL COMPUTE_LU_PRODUCT(MI, NJ, KW, LBLK, UBLK, PBLK)
    IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)

    IF (ABLK%IS_LOWRANK .AND. PBLK%IS_LOWRANK) THEN
      IF (PRESENT(T_RECOMP)) TT0 = OMP_GET_WTIME()
      IF (PRESENT(WS)) THEN
        CALL FAST_RECOMPRESS(MI, NJ, ABLK%RANK, ABLK%UFAC(1:MI,1:ABLK%RANK), ABLK%VFAC(1:ABLK%RANK,1:NJ), &
                              PBLK%RANK, PBLK%UFAC, PBLK%VFAC, TOL, ABLK%RMAX, NEWBLK, OK, WS)
      ELSE
        CALL FAST_RECOMPRESS(MI, NJ, ABLK%RANK, ABLK%UFAC(1:MI,1:ABLK%RANK), ABLK%VFAC(1:ABLK%RANK,1:NJ), &
                              PBLK%RANK, PBLK%UFAC, PBLK%VFAC, TOL, ABLK%RMAX, NEWBLK, OK)
      END IF
      IF (PRESENT(T_RECOMP)) T_RECOMP = T_RECOMP + (OMP_GET_WTIME() - TT0)
      IF (OK) THEN
        ABLK%IS_LOWRANK = .TRUE.
        ABLK%RANK = NEWBLK%RANK
        ABLK%UFAC(1:MI,1:NEWBLK%RANK) = NEWBLK%UFAC   ! ECRASEMENT en place, PAS de (re)ALLOCATE
        ABLK%VFAC(1:NEWBLK%RANK,1:NJ) = NEWBLK%VFAC
        RETURN
      END IF
      ! rang recompresse trop grand -> tombe dans le repli dense ci-dessous
    END IF

    IF (.NOT. ABLK%IS_LOWRANK) THEN
      ! deja dense (cas le plus frequent une fois la bascule faite) :
      ! mise a jour EN PLACE, ZERO allocation.
      IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
      CALL SUBTRACT_PRODUCT_DENSE(MI, NJ, PBLK, ABLK%DENSE)
      IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)
      RETURN
    END IF

    ! PREMIERE bascule bas-rang -> dense pour ce bloc (n'arrivera plus
    ! jamais ensuite, voir ci-dessus).
    IF (PRESENT(T_PRODUCT)) TT0 = OMP_GET_WTIME()
    ALLOCATE(ABLK%DENSE(MI,NJ))
    CALL EXPAND_BLOCK_FIXEDCAP(ABLK, ABLK%DENSE, MI)
    CALL SUBTRACT_PRODUCT_DENSE(MI, NJ, PBLK, ABLK%DENSE)
    ABLK%IS_LOWRANK = .FALSE.
    DEALLOCATE(ABLK%UFAC, ABLK%VFAC)
    IF (PRESENT(T_PRODUCT)) T_PRODUCT = T_PRODUCT + (OMP_GET_WTIME() - TT0)
  END SUBROUTINE SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP

END MODULE BLR_COMPRESS_MOD
