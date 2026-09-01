!===============================================================================
! OUT-OF-CORE (OOC) -- stockage disque des blocs BLR persistants pour la
! variante MPI_BLR_MATRIX_OOC_T (mpi_blr_matrix_mod.f90), demandee
! explicitement pour ne plus etre limite par la RAM disponible (ni pour la
! matrice assemblee, ni pour les facteurs LU apres remplissage de Schur, ni
! pour un second membre massif).
!
! PRINCIPE GENERAL : chaque processus MPI possede son PROPRE fichier scratch
! local (jamais partage -- la distribution 2D-cyclique existante garantit
! deja que chaque processus ne touche que ses propres blocs locaux, aucune
! coordination inter-processus n'est necessaire pour l'OOC lui-meme).
! Fichier ouvert en ACCES FLUX (ACCESS='STREAM').
!
! DEUX POLITIQUES D'ECRITURE ont coexiste au cours du developpement de ce
! module -- SEULE la seconde est utilisee desormais (DIAGBLK et OFFBLK) :
!   - APPEND-ONLY (OOC_WRITE_BLOCK_RAW) : un bloc reecrit ajoute un NOUVEL
!     enregistrement en fin de fichier plutot que d'ecraser en place --
!     l'ancien enregistrement devient orphelin. Encore utilisee pour le
!     second membre (OOC_RHS_T, prefixe "b") : chaque paquet est ecrit
!     UNE SEULE FOIS, jamais reecrit, l'aspect "orphelin" ne s'applique
!     donc jamais la -- append-only y est simplement le plus simple.
!   - EN PLACE, OFFSET FIXE (OOC_WRITE_BLOCK_FIXED) : reecrit TOUJOURS au
!     meme emplacement. Premier appel : reserve le PIRE CAS DENSE (M*N
!     elements, jamais la taille bas-rang courante) -- garantit qu'un bloc
!     dont la taille grossit a une reecriture ulterieure (mise a jour de
!     Schur : le rang n'est PAS garanti decroissant, ni pour un bloc
!     diagonal -- toujours dense, taille constante -- ni pour un bloc
!     bas-rang, ou COMPRESS_WORTHWHILE_RATIO garantit juste que la forme
!     bas-rang reste plus petite que le pire cas dense) tient toujours
!     dans la place reservee. Appels suivants : n'ecrit que les octets du
!     rang COURANT, jamais tout le pire cas reserve -- aucun octet
!     orphelin, contrairement a l'append-only. Sans cette politique, un
!     bloc touche par NBLK mises a jour de Schur gonflerait le fichier en
!     O(NBLK^2), soit O(N^2) -- mesure et documente dans l'historique du
!     projet, initialement seulement pour DIAGBLK puis etendu a OFFBLK.
!
! TROIS NIVEAUX AU-DESSUS DU FICHIER :
!   1) OOC_STORE_T / OOC_WRITE_BLOCK_RAW / OOC_WRITE_BLOCK_FIXED /
!      OOC_READ_BLOCK_RAW : le fichier lui-meme + la (des)serialisation
!      d'un BLR_BLOCK_T (les deux routines d'ecriture partagent le meme
!      format binaire et la meme lecture -- seul l'endroit ou l'on ecrit
!      differe).
!   2) OOC_CACHE_T : cache borne en octets (MAX_MB_PER_PROCESS), PAR
!      (li,lj), politique LRU simple -- pour les blocs hors-diagonale
!      (OFFBLK), ecrits EN PLACE (OOC_WRITE_BLOCK_FIXED) a l'eviction.
!   3) OOC_DIAG_CACHE_T : meme principe LRU, cle par un entier K -- pour
!      les blocs diagonaux (DIAGBLK), ecrits EN PLACE (OFFSET FIXE) a
!      l'eviction. Sans ce cache, un bloc diagonal serait charge/modifie/
!      reecrit a CHAQUE mise a jour de Schur qui le touche (une par etape
!      K<I) ; avec, il reste resident entre deux touches et n'est ecrit
!      qu'a son eviction (LRU ou flush final).
! Les deux caches sont orthogonaux au choix d'ecriture : ils decident QUAND
! toucher le disque, pas OU.
!===============================================================================
MODULE BLR_OOC_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE BLR_COMPRESS_MOD, ONLY: BLR_BLOCK_T, BLR_BLOCK_BYTES
  IMPLICIT NONE

  !---------------------------------------------------------------------------
  ! Descripteur PERSISTANT (toujours resident, tres petit) d'un bloc stocke
  ! sur disque -- remplace UFAC/VFAC/DENSE dans MPI_BLR_MATRIX_OOC_T%OFFBLK/
  ! DIAGBLK. OFFSET=-1 signifie "jamais encore ecrit".
  !---------------------------------------------------------------------------
  TYPE :: OOC_BLOCK_DESC_T
    LOGICAL      :: IS_LOWRANK = .FALSE.
    INTEGER      :: M = 0, N = 0, RANK = 0
    INTEGER(8)   :: OFFSET = -1_8
  END TYPE OOC_BLOCK_DESC_T

  !---------------------------------------------------------------------------
  ! Fichier scratch d'UN processus -- ACCESS='STREAM' (positionnement par
  ! octet via POS=), toujours ouvert en ecriture-ajout + lecture directe.
  !---------------------------------------------------------------------------
  TYPE :: OOC_STORE_T
    INTEGER              :: UNIT = -1
    CHARACTER(LEN=512)   :: PATH = ''
    INTEGER(8)           :: NEXT_OFFSET = 0_8
  END TYPE OOC_STORE_T

  INTEGER(8), PARAMETER :: OOC_BYTES_PER_ELEM = 8_8   ! COMPLEX(WP), WP=simple precision

  !---------------------------------------------------------------------------
  ! Un emplacement du cache -- bloc RESIDENT (BLK, forme ordinaire
  ! BLR_BLOCK_T) tant qu'il n'a pas ete evince. LI=LJ=0 signifie emplacement
  ! libre. DIRTY : le bloc a ete modifie depuis son dernier chargement/
  ! ecriture, doit etre reecrit avant d'etre evince ou en fin de vie.
  !---------------------------------------------------------------------------
  TYPE :: OOC_CACHE_SLOT_T
    INTEGER            :: LI = 0, LJ = 0
    LOGICAL            :: DIRTY = .FALSE.
    INTEGER(8)         :: LAST_USE = 0_8
    INTEGER(8)         :: NBYTES = 0_8
    TYPE(BLR_BLOCK_T)  :: BLK
  END TYPE OOC_CACHE_SLOT_T

  ! SLOT_OF(li,lj) : indice du slot residant pour ce bloc, 0 si absent --
  ! lookup DIRECT O(1), remplace la recherche LINEAIRE O(NSLOTS) d'origine
  ! (limite connue "v1" : meme un ACCES DEJA RESIDENT payait un scan complet
  ! de tous les emplacements a chaque appel, cout qui ne depend PAS du taux
  ! d'eviction -- mesure comme significatif meme a budget genereux, voir la
  ! memoire du projet). FREE_STACK/NFREE : pile des emplacements libres,
  ! evite AUSSI un scan lineaire pour trouver ou inserer un nouveau bloc.
  TYPE :: OOC_CACHE_T
    TYPE(OOC_CACHE_SLOT_T), ALLOCATABLE :: SLOTS(:)
    INTEGER, ALLOCATABLE :: SLOT_OF(:,:)
    INTEGER, ALLOCATABLE :: FREE_STACK(:)
    INTEGER      :: NFREE = 0
    INTEGER      :: NSLOTS = 0
    INTEGER(8)   :: BUDGET_BYTES = 0_8
    INTEGER(8)   :: USED_BYTES = 0_8
    INTEGER(8)   :: CLOCK = 0_8
  END TYPE OOC_CACHE_T

  !---------------------------------------------------------------------------
  ! Cache DEDIE aux blocs DIAGONAUX (DIAGBLK) -- meme principe LRU borne en
  ! octets que OOC_CACHE_T ci-dessus, mais cle par un SEUL entier K (indice
  ! de bloc global), pas par une paire (li,lj) : DIAGBLK est un tableau 1D,
  ! Fortran exige un type de descripteur EXACT (voir LOCAL_ROW_OFFSET_OOC
  ! dans mpi_blr_matrix_mod.f90 pour le meme genre de duplication deja
  ! acceptee dans ce module plutot que de forcer une forme generique).
  ! Ecrit a l'eviction via OOC_WRITE_BLOCK_FIXED (EN PLACE, voir plus haut)
  ! -- un bloc diagonal etant toujours dense et de taille fixe pour un K
  ! donne, aucun octet orphelin.
  !---------------------------------------------------------------------------
  TYPE :: OOC_DIAG_CACHE_SLOT_T
    INTEGER            :: K = 0
    LOGICAL            :: DIRTY = .FALSE.
    INTEGER(8)         :: LAST_USE = 0_8
    INTEGER(8)         :: NBYTES = 0_8
    TYPE(BLR_BLOCK_T)  :: BLK
  END TYPE OOC_DIAG_CACHE_SLOT_T

  ! SLOT_OF(K) : meme principe de lookup direct O(1) que OOC_CACHE_T ci-
  ! dessus, indexe par l'entier K plutot que par (li,lj).
  TYPE :: OOC_DIAG_CACHE_T
    TYPE(OOC_DIAG_CACHE_SLOT_T), ALLOCATABLE :: SLOTS(:)
    INTEGER, ALLOCATABLE :: SLOT_OF(:)
    INTEGER, ALLOCATABLE :: FREE_STACK(:)
    INTEGER      :: NFREE = 0
    INTEGER      :: NSLOTS = 0
    INTEGER(8)   :: BUDGET_BYTES = 0_8
    INTEGER(8)   :: USED_BYTES = 0_8
    INTEGER(8)   :: CLOCK = 0_8
  END TYPE OOC_DIAG_CACHE_T

CONTAINS

  !-----------------------------------------------------------------------
  ! Cree le sous-dossier scratch (si besoin) et ouvre le fichier local a CE
  ! processus (un fichier par rang MPI -- voir commentaire de tete).
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_STORE_INIT(STORE, SCRATCH_DIR, PREFIX, RANK)
    TYPE(OOC_STORE_T), INTENT(OUT) :: STORE
    CHARACTER(LEN=*), INTENT(IN)   :: SCRATCH_DIR, PREFIX
    INTEGER, INTENT(IN)            :: RANK
    CHARACTER(LEN=16) :: RANKSTR
    LOGICAL :: EXISTS

    INQUIRE(FILE=TRIM(SCRATCH_DIR), EXIST=EXISTS)
    IF (.NOT. EXISTS) CALL EXECUTE_COMMAND_LINE('mkdir -p "'//TRIM(SCRATCH_DIR)//'"')

    WRITE(RANKSTR,'(I6.6)') RANK
    STORE%PATH = TRIM(SCRATCH_DIR)//'/'//TRIM(PREFIX)//'_rank_'//TRIM(RANKSTR)//'.dat'
    STORE%NEXT_OFFSET = 0_8
    OPEN(NEWUNIT=STORE%UNIT, FILE=TRIM(STORE%PATH), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='READWRITE')
  END SUBROUTINE OOC_STORE_INIT

  SUBROUTINE OOC_STORE_FREE(STORE)
    TYPE(OOC_STORE_T), INTENT(INOUT) :: STORE
    IF (STORE%UNIT .GE. 0) CLOSE(STORE%UNIT, STATUS='DELETE')
    STORE%UNIT = -1
  END SUBROUTINE OOC_STORE_FREE

  !-----------------------------------------------------------------------
  ! Ecrit BLK (UFAC+VFAC concatenes si bas-rang, DENSE sinon) EN FIN DE
  ! FICHIER (append-only, voir commentaire de tete) et met a jour DESC.
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_WRITE_BLOCK_RAW(STORE, DESC, BLK)
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESC
    TYPE(BLR_BLOCK_T), INTENT(IN)         :: BLK
    INTEGER(8) :: NELEMS

    DESC%IS_LOWRANK = BLK%IS_LOWRANK
    DESC%M = BLK%M; DESC%N = BLK%N; DESC%RANK = BLK%RANK
    DESC%OFFSET = STORE%NEXT_OFFSET

    IF (BLK%IS_LOWRANK) THEN
      NELEMS = INT(BLK%M,8)*INT(BLK%RANK,8) + INT(BLK%RANK,8)*INT(BLK%N,8)
      IF (NELEMS .GT. 0_8) THEN
        WRITE(STORE%UNIT, POS=STORE%NEXT_OFFSET+1_8) BLK%UFAC, BLK%VFAC
      END IF
    ELSE
      NELEMS = INT(BLK%M,8)*INT(BLK%N,8)
      IF (NELEMS .GT. 0_8) WRITE(STORE%UNIT, POS=STORE%NEXT_OFFSET+1_8) BLK%DENSE
    END IF
    STORE%NEXT_OFFSET = STORE%NEXT_OFFSET + NELEMS*OOC_BYTES_PER_ELEM
  END SUBROUTINE OOC_WRITE_BLOCK_RAW

  !-----------------------------------------------------------------------
  ! Ecrit BLK EN PLACE, a un OFFSET FIXE reutilise a chaque appel. Premier
  ! appel (DESC%OFFSET=-1) : reserve TOUJOURS le PIRE CAS DENSE (M*N
  ! elements, cf. BLK%M/BLK%N -- la forme dense pleine du bloc, JAMAIS
  ! renseignes a la valeur du rang) en fin de fichier, meme si l'ecriture
  ! elle-meme est plus petite (bloc bas-rang) -- garantit qu'un rang qui
  ! grossit a une reecriture ulterieure (mise a jour de Schur : le rang
  ! n'est PAS garanti decroissant) tient toujours dans la place reservee,
  ! puisque COMPRESS_WORTHWHILE_RATIO (voir blr_compress_mod.f90, <=1/3 par
  ! defaut) garantit M*rang+rang*N < M*N des qu'un bloc reste bas-rang --
  ! le repli dense, lui, occupe EXACTEMENT M*N, jamais plus. Appels
  ! suivants : REECRIT AU MEME OFFSET, seulement les octets REELLEMENT
  ! utiles au rang COURANT (pas tout le pire cas reserve) -- aucun octet
  ! orphelin, par opposition a OOC_WRITE_BLOCK_RAW qui ajoute TOUJOURS un
  ! nouvel enregistrement. Utilise pour les blocs DIAGONAUX (toujours
  ! denses, M*N==NELEMS, reservation = ecriture reelle, comportement
  ! inchange) ET desormais pour les blocs hors-diagonale OFFBLK (voir
  ! OOC_CACHE_EVICT_SLOT plus bas).
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_WRITE_BLOCK_FIXED(STORE, DESC, BLK)
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESC
    TYPE(BLR_BLOCK_T), INTENT(IN)         :: BLK
    INTEGER(8) :: NELEMS, RESERVE_NELEMS, WRITE_POS

    DESC%IS_LOWRANK = BLK%IS_LOWRANK
    DESC%M = BLK%M; DESC%N = BLK%N; DESC%RANK = BLK%RANK

    IF (BLK%IS_LOWRANK) THEN
      NELEMS = INT(BLK%M,8)*INT(BLK%RANK,8) + INT(BLK%RANK,8)*INT(BLK%N,8)
    ELSE
      NELEMS = INT(BLK%M,8)*INT(BLK%N,8)
    END IF
    RESERVE_NELEMS = INT(BLK%M,8)*INT(BLK%N,8)   ! pire cas dense, toujours

    IF (DESC%OFFSET .LT. 0_8) THEN
      WRITE_POS = STORE%NEXT_OFFSET
      DESC%OFFSET = WRITE_POS
      STORE%NEXT_OFFSET = STORE%NEXT_OFFSET + RESERVE_NELEMS*OOC_BYTES_PER_ELEM
    ELSE
      WRITE_POS = DESC%OFFSET
    END IF

    IF (NELEMS .GT. 0_8) THEN
      IF (BLK%IS_LOWRANK) THEN
        WRITE(STORE%UNIT, POS=WRITE_POS+1_8) BLK%UFAC, BLK%VFAC
      ELSE
        WRITE(STORE%UNIT, POS=WRITE_POS+1_8) BLK%DENSE
      END IF
    END IF
  END SUBROUTINE OOC_WRITE_BLOCK_FIXED

  !-----------------------------------------------------------------------
  ! Relit un bloc depuis DESC%OFFSET -- alloue BLK a la forme decrite par
  ! DESC (M/N/RANK/IS_LOWRANK, deja connus SANS lire le fichier).
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_READ_BLOCK_RAW(STORE, DESC, BLK)
    TYPE(OOC_STORE_T), INTENT(INOUT)   :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(IN) :: DESC
    TYPE(BLR_BLOCK_T), INTENT(OUT)     :: BLK

    BLK%IS_LOWRANK = DESC%IS_LOWRANK
    BLK%M = DESC%M; BLK%N = DESC%N; BLK%RANK = DESC%RANK
    IF (DESC%IS_LOWRANK) THEN
      ALLOCATE(BLK%UFAC(DESC%M,DESC%RANK), BLK%VFAC(DESC%RANK,DESC%N))
      IF (DESC%M*DESC%RANK+DESC%RANK*DESC%N .GT. 0) THEN
        READ(STORE%UNIT, POS=DESC%OFFSET+1_8) BLK%UFAC, BLK%VFAC
      END IF
    ELSE
      ALLOCATE(BLK%DENSE(DESC%M,DESC%N))
      IF (DESC%M*DESC%N .GT. 0) READ(STORE%UNIT, POS=DESC%OFFSET+1_8) BLK%DENSE
    END IF
  END SUBROUTINE OOC_READ_BLOCK_RAW

  !-----------------------------------------------------------------------
  ! Cache borne en octets pour les blocs hors-diagonale (OFFBLK) -- NSLOTS
  ! dimensionne genereusement (le budget en octets est la VRAIE limite,
  ! NSLOTS n'est qu'un plafond structurel sur le nombre d'entrees
  ! simultanees, jamais atteint en pratique si BUDGET_MB est realiste face
  ! a la taille des blocs). NROWS/NCOLS : forme de DESCS(:,:) (donc de
  ! A_BLR%OFFBLK), necessaire pour dimensionner SLOT_OF.
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_CACHE_INIT(CACHE, BUDGET_MB, MAX_SLOTS, NROWS, NCOLS)
    TYPE(OOC_CACHE_T), INTENT(OUT) :: CACHE
    REAL(WP), INTENT(IN)           :: BUDGET_MB
    INTEGER, INTENT(IN)            :: MAX_SLOTS, NROWS, NCOLS
    INTEGER :: S
    CACHE%NSLOTS = MAX(MAX_SLOTS, 1)
    ALLOCATE(CACHE%SLOTS(CACHE%NSLOTS))
    ALLOCATE(CACHE%SLOT_OF(MAX(NROWS,1), MAX(NCOLS,1)))
    CACHE%SLOT_OF = 0
    ALLOCATE(CACHE%FREE_STACK(CACHE%NSLOTS))
    DO S = 1, CACHE%NSLOTS
      CACHE%FREE_STACK(S) = CACHE%NSLOTS - S + 1
    END DO
    CACHE%NFREE = CACHE%NSLOTS
    CACHE%BUDGET_BYTES = INT(BUDGET_MB*1024.0_WP*1024.0_WP, 8)
    CACHE%USED_BYTES = 0_8
    CACHE%CLOCK = 0_8
  END SUBROUTINE OOC_CACHE_INIT

  !-----------------------------------------------------------------------
  ! Ecrit (si DIRTY) et libere l'emplacement S -- DESCS(LI,LJ) est mis a
  ! jour par l'ecriture (nouvel offset si le bloc a change de taille).
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_CACHE_EVICT_SLOT(CACHE, STORE, DESCS, S)
    TYPE(OOC_CACHE_T), INTENT(INOUT)        :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)        :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT)   :: DESCS(:,:)
    INTEGER, INTENT(IN)                     :: S

    IF (CACHE%SLOTS(S)%LI .EQ. 0) RETURN
    IF (CACHE%SLOTS(S)%DIRTY) THEN
      CALL OOC_WRITE_BLOCK_FIXED(STORE, DESCS(CACHE%SLOTS(S)%LI,CACHE%SLOTS(S)%LJ), CACHE%SLOTS(S)%BLK)
    END IF
    IF (ALLOCATED(CACHE%SLOTS(S)%BLK%UFAC)) DEALLOCATE(CACHE%SLOTS(S)%BLK%UFAC)
    IF (ALLOCATED(CACHE%SLOTS(S)%BLK%VFAC)) DEALLOCATE(CACHE%SLOTS(S)%BLK%VFAC)
    IF (ALLOCATED(CACHE%SLOTS(S)%BLK%DENSE)) DEALLOCATE(CACHE%SLOTS(S)%BLK%DENSE)
    CACHE%USED_BYTES = CACHE%USED_BYTES - CACHE%SLOTS(S)%NBYTES
    CACHE%SLOT_OF(CACHE%SLOTS(S)%LI, CACHE%SLOTS(S)%LJ) = 0
    CACHE%SLOTS(S)%LI = 0; CACHE%SLOTS(S)%LJ = 0; CACHE%SLOTS(S)%DIRTY = .FALSE.; CACHE%SLOTS(S)%NBYTES = 0_8
    CACHE%NFREE = CACHE%NFREE + 1
    CACHE%FREE_STACK(CACHE%NFREE) = S
  END SUBROUTINE OOC_CACHE_EVICT_SLOT

  !-----------------------------------------------------------------------
  ! Libere de la place pour NEED_BYTES supplementaires -- evince le moins
  ! recemment utilise (plus petit LAST_USE parmi les emplacements occupes)
  ! jusqu'a ce que USED_BYTES+NEED_BYTES <= BUDGET_BYTES, ou qu'il n'y ait
  ! plus rien a evincer (accepte alors de depasser le budget plutot que de
  ! planter -- un budget trop petit pour UN SEUL bloc reste correct, juste
  ! plus lent que prevu).
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_CACHE_MAKE_ROOM(CACHE, STORE, DESCS, NEED_BYTES)
    TYPE(OOC_CACHE_T), INTENT(INOUT)      :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESCS(:,:)
    INTEGER(8), INTENT(IN)                :: NEED_BYTES
    INTEGER :: S, LRU_S
    INTEGER(8) :: LRU_VAL

    DO WHILE (CACHE%USED_BYTES + NEED_BYTES .GT. CACHE%BUDGET_BYTES)
      LRU_S = 0; LRU_VAL = HUGE(LRU_VAL)
      DO S = 1, CACHE%NSLOTS
        IF (CACHE%SLOTS(S)%LI .NE. 0 .AND. CACHE%SLOTS(S)%LAST_USE .LT. LRU_VAL) THEN
          LRU_VAL = CACHE%SLOTS(S)%LAST_USE; LRU_S = S
        END IF
      END DO
      IF (LRU_S .EQ. 0) RETURN   ! rien de plus a evincer
      CALL OOC_CACHE_EVICT_SLOT(CACHE, STORE, DESCS, LRU_S)
    END DO
  END SUBROUTINE OOC_CACHE_MAKE_ROOM

  !-----------------------------------------------------------------------
  ! Retourne l'indice de l'emplacement de cache contenant le bloc (LI,LJ),
  ! le chargeant depuis DESCS/STORE s'il n'est pas deja resident (evince
  ! au besoin pour respecter le budget AVANT de charger). Un bloc jamais
  ! encore ecrit (DESC%OFFSET=-1, ex. juste apres COMPRESS_BLOCK sans
  ! passage disque) doit etre INSERE directement via OOC_CACHE_PUT_NEW,
  ! pas via cette fonction.
  !-----------------------------------------------------------------------
  FUNCTION OOC_CACHE_GET(CACHE, STORE, DESCS, LI, LJ) RESULT(S)
    TYPE(OOC_CACHE_T), INTENT(INOUT)      :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESCS(:,:)
    INTEGER, INTENT(IN)                   :: LI, LJ
    INTEGER :: S, FREE_S
    INTEGER(8) :: NBYTES

    S = CACHE%SLOT_OF(LI,LJ)
    IF (S .NE. 0) THEN
      CACHE%CLOCK = CACHE%CLOCK + 1_8
      CACHE%SLOTS(S)%LAST_USE = CACHE%CLOCK
      RETURN
    END IF

    ! non resident -- charge depuis le disque
    NBYTES = 8_8*(INT(DESCS(LI,LJ)%M,8)*INT(MAX(DESCS(LI,LJ)%RANK,DESCS(LI,LJ)%N),8) + &
                  INT(DESCS(LI,LJ)%RANK,8)*INT(DESCS(LI,LJ)%N,8))
    IF (.NOT. DESCS(LI,LJ)%IS_LOWRANK) NBYTES = 8_8*INT(DESCS(LI,LJ)%M,8)*INT(DESCS(LI,LJ)%N,8)
    CALL OOC_CACHE_MAKE_ROOM(CACHE, STORE, DESCS, NBYTES)

    IF (CACHE%NFREE .EQ. 0) THEN
      WRITE(6,'(A)') 'ERREUR OOC_CACHE_GET : plus d''emplacement libre (augmenter MAX_SLOTS).'
      STOP 1
    END IF
    FREE_S = CACHE%FREE_STACK(CACHE%NFREE)
    CACHE%NFREE = CACHE%NFREE - 1

    CALL OOC_READ_BLOCK_RAW(STORE, DESCS(LI,LJ), CACHE%SLOTS(FREE_S)%BLK)
    CACHE%SLOTS(FREE_S)%LI = LI; CACHE%SLOTS(FREE_S)%LJ = LJ; CACHE%SLOTS(FREE_S)%DIRTY = .FALSE.
    CACHE%SLOTS(FREE_S)%NBYTES = NBYTES
    CACHE%SLOT_OF(LI,LJ) = FREE_S
    CACHE%USED_BYTES = CACHE%USED_BYTES + NBYTES
    CACHE%CLOCK = CACHE%CLOCK + 1_8
    CACHE%SLOTS(FREE_S)%LAST_USE = CACHE%CLOCK
    S = FREE_S
  END FUNCTION OOC_CACHE_GET

  !-----------------------------------------------------------------------
  ! Insere un bloc FRAICHEMENT calcule (jamais encore sur disque) dans le
  ! cache -- MOVE_ALLOC transfere la propriete des tableaux de BLK (zero
  ! copie), BLK ressort vide (comme apres un MOVE_ALLOC standard).
  !-----------------------------------------------------------------------
  FUNCTION OOC_CACHE_PUT_NEW(CACHE, STORE, DESCS, LI, LJ, BLK) RESULT(S)
    TYPE(OOC_CACHE_T), INTENT(INOUT)      :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESCS(:,:)
    INTEGER, INTENT(IN)                   :: LI, LJ
    TYPE(BLR_BLOCK_T), INTENT(INOUT)      :: BLK
    INTEGER :: S, FREE_S
    INTEGER(8) :: NBYTES

    NBYTES = BLR_BLOCK_BYTES(BLK)
    CALL OOC_CACHE_MAKE_ROOM(CACHE, STORE, DESCS, NBYTES)
    IF (CACHE%NFREE .EQ. 0) THEN
      WRITE(6,'(A)') 'ERREUR OOC_CACHE_PUT_NEW : plus d''emplacement libre (augmenter MAX_SLOTS).'
      STOP 1
    END IF
    FREE_S = CACHE%FREE_STACK(CACHE%NFREE)
    CACHE%NFREE = CACHE%NFREE - 1

    CACHE%SLOTS(FREE_S)%BLK%IS_LOWRANK = BLK%IS_LOWRANK
    CACHE%SLOTS(FREE_S)%BLK%M = BLK%M; CACHE%SLOTS(FREE_S)%BLK%N = BLK%N; CACHE%SLOTS(FREE_S)%BLK%RANK = BLK%RANK
    IF (BLK%IS_LOWRANK) THEN
      CALL MOVE_ALLOC(BLK%UFAC, CACHE%SLOTS(FREE_S)%BLK%UFAC)
      CALL MOVE_ALLOC(BLK%VFAC, CACHE%SLOTS(FREE_S)%BLK%VFAC)
    ELSE
      CALL MOVE_ALLOC(BLK%DENSE, CACHE%SLOTS(FREE_S)%BLK%DENSE)
    END IF
    CACHE%SLOTS(FREE_S)%LI = LI; CACHE%SLOTS(FREE_S)%LJ = LJ
    CACHE%SLOTS(FREE_S)%DIRTY = .TRUE.   ! jamais ecrit sur disque -- DESC%OFFSET encore -1
    CACHE%SLOTS(FREE_S)%NBYTES = NBYTES
    CACHE%SLOT_OF(LI,LJ) = FREE_S
    CACHE%USED_BYTES = CACHE%USED_BYTES + NBYTES
    CACHE%CLOCK = CACHE%CLOCK + 1_8
    CACHE%SLOTS(FREE_S)%LAST_USE = CACHE%CLOCK
    S = FREE_S
  END FUNCTION OOC_CACHE_PUT_NEW

  SUBROUTINE OOC_CACHE_MARK_DIRTY(CACHE, S)
    TYPE(OOC_CACHE_T), INTENT(INOUT) :: CACHE
    INTEGER, INTENT(IN)              :: S
    CACHE%SLOTS(S)%DIRTY = .TRUE.
  END SUBROUTINE OOC_CACHE_MARK_DIRTY

  !-----------------------------------------------------------------------
  ! Force l'ecriture+liberation de TOUS les emplacements occupes -- appele
  ! en fin de factorisation (tous les blocs restant residents doivent
  ! atterrir sur disque pour etre relus a la resolution).
  !-----------------------------------------------------------------------
  SUBROUTINE OOC_CACHE_FLUSH_ALL(CACHE, STORE, DESCS)
    TYPE(OOC_CACHE_T), INTENT(INOUT)      :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESCS(:,:)
    INTEGER :: S
    DO S = 1, CACHE%NSLOTS
      IF (CACHE%SLOTS(S)%LI .NE. 0) CALL OOC_CACHE_EVICT_SLOT(CACHE, STORE, DESCS, S)
    END DO
  END SUBROUTINE OOC_CACHE_FLUSH_ALL

  SUBROUTINE OOC_CACHE_FREE(CACHE)
    TYPE(OOC_CACHE_T), INTENT(INOUT) :: CACHE
    IF (ALLOCATED(CACHE%SLOTS)) DEALLOCATE(CACHE%SLOTS)
    IF (ALLOCATED(CACHE%SLOT_OF)) DEALLOCATE(CACHE%SLOT_OF)
    IF (ALLOCATED(CACHE%FREE_STACK)) DEALLOCATE(CACHE%FREE_STACK)
    CACHE%NSLOTS = 0
  END SUBROUTINE OOC_CACHE_FREE

  !===========================================================================
  ! CACHE DEDIE AUX BLOCS DIAGONAUX (DIAGBLK) -- voir les types OOC_DIAG_
  ! CACHE_SLOT_T/OOC_DIAG_CACHE_T declares plus haut (avant CONTAINS).
  !===========================================================================
  SUBROUTINE OOC_DIAG_CACHE_INIT(CACHE, BUDGET_MB, MAX_SLOTS, NMAX)
    TYPE(OOC_DIAG_CACHE_T), INTENT(OUT) :: CACHE
    REAL(WP), INTENT(IN)                :: BUDGET_MB
    INTEGER, INTENT(IN)                 :: MAX_SLOTS, NMAX
    INTEGER :: S
    CACHE%NSLOTS = MAX(MAX_SLOTS, 1)
    ALLOCATE(CACHE%SLOTS(CACHE%NSLOTS))
    ALLOCATE(CACHE%SLOT_OF(MAX(NMAX,1)))
    CACHE%SLOT_OF = 0
    ALLOCATE(CACHE%FREE_STACK(CACHE%NSLOTS))
    DO S = 1, CACHE%NSLOTS
      CACHE%FREE_STACK(S) = CACHE%NSLOTS - S + 1
    END DO
    CACHE%NFREE = CACHE%NSLOTS
    CACHE%BUDGET_BYTES = INT(BUDGET_MB*1024.0_WP*1024.0_WP, 8)
    CACHE%USED_BYTES = 0_8
    CACHE%CLOCK = 0_8
  END SUBROUTINE OOC_DIAG_CACHE_INIT

  SUBROUTINE OOC_DIAG_CACHE_EVICT_SLOT(CACHE, STORE, DESCS, S)
    TYPE(OOC_DIAG_CACHE_T), INTENT(INOUT)   :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)        :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT)   :: DESCS(:)
    INTEGER, INTENT(IN)                     :: S

    IF (CACHE%SLOTS(S)%K .EQ. 0) RETURN
    IF (CACHE%SLOTS(S)%DIRTY) THEN
      CALL OOC_WRITE_BLOCK_FIXED(STORE, DESCS(CACHE%SLOTS(S)%K), CACHE%SLOTS(S)%BLK)
    END IF
    IF (ALLOCATED(CACHE%SLOTS(S)%BLK%DENSE)) DEALLOCATE(CACHE%SLOTS(S)%BLK%DENSE)
    IF (ALLOCATED(CACHE%SLOTS(S)%BLK%UFAC)) DEALLOCATE(CACHE%SLOTS(S)%BLK%UFAC)
    IF (ALLOCATED(CACHE%SLOTS(S)%BLK%VFAC)) DEALLOCATE(CACHE%SLOTS(S)%BLK%VFAC)
    CACHE%USED_BYTES = CACHE%USED_BYTES - CACHE%SLOTS(S)%NBYTES
    CACHE%SLOT_OF(CACHE%SLOTS(S)%K) = 0
    CACHE%SLOTS(S)%K = 0; CACHE%SLOTS(S)%DIRTY = .FALSE.; CACHE%SLOTS(S)%NBYTES = 0_8
    CACHE%NFREE = CACHE%NFREE + 1
    CACHE%FREE_STACK(CACHE%NFREE) = S
  END SUBROUTINE OOC_DIAG_CACHE_EVICT_SLOT

  SUBROUTINE OOC_DIAG_CACHE_MAKE_ROOM(CACHE, STORE, DESCS, NEED_BYTES)
    TYPE(OOC_DIAG_CACHE_T), INTENT(INOUT) :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)      :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT) :: DESCS(:)
    INTEGER(8), INTENT(IN)                :: NEED_BYTES
    INTEGER :: S, LRU_S
    INTEGER(8) :: LRU_VAL

    DO WHILE (CACHE%USED_BYTES + NEED_BYTES .GT. CACHE%BUDGET_BYTES)
      LRU_S = 0; LRU_VAL = HUGE(LRU_VAL)
      DO S = 1, CACHE%NSLOTS
        IF (CACHE%SLOTS(S)%K .NE. 0 .AND. CACHE%SLOTS(S)%LAST_USE .LT. LRU_VAL) THEN
          LRU_VAL = CACHE%SLOTS(S)%LAST_USE; LRU_S = S
        END IF
      END DO
      IF (LRU_S .EQ. 0) RETURN
      CALL OOC_DIAG_CACHE_EVICT_SLOT(CACHE, STORE, DESCS, LRU_S)
    END DO
  END SUBROUTINE OOC_DIAG_CACHE_MAKE_ROOM

  ! Retourne l'indice de l'emplacement contenant le bloc diagonal K,
  ! le chargeant depuis DESCS(K)/STORE s'il n'est pas deja resident.
  FUNCTION OOC_DIAG_CACHE_GET(CACHE, STORE, DESCS, K) RESULT(S)
    TYPE(OOC_DIAG_CACHE_T), INTENT(INOUT)  :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)       :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT)  :: DESCS(:)
    INTEGER, INTENT(IN)                    :: K
    INTEGER :: S, FREE_S
    INTEGER(8) :: NBYTES

    S = CACHE%SLOT_OF(K)
    IF (S .NE. 0) THEN
      CACHE%CLOCK = CACHE%CLOCK + 1_8
      CACHE%SLOTS(S)%LAST_USE = CACHE%CLOCK
      RETURN
    END IF

    NBYTES = 8_8*INT(DESCS(K)%M,8)*INT(DESCS(K)%N,8)
    CALL OOC_DIAG_CACHE_MAKE_ROOM(CACHE, STORE, DESCS, NBYTES)

    IF (CACHE%NFREE .EQ. 0) THEN
      WRITE(6,'(A)') 'ERREUR OOC_DIAG_CACHE_GET : plus d''emplacement libre (augmenter MAX_SLOTS).'
      STOP 1
    END IF
    FREE_S = CACHE%FREE_STACK(CACHE%NFREE)
    CACHE%NFREE = CACHE%NFREE - 1

    CALL OOC_READ_BLOCK_RAW(STORE, DESCS(K), CACHE%SLOTS(FREE_S)%BLK)
    CACHE%SLOTS(FREE_S)%K = K; CACHE%SLOTS(FREE_S)%DIRTY = .FALSE.
    CACHE%SLOTS(FREE_S)%NBYTES = NBYTES
    CACHE%SLOT_OF(K) = FREE_S
    CACHE%USED_BYTES = CACHE%USED_BYTES + NBYTES
    CACHE%CLOCK = CACHE%CLOCK + 1_8
    CACHE%SLOTS(FREE_S)%LAST_USE = CACHE%CLOCK
    S = FREE_S
  END FUNCTION OOC_DIAG_CACHE_GET

  SUBROUTINE OOC_DIAG_CACHE_MARK_DIRTY(CACHE, S)
    TYPE(OOC_DIAG_CACHE_T), INTENT(INOUT) :: CACHE
    INTEGER, INTENT(IN)                   :: S
    CACHE%SLOTS(S)%DIRTY = .TRUE.
  END SUBROUTINE OOC_DIAG_CACHE_MARK_DIRTY

  SUBROUTINE OOC_DIAG_CACHE_FLUSH_ALL(CACHE, STORE, DESCS)
    TYPE(OOC_DIAG_CACHE_T), INTENT(INOUT)  :: CACHE
    TYPE(OOC_STORE_T), INTENT(INOUT)       :: STORE
    TYPE(OOC_BLOCK_DESC_T), INTENT(INOUT)  :: DESCS(:)
    INTEGER :: S
    DO S = 1, CACHE%NSLOTS
      IF (CACHE%SLOTS(S)%K .NE. 0) CALL OOC_DIAG_CACHE_EVICT_SLOT(CACHE, STORE, DESCS, S)
    END DO
  END SUBROUTINE OOC_DIAG_CACHE_FLUSH_ALL

  SUBROUTINE OOC_DIAG_CACHE_FREE(CACHE)
    TYPE(OOC_DIAG_CACHE_T), INTENT(INOUT) :: CACHE
    IF (ALLOCATED(CACHE%SLOTS)) DEALLOCATE(CACHE%SLOTS)
    IF (ALLOCATED(CACHE%SLOT_OF)) DEALLOCATE(CACHE%SLOT_OF)
    IF (ALLOCATED(CACHE%FREE_STACK)) DEALLOCATE(CACHE%FREE_STACK)
    CACHE%NSLOTS = 0
  END SUBROUTINE OOC_DIAG_CACHE_FREE

END MODULE BLR_OOC_MOD
