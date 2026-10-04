!===============================================================================
! Solveur LU MPI+BLR : le stockage local de chaque processus est lui-meme
! compresse, pas seulement la mise a jour de Schur ou la diffusion des
! panneaux.
!
! DIFFERENCE avec les autres solveurs MPI de ce projet : ces modules
! stockent toujours A_LOCAL(MLOC,NLOC) comme un tableau dense (la
! compression n'accelere que la mise a jour de Schur, ou la diffusion des
! panneaux -- jamais le stockage local lui-meme). Ici, le stockage local
! de chaque processus est un tableau de BLR_BLOCK_T (compresse des
! l'assemblage, comme blr_matrix_mod.f90 en sequentiel) -- le passage au
! bloc plein est evite au maximum.
!
! STRUCTURE DE DONNEES : au lieu de l'indexation element-level habituelle
! (NUMROC/INDXL2G sur des lignes/colonnes individuelles), l'indexation se
! fait ici au niveau bloc : chaque processus possede NBLKR_LOC blocs-lignes
! locaux et NBLKC_LOC blocs-colonnes locaux (BLOCK_NUMROC, meme principe
! que NUMROC mais un "element" = un bloc NB x NB entier, cycle = 1 bloc
! par processus). Bloc-ligne local li <-> bloc-ligne global I =
! (li-1)*PR + MYROW + 1. Blocs-colonnes : cyclique DECALE par defaut
! (CMAP=1, BLR_COL_OWNER/BLR_COL_LOCAL/BLR_COL_GLOBAL/BLR_NBLKC_LOC, option
! DIAG_SPREAD des assemblages) pour repartir les blocs diagonaux sur tous
! les processus ; CMAP=0 = cyclique pur, symetrique des lignes -- chaque
! processus d'une meme ligne de grille (meme MYROW) partage exactement le
! meme ensemble de lignes-bloc locales (la formule ne depend que de
! MYROW/PR), ce qui est la cle de l'algorithme de diffusion ci-dessous.
!
! MESSAGES MPI MINIMISES : les panneaux L21/U12 sont diffuses deja
! compresses (BCAST_BLR_BLOCK, mpi_blr_lu_pivot_bcast_mod.f90) -- jamais
! dense-puis-compression, ni meme dense du tout ici. Comme chaque bloc du
! panneau garde sa propre representation (dense ou UFAC/VFAC de rang tres
! inferieur a NB), le volume envoye est proportionnel au rang reel de la
! matrice, pas a NB.
!
! MISE A JOUR DE SCHUR : reste compressee aussi souvent que possible via
! SCHUR_UPDATE_BLR_BLOCK (blr_compress_mod.f90) -- meme mecanisme,
! applique ici bloc local par bloc local.
!
! RESOLUTION : reprend l'algorithme 2D complet de mpi_lu_solve_2d_mod.f90
! (ALLGATHERV+diffusion+SCATTERV sur B, aucune replication residuelle) mais
! remplace le CGEMM dense par BLR_BLOCK_APPLY, bloc local par bloc local --
! exploite directement la compression des facteurs stockes.
!
! LES VARIANTES NOPIV DE CE MODULE RESTENT SANS PIVOT (matrices a diagonale
! dominante uniquement) ; aucune ne s'appuie sur la memoire partagee MPI
! (MPI_Win/RMA) -- l'architecture ne l'empeche pas d'etre ajoutee plus
! tard (chaque processus garde une frontiere claire entre ses donnees
! locales et les messages recus), mais ce n'est pas necessaire ici.
!===============================================================================
MODULE MPI_BLR_MATRIX_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE BLR_OOC_MOD
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE, DENSE_PANEL_T
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_MILESTONE, PRINT_PROGRESS_BAR
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB, BLR_TRIM_HEAP
  USE MKL_INTERFACES_MOD, ONLY: CGETRF, CSWAP, CTRSM
  USE MPI
  IMPLICIT NONE

  ! Hook optionnel pour un suivi fin (ex: RSS/CPU) de l'avancement DANS
  ! l'assemblage/la factorisation/la resolution -- ces boucles sont
  ! internes a ce module, invisibles depuis l'appelant, donc impossibles a
  ! echantillonner a une fraction donnee sans un tel callback. Appele au
  ! franchissement de chaque palier de 5% (voir STAGE_PROGRESS_SHOULD_FIRE
  ! ci-dessous), jamais a chaque iteration -- NTOTAL peut valoir plusieurs
  ! milliers. NDONE/NTOTAL ont le meme sens que pour PRINT_PROGRESS_BAR.
  ABSTRACT INTERFACE
    SUBROUTINE STAGE_PROGRESS_CALLBACK_IFACE(NDONE, NTOTAL)
      INTEGER, INTENT(IN) :: NDONE, NTOTAL
    END SUBROUTINE STAGE_PROGRESS_CALLBACK_IFACE
  END INTERFACE

  TYPE :: MPI_BLR_MATRIX_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (BLR_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)           ! decoupage GLOBAL, partage
    TYPE(DENSE_PANEL_T), ALLOCATABLE :: DIAGBLK(:)       ! indexe par bloc GLOBAL K -- alloue SEULEMENT si proprietaire local
    TYPE(BLR_BLOCK_T), ALLOCATABLE   :: OFFBLK(:,:)      ! indexe par (bloc-ligne LOCAL li, bloc-colonne LOCAL lj)
    INTEGER, ALLOCATABLE :: IPIV(:,:)                    ! (NB,NBLK) -- pivot CONFINE au bloc diagonal K (convention
                                                          ! CLASWP), rempli par MPI_BLR_MATRIX_LU_FACTOR_PIVOT SEULEMENT
                                                          ! (diffuse a TOUS les processus -- petit, cf. commentaire de la
                                                          ! routine), non alloue/non pertinent pour la variante NOPIV.
  END TYPE MPI_BLR_MATRIX_T

  !-----------------------------------------------------------------------
  ! Variante ARENA de MPI_BLR_MATRIX_T (voir "Peux tu mettre en place le
  ! seul vrai correctif, arena sur toute la factorisation, tout en gardant
  ! les executables actuels ?" -- correctif de l'ecart RSS-vs-theorique
  ! documente dans blr_arena_mod.f90/examples/README.md). SEULE difference
  ! structurelle : OFFBLK est un tableau de BLR_BLOCK_ARENA_T (composantes
  ! POINTER, vues carvees dans ARENA) au lieu de BLR_BLOCK_T (composantes
  ! ALLOCATABLE, un ALLOCATE/DEALLOCATE par bloc et par recompression --
  ! source de la fragmentation de tas mesuree). DIAGBLK reste ALLOCATABLE
  ! (dense, jamais recompresse, pas concerne par la fragmentation).
  !
  ! Type ENTIEREMENT NOUVEAU : MPI_BLR_MATRIX_T et toutes les routines qui
  ! l'utilisent (MPI_ASSEMBLE_BLR_MATRIX, MPI_BLR_MATRIX_LU_FACTOR_PIVOT,
  ! MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D, etc.) restent INCHANGEES -- tous les
  ! executables existants continuent d'utiliser MPI_BLR_MATRIX_T tel quel,
  ! permettant une comparaison A/B directe avec les nouvelles routines
  ! *_ARENA ci-dessous.
  !-----------------------------------------------------------------------
  TYPE :: MPI_BLR_MATRIX_ARENA_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (BLR_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(DENSE_PANEL_T), ALLOCATABLE      :: DIAGBLK(:)
    TYPE(BLR_BLOCK_ARENA_T), ALLOCATABLE  :: OFFBLK(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:,:)
    TYPE(BLR_ARENA_T) :: ARENA
  END TYPE MPI_BLR_MATRIX_ARENA_T

  !-----------------------------------------------------------------------
  ! Variante "CAPACITE FIXE" : les blocs bas-rang stockent en permanence
  ! le rang maximal que COMPRESS_WORTHWHILE_RATIO autorise ; les echanges
  ! MPI et les operations associees ne portent, eux, que sur le rang reel.
  ! Voir BLR_BLOCK_FIXEDCAP_T (blr_compress_mod.f90) pour le
  ! principe complet. Contrairement a MPI_BLR_MATRIX_ARENA_T, AUCUN champ
  ! ARENA ici -- chaque bloc gere sa propre capacite RMAX via de simples
  ! composantes ALLOCATABLE, aucun allocateur/pool partage, aucun attribut
  ! TARGET necessaire nulle part dans cette variante.
  !-----------------------------------------------------------------------
  TYPE :: MPI_BLR_MATRIX_FIXEDCAP_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (BLR_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(DENSE_PANEL_T), ALLOCATABLE         :: DIAGBLK(:)
    TYPE(BLR_BLOCK_FIXEDCAP_T), ALLOCATABLE  :: OFFBLK(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:,:)
  END TYPE MPI_BLR_MATRIX_FIXEDCAP_T

  !-----------------------------------------------------------------------
  ! Variante OUT-OF-CORE : ne plus etre limite par la RAM disponible, ni
  ! pour la matrice assemblee, ni pour les facteurs LU apres remplissage
  ! de Schur, ni pour un second membre massif -- basee sur l'algorithme
  ! classique (pas l'arena : le cycle allocate/deallocate individuel de
  ! BLR_BLOCK_T se prete naturellement a l'eviction disque, contrairement
  ! a l'arena qui ne libere jamais un bloc individuellement).
  !
  ! DIAGBLK/OFFBLK ne contiennent plus les donnees directement (UFAC/VFAC/
  ! DENSE) mais un DESCRIPTEUR LEGER (OOC_BLOCK_DESC_T, blr_ooc_mod.f90 --
  ! M/N/RANG/IS_LOWRANK + position sur le fichier scratch) : toujours
  ! resident, negligeable en RAM quel que soit N. Les donnees REELLES
  ! vivent sur STORE_A (fichier scratch local a CE processus, jamais
  ! partage -- chaque processus ne possede deja que ses propres blocs
  ! locaux dans la distribution 2D-cyclique existante). CACHE_A borne la
  ! RAM consommee par les blocs EFFECTIVEMENT residents a un instant donne,
  ! bornee par MAX_MB_PER_PROCESS (parametre utilisateur) -- CACHE_A pour
  ! les blocs hors-diagonale (OFFBLK), CACHE_D pour les blocs diagonaux
  ! (DIAGBLK). Les deux n'ecrivent sur disque qu'a l'eviction, et TOUJOURS
  ! EN PLACE (OOC_WRITE_BLOCK_FIXED, voir blr_ooc_mod.f90) : premier appel
  ! reserve le PIRE CAS DENSE (M*N, jamais la taille bas-rang courante),
  ! garantissant qu'un rang qui grossit plus tard (mise a jour de Schur :
  ! le rang n'est PAS garanti decroissant) tient toujours dans la place
  ! reservee -- aucun octet orphelin, ni pour DIAGBLK (toujours dense, M*N
  ! est deja sa taille reelle) ni pour OFFBLK (bas-rang, mais toujours
  ! strictement plus petit que M*N par construction de COMPRESS_WORTHWHILE_
  ! RATIO, voir blr_compress_mod.f90).
  !
  ! DIMENSIONNER MAX_MB_PER_PROCESS (mesure sur N=20000, NB=400, grille
  ! 2x2, voir bench/bench_blr_mpi_pivot_ooc_compare.f90) : un budget
  ! INSUFFISANT face a l'ensemble de travail REEL du processus ne degrade
  ! pas gracieusement -- il fait THRASHER le cache LRU, avec un cout
  ! disproportionne par rapport a l'ecart de RAM economise. Cause : un bloc
  ! hors-diagonale n'est retouche par une mise a jour de Schur qu'une fois
  ! tous les ~NBLK pas de la factorisation, et la TOTALITE de l'ensemble de
  ! travail local est parcourue entre deux touches du MEME bloc -- un
  ! schema d'acces quasi-SEQUENTIEL ou AUCUNE taille de cache inferieure a
  ! l'ensemble de travail complet n'aide, quelle que soit la politique
  ! d'eviction (LRU ou autre). Mesure concrete : a budget=200 Mo (RSS
  ! reelle observee ~127 Mo, donc en apparence "assez"), le temps de
  ! factorisation etait domine a 68% par les acces cache (10.1s sur 14.9s,
  ! CGETRF/TRSM/Schur eux-memes ne prenant que 0.97s, deja comparable au
  ! classique) ; passer a budget=800 Mo (RSS reelle quasi inchangee, ~128
  ! Mo -- l'ensemble de travail REEL etait deja atteint, PAS le budget)
  ! a fait chuter ce cout cache de 10.1s a 0.75s (/13.5), et le temps de
  ! factorisation total de 14.9s a 3.1s (/4.9), pour une RAM quasi
  ! identique. CONCLUSION PRATIQUE : ne pas choisir MAX_MB_PER_PROCESS "le
  ! plus petit possible" pour economiser de la RAM -- le choisir au moins
  ! egal a l'ensemble de travail reel du processus (verifiable via
  ! MPI_BLR_MATRIX_OOC_RESIDENT_BYTES a un budget genereux, ou en
  ! augmentant le budget jusqu'a ce que la RSS reelle cesse de croitre) ;
  ! l'interet de l'OOC n'est PAS de tourner avec le budget minimal, mais de
  ! rester CORRECT et de degrader la RAM au lieu de planter quand le
  ! probleme est trop gros pour tenir en RAM du tout -- en-dessous de ce
  ! seuil "ensemble de travail", le compromis temps/memoire devient tres
  ! defavorable et croit sans borne a mesure que N grandit.
  !-----------------------------------------------------------------------
  TYPE :: MPI_BLR_MATRIX_OOC_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER :: CMAP = 1     ! distribution des blocs-colonnes (BLR_COL_OWNER) : 0 cyclique, 1 decalee
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(OOC_BLOCK_DESC_T), ALLOCATABLE :: DIAGBLK(:)
    TYPE(OOC_BLOCK_DESC_T), ALLOCATABLE :: OFFBLK(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:,:)
    TYPE(OOC_STORE_T) :: STORE_A
    TYPE(OOC_CACHE_T) :: CACHE_A
    TYPE(OOC_DIAG_CACHE_T) :: CACHE_D
    REAL(WP) :: MAX_MB_PER_PROCESS = 0.0_WP
  END TYPE MPI_BLR_MATRIX_OOC_T

  !-----------------------------------------------------------------------
  ! Second membre OUT-OF-CORE : assemble et stocke SEPAREMENT de la
  ! matrice (fichier scratch propre, prefixe "b" -- voir
  ! MPI_ASSEMBLE_RHS_OOC), par PAQUETS d'au plus NRHS_BATCH colonnes
  ! GLOBALES -- pour borner la memoire meme dans un cas massivement
  ! second membre. Chaque paquet est un bloc dense ordinaire
  ! (BLR_BLOCK_T, IS_LOWRANK=.FALSE.) au sens de blr_ooc_mod.f90 -- aucun
  ! nouveau mecanisme d'(e)serialisation necessaire, juste reutilise. Pas
  ! de cache ici (contrairement a OFFBLK) : chaque paquet n'est relu
  ! qu'UNE FOIS par MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D, un acces direct
  ! disque suffit.
  !-----------------------------------------------------------------------
  TYPE :: OOC_RHS_T
    TYPE(OOC_STORE_T) :: STORE_B
    TYPE(OOC_BLOCK_DESC_T), ALLOCATABLE :: BATCH(:)
    INTEGER :: NRHS = 0, NRHS_BATCH = 0, NBATCHES = 0, MLOC = 0
  END TYPE OOC_RHS_T

  ABSTRACT INTERFACE
    SUBROUTINE FILL_RHS_IFACE(ISTART, IEND, RSTART, REND, BLK)
      IMPORT :: WP
      INTEGER, INTENT(IN)       :: ISTART, IEND, RSTART, REND
      COMPLEX(WP), INTENT(OUT)  :: BLK(IEND-ISTART+1, REND-RSTART+1)
    END SUBROUTINE FILL_RHS_IFACE

    ! Rappelee UNE FOIS PAR PAQUET par MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D,
    ! avec la solution LOCALE (lignes locales, colonnes GSTART:GEND
    ! GLOBALES) de CE paquet -- a l'appelant de decider quoi en faire
    ! (rassembler dans un tableau complet si NRHS modeste comme dans les
    ! exemples/tests, ou ecrire directement sur disque pour un cas
    ! massivement second membre reel). Aucune accumulation NRHS-large
    ! n'est faite par la bibliotheque elle-meme -- c'est precisement le
    ! but de cette interface : borner la memoire de la resolution
    ! independamment de NRHS.
    ! BATCH_RSTART/BATCH_REND : plage de colonnes GLOBALE du paquet ENTIER
    ! (identique sur tous les processus -- necessaire si l'appelant doit
    ! rassembler entre processus, ex. via GATHER_2D_RHS avec
    ! NRHS:=BATCH_REND-BATCH_RSTART+1). GSTART/GEND : sous-plage GLOBALE
    ! possedee par CE processus (coherent avec X_LOCAL/NRHS_LOC_B).
    SUBROUTINE CONSUME_X_BATCH_IFACE(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
      IMPORT :: WP
      INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
      COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    END SUBROUTINE CONSUME_X_BATCH_IFACE
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------
  ! .TRUE. exactement quand NDONE fait franchir un nouveau palier de 5%
  ! depuis NDONE-1 (comparaison par division entiere, meme technique
  ! qu'un redessin de barre de progression), ou a la toute derniere
  ! iteration -- garantit un dernier appel a 100% meme si NTOTAL n'est
  ! pas un multiple de 20.
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION STAGE_PROGRESS_SHOULD_FIRE(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    INTEGER :: NT
    NT = MAX(NTOTAL, 1)
    STAGE_PROGRESS_SHOULD_FIRE = (NDONE .EQ. NTOTAL) .OR. &
         ((NDONE*20)/NT .NE. (MAX(NDONE-1,0)*20)/NT)
  END FUNCTION STAGE_PROGRESS_SHOULD_FIRE

  !-----------------------------------------------------------------------
  ! Repartition round-robin de NBLK blocs sur P proprietaires (1 bloc par
  ! cycle) -- meme principe que NUMROC mais a la granularite BLOC, pas
  ! ELEMENT (donc pas de dependance a NB ici).
  !-----------------------------------------------------------------------
  FUNCTION BLOCK_NUMROC(NBLK, P, MYP) RESULT(NLOC)
    INTEGER, INTENT(IN) :: NBLK, P, MYP
    INTEGER :: NLOC
    IF (MYP .LT. MOD(NBLK,P)) THEN
      NLOC = NBLK/P + 1
    ELSE
      NLOC = NBLK/P
    END IF
  END FUNCTION BLOCK_NUMROC

  !-----------------------------------------------------------------------
  ! Distribution des blocs-COLONNES sur la grille (LU et LDL^T). Les
  ! blocs-lignes restent cycliques : I -> ligne de grille MOD(I-1,PR),
  ! local li=(I-1)/PR+1.
  !   CMAP = 0 : cyclique pur, J -> MOD(J-1,PC). Le bloc diagonal K tombe
  !     sur (MOD(K-1,PR), MOD(K-1,PC)) : si g = PGCD(PR,PC) > 1, seuls
  !     PR*PC/g processus possedent des blocs diagonaux (2 sur 4 en 2x2,
  !     4 sur 8 en 2x4) et concentrent leur factorisation, leur Schur dense
  !     et leur memoire.
  !   CMAP = 1 (defaut, champ CMAP des types, option DIAG_SPREAD des
  !     assemblages) : cyclique DECALE d'une colonne a chaque periode de
  !     L = PPCM(PR,PC) blocs, J -> MOD(J-1 + (J-1)/L, PC). Sur g periodes
  !     (g*L = PR*PC blocs), les blocs diagonaux visitent chacun des PR*PC
  !     processus exactement une fois ; chaque colonne de grille recoit
  !     toujours L/PC blocs par periode. Si g = 1 (PR et PC premiers entre
  !     eux, dont PR=1 ou PC=1), la distribution cyclique couvre deja toute
  !     la grille : CMAP=1 y est identique a CMAP=0.
  ! Indice local lj dans la colonne proprietaire : rang de J parmi les
  ! blocs-colonnes de cette colonne (croissant avec J dans les deux cas).
  !-----------------------------------------------------------------------
  INTEGER FUNCTION BLR_CMAP_PERIOD(CMAP, GRID)
    INTEGER, INTENT(IN)      :: CMAP
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER :: A1, B1, T
    BLR_CMAP_PERIOD = 0                       ! 0 : pas de decalage
    IF (CMAP .EQ. 0) RETURN
    A1 = GRID%PR; B1 = GRID%PC
    DO WHILE (B1 .NE. 0)
      T = MOD(A1, B1); A1 = B1; B1 = T
    END DO
    IF (A1 .GT. 1) BLR_CMAP_PERIOD = (GRID%PR/A1)*GRID%PC   ! PPCM
  END FUNCTION BLR_CMAP_PERIOD

  ! Colonne de grille proprietaire du bloc-colonne global J.
  INTEGER FUNCTION BLR_COL_OWNER(CMAP, J, GRID)
    INTEGER, INTENT(IN)      :: CMAP, J
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER :: L
    L = BLR_CMAP_PERIOD(CMAP, GRID)
    IF (L .EQ. 0) THEN
      BLR_COL_OWNER = MOD(J-1, GRID%PC)
    ELSE
      BLR_COL_OWNER = MOD(J-1 + (J-1)/L, GRID%PC)
    END IF
  END FUNCTION BLR_COL_OWNER

  ! Indice local du bloc-colonne global J dans sa colonne proprietaire.
  INTEGER FUNCTION BLR_COL_LOCAL(CMAP, J, GRID)
    INTEGER, INTENT(IN)      :: CMAP, J
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER :: L
    L = BLR_CMAP_PERIOD(CMAP, GRID)
    IF (L .EQ. 0) THEN
      BLR_COL_LOCAL = (J-1)/GRID%PC + 1
    ELSE
      BLR_COL_LOCAL = ((J-1)/L)*(L/GRID%PC) + MOD(J-1, L)/GRID%PC + 1
    END IF
  END FUNCTION BLR_COL_LOCAL

  ! Bloc-colonne global du bloc-colonne local LJ de MA colonne de grille.
  INTEGER FUNCTION BLR_COL_GLOBAL(CMAP, LJ, GRID)
    INTEGER, INTENT(IN)      :: CMAP, LJ
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER :: L, P, Q
    L = BLR_CMAP_PERIOD(CMAP, GRID)
    IF (L .EQ. 0) THEN
      BLR_COL_GLOBAL = (LJ-1)*GRID%PC + GRID%MYCOL + 1
    ELSE
      P = (LJ-1)/(L/GRID%PC)
      Q = MOD(LJ-1, L/GRID%PC)
      BLR_COL_GLOBAL = P*L + Q*GRID%PC + MODULO(GRID%MYCOL - P, GRID%PC) + 1
    END IF
  END FUNCTION BLR_COL_GLOBAL

  ! Nombre de blocs-colonnes de MA colonne de grille.
  INTEGER FUNCTION BLR_NBLKC_LOC(CMAP, NBLK, GRID)
    INTEGER, INTENT(IN)      :: CMAP, NBLK
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER :: J
    BLR_NBLKC_LOC = 0
    DO J = 1, NBLK
      IF (BLR_COL_OWNER(CMAP, J, GRID) .EQ. GRID%MYCOL) BLR_NBLKC_LOC = BLR_NBLKC_LOC + 1
    END DO
  END FUNCTION BLR_NBLKC_LOC

  !-----------------------------------------------------------------------
  ! Assemblage bloc par bloc, avec compression immediate des blocs hors-
  ! diagonale locaux (meme principe que ASSEMBLE_BLR_MATRIX/blr_matrix_
  ! mod.f90, mais chaque processus n'assemble que sa part -- GET_BLOCK
  ! recoit des indices globaux, comme en sequentiel).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, &
                                      T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                      SHOW_PROGRESS, T_GET_BLOCK, T_COMPRESS, T_GET_BLOCK_MIN, T_COMPRESS_MIN, &
                                      USE_RSVD, RSVD_OVERSAMPLE, RSVD_NPOWER, PROGRESS_CALLBACK, DIAG_SPREAD)
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir BLR_COL_OWNER) -- les blocs diagonaux visitent tous les
    ! processus ; .FALSE. : cyclique pur (ancien comportement).
    LOGICAL, INTENT(IN), OPTIONAL :: DIAG_SPREAD
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_T), INTENT(OUT) :: A_BLR
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS
    ! USE_RSVD (defaut .TRUE. depuis 2026-10-03 -- passer .FALSE. pour
    ! revenir a la SVD complete CGESVD, rang optimal garanti) :
    ! remplace la SVD complete (COMPRESS_BLOCK, CGESVD) par la SVD
    ! randomisee ADAPTATIVE (COMPRESS_BLOCK_RSVD_ADAPTIVE, blr_compress_
    ! mod.f90 -- rang decouvert par doublement, TOL seul, sans connaissance
    ! prealable) pour compresser chaque bloc hors-diagonale. MESURE (voir
    ! test/bench_rsvd_compress.f90) : x34-54 plus rapide que CGESVD sur un
    ! bloc bien compressible, erreur de reconstruction du meme ordre
    ! (ratio 1.00) une fois le bogue d'orthonormalisation corrige -- mais
    ! GARANTIE HEURISTIQUE, pas une preuve formelle (voir la documentation
    ! de COMPRESS_BLOCK_RSVD_ADAPTIVE), et une ZONE MORTE CARACTERISEE pres
    ! du plafond COMPRESS_WORTHWHILE_RATIO ou elle replie en dense la ou
    ! COMPRESS_BLOCK compresserait (comportement accepte, voir son
    ! commentaire). RSVD_OVERSAMPLE/RSVD_NPOWER (optionnels, defauts 10/1
    ! -- memes valeurs par defaut que COMPRESS_BLOCK_RSVD_ADAPTIVE) :
    ! transmis tels quels, voir son commentaire pour leur role.
    LOGICAL, INTENT(IN), OPTIONAL   :: USE_RSVD
    INTEGER, INTENT(IN), OPTIONAL   :: RSVD_OVERSAMPLE, RSVD_NPOWER
    ! Appele par CHAQUE processus au franchissement d'un palier de 5% de
    ! SA PROPRE progression LOCALE (NDONE/NTOTAL=NBLKR_LOC*NBLKC_LOC) --
    ! l'assemblage n'etant pas synchronise entre processus (voir plus
    ! bas), les instants de franchissement different d'un processus a
    ! l'autre, ce qui est attendu (voir example_blr_mpi_pivot.f90,
    ! PROFILE_USAGE). Absent par defaut, aucun surcout si non fourni.
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK
    ! DIAGNOSTIC (utile pour evaluer si une SVD randomisee vaut le coup a
    ! la place de CGESVD dans COMPRESS_BLOCK) : T_GET_BLOCK cumule le
    ! temps dans tous les appels GET_BLOCK (evaluation du noyau, diagonal +
    ! hors-diagonale) ; T_COMPRESS cumule le temps dans tous les appels
    ! COMPRESS_BLOCK (SVD complete, hors-diagonale uniquement). Optionnels,
    ! aucun surcout si absents.
    !
    ! T_ASSEMBLE_COMPRESS/T_GET_BLOCK/T_COMPRESS sont un vrai MAX sur
    ! GRID%COMM, jamais la valeur locale du rang appelant -- important ici
    ! car cette boucle d'assemblage est "embarrassingly parallel" (aucune
    ! synchronisation entre processus, voir commentaire plus bas) : sans
    ! cette reduction, un rang plus lent ne ferait pas attendre les autres
    ! avant que leur propre chronometre local ne s'arrete, et le temps
    ! qu'il lui reste a passer fuirait silencieusement dans la phase
    ! suivante (premiere communication collective de la factorisation, qui
    ! l'attend). Nouveaux _MIN (optionnels) : ecart MAX/MIN direct =
    ! indicateur de desequilibre de charge entre processus (une
    ! distribution bloc-cyclique devrait rester equilibree, un ecart
    ! important signalerait un noyau GET_BLOCK dont le cout varie fortement
    ! selon les blocs, ou un decoupage NBLK non multiple de PR*PC).
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_GET_BLOCK, T_COMPRESS, T_GET_BLOCK_MIN, T_COMPRESS_MIN

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, PBAR_STATE
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0, T_GETBLK_L, T_COMPRESS_L, T_ASM_L
    LOGICAL :: PROGRESS, RSVD_FLAG
    INTEGER :: RSVD_OVS, RSVD_NPOW

    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD
    RSVD_OVS = 10; IF (PRESENT(RSVD_OVERSAMPLE)) RSVD_OVS = RSVD_OVERSAMPLE
    RSVD_NPOW = 1; IF (PRESENT(RSVD_NPOWER)) RSVD_NPOW = RSVD_NPOWER

    ! La barre de progression est affichee UNIQUEMENT par le processus 0,
    ! sur SA PROPRE part locale (NBLKR_LOC*NBLKC_LOC) -- pas de
    ! synchronisation entre processus dans cette boucle (assemblage
    ! "embarrassingly parallel", chaque processus assemble ses blocs sans
    ! communication), donc pas de pourcentage GLOBAL exact possible sans
    ! communication supplementaire. Approximation raisonnable : le
    ! decoupage bloc-cyclique (BLOCK_NUMROC) equilibre la charge entre
    ! processus, la progression du processus 0 suit donc de pres celle
    ! des autres.
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NDONE = 0
    PBAR_STATE = -1
    T_GETBLK_L = 0.0D0; T_COMPRESS_L = 0.0D0

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A_BLR%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A_BLR%NBLKC_LOC = BLR_NBLKC_LOC(A_BLR%CMAP, A_BLR%NBLK, GRID)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(li,lj), RSVD_OVS, RSVD_NPOW)
          ELSE
            CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(li,lj))
          END IF
          T_COMPRESS_L = T_COMPRESS_L + (MPI_WTIME() - TT0)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(MI,NJ)
          END IF
          DEALLOCATE(TMP)
        END IF
        NDONE = NDONE + 1
        ! Barre qui se remplit progressivement (voir progress_bar_mod.f90)
        ! -- robuste sous MPI puisqu'elle n'emet aucun caractere de
        ! controle. Seul le rang 0 ecrit ici (PROGRESS deja gate sur
        ! GRID%MYID.EQ.0).
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('1) Assemblage+compression (rang 0)', NDONE, NTOTAL, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(NDONE, NTOTAL)) CALL PROGRESS_CALLBACK(NDONE, NTOTAL)
        END IF
      END DO
    END DO
    T_ASM_L = MPI_WTIME() - T0
    CALL MPI_ALLREDUCE(T_ASM_L, T_ASSEMBLE_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_GET_BLOCK)) CALL MPI_ALLREDUCE(T_GETBLK_L,   T_GET_BLOCK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_COMPRESS))  CALL MPI_ALLREDUCE(T_COMPRESS_L, T_COMPRESS,  1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_GET_BLOCK_MIN)) CALL MPI_ALLREDUCE(T_GETBLK_L,   T_GET_BLOCK_MIN, 1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (PRESENT(T_COMPRESS_MIN))  CALL MPI_ALLREDUCE(T_COMPRESS_L, T_COMPRESS_MIN,  1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX

  !-----------------------------------------------------------------------
  ! Empreinte memoire LOCALE (octets) de la part de A_BLR possedee par CE
  ! processus -- meme esprit que BLR_MATRIX_MEMORY_BYTES ("A" sequentiel),
  ! l'appelant reduit (MAX ou SUM) selon ce qu'il veut mesurer.
  !-----------------------------------------------------------------------
  FUNCTION MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR) RESULT(NBYTES)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER(8) :: NBYTES
    INTEGER :: li, lj, K
    NBYTES = 0_8
    DO K = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(K)%D)) NBYTES = NBYTES + 8_8*INT(SIZE(A_BLR%DIAGBLK(K)%D),8)
    END DO
    DO lj = 1, A_BLR%NBLKC_LOC
      DO li = 1, A_BLR%NBLKR_LOC
        NBYTES = NBYTES + BLR_BLOCK_BYTES(A_BLR%OFFBLK(li,lj))
      END DO
    END DO
  END FUNCTION MPI_BLR_MATRIX_MEMORY_BYTES

  !-----------------------------------------------------------------------
  ! Libere la RAM locale de A_BLR (DIAGBLK/OFFBLK -- l'essentiel de la
  ! memoire, cf. MPI_BLR_MATRIX_MEMORY_BYTES ci-dessus -- puis BSTART/BEND/
  ! IPIV, petits). A_BLR n'est pas ALLOCATABLE (ses composantes le sont) :
  ! desalouer chaque tableau derive (DIAGBLK, OFFBLK) desaloue AUSSI leurs
  ! composantes allouables (UFAC/VFAC/DENSE par element), automatiquement,
  ! norme Fortran. A appeler quand A_BLR n'est plus necessaire (ex. juste
  ! apres la resolution) pour reduire le pic RSS avant un post-traitement.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_FREE(A_BLR)
    TYPE(MPI_BLR_MATRIX_T), INTENT(INOUT) :: A_BLR
    IF (ALLOCATED(A_BLR%DIAGBLK)) DEALLOCATE(A_BLR%DIAGBLK)
    IF (ALLOCATED(A_BLR%OFFBLK))  DEALLOCATE(A_BLR%OFFBLK)
    IF (ALLOCATED(A_BLR%BSTART))  DEALLOCATE(A_BLR%BSTART)
    IF (ALLOCATED(A_BLR%BEND))    DEALLOCATE(A_BLR%BEND)
    IF (ALLOCATED(A_BLR%IPIV))    DEALLOCATE(A_BLR%IPIV)
    A_BLR%N = 0; A_BLR%NB = 0; A_BLR%NBLK = 0; A_BLR%NBLKR_LOC = 0; A_BLR%NBLKC_LOC = 0
  END SUBROUTINE MPI_BLR_MATRIX_FREE

  !-----------------------------------------------------------------------
  ! Factorisation LU en place, sans pivot, sur A_BLR -- reste compressee
  ! tout du long (voir commentaire de tete de fichier pour l'algorithme
  ! complet). N_COMPRESSED/N_OFFDIAG/RANK_SUM caracterisent les facteurs
  ! L/U finaux (sommes sur tous les processus).
  !-----------------------------------------------------------------------
  ! T_DIAG/T_BCAST_DKK/T_TRSM/T_BCAST_PANEL/T_ALLOC/T_MPICALL/T_SCHUR_PRODUCT/
  ! T_RECOMPRESS/T_SCHUR_DIAG (tous optionnels, INTENT(OUT), max/processus) :
  ! decoupage fin du temps de factorisation, categorie par categorie --
  ! T_BCAST_PANEL est le temps total de l'etape 5 (diffusion des panneaux
  ! L/U compresses) et INCLUT T_ALLOC/T_MPICALL comme sous-parties (le
  ! reste, non isole separement, est le pack/unpack RESHAPE dans
  ! BCAST_BLR_BLOCK) ; T_SCHUR_PRODUCT = temps de calcul du produit L*U
  ! PLUS le repli dense eventuel ("Schur SANS la recompression") ;
  ! T_RECOMPRESS = FAST_RECOMPRESS seul ; T_SCHUR_DIAG = mise a jour des
  ! blocs diagonaux du Schur (SCHUR_UPDATE_BLOCK, toujours dense, sous-cas
  ! separe de SCHUR_UPDATE_BLR_BLOCK). Aucun surcout si non fournis.
  SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_NOPIV(GRID, A_BLR, TOL, INFO, &
                                             N_COMPRESSED, N_OFFDIAG, RANK_SUM, T_FACTOR, &
                                             T_DIAG, T_BCAST_DKK, T_TRSM, T_BCAST_PANEL, &
                                             T_ALLOC, T_MPICALL, T_SCHUR_PRODUCT, T_RECOMPRESS, T_SCHUR_DIAG, &
                                             SHOW_PROGRESS)
    TYPE(GRID_T), INTENT(IN)          :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(INOUT) :: A_BLR
    REAL(WP), INTENT(IN)              :: TOL
    INTEGER, INTENT(OUT)              :: INFO
    INTEGER, INTENT(OUT)              :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)           :: RANK_SUM
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_DIAG, T_BCAST_DKK, T_TRSM, T_BCAST_PANEL
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_ALLOC, T_MPICALL, T_SCHUR_PRODUCT, T_RECOMPRESS, T_SCHUR_DIAG
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0
    DOUBLE PRECISION :: T_DIAG_L, T_BCASTDKK_L, T_TRSM_L, T_BCASTPANEL_L
    DOUBLE PRECISION :: T_ALLOC_L, T_MPICALL_L, T_PRODUCT_L, T_RECOMP_L, T_SCHURDIAG_L
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    LAST_PCT = -5

    T0 = MPI_WTIME()
    T_DIAG_L = 0.0D0; T_BCASTDKK_L = 0.0D0; T_TRSM_L = 0.0D0; T_BCASTPANEL_L = 0.0D0
    T_ALLOC_L = 0.0D0; T_MPICALL_L = 0.0D0; T_PRODUCT_L = 0.0D0; T_RECOMP_L = 0.0D0; T_SCHURDIAG_L = 0.0D0
    NBLK = A_BLR%NBLK
    INFO = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8
    ALLOCATE(LRECV(MAX(A_BLR%NBLKR_LOC,1)), URECV(MAX(A_BLR%NBLKC_LOC,1)))
    CALL BLR_WORKSPACE_ALLOC(WS, A_BLR%NB)

    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)
      INFO_STEP = 0

      ! 1) facteur diagonal, local au proprietaire
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        TT0 = MPI_WTIME()
        CALL BLR_LU_NOPIV_BLOCKED(KW, A_BLR%DIAGBLK(K)%D, KW, IINFO)
        T_DIAG_L = T_DIAG_L + (MPI_WTIME() - TT0)
        IF (IINFO .NE. 0) INFO_STEP = A_BLR%BSTART(K) - 1 + IINFO
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(LRECV, URECV)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF

      ! 2) DKK diffuse le long de la ligne PROW_K (pour U12) ET de la
      !    colonne PCOL_K (pour L21) -- toujours dense, KW<=NB, jamais
      !    interessant a compresser.
      ALLOCATE(DKK(KW,KW))
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
      TT0 = MPI_WTIME()
      IF (GRID%MYROW .EQ. PROW_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      IF (GRID%MYCOL .EQ. PCOL_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      T_BCASTDKK_L = T_BCASTDKK_L + (MPI_WTIME() - TT0)

      ! 3) U12 : la ligne PROW_K applique le TRSM sur SES blocs locaux
      !    deja compresses -- aucune communication supplementaire.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        TT0 = MPI_WTIME()
        li = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BLR_TRSM_LEFT_LOWER_UNIT(DKK, KW, A_BLR%OFFBLK(li,lj))
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
        T_TRSM_L = T_TRSM_L + (MPI_WTIME() - TT0)
      END IF

      ! 4) L21 : la colonne PCOL_K applique le TRSM sur SES blocs locaux
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        TT0 = MPI_WTIME()
        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BLR_TRSM_RIGHT_UPPER_NONUNIT(A_BLR%OFFBLK(li,lj), DKK, KW)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
        T_TRSM_L = T_TRSM_L + (MPI_WTIME() - TT0)
      END IF
      DEALLOCATE(DKK)

      ! 5) diffusion des panneaux L/U deja compresses (reduit au maximum
      !    le volume envoye, voir BCAST_BLR_BLOCK) : le
      !    panneau L(:,K) local a la colonne PCOL_K est diffuse le long de
      !    ROW_COMM (touche toute la ligne MYROW -- MEME ensemble de
      !    lignes-bloc locales pour tous les processus de cette ligne).
      !    Symetriquement pour U(K,:) le long de COL_COMM.
      TT0 = MPI_WTIME()
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, &
                                A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)), T_ALLOC_L, T_MPICALL_L)
          LRECV(li) = A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID))
        END DO
      ELSE
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li), &
                                T_ALLOC_L, T_MPICALL_L)
        END DO
      END IF

      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, &
                                A_BLR%OFFBLK((K-1)/GRID%PR+1,lj), T_ALLOC_L, T_MPICALL_L)
          URECV(lj) = A_BLR%OFFBLK((K-1)/GRID%PR+1,lj)
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj), &
                                T_ALLOC_L, T_MPICALL_L)
        END DO
      END IF
      T_BCASTPANEL_L = T_BCASTPANEL_L + (MPI_WTIME() - TT0)

      ! 6) mise a jour de Schur locale, bloc par bloc, en gardant chaque
      !    bloc hors-diagonale compresse aussi souvent que possible.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          IF (I .EQ. J) THEN
            TT0 = MPI_WTIME()
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LRECV(li), URECV(lj), A_BLR%DIAGBLK(I)%D, MI)
            T_SCHURDIAG_L = T_SCHURDIAG_L + (MPI_WTIME() - TT0)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK(A_BLR%OFFBLK(li,lj), MI, NJ, KW, LRECV(li), URECV(lj), TOL, WS, &
                                         T_PRODUCT_L, T_RECOMP_L)
          END IF
        END DO
      END DO

      CALL BLR_TRIM_HEAP()        ! fin d'etape : memoire liberee rendue au systeme
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('Factorisation LU (sans pivot)', K, NBLK, LAST_PCT)
    END DO

    DEALLOCATE(LRECV, URECV)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) T_FACTOR = MPI_WTIME() - T0

    IF (PRESENT(T_DIAG))          CALL MPI_ALLREDUCE(T_DIAG_L,      T_DIAG,          1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_BCAST_DKK))     CALL MPI_ALLREDUCE(T_BCASTDKK_L,  T_BCAST_DKK,     1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_TRSM))          CALL MPI_ALLREDUCE(T_TRSM_L,      T_TRSM,          1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_BCAST_PANEL))   CALL MPI_ALLREDUCE(T_BCASTPANEL_L,T_BCAST_PANEL,   1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_ALLOC))         CALL MPI_ALLREDUCE(T_ALLOC_L,     T_ALLOC,         1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_MPICALL))       CALL MPI_ALLREDUCE(T_MPICALL_L,   T_MPICALL,       1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_SCHUR_PRODUCT)) CALL MPI_ALLREDUCE(T_PRODUCT_L,   T_SCHUR_PRODUCT, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_RECOMPRESS))    CALL MPI_ALLREDUCE(T_RECOMP_L,    T_RECOMPRESS,    1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_SCHUR_DIAG))    CALL MPI_ALLREDUCE(T_SCHURDIAG_L, T_SCHUR_DIAG,    1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_NOPIV

  !-----------------------------------------------------------------------
  ! Factorisation LU en place, avec pivot partiel confine au bloc diagonal
  ! (meme philosophie que block_lu_mod.f90 et mpi_lu_pivot_mod.f90 -- pivot
  ! recherche uniquement dans les KW lignes du bloc diagonal courant,
  ! jamais globalement). Choisi plutot qu'une methode "ne decompresse
  ! jamais" (plus risquee) car permuter des lignes d'un bloc bas-rang ne
  ! touche que UFAC (cout O(rang), jamais de decompression) -- voir
  ! BLR_BLOCK_ROW_SWAP.
  !
  ! DIAGONALE : CGETRF directement (plus de conflit avec un invariant
  ! "sans pivot" puisqu'on l'assouplit ICI expres) -- recupere gratuitement
  ! le blocage/threading MKL deja explore pour BLR_LU_NOPIV_BLOCKED, ET les
  ! indices de pivot en un seul appel.
  !
  ! PROPAGATION : la sequence de pivots (KW entiers, convention CLASWP --
  ! swaps ORDONNES, PAS un vecteur de permutation direct) est diffusee a
  ! TOUT le communicateur GLOBAL (pas seulement ROW_COMM) et stockee dans
  ! A_BLR%IPIV(:,K) -- volontairement un peu plus large que necessaire a la
  ! factorisation (ROW_COMM suffirait ici) pour que la resolution puisse la
  ! relire directement SANS communication supplementaire (elle en a besoin
  ! sur TOUTE la grille, cf. MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D) -- message
  ! minuscule (KW entiers), le surcout est negligeable.
  !
  ! Applique ensuite a TOUS les blocs de la MEME ligne-bloc K deja stockes
  ! localement (OFFBLK(li_K,lj) pour tout J.NE.K, aussi bien J<K -- facteurs
  ! L deja finalises aux etapes precedentes, meme raisonnement que CLASWP
  ! appliquee aux colonnes deja traitees dans LAPACK -- que J>K, blocs bruts
  ! pas encore traites). RIEN a faire sur la colonne K (blocs I>K) : leurs
  ! lignes ne sont pas celles du bloc diagonal K, non concernees.
  !
  ! N_SWAPS (diagnostic, optionnel) : nombre TOTAL de swaps reellement
  ! effectues (somme sur tous les K) -- 0 signifie qu'aucun pivot n'a
  ! jamais ete necessaire (matrice testee pas assez agressive), utile pour
  ! verifier qu'un test de validation exerce reellement le code de pivot.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                             N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                             SHOW_PROGRESS, SHOW_RSS_TRACE, USE_MOVE_ALLOC, PROGRESS_CALLBACK)
    TYPE(GRID_T), INTENT(IN)          :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(INOUT) :: A_BLR
    REAL(WP), INTENT(IN)              :: TOL
    INTEGER, INTENT(OUT)              :: INFO
    INTEGER, INTENT(OUT)              :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)           :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL    :: N_SWAPS
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL     :: SHOW_PROGRESS
    ! Diagnostic memoire : affiche la RSS maximale du processus 0 observee
    ! pendant la factorisation (voir plus bas). N'affecte aucun resultat --
    ! uniquement un WRITE en fin de boucle.
    LOGICAL, INTENT(IN), OPTIONAL     :: SHOW_RSS_TRACE
    ! USE_MOVE_ALLOC (optionnel, defaut .FALSE.) : transmis tel quel a
    ! SCHUR_UPDATE_BLR_BLOCK (voir blr_compress_mod.f90) -- evite la
    ! reallocation implicite de "ABLK = NEWBLK" au profit d'un transfert
    ! de propriete (MOVE_ALLOC). Resultat numerique IDENTIQUE.
    LOGICAL, INTENT(IN), OPTIONAL     :: USE_MOVE_ALLOC
    ! Appele par TOUS les processus au franchissement d'un palier de 5%
    ! de K/NBLK -- contrairement a l'assemblage, cette boucle est
    ! GLOBALEMENT synchronisee (voir plus bas), donc tous les processus
    ! franchissent chaque palier au meme K, pas d'ecart entre eux comme a
    ! l'assemblage. Absent par defaut, aucun surcout si non fourni.
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, NS_LOCAL, li_K, PBAR_STATE
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, T_FACTOR_L, RSS_MAX_MB
    LOGICAL :: PROGRESS, RSS_TRACE

    ! Boucle K=1..NBLK GLOBALEMENT synchronisee (MPI_BCAST a chaque etape,
    ! voir plus bas) -- contrairement a l'assemblage, K/NBLK est ici un
    ! pourcentage d'avancement EXACT (tous les processus traitent le meme K
    ! au meme moment), pas une approximation. Affiche uniquement par le
    ! processus 0 pour eviter des lignes entrelacees.
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    PBAR_STATE = -1
    RSS_TRACE = .FALSE.
    IF (PRESENT(SHOW_RSS_TRACE)) RSS_TRACE = SHOW_RSS_TRACE .AND. (GRID%MYID.EQ.0)
    RSS_MAX_MB = 0.0D0

    T0 = MPI_WTIME()
    NBLK = A_BLR%NBLK
    INFO = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8; NS_LOCAL = 0
    ALLOCATE(LRECV(MAX(A_BLR%NBLKR_LOC,1)), URECV(MAX(A_BLR%NBLKC_LOC,1)))
    ALLOCATE(A_BLR%IPIV(A_BLR%NB, NBLK))
    A_BLR%IPIV = 0
    CALL BLR_WORKSPACE_ALLOC(WS, A_BLR%NB)

    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW))
      IPIV_K = 0

      ! 1) facteur diagonal AVEC PIVOT (CGETRF), local au proprietaire
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        CALL CGETRF(KW, KW, A_BLR%DIAGBLK(K)%D, KW, IPIV_K, IINFO)
        IF (IINFO .NE. 0) INFO_STEP = A_BLR%BSTART(K) - 1 + IINFO
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(LRECV, URECV, IPIV_K)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF

      ! Diffusion GLOBALE (voir commentaire de tete) -- petite (KW entiers)
      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      A_BLR%IPIV(1:KW,K) = IPIV_K

      ! 2) propage le pivot aux AUTRES blocs de la MEME ligne-bloc K, deja
      !    stockes localement (J<K finalises, J>K encore bruts) -- SEUL le
      !    proprietaire de la ligne PROW_K en possede.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li_K = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .EQ. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP(A_BLR%OFFBLK(li_K,lj), KW, IPIV_K)
        END DO
      END IF
      DO I = 1, KW
        IF (IPIV_K(I) .NE. I) NS_LOCAL = NS_LOCAL + 1
      END DO
      DEALLOCATE(IPIV_K)

      ! 3) DKK diffuse le long de la ligne PROW_K (pour U12) ET de la
      !    colonne PCOL_K (pour L21) -- inchange (DKK contient maintenant
      !    L+U de CGETRF, meme convention que BLR_LU_NOPIV_BLOCKED : L
      !    unitaire sous la diagonale, U non-unitaire sur/au-dessus).
      ALLOCATE(DKK(KW,KW))
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
      IF (GRID%MYROW .EQ. PROW_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      IF (GRID%MYCOL .EQ. PCOL_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      ! 4) U12 : la ligne PROW_K applique le TRSM sur SES blocs locaux
      !    (deja pivotes a l'etape 2 ci-dessus si J>K).
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BLR_TRSM_LEFT_LOWER_UNIT(DKK, KW, A_BLR%OFFBLK(li,lj))
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF

      ! 5) L21 : la colonne PCOL_K applique le TRSM sur SES blocs locaux
      !    (colonne K, jamais concernee par le pivot de l'etape 2).
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BLR_TRSM_RIGHT_UPPER_NONUNIT(A_BLR%OFFBLK(li,lj), DKK, KW)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF
      DEALLOCATE(DKK)

      ! 6) diffusion des panneaux L/U DEJA COMPRESSES -- inchange.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, &
                                A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)))
          LRECV(li) = A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID))
        END DO
      ELSE
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li))
        END DO
      END IF

      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, &
                                A_BLR%OFFBLK((K-1)/GRID%PR+1,lj))
          URECV(lj) = A_BLR%OFFBLK((K-1)/GRID%PR+1,lj)
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      END IF

      ! 7) mise a jour de Schur locale -- inchange.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          IF (I .EQ. J) THEN
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LRECV(li), URECV(lj), A_BLR%DIAGBLK(I)%D, MI)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK(A_BLR%OFFBLK(li,lj), MI, NJ, KW, LRECV(li), URECV(lj), TOL, WS, &
                                         USE_MOVE_ALLOC=USE_MOVE_ALLOC)
          END IF
        END DO
      END DO

      CALL BLR_TRIM_HEAP()        ! fin d'etape : memoire liberee rendue au systeme
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('2) Factorisation LU (pivot)', K, NBLK, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(K, NBLK)) CALL PROGRESS_CALLBACK(K, NBLK)
      END IF
      ! N'affiche que le maximum observe, pas une ligne par K -- le detail
      ! par etape sert au diagnostic fin de l'ecart RSS-vs-theorique (voir
      ! examples/README.md), le maximum seul suffit en usage courant.
      ! Toujours GET_RSS_KB() a chaque K (cout negligeable), juste pas de
      ! WRITE intermediaire.
      IF (RSS_TRACE) RSS_MAX_MB = MAX(RSS_MAX_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
    END DO
    IF (RSS_TRACE) WRITE(6,'(A,F8.1,A)') '     [RSS trace] MAX observe pendant la factorisation = ', RSS_MAX_MB, ' Mo'

    DEALLOCATE(LRECV, URECV)
    CALL BLR_WORKSPACE_FREE(WS)
    ! T_FACTOR est un vrai MAX sur GRID%COMM, jamais la valeur locale du
    ! rang appelant -- coherent avec T_DIAG/T_TRSM/etc. de MPI_BLR_MATRIX_
    ! LU_FACTOR_NOPIV (voir plus haut). Necessaire en toute rigueur meme si
    ! la synchronisation implicite des broadcasts par etape maintient
    ! generalement les rangs au coude-a-coude : un rang plus charge
    ! pourrait rattraper son retard pendant l'attente sur un broadcast
    ! sans que cela n'allonge le T0-a-T_WTIME() local du rang 0, faussant
    ! une simple lecture de la copie du rang 0.
    IF (PRESENT(T_FACTOR)) THEN
      T_FACTOR_L = MPI_WTIME() - T0
      CALL MPI_ALLREDUCE(T_FACTOR_L, T_FACTOR, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    ! NS_LOCAL est deja IDENTIQUE sur tous les processus (chacun le calcule
    ! a partir de SA PROPRE copie de IPIV_K, diffusee identique a tous via
    ! GRID%COMM ci-dessus) -- aucune reduction necessaire, juste le lire.
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_LOCAL
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT

  !-----------------------------------------------------------------------
  ! Resolution A*X=B, distribution 2D complete de B (comme mpi_lu_solve_
  ! 2d_mod.f90 -- ALLGATHERV+diffusion+SCATTERV, RSTARTS/RLOCS calcules
  ! localement sans communication), mais la contribution MI x NRHS de
  ! chaque bloc trainant local est calculee via BLR_BLOCK_APPLY (exploite
  ! directement la compression), pas un CGEMM dense sur tout le panneau --
  ! boucle bloc local par bloc local.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS)
    TYPE(GRID_T), INTENT(IN)           :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)         :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, RS, RE, NROWTRAIL, C, IERR
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, LAST_PCT
    LOGICAL :: PROGRESS

    NBLK = A_BLR%NBLK
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    STEP = 0
    NSTEPS = 2 * NBLK
    LAST_PCT = -5
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    ! === 0) P*B : identite ici (SANS PIVOT) -- rien a faire. ===

    ! === 1) L*Y = P*B (substitution avant) ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ! DKK n'est ALLOUE que chez le proprietaire EXACT (MYROW==PROW_K
        ! ET MYCOL==PCOL_K) -- toutes les AUTRES colonnes de cette ligne
        ! doivent recevoir une copie via diffusion avant de l'utiliser
        ! (meme motif que MPI_BLR_MATRIX_LU_FACTOR_NOPIV -- bug trouve en
        ! debuggant : le code utilisait DIAGBLK(K)%D directement, non
        ! alloue chez les non-proprietaires).
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      ! NROWTRAIL ne depend que de MYROW (via NBLKR_LOC/BSTART/BEND), donc
      ! est IDENTIQUE pour tous les processus de cette ligne -- calcule ICI
      ! par TOUS (pas seulement le proprietaire de la colonne K), condition
      ! NECESSAIRE pour que MPI_SCATTERV (appel COLLECTIF) soit appele par
      ! TOUS les membres de ROW_COMM et non par le seul proprietaire (bug
      ! trouve en debuggant : la version precedente imbriquait le
      ! SCATTERV DANS le "IF (MYCOL.EQ.PCOL_K)", violant la semantique
      ! collective de MPI_SCATTERV).
      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .GT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          ! SEUL le proprietaire de la colonne K possede les blocs OFFBLK
          ! necessaires -- lui seul remplit CONTRIB.
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .LE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          ! RS avance pour TOUS les blocs-lignes locaux (I<=K inclus,
          ! puisque B_LOCAL les concatene tous dans l'ordre li=1..),
          ! RE avance seulement pour les blocs TRAINANTS (ceux
          ! effectivement empaquetes dans CONTRIB_LOC) -- eviter de
          ! rappeler LOCAL_ROW_OFFSET (cout O(li) par appel) dans cette
          ! boucle.
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('Resolution (sans pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    ! === 2) U*X = Y (substitution arriere) -- meme principe, tetes ===
    DO K = NBLK, 1, -1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      ! Meme correction qu'en substitution avant : NROWTRAIL calcule par
      ! TOUS (ne depend que de MYROW), CONTRIB rempli par le seul
      ! proprietaire mais SCATTERV appele par TOUTE la ligne.
      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .LT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .GE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            IF (I .LT. K) THEN
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('Resolution (sans pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_2D

  !-----------------------------------------------------------------------
  ! Identique a MPI_BLR_MATRIX_LU_SOLVE_2D, pour une matrice factorisee par
  ! MPI_BLR_MATRIX_LU_FACTOR_PIVOT -- SEULE difference : "P*B" n'est plus
  ! l'identite (etape 0 du commentaire de la version sans pivot). Applique,
  ! pour chaque bloc K, la MEME sequence de swaps que celle utilisee a la
  ! factorisation (A_BLR%IPIV(:,K), deja diffusee a TOUTE la grille -- donc
  ! disponible ICI SANS communication supplementaire) sur B_LOCAL, AVANT le
  ! CTRSM de substitution avant. La substitution arriere n'a besoin d'AUCUN
  ! pivot (U*X=Y ne depend pas de P -- convention standard P*A=L*U).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS, &
                                               PROGRESS_CALLBACK)
    TYPE(GRID_T), INTENT(IN)           :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)         :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS
    ! Appele par TOUS les processus au franchissement d'un palier de 5%
    ! de STEP/NSTEPS -- boucle GLOBALEMENT synchronisee (voir plus bas),
    ! meme remarque que pour MPI_BLR_MATRIX_LU_FACTOR_PIVOT. Absent par
    ! defaut, aucun surcout si non fourni.
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, RS, RE, NROWTRAIL, C, IERR
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, PBAR_STATE
    LOGICAL :: PROGRESS

    ! Une seule barre pour toute la resolution (substitution avant PUIS
    ! arriere, 2*NBLK etapes au total) -- comme la factorisation, chaque K
    ! est GLOBALEMENT synchronise (MPI_BCAST a chaque etape), K/NBLK (ou
    ! ici STEP/NSTEPS) est donc un pourcentage exact, pas une approximation.
    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    STEP = 0
    PBAR_STATE = -1

    NBLK = A_BLR%NBLK
    NSTEPS = 2 * NBLK
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    ! === 0) P*B : applique le swap de CHAQUE bloc a B_LOCAL en une fois,
    !     avant toute substitution -- important : ceci doit preceder toute
    !     soustraction croisee entre blocs (etape 1 ci-dessous), car P_K ne
    !     doit s'appliquer qu'au B original, jamais a une valeur
    !     partiellement mise a jour (appliquer le swap de l'etape K "en
    !     ligne" au debut de l'iteration K permuterait aussi les
    !     contributions deja soustraites par les etapes m<K, ce qui serait
    !     faux). Chaque processus traite uniquement les blocs dont il
    !     possede des lignes (MYROW==PROW_K) ; IPIV deja diffuse a toute la
    !     grille, aucune communication necessaire ici. ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      IF (GRID%MYROW.NE.PROW_K .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET(A_BLR, GRID, K)
      DO I = 1, KW
        IF (A_BLR%IPIV(I,K) .NE. I) THEN
          CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+A_BLR%IPIV(I,K)-1,1), LDB)
        END IF
      END DO
    END DO

    ! === 1) L*Y = P*B (substitution avant) ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
      END IF
      ! Seule la colonne PCOL_K possede les blocs (I,K) necessaires a la
      ! correction des lignes en aval (voir BLR_BLOCK_APPLY plus bas,
      ! toujours gate sur MYCOL==PCOL_K) -- inutile de reconstituer puis
      ! diffuser le panneau X_K en LARGEUR COMPLETE (NRHS) sur TOUTES les
      ! colonnes, les autres ne l'auraient jamais lu. Un GATHERV (pas un
      ! ALLGATHERV) suffit, avec un tampon complet SEULEMENT chez la
      ! racine PCOL_K ; ailleurs, un tampon factice de taille 1 (ignore
      ! par MPI pour un appelant non-racine).
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        ALLOCATE(XK_FULL(KW,NRHS))
      ELSE
        ALLOCATE(XK_FULL(1,1))
      END IF
      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_GATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                          XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      ! Diffusion DANS la colonne PCOL_K uniquement (jamais vers les
      ! autres colonnes, qui ne l'utiliseraient pas) -- COL_COMM ne
      ! regroupe de toute facon que les processus de MEME MYCOL, donc cet
      ! IF ne desynchronise rien : tous les membres de CE COL_COMM (ceux
      ! a MYCOL==PCOL_K) l'executent identiquement.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      END IF

      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .GT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .LE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        ELSE
          ALLOCATE(CONTRIB(1,1))
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
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
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution (pivot)', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO

    ! === 2) U*X = Y (substitution arriere) -- AUCUN pivot necessaire ===
    DO K = NBLK, 1, -1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
      END IF
      ! Voir le meme changement dans la substitution avant ci-dessus :
      ! seule la colonne PCOL_K possede les blocs (I,K) necessaires ici.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        ALLOCATE(XK_FULL(KW,NRHS))
      ELSE
        ALLOCATE(XK_FULL(1,1))
      END IF
      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_GATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                          XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      END IF

      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .LT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .GE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        ELSE
          ALLOCATE(CONTRIB(1,1))
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            IF (I .LT. K) THEN
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
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution (pivot)', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D

  !-----------------------------------------------------------------------
  ! Decoupage simple par blocs de NRHS sur PC (identique a RHS_LOCAL_RANGE
  ! de mpi_lu_solve_2d_mod.f90 -- duplique ici pour eviter de dependre du
  ! module de resolution DENSE dans ce module BLR ; meme formule).
  !-----------------------------------------------------------------------
  SUBROUTINE RHS_LOCAL_RANGE_MBM(NRHS, PC, C, RSTART, RLOC)
    INTEGER, INTENT(IN)  :: NRHS, PC, C
    INTEGER, INTENT(OUT) :: RSTART, RLOC
    INTEGER :: BASE, REM
    BASE = NRHS / PC
    REM  = MOD(NRHS, PC)
    IF (C .LT. REM) THEN
      RLOC = BASE + 1; RSTART = C*(BASE+1) + 1
    ELSE
      RLOC = BASE; RSTART = REM*(BASE+1) + (C-REM)*BASE + 1
    END IF
  END SUBROUTINE RHS_LOCAL_RANGE_MBM

  FUNCTION RHS_NUMROC_MBM(NRHS, PC, MYCOL) RESULT(NRHS_LOC)
    INTEGER, INTENT(IN) :: NRHS, PC, MYCOL
    INTEGER :: NRHS_LOC, DUMMY
    CALL RHS_LOCAL_RANGE_MBM(NRHS, PC, MYCOL, DUMMY, NRHS_LOC)
  END FUNCTION RHS_NUMROC_MBM

  !-----------------------------------------------------------------------
  ! Decalage LOCAL (1-based) du premier element de B_LOCAL correspondant
  ! au bloc GLOBAL I (I doit etre possede localement en ligne, I.e.
  ! MOD(I-1,PR)==MYROW) -- B_LOCAL concatene les blocs-lignes locaux DANS
  ! L'ORDRE de li=1..NBLKR_LOC, chacun de largeur BEND(I)-BSTART(I)+1.
  !-----------------------------------------------------------------------
  FUNCTION LOCAL_ROW_OFFSET(A_BLR, GRID, IGLOBAL) RESULT(RS)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    TYPE(GRID_T), INTENT(IN)           :: GRID
    INTEGER, INTENT(IN)                :: IGLOBAL
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (IGLOBAL-1-GRID%MYROW)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
    END DO
  END FUNCTION LOCAL_ROW_OFFSET

  !=========================================================================
  ! VARIANTES ARENA (MPI_BLR_MATRIX_ARENA_T) -- voir le type ci-dessus et
  ! blr_arena_mod.f90/blr_compress_mod.f90 pour la justification complete.
  ! Chaque routine ci-dessous est le miroir EXACT de sa version non-arena
  ! (meme algorithme, meme sequence de communications MPI) -- seule la
  ! mecanique de stockage des blocs hors-diagonale change (carve dans
  ! ARENA au lieu d'ALLOCATE/DEALLOCATE individuels).
  !=========================================================================

  !-----------------------------------------------------------------------
  ! Variante ARENA de l'etape 1 (voir MPI_ASSEMBLE_BLR_MATRIX) --
  ! COMPRESS_BLOCK_ARENA au lieu de COMPRESS_BLOCK, sinon identique.
  !
  ! TARGET sur A_BLR (dummy) : NECESSAIRE, pas juste prudent -- par les
  ! regles Fortran d'association d'argument, un pointeur carve par
  ! ARENA_CARVE (au fond de la chaine d'appel, dans COMPRESS_BLOCK_ARENA)
  ! ne reste valide APRES le retour de cet appel QUE SI l'argument effectif
  ! passe a ce niveau (ici A_BLR%ARENA) a lui-meme l'attribut TARGET a CE
  ! niveau d'appel -- ce qui remonte, de proche en proche, jusqu'a exiger
  ! que A_BLR (dummy ICI) soit TARGET, ET que la variable A_BLR de
  ! l'appelant (l'exemple) le soit aussi (voir
  ! example_blr_mpi_pivot_arena.f90). Sans cela, les pointeurs UFAC/VFAC/
  ! DENSE des blocs deviendraient indefinis des la sortie de cette routine.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX_ARENA(GRID, N, NB, TOL, GET_BLOCK, A_BLR, &
                                            T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                            SHOW_PROGRESS, USE_RSVD, DIAG_SPREAD)
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir BLR_COL_OWNER) -- les blocs diagonaux visitent tous les
    ! processus ; .FALSE. : cyclique pur (ancien comportement).
    LOGICAL, INTENT(IN), OPTIONAL :: DIAG_SPREAD
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(OUT), TARGET :: A_BLR
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS
    ! USE_RSVD (defaut .TRUE.) : SVD randomisee adaptative dans un bloc
    ! ordinaire, puis copie dans l'arena (COPY_REGULAR_TO_ARENA) ; .FALSE. =
    ! CGESVD directement dans l'arena (COMPRESS_BLOCK_ARENA).
    LOGICAL, INTENT(IN), OPTIONAL   :: USE_RSVD

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    TYPE(BLR_BLOCK_T) :: TMPBLK
    LOGICAL :: RSVD_FLAG
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NDONE = 0
    LAST_PCT = -5
    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A_BLR%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A_BLR%NBLKC_LOC = BLR_NBLKC_LOC(A_BLR%CMAP, A_BLR%NBLK, GRID)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, TMPBLK)
            CALL COPY_REGULAR_TO_ARENA(TMPBLK, A_BLR%OFFBLK(li,lj), A_BLR%ARENA)
          ELSE
            CALL COMPRESS_BLOCK_ARENA(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(li,lj), A_BLR%ARENA)
          END IF
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(MI,NJ)
          END IF
          DEALLOCATE(TMP)
        END IF
        NDONE = NDONE + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('1) Assemblage+compression arena (rang 0)', NDONE, NTOTAL, LAST_PCT)
      END DO
    END DO
    T_ASSEMBLE_COMPRESS = MPI_WTIME() - T0

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX_ARENA

  !-----------------------------------------------------------------------
  ! Empreinte memoire ANALYTIQUE de la part LOCALE de A_BLR (octets) --
  ! diagonale (ALLOCATABLE, inchangee) + ARENA_BYTES (voir blr_arena_mod.f90
  ! -- INCLUT l'espace gaspille par les anciennes versions de blocs jamais
  ! recuperees : chiffre HONNETE, pas seulement la donnee utile).
  !-----------------------------------------------------------------------
  FUNCTION MPI_BLR_MATRIX_ARENA_MEMORY_BYTES(A_BLR) RESULT(NBYTES)
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: A_BLR
    INTEGER(8) :: NBYTES
    INTEGER :: K
    NBYTES = 0_8
    DO K = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(K)%D)) NBYTES = NBYTES + 8_8*INT(SIZE(A_BLR%DIAGBLK(K)%D),8)
    END DO
    NBYTES = NBYTES + ARENA_BYTES(A_BLR%ARENA)
  END FUNCTION MPI_BLR_MATRIX_ARENA_MEMORY_BYTES

  !-----------------------------------------------------------------------
  ! Diagnostic memoire de l'arena : decompose l'empreinte OFFBLK locale
  ! (ARENA_BYTES(A_BLR%ARENA)) en 3 categories :
  !
  !   USED_BYTES   : taille LOGIQUE actuelle (ce que BLR_BLOCK_BYTES
  !                  donnerait) -- la donnee VRAIMENT utile.
  !   SLACK_BYTES  : capacite haute-marque INUTILISEE des blocs ENCORE
  !                  bas-rang (UFAC_FLAT/VFAC_FLAT plus grands que le rang
  !                  ACTUEL, souvenir d'un rang plus grand a une etape
  !                  passee -- jamais retreci, voir ARENA_CARVE_REUSE).
  !   ORPHAN_BYTES : buffers UFAC_FLAT/VFAC_FLAT ABANDONNES d'un bloc
  !                  DEVENU dense entre-temps -- la bascule bas-rang->dense
  !                  est DEFINITIVE (voir SCHUR_UPDATE_BLR_BLOCK_ARENA),
  !                  donc ces buffers ne seront JAMAIS reutilises, mais
  !                  restent comptes dans ARENA_BYTES pour toujours.
  !
  ! USED+SLACK+ORPHAN doit etre <= ARENA_BYTES(A_BLR%ARENA) ; l'ecart
  ! residuel (s'il existe) est la fragmentation INTERNE aux chunks de
  ! l'arena elle-meme (fin de chunk non comblee avant un changement de
  ! chunk, voir blr_arena_mod.f90/ARENA_CARVE) -- generalement negligeable
  ! face aux deux categories ci-dessus vu la taille des chunks (16 Mo).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_ARENA_MEMORY_DIAG(A_BLR, USED_BYTES, SLACK_BYTES, ORPHAN_BYTES)
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: A_BLR
    INTEGER(8), INTENT(OUT) :: USED_BYTES, SLACK_BYTES, ORPHAN_BYTES
    INTEGER :: li, lj

    USED_BYTES = 0_8; SLACK_BYTES = 0_8; ORPHAN_BYTES = 0_8
    DO lj = 1, A_BLR%NBLKC_LOC
      DO li = 1, A_BLR%NBLKR_LOC
        ASSOCIATE (B => A_BLR%OFFBLK(li,lj))
          IF (B%IS_LOWRANK) THEN
            USED_BYTES = USED_BYTES + 8_8*(INT(B%M,8)*INT(B%RANK,8) + INT(B%RANK,8)*INT(B%N,8))
            IF (ASSOCIATED(B%UFAC_FLAT)) &
                 SLACK_BYTES = SLACK_BYTES + 8_8*(SIZE(B%UFAC_FLAT,KIND=8) - INT(B%M,8)*INT(B%RANK,8))
            IF (ASSOCIATED(B%VFAC_FLAT)) &
                 SLACK_BYTES = SLACK_BYTES + 8_8*(SIZE(B%VFAC_FLAT,KIND=8) - INT(B%RANK,8)*INT(B%N,8))
            IF (ASSOCIATED(B%DENSE_FLAT)) ORPHAN_BYTES = ORPHAN_BYTES + 8_8*SIZE(B%DENSE_FLAT,KIND=8)
          ELSE
            USED_BYTES = USED_BYTES + 8_8*INT(B%M,8)*INT(B%N,8)
            IF (ASSOCIATED(B%UFAC_FLAT)) ORPHAN_BYTES = ORPHAN_BYTES + 8_8*SIZE(B%UFAC_FLAT,KIND=8)
            IF (ASSOCIATED(B%VFAC_FLAT)) ORPHAN_BYTES = ORPHAN_BYTES + 8_8*SIZE(B%VFAC_FLAT,KIND=8)
          END IF
        END ASSOCIATE
      END DO
    END DO
  END SUBROUTINE MPI_BLR_MATRIX_ARENA_MEMORY_DIAG

  !-----------------------------------------------------------------------
  ! Variante ARENA de MPI_BLR_MATRIX_LU_FACTOR_PIVOT -- voir ce dernier
  ! pour le commentaire algorithmique complet (pivot confine au bloc
  ! diagonal, CLASWP, etc.), INCHANGE ici. SEULES differences :
  !   - etape 2 : BLR_BLOCK_ROW_SWAP_ARENA (idem, touche seulement UFAC)
  !   - etapes 4/5 : BLR_TRSM_*_ARENA (idem, corps strictement identique)
  !   - etape 6 : le proprietaire fait D'ABORD COPY_ARENA_TO_REGULAR (copie
  !     profonde du bloc arena vers LRECV/URECV, BLR_BLOCK_T ordinaire),
  !     PUIS diffuse DEPUIS cette copie via BCAST_BLR_BLOCK -- INCHANGEE,
  !     jamais modifiee (elle ne sait pas lire un BLR_BLOCK_ARENA_T, type
  !     different) : reordonner "copie-puis-diffuse" au lieu de
  !     "diffuse-puis-copie" (ordre de la version non-arena) donne un
  !     resultat RIGOUREUSEMENT identique (BCAST_BLR_BLOCK, cote racine,
  !     ne lit jamais BLK%M/N -- seulement IS_LOWRANK/RANK/UFAC/VFAC/DENSE,
  !     tous deja corrects dans la copie) sans toucher a une routine
  !     partagee avec les executables existants.
  !   - etape 7 : SCHUR_UPDATE_BLR_BLOCK_ARENA (carve dans A_BLR%ARENA) --
  !     pas de USE_MOVE_ALLOC ici, non pertinent : la reaffectation de
  !     pointeur DANS SCHUR_UPDATE_BLR_BLOCK_ARENA EST deja l'operation
  !     "zero-copie" que MOVE_ALLOC approchait cote ALLOCATABLE.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA(GRID, A_BLR, TOL, INFO, &
                                                   N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                                   SHOW_PROGRESS, SHOW_RSS_TRACE)
    TYPE(GRID_T), INTENT(IN)          :: GRID
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(INOUT), TARGET :: A_BLR
    REAL(WP), INTENT(IN)              :: TOL
    INTEGER, INTENT(OUT)              :: INFO
    INTEGER, INTENT(OUT)              :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)           :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL    :: N_SWAPS
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS
    ! Meme role que dans MPI_BLR_MATRIX_LU_FACTOR_PIVOT -- affiche la RSS
    ! du processus 0 a chaque etape K, cette fois pour MESURER si l'arena
    ! reduit effectivement l'ecart RSS-vs-theorique (voir memoire projet).
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_RSS_TRACE

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, NS_LOCAL, li_K, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS, RSS_TRACE

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    LAST_PCT = -5
    RSS_TRACE = .FALSE.
    IF (PRESENT(SHOW_RSS_TRACE)) RSS_TRACE = SHOW_RSS_TRACE .AND. (GRID%MYID.EQ.0)

    T0 = MPI_WTIME()
    NBLK = A_BLR%NBLK
    INFO = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8; NS_LOCAL = 0
    ALLOCATE(LRECV(MAX(A_BLR%NBLKR_LOC,1)), URECV(MAX(A_BLR%NBLKC_LOC,1)))
    ALLOCATE(A_BLR%IPIV(A_BLR%NB, NBLK))
    A_BLR%IPIV = 0
    CALL BLR_WORKSPACE_ALLOC(WS, A_BLR%NB)

    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW))
      IPIV_K = 0

      ! 1) facteur diagonal AVEC PIVOT (CGETRF) -- inchange (dense).
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        CALL CGETRF(KW, KW, A_BLR%DIAGBLK(K)%D, KW, IPIV_K, IINFO)
        IF (IINFO .NE. 0) INFO_STEP = A_BLR%BSTART(K) - 1 + IINFO
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(LRECV, URECV, IPIV_K)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF

      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      A_BLR%IPIV(1:KW,K) = IPIV_K

      ! 2) propage le pivot -- BLR_BLOCK_ROW_SWAP_ARENA.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li_K = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .EQ. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP_ARENA(A_BLR%OFFBLK(li_K,lj), KW, IPIV_K)
        END DO
      END IF
      DO I = 1, KW
        IF (IPIV_K(I) .NE. I) NS_LOCAL = NS_LOCAL + 1
      END DO
      DEALLOCATE(IPIV_K)

      ! 3) DKK -- inchange.
      ALLOCATE(DKK(KW,KW))
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
      IF (GRID%MYROW .EQ. PROW_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      IF (GRID%MYCOL .EQ. PCOL_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      ! 4) U12 -- BLR_TRSM_LEFT_LOWER_UNIT_ARENA.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BLR_TRSM_LEFT_LOWER_UNIT_ARENA(DKK, KW, A_BLR%OFFBLK(li,lj))
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF

      ! 5) L21 -- BLR_TRSM_RIGHT_UPPER_NONUNIT_ARENA.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BLR_TRSM_RIGHT_UPPER_NONUNIT_ARENA(A_BLR%OFFBLK(li,lj), DKK, KW)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF
      DEALLOCATE(DKK)

      ! 6) diffusion des panneaux L/U DEJA COMPRESSES -- copie-puis-diffuse
      !    (voir commentaire de tete), BCAST_BLR_BLOCK INCHANGEE.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL COPY_ARENA_TO_REGULAR(A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)), LRECV(li))
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li))
        END DO
      ELSE
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li))
        END DO
      END IF

      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL COPY_ARENA_TO_REGULAR(A_BLR%OFFBLK((K-1)/GRID%PR+1,lj), URECV(lj))
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      END IF

      ! 7) mise a jour de Schur -- SCHUR_UPDATE_BLR_BLOCK_ARENA.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          IF (I .EQ. J) THEN
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LRECV(li), URECV(lj), A_BLR%DIAGBLK(I)%D, MI)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK_ARENA(A_BLR%OFFBLK(li,lj), MI, NJ, KW, LRECV(li), URECV(lj), TOL, &
                                               A_BLR%ARENA, WS)
          END IF
        END DO
      END DO

      CALL BLR_TRIM_HEAP()        ! fin d'etape : memoire liberee rendue au systeme
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('3) Factorisation LU arena (pivot)', K, NBLK, LAST_PCT)
      IF (RSS_TRACE) WRITE(6,'(A,I4,A,F8.1,A)') '     [RSS trace arena] K=', K, ' : RSS = ', &
           REAL(GET_RSS_KB(),8)/1024.0D0, ' Mo'
    END DO

    DEALLOCATE(LRECV, URECV)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) T_FACTOR = MPI_WTIME() - T0

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_LOCAL
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA

  !-----------------------------------------------------------------------
  ! Variante ARENA de LOCAL_ROW_OFFSET -- corps STRICTEMENT identique,
  ! dupliquee car Fortran exige un type de dummy exact (pas de conversion
  ! implicite MPI_BLR_MATRIX_T <-> MPI_BLR_MATRIX_ARENA_T).
  !-----------------------------------------------------------------------
  FUNCTION LOCAL_ROW_OFFSET_ARENA(A_BLR, GRID, IGLOBAL) RESULT(RS)
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: A_BLR
    TYPE(GRID_T), INTENT(IN)                 :: GRID
    INTEGER, INTENT(IN)                      :: IGLOBAL
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (IGLOBAL-1-GRID%MYROW)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
    END DO
  END FUNCTION LOCAL_ROW_OFFSET_ARENA

  !-----------------------------------------------------------------------
  ! Variante ARENA de MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D -- IDENTIQUE (voir
  ! ce dernier pour le commentaire algorithmique complet), seule
  ! BLR_BLOCK_APPLY -> BLR_BLOCK_APPLY_ARENA change (lecture SEULE des
  ! blocs, aucune creation de pointeur persistant -- TARGET non necessaire
  ! sur A_BLR ici, contrairement a l'assemblage/factorisation).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_ARENA_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS)
    TYPE(GRID_T), INTENT(IN)                 :: GRID
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                      :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)               :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL            :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, RS, RE, NROWTRAIL, C, IERR
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, LAST_PCT
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    STEP = 0
    LAST_PCT = -5

    NBLK = A_BLR%NBLK
    NSTEPS = 2 * NBLK
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    ! === 0) P*B -- inchange ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      IF (GRID%MYROW.NE.PROW_K .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_ARENA(A_BLR, GRID, K)
      DO I = 1, KW
        IF (A_BLR%IPIV(I,K) .NE. I) THEN
          CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+A_BLR%IPIV(I,K)-1,1), LDB)
        END IF
      END DO
    END DO

    ! === 1) L*Y = P*B (substitution avant) ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET_ARENA(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .GT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .LE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY_ARENA(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('4) Resolution arena (pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    ! === 2) U*X = Y (substitution arriere) -- AUCUN pivot necessaire ===
    DO K = NBLK, 1, -1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET_ARENA(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .LT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .GE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY_ARENA(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            IF (I .LT. K) THEN
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('4) Resolution arena (pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_ARENA_2D

  !=========================================================================
  ! VARIANTES FIXEDCAP (MPI_BLR_MATRIX_FIXEDCAP_T) -- voir le type
  ! ci-dessus et BLR_BLOCK_FIXEDCAP_T (blr_compress_mod.f90) pour la
  ! justification complete. Chaque routine ci-dessous est le miroir
  ! algorithmique EXACT de sa version non-FIXEDCAP -- seule la mecanique
  ! de stockage des blocs hors-diagonale change (UFAC/VFAC alloues UNE
  ! SEULE FOIS a la capacite RMAX, jamais realloues ensuite).
  !=========================================================================

  !-----------------------------------------------------------------------
  ! Variante FIXEDCAP de l'etape 1 (voir MPI_ASSEMBLE_BLR_MATRIX) --
  ! COMPRESS_BLOCK (INCHANGEE) produit un BLR_BLOCK_T ordinaire (taille
  ! exacte) dans un temporaire, immediatement copie (COPY_REGULAR_TO_
  ! FIXEDCAP) vers la capacite RMAX definitive de A_BLR%OFFBLK(li,lj) --
  ! AUCUN attribut TARGET necessaire (composantes ALLOCATABLE, pas
  ! POINTER).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX_FIXEDCAP(GRID, N, NB, TOL, GET_BLOCK, A_BLR, &
                                               T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                               SHOW_PROGRESS, USE_RSVD, DIAG_SPREAD)
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir BLR_COL_OWNER) -- les blocs diagonaux visitent tous les
    ! processus ; .FALSE. : cyclique pur (ancien comportement).
    LOGICAL, INTENT(IN), OPTIONAL :: DIAG_SPREAD
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(OUT) :: A_BLR
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS
    LOGICAL, INTENT(IN), OPTIONAL   :: USE_RSVD   ! defaut .TRUE. (SVD randomisee adaptative)
    LOGICAL :: RSVD_FLAG

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    TYPE(BLR_BLOCK_T) :: TMPBLK
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NDONE = 0
    LAST_PCT = -5
    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A_BLR%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A_BLR%NBLKC_LOC = BLR_NBLKC_LOC(A_BLR%CMAP, A_BLR%NBLK, GRID)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, TMPBLK)
          ELSE
            CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, TMPBLK)
          END IF
          CALL COPY_REGULAR_TO_FIXEDCAP(TMPBLK, A_BLR%OFFBLK(li,lj))
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(MI,NJ)
          END IF
          DEALLOCATE(TMP)
        END IF
        NDONE = NDONE + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('1) Assemblage+compression (rang 0)', NDONE, NTOTAL, LAST_PCT)
      END DO
    END DO
    T_ASSEMBLE_COMPRESS = MPI_WTIME() - T0

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Empreinte memoire ANALYTIQUE de la part LOCALE de A_BLR (octets) --
  ! CAPACITE REELLEMENT ALLOUEE (M*RMAX+RMAX*N pour un bloc bas-rang, PAS
  ! M*RANK+RANK*N) : c'est la grandeur comparable a la RSS mesuree, le but
  ! explicite de cette variante etant d'accepter un peu plus de memoire
  ! ALLOUEE en echange de ZERO reallocation (donc, en principe, une RSS
  ! bien plus proche de cette valeur theorique que l'existant).
  !-----------------------------------------------------------------------
  FUNCTION MPI_BLR_MATRIX_FIXEDCAP_MEMORY_BYTES(A_BLR) RESULT(NBYTES)
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(IN) :: A_BLR
    INTEGER(8) :: NBYTES
    INTEGER :: li, lj, K
    NBYTES = 0_8
    DO K = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(K)%D)) NBYTES = NBYTES + 8_8*INT(SIZE(A_BLR%DIAGBLK(K)%D),8)
    END DO
    DO lj = 1, A_BLR%NBLKC_LOC
      DO li = 1, A_BLR%NBLKR_LOC
        ASSOCIATE (B => A_BLR%OFFBLK(li,lj))
          IF (B%IS_LOWRANK) THEN
            NBYTES = NBYTES + 8_8*(INT(B%M,8)*INT(B%RMAX,8) + INT(B%RMAX,8)*INT(B%N,8))
          ELSE IF (ALLOCATED(B%DENSE)) THEN
            NBYTES = NBYTES + 8_8*INT(B%M,8)*INT(B%N,8)
          END IF
        END ASSOCIATE
      END DO
    END DO
  END FUNCTION MPI_BLR_MATRIX_FIXEDCAP_MEMORY_BYTES

  !-----------------------------------------------------------------------
  ! DIAGNOSTIC : donnee VRAIMENT utile (rang REEL, comme BLR_BLOCK_BYTES),
  ! par opposition a la capacite ALLOUEE ci-dessus -- la difference est le
  ! "gaspillage RMAX-vs-RANK" assume par cette variante des la premiere
  ! compression de chaque bloc (PAS un gaspillage accumule au fil du
  ! temps comme l'ancien arena carve-systematique : ici il est CONSTANT,
  ! connu et borne des le depart, jamais aggrave par des reallocations).
  !-----------------------------------------------------------------------
  FUNCTION MPI_BLR_MATRIX_FIXEDCAP_USED_BYTES(A_BLR) RESULT(NBYTES)
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(IN) :: A_BLR
    INTEGER(8) :: NBYTES
    INTEGER :: li, lj, K
    NBYTES = 0_8
    DO K = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(K)%D)) NBYTES = NBYTES + 8_8*INT(SIZE(A_BLR%DIAGBLK(K)%D),8)
    END DO
    DO lj = 1, A_BLR%NBLKC_LOC
      DO li = 1, A_BLR%NBLKR_LOC
        ASSOCIATE (B => A_BLR%OFFBLK(li,lj))
          IF (B%IS_LOWRANK) THEN
            NBYTES = NBYTES + 8_8*(INT(B%M,8)*INT(B%RANK,8) + INT(B%RANK,8)*INT(B%N,8))
          ELSE IF (ALLOCATED(B%DENSE)) THEN
            NBYTES = NBYTES + 8_8*INT(B%M,8)*INT(B%N,8)
          END IF
        END ASSOCIATE
      END DO
    END DO
  END FUNCTION MPI_BLR_MATRIX_FIXEDCAP_USED_BYTES

  !-----------------------------------------------------------------------
  ! Variante FIXEDCAP de MPI_BLR_MATRIX_LU_FACTOR_PIVOT -- voir ce dernier
  ! pour le commentaire algorithmique complet (pivot confine au bloc
  ! diagonal, CLASWP, etc.), INCHANGE ici. SEULES differences :
  !   - etape 2 : BLR_BLOCK_ROW_SWAP_FIXEDCAP (identique, opere sur les
  !     lignes de UFAC, jamais paddees)
  !   - etapes 4/5 : BLR_TRSM_*_FIXEDCAP (identiques cote UFAC, vue
  !     tranchee cote VFAC -- voir blr_compress_mod.f90)
  !   - etape 6 : copie-puis-diffuse (COPY_FIXEDCAP_TO_REGULAR, MEME
  !     pattern que l'arena) vers LRECV/URECV (BLR_BLOCK_T ordinaires,
  !     EXACTEMENT dimensionnes au rang reel), PUIS BCAST_BLR_BLOCK
  !     INCHANGEE -- le volume MPI ne depend que du rang reel, jamais de
  !     RMAX.
  !   - etape 7 : SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP (ecrase en place, ZERO
  !     reallocation apres la premiere compression du bloc) -- pas de
  !     USE_MOVE_ALLOC ici, non pertinent (aucune reallocation a eviter).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_FIXEDCAP(GRID, A_BLR, TOL, INFO, &
                                                      N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                                      SHOW_PROGRESS, SHOW_RSS_TRACE)
    TYPE(GRID_T), INTENT(IN)          :: GRID
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(INOUT) :: A_BLR
    REAL(WP), INTENT(IN)              :: TOL
    INTEGER, INTENT(OUT)              :: INFO
    INTEGER, INTENT(OUT)              :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)           :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL    :: N_SWAPS
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_RSS_TRACE

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, NS_LOCAL, li_K, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS, RSS_TRACE

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    LAST_PCT = -5
    RSS_TRACE = .FALSE.
    IF (PRESENT(SHOW_RSS_TRACE)) RSS_TRACE = SHOW_RSS_TRACE .AND. (GRID%MYID.EQ.0)

    T0 = MPI_WTIME()
    NBLK = A_BLR%NBLK
    INFO = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8; NS_LOCAL = 0
    ALLOCATE(LRECV(MAX(A_BLR%NBLKR_LOC,1)), URECV(MAX(A_BLR%NBLKC_LOC,1)))
    ALLOCATE(A_BLR%IPIV(A_BLR%NB, NBLK))
    A_BLR%IPIV = 0
    CALL BLR_WORKSPACE_ALLOC(WS, A_BLR%NB)

    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW))
      IPIV_K = 0

      ! 1) facteur diagonal AVEC PIVOT (CGETRF) -- inchange (dense).
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        CALL CGETRF(KW, KW, A_BLR%DIAGBLK(K)%D, KW, IPIV_K, IINFO)
        IF (IINFO .NE. 0) INFO_STEP = A_BLR%BSTART(K) - 1 + IINFO
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(LRECV, URECV, IPIV_K)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF

      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      A_BLR%IPIV(1:KW,K) = IPIV_K

      ! 2) propage le pivot -- BLR_BLOCK_ROW_SWAP_FIXEDCAP.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li_K = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .EQ. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP_FIXEDCAP(A_BLR%OFFBLK(li_K,lj), KW, IPIV_K)
        END DO
      END IF
      DO I = 1, KW
        IF (IPIV_K(I) .NE. I) NS_LOCAL = NS_LOCAL + 1
      END DO
      DEALLOCATE(IPIV_K)

      ! 3) DKK -- inchange.
      ALLOCATE(DKK(KW,KW))
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
      IF (GRID%MYROW .EQ. PROW_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      IF (GRID%MYCOL .EQ. PCOL_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      ! 4) U12 -- BLR_TRSM_LEFT_LOWER_UNIT_FIXEDCAP.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BLR_TRSM_LEFT_LOWER_UNIT_FIXEDCAP(DKK, KW, A_BLR%OFFBLK(li,lj))
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF

      ! 5) L21 -- BLR_TRSM_RIGHT_UPPER_NONUNIT_FIXEDCAP.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BLR_TRSM_RIGHT_UPPER_NONUNIT_FIXEDCAP(A_BLR%OFFBLK(li,lj), DKK, KW)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF
      DEALLOCATE(DKK)

      ! 6) diffusion des panneaux L/U DEJA COMPRESSES -- copie-puis-diffuse
      !    (voir commentaire de tete), BCAST_BLR_BLOCK INCHANGEE.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL COPY_FIXEDCAP_TO_REGULAR(A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)), LRECV(li))
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li))
        END DO
      ELSE
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li))
        END DO
      END IF

      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL COPY_FIXEDCAP_TO_REGULAR(A_BLR%OFFBLK((K-1)/GRID%PR+1,lj), URECV(lj))
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      END IF

      ! 7) mise a jour de Schur -- SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          IF (I .EQ. J) THEN
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LRECV(li), URECV(lj), A_BLR%DIAGBLK(I)%D, MI)
          ELSE
            CALL SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP(A_BLR%OFFBLK(li,lj), MI, NJ, KW, LRECV(li), URECV(lj), TOL, WS)
          END IF
        END DO
      END DO

      CALL BLR_TRIM_HEAP()        ! fin d'etape : memoire liberee rendue au systeme
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('3) Factorisation LU fixedcap (pivot)', K, NBLK, LAST_PCT)
      IF (RSS_TRACE) WRITE(6,'(A,I4,A,F8.1,A)') '     [RSS trace fixedcap] K=', K, ' : RSS = ', &
           REAL(GET_RSS_KB(),8)/1024.0D0, ' Mo'
    END DO

    DEALLOCATE(LRECV, URECV)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) T_FACTOR = MPI_WTIME() - T0

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_LOCAL
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Variante FIXEDCAP de LOCAL_ROW_OFFSET -- corps STRICTEMENT identique,
  ! dupliquee car Fortran exige un type de dummy exact.
  !-----------------------------------------------------------------------
  FUNCTION LOCAL_ROW_OFFSET_FIXEDCAP(A_BLR, GRID, IGLOBAL) RESULT(RS)
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(IN) :: A_BLR
    TYPE(GRID_T), INTENT(IN)                    :: GRID
    INTEGER, INTENT(IN)                         :: IGLOBAL
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (IGLOBAL-1-GRID%MYROW)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
    END DO
  END FUNCTION LOCAL_ROW_OFFSET_FIXEDCAP

  !-----------------------------------------------------------------------
  ! Variante FIXEDCAP de MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D -- IDENTIQUE
  ! (voir ce dernier pour le commentaire algorithmique complet), seule
  ! BLR_BLOCK_APPLY -> BLR_BLOCK_APPLY_FIXEDCAP change (lecture SEULE,
  ! rang REEL uniquement -- voir blr_compress_mod.f90).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_FIXEDCAP_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS)
    TYPE(GRID_T), INTENT(IN)                    :: GRID
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                         :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)                  :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL               :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, RS, RE, NROWTRAIL, C, IERR
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, LAST_PCT
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    STEP = 0
    LAST_PCT = -5

    NBLK = A_BLR%NBLK
    NSTEPS = 2 * NBLK
    ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
    DO C = 0, GRID%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    ! === 0) P*B -- inchange ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      IF (GRID%MYROW.NE.PROW_K .OR. NRHS_LOC.LE.0) CYCLE
      RS = LOCAL_ROW_OFFSET_FIXEDCAP(A_BLR, GRID, K)
      DO I = 1, KW
        IF (A_BLR%IPIV(I,K) .NE. I) THEN
          CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+A_BLR%IPIV(I,K)-1,1), LDB)
        END IF
      END DO
    END DO

    ! === 1) L*Y = P*B (substitution avant) ===
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET_FIXEDCAP(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .GT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .LE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY_FIXEDCAP(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('4) Resolution fixedcap (pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    ! === 2) U*X = Y (substitution arriere) -- AUCUN pivot necessaire ===
    DO K = NBLK, 1, -1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

      ALLOCATE(XK_FULL(KW,NRHS))
      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) DKK = A_BLR%DIAGBLK(K)%D
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

        RS = LOCAL_ROW_OFFSET_FIXEDCAP(A_BLR, GRID, K)
        IF (NRHS_LOC .GT. 0) THEN
          CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS_LOC, &
                      CONE, DKK, KW, B_LOCAL(RS,1), LDB)
        END IF
        DEALLOCATE(DKK)
        ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
        IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
      NROWTRAIL = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .LT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END DO
      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CONTRIB = CZERO
          RE = 0
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            IF (I .GE. K) CYCLE
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            CALL BLR_BLOCK_APPLY_FIXEDCAP(A_BLR%OFFBLK(li,lj), NRHS, CONE, XK_FULL, KW, CONE, CONTRIB(RE+1,1), NROWTRAIL)
            RE = RE + MI
          END DO
        END IF
        ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
        END DO
        CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                           CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NRHS_LOC .GT. 0) THEN
          RE = 0; RS = 1
          DO li = 1, A_BLR%NBLKR_LOC
            I = (li-1)*GRID%PR + GRID%MYROW + 1
            MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
            IF (I .LT. K) THEN
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('4) Resolution fixedcap (pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_FIXEDCAP_2D

  !=========================================================================
  ! VARIANTES OUT-OF-CORE (MPI_BLR_MATRIX_OOC_T) -- voir le type ci-dessus
  ! et blr_ooc_mod.f90 pour le principe complet. Basees sur l'algorithme
  ! CLASSIQUE (MPI_BLR_MATRIX_LU_FACTOR_PIVOT), jamais l'arena.
  !=========================================================================

  !-----------------------------------------------------------------------
  ! Assemblage + compression, identique algorithmiquement a
  ! MPI_ASSEMBLE_BLR_MATRIX -- seule difference : les blocs diagonaux sont
  ! ecrits directement sur STORE_A (jamais residents en dehors de cette
  ! iteration), les blocs hors-diagonale passent par CACHE_A (borne a
  ! MAX_MB_PER_PROCESS -- l'eviction naturelle du cache pendant la boucle
  ! ecrit deja la plupart des blocs sur disque au fil de l'eau ;
  ! OOC_CACHE_FLUSH_ALL en fin de routine garantit que TOUT est sur disque
  ! avant que l'appelant ne passe a l'etape suivante, meme les derniers
  ! blocs encore residents faute d'avoir ete evinces).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                          T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                          SHOW_PROGRESS, T_GET_BLOCK, T_COMPRESS, T_IO, &
                                          T_GET_BLOCK_MIN, T_COMPRESS_MIN, T_IO_MIN, &
                                          USE_RSVD, RSVD_OVERSAMPLE, RSVD_NPOWER, PROGRESS_CALLBACK, DIAG_SPREAD)
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    ! DIAG_SPREAD (defaut .TRUE.) : distribution decalee des blocs-colonnes
    ! (CMAP=1, voir BLR_COL_OWNER) -- les blocs diagonaux visitent tous les
    ! processus ; .FALSE. : cyclique pur (ancien comportement).
    LOGICAL, INTENT(IN), OPTIONAL :: DIAG_SPREAD
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(OUT) :: A_BLR
    CHARACTER(LEN=*), INTENT(IN)    :: SCRATCH_DIR
    REAL(WP), INTENT(IN)            :: MAX_MB_PER_PROCESS
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS
    ! DIAGNOSTIC (utile pour evaluer si une SVD randomisee vaut le coup a
    ! la place de CGESVD dans COMPRESS_BLOCK) : T_GET_BLOCK cumule le
    ! temps dans tous les appels GET_BLOCK ; T_COMPRESS cumule le temps dans
    ! tous les appels COMPRESS_BLOCK (SVD complete) ; T_IO cumule le temps
    ! d'ecriture disque (bloc diagonal, insertion dans CACHE_A -- qui peut
    ! evincer donc ecrire d'autres blocs -- et le flush final des deux
    ! caches). Optionnels, aucun surcout si absents.
    !
    ! T_ASSEMBLE_COMPRESS/T_GET_BLOCK/T_COMPRESS/T_IO sont un vrai MAX sur
    ! GRID%COMM (voir le commentaire equivalent sur MPI_ASSEMBLE_BLR_MATRIX
    ! ci-dessus pour la justification) ; _MIN optionnels en plus pour
    ! reveler un desequilibre de charge entre processus.
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_GET_BLOCK, T_COMPRESS, T_IO
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_GET_BLOCK_MIN, T_COMPRESS_MIN, T_IO_MIN
    ! Voir MPI_ASSEMBLE_BLR_MATRIX pour le detail (USE_RSVD/RSVD_OVERSAMPLE/
    ! RSVD_NPOWER, meme convention).
    LOGICAL, INTENT(IN), OPTIONAL   :: USE_RSVD
    INTEGER, INTENT(IN), OPTIONAL   :: RSVD_OVERSAMPLE, RSVD_NPOWER
    ! Voir MPI_ASSEMBLE_BLR_MATRIX pour le detail (meme convention -- appele
    ! par CHAQUE processus au franchissement d'un palier de 5% de SA PROPRE
    ! progression LOCALE). Absent par defaut, aucun surcout si non fourni.
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:), TMPD(:,:)
    TYPE(BLR_BLOCK_T) :: TMPBLK, DBLK
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, IDX, PBAR_STATE
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0, T_GETBLK_L, T_COMPRESS_L, T_IO_L, T_ASM_L
    LOGICAL :: PROGRESS, RSVD_FLAG
    INTEGER :: RSVD_OVS, RSVD_NPOW

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NDONE = 0
    PBAR_STATE = -1
    T_GETBLK_L = 0.0D0; T_COMPRESS_L = 0.0D0; T_IO_L = 0.0D0
    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD
    RSVD_OVS = 10; IF (PRESENT(RSVD_OVERSAMPLE)) RSVD_OVS = RSVD_OVERSAMPLE
    RSVD_NPOW = 1; IF (PRESENT(RSVD_NPOWER)) RSVD_NPOW = RSVD_NPOWER

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%CMAP = 1
    IF (PRESENT(DIAG_SPREAD)) A_BLR%CMAP = MERGE(1, 0, DIAG_SPREAD)
    A_BLR%NBLKC_LOC = BLR_NBLKC_LOC(A_BLR%CMAP, A_BLR%NBLK, GRID)
    A_BLR%MAX_MB_PER_PROCESS = MAX_MB_PER_PROCESS
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    CALL OOC_STORE_INIT(A_BLR%STORE_A, SCRATCH_DIR, 'a', GRID%MYID)
    CALL OOC_CACHE_INIT(A_BLR%CACHE_A, MAX_MB_PER_PROCESS, NTOTAL+4, &
                         MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1))
    CALL OOC_DIAG_CACHE_INIT(A_BLR%CACHE_D, MAX_MB_PER_PROCESS, A_BLR%NBLKR_LOC+2, A_BLR%NBLK)

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(TMPD(MI,MI))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMPD)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
          DBLK%IS_LOWRANK = .FALSE.
          DBLK%M = MI; DBLK%N = MI
          CALL MOVE_ALLOC(TMPD, DBLK%DENSE)
          TT0 = MPI_WTIME()
          CALL OOC_WRITE_BLOCK_FIXED(A_BLR%STORE_A, A_BLR%DIAGBLK(I), DBLK)
          T_IO_L = T_IO_L + (MPI_WTIME() - TT0)
          DEALLOCATE(DBLK%DENSE)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          TT0 = MPI_WTIME()
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          T_GETBLK_L = T_GETBLK_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, TMPBLK, RSVD_OVS, RSVD_NPOW)
          ELSE
            CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, TMPBLK)
          END IF
          T_COMPRESS_L = T_COMPRESS_L + (MPI_WTIME() - TT0)
          DEALLOCATE(TMP)
          NO_LOCAL = NO_LOCAL + 1
          IF (TMPBLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + TMPBLK%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(MI,NJ)
          END IF
          A_BLR%OFFBLK(li,lj)%IS_LOWRANK = TMPBLK%IS_LOWRANK
          A_BLR%OFFBLK(li,lj)%M = TMPBLK%M; A_BLR%OFFBLK(li,lj)%N = TMPBLK%N; A_BLR%OFFBLK(li,lj)%RANK = TMPBLK%RANK
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_PUT_NEW(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, lj, TMPBLK)
          T_IO_L = T_IO_L + (MPI_WTIME() - TT0)
        END IF
        NDONE = NDONE + 1
        ! Voir MPI_ASSEMBLE_BLR_MATRIX pour le detail du choix barre vs paliers.
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('1) Assemblage+compression ooc (rang 0)', NDONE, NTOTAL, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(NDONE, NTOTAL)) CALL PROGRESS_CALLBACK(NDONE, NTOTAL)
        END IF
      END DO
    END DO

    ! Garantit que TOUT est sur disque avant l'etape suivante (voir
    ! commentaire de tete) -- les descripteurs OFFBLK/DIAGBLK restent
    ! valides quel que soit l'etat du cache apres cet appel.
    TT0 = MPI_WTIME()
    CALL OOC_CACHE_FLUSH_ALL(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK)
    T_IO_L = T_IO_L + (MPI_WTIME() - TT0)
    T_ASM_L = MPI_WTIME() - T0
    CALL MPI_ALLREDUCE(T_ASM_L, T_ASSEMBLE_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_GET_BLOCK)) CALL MPI_ALLREDUCE(T_GETBLK_L,   T_GET_BLOCK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_COMPRESS))  CALL MPI_ALLREDUCE(T_COMPRESS_L, T_COMPRESS,  1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_IO))        CALL MPI_ALLREDUCE(T_IO_L,       T_IO,        1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_GET_BLOCK_MIN)) CALL MPI_ALLREDUCE(T_GETBLK_L,   T_GET_BLOCK_MIN, 1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (PRESENT(T_COMPRESS_MIN))  CALL MPI_ALLREDUCE(T_COMPRESS_L, T_COMPRESS_MIN,  1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (PRESENT(T_IO_MIN))        CALL MPI_ALLREDUCE(T_IO_L,       T_IO_MIN,        1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
  END SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX_OOC

  !-----------------------------------------------------------------------
  ! Assemblage du second membre, separement de la matrice (fichier
  ! scratch propre, prefixe "b"), par paquets d'au plus NRHS_BATCH
  ! colonnes globales -- chaque paquet assemble est ecrit sur disque puis
  ! libere avant de passer au suivant (voir OOC_RHS_T) : la memoire a
  ! cette etape est bornee par NRHS_BATCH, jamais par NRHS total.
  ! L'assemblage reste une etape separee de la resolution.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
    TYPE(GRID_T), INTENT(IN)                :: GRID
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(IN)  :: A_BLR
    INTEGER, INTENT(IN)                     :: NRHS, NRHS_BATCH
    PROCEDURE(FILL_RHS_IFACE)               :: FILL_RHS
    CHARACTER(LEN=*), INTENT(IN)            :: SCRATCH_DIR
    TYPE(OOC_RHS_T), INTENT(OUT)            :: RHS

    TYPE(BLR_BLOCK_T) :: BBLK
    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, I, RS, WIDTH, B, RSTART, REND

    RHS%NRHS = NRHS
    RHS%NRHS_BATCH = MAX(NRHS_BATCH,1)
    RHS%NBATCHES = (NRHS + RHS%NRHS_BATCH - 1) / RHS%NRHS_BATCH

    RHS%MLOC = 0
    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RHS%MLOC = RHS%MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
    END DO

    CALL OOC_STORE_INIT(RHS%STORE_B, SCRATCH_DIR, 'b', GRID%MYID)
    ALLOCATE(RHS%BATCH(RHS%NBATCHES))

    ! Pour CHAQUE paquet, seule la PART LOCALE de CE processus (lignes ET
    ! colonnes -- meme distribution 2D que B_LOCAL dans la resolution
    ! classique, RHS_LOCAL_RANGE_MBM sur la LARGEUR du paquet) est
    ! assemblee et ecrite : la resolution OOC (ci-dessous) doit retrouver
    ! EXACTEMENT ce que B_LOCAL contiendrait pour ce paquet dans
    ! l'algorithme non-OOC, aucune redistribution supplementaire.
    DO B = 1, RHS%NBATCHES
      RSTART = (B-1)*RHS%NRHS_BATCH + 1
      REND = MIN(B*RHS%NRHS_BATCH, NRHS)
      WIDTH = REND - RSTART + 1
      BLOCK
        INTEGER :: LOCAL_OFFSET, NRHS_LOC_B, GSTART, GEND
        CALL RHS_LOCAL_RANGE_MBM(WIDTH, GRID%PC, GRID%MYCOL, LOCAL_OFFSET, NRHS_LOC_B)
        GSTART = RSTART + LOCAL_OFFSET - 1
        GEND = GSTART + NRHS_LOC_B - 1
        ALLOCATE(TMP(MAX(RHS%MLOC,1),MAX(NRHS_LOC_B,1)))
        RS = 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (NRHS_LOC_B .GT. 0) THEN
            CALL FILL_RHS(A_BLR%BSTART(I), A_BLR%BEND(I), GSTART, GEND, &
                           TMP(RS:RS+A_BLR%BEND(I)-A_BLR%BSTART(I),1:NRHS_LOC_B))
          END IF
          RS = RS + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
        END DO
        BBLK%IS_LOWRANK = .FALSE.
        BBLK%M = MAX(RHS%MLOC,1); BBLK%N = MAX(NRHS_LOC_B,1)
        CALL MOVE_ALLOC(TMP, BBLK%DENSE)
        CALL OOC_WRITE_BLOCK_RAW(RHS%STORE_B, RHS%BATCH(B), BBLK)
        DEALLOCATE(BBLK%DENSE)
      END BLOCK
    END DO
  END SUBROUTINE MPI_ASSEMBLE_RHS_OOC

  SUBROUTINE OOC_RHS_FREE(RHS)
    TYPE(OOC_RHS_T), INTENT(INOUT) :: RHS
    CALL OOC_STORE_FREE(RHS%STORE_B)
    IF (ALLOCATED(RHS%BATCH)) DEALLOCATE(RHS%BATCH)
  END SUBROUTINE OOC_RHS_FREE

  SUBROUTINE MPI_BLR_MATRIX_OOC_FREE(A_BLR)
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(INOUT) :: A_BLR
    CALL OOC_CACHE_FREE(A_BLR%CACHE_A)
    CALL OOC_DIAG_CACHE_FREE(A_BLR%CACHE_D)
    CALL OOC_STORE_FREE(A_BLR%STORE_A)
  END SUBROUTINE MPI_BLR_MATRIX_OOC_FREE

  !-----------------------------------------------------------------------
  ! DIAGNOSTIC : octets REELLEMENT residents en RAM pour CE processus a
  ! l'instant de l'appel (CACHE_A%USED_BYTES + CACHE_D%USED_BYTES, chacun
  ! <= son budget MAX_MB_PER_PROCESS converti en octets, a l'exception du
  ! cas degenere ou un SEUL bloc depasse deja le budget a lui seul -- voir
  ! OOC_CACHE_MAKE_ROOM/OOC_DIAG_CACHE_MAKE_ROOM). Sert a VERIFIER que le
  ! budget est bien respecte, par opposition a la taille TOTALE sur disque
  ! (non bornee par design).
  !-----------------------------------------------------------------------
  FUNCTION MPI_BLR_MATRIX_OOC_RESIDENT_BYTES(A_BLR) RESULT(NBYTES)
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(IN) :: A_BLR
    INTEGER(8) :: NBYTES
    NBYTES = A_BLR%CACHE_A%USED_BYTES + A_BLR%CACHE_D%USED_BYTES
  END FUNCTION MPI_BLR_MATRIX_OOC_RESIDENT_BYTES

  !-----------------------------------------------------------------------
  ! Factorisation LU en place, avec pivot -- meme algorithme exact que
  ! MPI_BLR_MATRIX_LU_FACTOR_PIVOT (voir son commentaire de tete pour le
  ! detail complet), avec des points de chargement/ecriture explicites la
  ! ou l'original lit/ecrit
  ! A_BLR%DIAGBLK(K)%D ou A_BLR%OFFBLK(li,lj) directement :
  !   - DIAGBLK : passe par CACHE_D (voir blr_ooc_mod.f90), meme principe
  !     LRU que CACHE_A -- un bloc diagonal reste resident entre deux
  !     touches successives (mises a jour de Schur, une par etape K<I) et
  !     n'est ecrit sur disque qu'a son eviction, EN PLACE via
  !     OOC_WRITE_BLOCK_FIXED (aucun octet orphelin, contrairement a
  !     l'append-only : un bloc diagonal est toujours dense et de taille
  !     fixe pour un K donne).
  !   - OFFBLK : passe par CACHE_A (voir blr_ooc_mod.f90), ecrit EN PLACE a
  !     l'eviction (OOC_WRITE_BLOCK_FIXED, meme routine que DIAGBLK -- le
  !     rang peut varier d'une reecriture a l'autre, mais le premier appel
  !     reserve toujours le pire cas dense, donc aucun octet orphelin
  !     malgre cette variation) -- des qu'un bloc panneau (etape 6) a fini
  !     son role (copie dans LRECV/URECV), il est explicitement EVINCE
  !     (ecrit + libere) : il ne sera plus jamais touche avant la
  !     resolution. Les blocs de la sous-matrice traineante (etape 7)
  !     restent dans le cache jusqu'a ce que le budget MAX_MB_PER_PROCESS
  !     force une eviction LRU -- la boucle naturelle sur (li,lj) suffit a
  !     borner la RAM, aucune logique de "streaming par colonne" explicite
  !     necessaire.
  !
  ! OOC_CACHE_GET/OOC_DIAG_CACHE_GET utilisent un lookup direct O(1)
  ! (SLOT_OF(li,lj)/SLOT_OF(K), voir blr_ooc_mod.f90) + une pile
  ! d'emplacements libres (FREE_STACK) pour l'insertion, plutot qu'une
  ! recherche lineaire parmi les emplacements du cache a chaque acces --
  ! seule la recherche du bloc LRU a evincer (OOC_CACHE_MAKE_ROOM/
  ! OOC_DIAG_CACHE_MAKE_ROOM) reste lineaire, mais seulement lors d'un
  ! miss necessitant une eviction, pas a chaque acces.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                                 N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                                 SHOW_PROGRESS, SHOW_RSS_TRACE, &
                                                 T_CACHE_GET, T_CACHE_EVICT, T_COMPUTE, T_MPICALL, &
                                                 T_CACHE_GET_MIN, T_CACHE_EVICT_MIN, T_COMPUTE_MIN, T_MPICALL_MIN, &
                                                 PROGRESS_CALLBACK)
    TYPE(GRID_T), INTENT(IN)              :: GRID
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(INOUT) :: A_BLR
    REAL(WP), INTENT(IN)                  :: TOL
    INTEGER, INTENT(OUT)                  :: INFO
    INTEGER, INTENT(OUT)                  :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)               :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL        :: N_SWAPS
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL         :: SHOW_PROGRESS
    LOGICAL, INTENT(IN), OPTIONAL         :: SHOW_RSS_TRACE
    ! DIAGNOSTIC (utile pour isoler ou part le temps de la
    ! factorisation OOC face au classique) : T_CACHE_GET cumule le temps
    ! passe dans TOUS les appels OOC_CACHE_GET/OOC_DIAG_CACHE_GET (lookup
    ! O(1) + lecture disque quand le bloc n'est pas deja resident) ;
    ! T_CACHE_EVICT cumule le temps des evictions de panneau (etape 6,
    ! TOUJOURS reelles, ecriture disque inconditionnelle quel que soit le
    ! budget) ; T_COMPUTE cumule le calcul pur (CGETRF, TRSM, mises a jour
    ! de Schur, swap de ligne) ; T_MPICALL cumule les diffusions MPI
    ! (BCAST de IPIV_K/DKK et BCAST_BLR_BLOCK des panneaux). Ces 4
    ! categories couvrent la totalite de T_FACTOR (a la precision des
    ! quelques operations non chronometrees pres, negligeables : ALLOCATE/
    ! DEALLOCATE/MOVE_ALLOC). Tous optionnels, aucun surcout si absents.
    ! Si T_CACHE_GET domine T_FACTOR alors que T_COMPUTE reste comparable
    ! au chemin classique, ce n'est PAS un defaut du cache -- voir le
    ! commentaire "DIMENSIONNER MAX_MB_PER_PROCESS" sur MPI_BLR_MATRIX_
    ! OOC_T plus haut : c'est le signe que le budget est sous l'ensemble
    ! de travail reel, augmenter MAX_MB_PER_PROCESS est le seul vrai levier.
    !
    ! T_CACHE_GET/T_CACHE_EVICT/T_COMPUTE/T_MPICALL sont un vrai MAX sur
    ! GRID%COMM, jamais la valeur locale du rang appelant -- coherent avec
    ! T_DIAG/T_TRSM/etc. de MPI_BLR_MATRIX_LU_FACTOR_NOPIV. Nouveaux
    ! _MIN (optionnels, meme cout : un MPI_ALLREDUCE de plus par
    ! categorie) donnent en prime un INDICATEUR DE DESEQUILIBRE DE CHARGE
    ! direct : un ecart MAX/MIN important sur T_CACHE_GET par exemple
    ! signale qu'un sous-ensemble de processus subit beaucoup plus
    ! d'evictions/lectures disque que les autres (distribution du travail
    ! non uniforme), invisible tant qu'on ne voit que la valeur MAX seule.
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_CACHE_GET, T_CACHE_EVICT, T_COMPUTE, T_MPICALL
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_CACHE_GET_MIN, T_CACHE_EVICT_MIN, T_COMPUTE_MIN, T_MPICALL_MIN
    ! Voir MPI_BLR_MATRIX_LU_FACTOR_PIVOT pour le detail (meme convention).
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), DIAGTMP(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP, IDX, IDX_D
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, NS_LOCAL, li_K, PBAR_STATE
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, TT0, T_FACTOR_L, RSS_MAX_MB, CACHE_MAX_MB
    DOUBLE PRECISION :: T_GET_L, T_EVICT_L, T_COMPUTE_L, T_MPI_L
    LOGICAL :: PROGRESS, RSS_TRACE

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    PBAR_STATE = -1
    RSS_TRACE = .FALSE.
    IF (PRESENT(SHOW_RSS_TRACE)) RSS_TRACE = SHOW_RSS_TRACE .AND. (GRID%MYID.EQ.0)
    RSS_MAX_MB = 0.0D0; CACHE_MAX_MB = 0.0D0

    T0 = MPI_WTIME()
    T_GET_L = 0.0D0; T_EVICT_L = 0.0D0; T_COMPUTE_L = 0.0D0; T_MPI_L = 0.0D0
    NBLK = A_BLR%NBLK
    INFO = 0
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8; NS_LOCAL = 0
    ALLOCATE(LRECV(MAX(A_BLR%NBLKR_LOC,1)), URECV(MAX(A_BLR%NBLKC_LOC,1)))
    ALLOCATE(A_BLR%IPIV(A_BLR%NB, NBLK))
    A_BLR%IPIV = 0
    CALL BLR_WORKSPACE_ALLOC(WS, A_BLR%NB)

    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW))
      IPIV_K = 0

      ! 1) facteur diagonal AVEC PIVOT (CGETRF) -- CHARGE DIAGBLK(K) via
      !    CACHE_D (proprietaire seulement, comme l'original), factorise
      !    dans DKK directement (evite un aller-retour supplementaire :
      !    DKK est de toute facon diffuse juste apres, voir etape 3). Le
      !    slot de cache reste "creux" (DENSE deplace vers DKK) jusqu'a
      !    l'etape 3, qui le repeuple -- rien d'autre ne touche CACHE_D
      !    pour ce K entre les deux (voir etape 2, uniquement CACHE_A).
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        TT0 = MPI_WTIME()
        IDX_D = OOC_DIAG_CACHE_GET(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK, K)
        T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
        TT0 = MPI_WTIME()
        CALL MOVE_ALLOC(A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%DENSE, DKK)
        CALL CGETRF(KW, KW, DKK, KW, IPIV_K, IINFO)
        T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
        IF (IINFO .NE. 0) INFO_STEP = A_BLR%BSTART(K) - 1 + IINFO
      END IF
      TT0 = MPI_WTIME()
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(LRECV, URECV)
        IF (ALLOCATED(IPIV_K)) DEALLOCATE(IPIV_K)
        IF (ALLOCATED(DKK)) DEALLOCATE(DKK)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF

      TT0 = MPI_WTIME()
      CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
      T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)
      A_BLR%IPIV(1:KW,K) = IPIV_K

      ! 2) propage le pivot -- charge/modifie/marque DIRTY via le cache.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li_K = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .EQ. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li_K, lj)
          T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BLR_BLOCK_ROW_SWAP(A_BLR%CACHE_A%SLOTS(IDX)%BLK, KW, IPIV_K)
          T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
          CALL OOC_CACHE_MARK_DIRTY(A_BLR%CACHE_A, IDX)
        END DO
      END IF
      DO I = 1, KW
        IF (IPIV_K(I) .NE. I) NS_LOCAL = NS_LOCAL + 1
      END DO
      DEALLOCATE(IPIV_K)

      ! 3) DKK diffuse -- deja factorise (etape 1) chez le proprietaire,
      !    ALLOUE ICI pour les non-proprietaires (l'original le fait tout
      !    aussi, DKK n'existe pas encore chez eux a ce point).
      IF (.NOT. ALLOCATED(DKK)) ALLOCATE(DKK(KW,KW))
      TT0 = MPI_WTIME()
      IF (GRID%MYROW .EQ. PROW_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      IF (GRID%MYCOL .EQ. PCOL_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)

      ! DIAGBLK(K) a fini son role dans la factorisation elle-meme (jamais
      ! retouche avant la resolution EN TANT QUE PIVOT -- mais peut encore
      ! recevoir des mises a jour de Schur si K est aussi un indice I>K'
      ! d'une etape ulterieure... non, K est fixe ici, ce cas ne se produit
      ! pas pour SON PROPRE indice K) -- le proprietaire repeuple le slot
      ! de cache (DIRTY), aucune ecriture disque immediate : elle n'aura
      ! lieu qu'a l'eviction (LRU ou flush final).
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        TT0 = MPI_WTIME()
        ALLOCATE(DIAGTMP(KW,KW))
        DIAGTMP = DKK
        A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%IS_LOWRANK = .FALSE.
        A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%M = KW; A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%N = KW
        CALL MOVE_ALLOC(DIAGTMP, A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%DENSE)
        CALL OOC_DIAG_CACHE_MARK_DIRTY(A_BLR%CACHE_D, IDX_D)
        T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
      END IF

      ! 4) U12 -- charge/modifie/marque DIRTY via le cache.
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, lj)
          T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BLR_TRSM_LEFT_LOWER_UNIT(DKK, KW, A_BLR%CACHE_A%SLOTS(IDX)%BLK)
          T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
          CALL OOC_CACHE_MARK_DIRTY(A_BLR%CACHE_A, IDX)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF

      ! 5) L21 -- idem.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, lj)
          T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BLR_TRSM_RIGHT_UPPER_NONUNIT(A_BLR%CACHE_A%SLOTS(IDX)%BLK, DKK, KW)
          T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
          CALL OOC_CACHE_MARK_DIRTY(A_BLR%CACHE_A, IDX)
          NO_LOCAL = NO_LOCAL + 1
          IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li,lj)%RANK
          ELSE
            RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li,lj)%M, A_BLR%OFFBLK(li,lj)%N)
          END IF
        END DO
      END IF
      DEALLOCATE(DKK)

      ! 6) diffusion des panneaux L/U DEJA COMPRESSES, copie dans
      !    LRECV/URECV -- le proprietaire EVINCE (ecrit+libere) le bloc
      !    panneau juste apres la copie, il est definitivement fini.
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, BLR_COL_LOCAL(A_BLR%CMAP, K, GRID))
          T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, &
                                A_BLR%CACHE_A%SLOTS(IDX)%BLK)
          T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)
          LRECV(li) = A_BLR%CACHE_A%SLOTS(IDX)%BLK
          TT0 = MPI_WTIME()
          CALL OOC_CACHE_EVICT_SLOT(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, IDX)
          T_EVICT_L = T_EVICT_L + (MPI_WTIME() - TT0)
        END DO
      ELSE
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, A_BLR%BEND(I)-A_BLR%BSTART(I)+1, KW, LRECV(li))
          T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)
        END DO
      END IF

      IF (GRID%MYROW .EQ. PROW_K) THEN
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          TT0 = MPI_WTIME()
          IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, (K-1)/GRID%PR+1, lj)
          T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, &
                                A_BLR%CACHE_A%SLOTS(IDX)%BLK)
          T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)
          URECV(lj) = A_BLR%CACHE_A%SLOTS(IDX)%BLK
          TT0 = MPI_WTIME()
          CALL OOC_CACHE_EVICT_SLOT(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, IDX)
          T_EVICT_L = T_EVICT_L + (MPI_WTIME() - TT0)
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .LE. K) CYCLE
          TT0 = MPI_WTIME()
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
          T_MPI_L = T_MPI_L + (MPI_WTIME() - TT0)
        END DO
      END IF

      ! 7) mise a jour de Schur locale -- diagonale : passe par CACHE_D ;
      !    hors-diagonale : passe par CACHE_A, DESCRIPTEUR resynchronise
      !    apres coup (RANG/IS_LOWRANK peuvent changer).
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          IF (I .EQ. J) THEN
            TT0 = MPI_WTIME()
            IDX_D = OOC_DIAG_CACHE_GET(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK, I)
            T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
            TT0 = MPI_WTIME()
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LRECV(li), URECV(lj), A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%DENSE, MI)
            T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
            CALL OOC_DIAG_CACHE_MARK_DIRTY(A_BLR%CACHE_D, IDX_D)
          ELSE
            TT0 = MPI_WTIME()
            IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, lj)
            T_GET_L = T_GET_L + (MPI_WTIME() - TT0)
            TT0 = MPI_WTIME()
            CALL SCHUR_UPDATE_BLR_BLOCK(A_BLR%CACHE_A%SLOTS(IDX)%BLK, MI, NJ, KW, LRECV(li), URECV(lj), TOL, WS)
            T_COMPUTE_L = T_COMPUTE_L + (MPI_WTIME() - TT0)
            A_BLR%OFFBLK(li,lj)%IS_LOWRANK = A_BLR%CACHE_A%SLOTS(IDX)%BLK%IS_LOWRANK
            A_BLR%OFFBLK(li,lj)%RANK = A_BLR%CACHE_A%SLOTS(IDX)%BLK%RANK
            CALL OOC_CACHE_MARK_DIRTY(A_BLR%CACHE_A, IDX)
          END IF
        END DO
      END DO

      CALL BLR_TRIM_HEAP()        ! fin d'etape : memoire liberee rendue au systeme
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Factorisation LU ooc (pivot)', K, NBLK, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(K, NBLK)) CALL PROGRESS_CALLBACK(K, NBLK)
      END IF
      ! N'affiche que le maximum observe (RSS et cache resident), pas une
      ! ligne par K -- voir le commentaire equivalent sur MPI_BLR_MATRIX_
      ! LU_FACTOR_PIVOT (version classique).
      IF (RSS_TRACE) THEN
        RSS_MAX_MB = MAX(RSS_MAX_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
        CACHE_MAX_MB = MAX(CACHE_MAX_MB, REAL(A_BLR%CACHE_A%USED_BYTES+A_BLR%CACHE_D%USED_BYTES,8)/1024.0D0/1024.0D0)
      END IF
    END DO
    IF (RSS_TRACE) WRITE(6,'(A,F8.1,A,F8.1,A)') '     [RSS trace ooc] MAX observe pendant la factorisation : RSS = ', &
         RSS_MAX_MB, ' Mo , cache resident = ', CACHE_MAX_MB, ' Mo'

    ! Fin de factorisation : tout ce qui reste resident (blocs traineants
    ! jamais evinces faute d'avoir atteint le budget) part sur disque --
    ! la resolution doit pouvoir tout relire depuis STORE_A. Compte comme
    ! T_CACHE_EVICT (ecriture disque, meme categorie que les evictions de
    ! panneau de l'etape 6).
    TT0 = MPI_WTIME()
    CALL OOC_CACHE_FLUSH_ALL(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK)
    T_EVICT_L = T_EVICT_L + (MPI_WTIME() - TT0)

    DEALLOCATE(LRECV, URECV)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) THEN
      T_FACTOR_L = MPI_WTIME() - T0
      CALL MPI_ALLREDUCE(T_FACTOR_L, T_FACTOR, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
    IF (PRESENT(T_CACHE_GET))   CALL MPI_ALLREDUCE(T_GET_L,     T_CACHE_GET,   1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_CACHE_EVICT)) CALL MPI_ALLREDUCE(T_EVICT_L,   T_CACHE_EVICT, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_COMPUTE))     CALL MPI_ALLREDUCE(T_COMPUTE_L, T_COMPUTE,     1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_MPICALL))     CALL MPI_ALLREDUCE(T_MPI_L,     T_MPICALL,     1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    IF (PRESENT(T_CACHE_GET_MIN))   CALL MPI_ALLREDUCE(T_GET_L,     T_CACHE_GET_MIN,   1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (PRESENT(T_CACHE_EVICT_MIN)) CALL MPI_ALLREDUCE(T_EVICT_L,   T_CACHE_EVICT_MIN, 1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (PRESENT(T_COMPUTE_MIN))     CALL MPI_ALLREDUCE(T_COMPUTE_L, T_COMPUTE_MIN,     1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (PRESENT(T_MPICALL_MIN))     CALL MPI_ALLREDUCE(T_MPI_L,     T_MPICALL_MIN,     1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_LOCAL
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC

  !-----------------------------------------------------------------------
  ! Variante OOC de LOCAL_ROW_OFFSET -- corps STRICTEMENT identique.
  !-----------------------------------------------------------------------
  FUNCTION LOCAL_ROW_OFFSET_OOC(A_BLR, GRID, IGLOBAL) RESULT(RS)
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(IN) :: A_BLR
    TYPE(GRID_T), INTENT(IN)               :: GRID
    INTEGER, INTENT(IN)                    :: IGLOBAL
    INTEGER :: RS, li, I
    RS = 1
    DO li = 1, (IGLOBAL-1-GRID%MYROW)/GRID%PR
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      RS = RS + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
    END DO
  END FUNCTION LOCAL_ROW_OFFSET_OOC

  !-----------------------------------------------------------------------
  ! Resolution, par paquets d'au plus NRHS_BATCH colonnes globales (RHS,
  ! deja assemble sur disque par MPI_ASSEMBLE_RHS_OOC). Pour chaque
  ! paquet : lit B depuis RHS%STORE_B (append-only, jamais recycle --
  ! chaque paquet n'est lu qu'une fois, jamais reutilise entre paquets, un
  ! cache n'apporterait rien ici), execute exactement l'algorithme de
  ! MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D (substitution avant puis arriere) avec
  ! WIDTH=largeur du paquet en guise de "NRHS", puis appelle
  ! CONSUME_X_BATCH avec la solution de ce paquet avant de passer au
  ! suivant -- la memoire de ce flux (B, X) est bornee par NRHS_BATCH,
  ! jamais par NRHS total.
  !
  ! DIAGBLK/OFFBLK, EUX, sont LUS VIA CACHE_A/CACHE_D (memes caches que la
  ! factorisation, reutilises TELS QUELS -- ils sont vides apres le flush
  ! de fin de factorisation, mais restent structurellement valides). Point
  ! important : contrairement a la factorisation, un bloc L/U donne est
  ! touche EXACTEMENT UNE FOIS par phase (avant OU arriere) PAR PAQUET --
  ! a l'interieur d'un seul paquet, un cache n'apporte donc rien (chaque
  ! acces est un premier acces). Le gain vient d'ENTRE deux paquets : les
  ! MEMES blocs L/U sont necessaires a CHAQUE paquet (ils ne dependent pas
  ! des colonnes de second membre en cours), donc un cache qui SURVIT au
  ! paquet en cours evite de relire NBATCHES fois les memes donnees --
  ! MESURE (N=20000, NB=400, blocs quasi tous denses -- 1% compresses --
  ! MAX_MB_PER_PROCESS=900, NRHS_BATCH=2 donc NBATCHES=3, voir test/bench_
  ! blr_mpi_pivot_ooc_compare.f90) : temps de resolution divise par 2.76
  ! (36.5s -> 13.2s), RSS quasi inchangee. Gain plus modeste que l'espoir
  ! initial (facteur NBATCHES=3 aurait suggere jusqu'a /3) -- coherent avec
  ! la mise en garde suivante : MAX_MB_PER_PROCESS=900 est proche mais pas
  ! tres au-dessus du besoin reel observe (~915-920 Mo), une partie du
  ! gain se perd donc encore en eviction. ATTENTION : ce gain suppose que
  ! MAX_MB_PER_PROCESS
  ! couvre l'ensemble complet des donnees L/U locales -- sinon, le meme
  ! phenomene de thrashing que la factorisation (acces quasi-sequentiel,
  ! l'integralite de l'ensemble de travail est parcourue entre deux
  ! passages sur un meme bloc) annule le benefice (voir le commentaire
  ! "DIMENSIONNER MAX_MB_PER_PROCESS" sur MPI_BLR_MATRIX_OOC_T).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH, SHOW_PROGRESS, &
                                                   PROGRESS_CALLBACK)
    TYPE(GRID_T), INTENT(IN)                  :: GRID
    TYPE(MPI_BLR_MATRIX_OOC_T), INTENT(INOUT) :: A_BLR
    TYPE(OOC_RHS_T), INTENT(INOUT)            :: RHS
    PROCEDURE(CONSUME_X_BATCH_IFACE)          :: CONSUME_X_BATCH
    LOGICAL, INTENT(IN), OPTIONAL             :: SHOW_PROGRESS
    ! Voir MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D pour le detail (meme
    ! convention) -- NSTEPS/STEP redemarrent a chaque PAQUET de second
    ! membre (RHS%NBATCHES), donc les paliers de 5% le font aussi : chaque
    ! paquet est traite comme sa propre etape du point de vue du callback.
    PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE), OPTIONAL :: PROGRESS_CALLBACK

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:)
    TYPE(BLR_BLOCK_T) :: DTMP
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, RS, RE, NROWTRAIL, C, IERR, LDB
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, IDX, IDX_D, PBAR_STATE
    INTEGER :: B, RSTART, REND, WIDTH, NRHS, NRHS_LOC, LOCAL_OFFSET, GSTART, GEND
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NBLK = A_BLR%NBLK

    DO B = 1, RHS%NBATCHES
      RSTART = (B-1)*RHS%NRHS_BATCH + 1
      REND = MIN(B*RHS%NRHS_BATCH, RHS%NRHS)
      WIDTH = REND - RSTART + 1
      NRHS = WIDTH
      CALL RHS_LOCAL_RANGE_MBM(WIDTH, GRID%PC, GRID%MYCOL, LOCAL_OFFSET, NRHS_LOC)
      GSTART = RSTART + LOCAL_OFFSET - 1
      GEND = GSTART + NRHS_LOC - 1

      CALL OOC_READ_BLOCK_RAW(RHS%STORE_B, RHS%BATCH(B), DTMP)
      LDB = MAX(RHS%MLOC,1)
      ALLOCATE(B_LOCAL(LDB,MAX(NRHS_LOC,1)))
      CALL MOVE_ALLOC(DTMP%DENSE, B_LOCAL)

      STEP = 0
      PBAR_STATE = -1
      NSTEPS = 2 * NBLK
      ALLOCATE(RSTARTS(0:GRID%PC-1), RLOCS(0:GRID%PC-1), RCOUNTS(0:GRID%PC-1), RDISPLS(0:GRID%PC-1))
      DO C = 0, GRID%PC-1
        CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, C, RSTARTS(C), RLOCS(C))
      END DO

      ! === 0) P*B -- inchange (IPIV deja diffuse a la factorisation) ===
      DO K = 1, NBLK
        KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
        PROW_K = MOD(K-1, GRID%PR)
        IF (GRID%MYROW.NE.PROW_K .OR. NRHS_LOC.LE.0) CYCLE
        RS = LOCAL_ROW_OFFSET_OOC(A_BLR, GRID, K)
        DO I = 1, KW
          IF (A_BLR%IPIV(I,K) .NE. I) THEN
            CALL CSWAP(NRHS_LOC, B_LOCAL(RS+I-1,1), LDB, B_LOCAL(RS+A_BLR%IPIV(I,K)-1,1), LDB)
          END IF
        END DO
      END DO

      ! === 1) L*Y = P*B (substitution avant) -- DIAGBLK/OFFBLK LUS via
      !     CACHE_D/CACHE_A (voir commentaire de tete : cache PERSISTANT
      !     entre paquets, en LECTURE SEULE). ===
      DO K = 1, NBLK
        KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
        PROW_K = MOD(K-1, GRID%PR)
        PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

        ALLOCATE(XK_FULL(KW,NRHS))
        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(DKK(KW,KW))
          IF (GRID%MYCOL .EQ. PCOL_K) THEN
            IDX_D = OOC_DIAG_CACHE_GET(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK, K)
            DKK = A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%DENSE
          END IF
          CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

          RS = LOCAL_ROW_OFFSET_OOC(A_BLR, GRID, K)
          IF (NRHS_LOC .GT. 0) THEN
            CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS_LOC, &
                        CONE, DKK, KW, B_LOCAL(RS,1), LDB)
          END IF
          DEALLOCATE(DKK)
          ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
          IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
          END DO
          CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                               XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
          DEALLOCATE(XK_LOC)
        END IF
        CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        NROWTRAIL = 0
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .GT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
        END DO
        IF (NROWTRAIL .GT. 0) THEN
          ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
          IF (GRID%MYCOL .EQ. PCOL_K) THEN
            CONTRIB = CZERO
            RE = 0
            DO li = 1, A_BLR%NBLKR_LOC
              I = (li-1)*GRID%PR + GRID%MYROW + 1
              IF (I .LE. K) CYCLE
              MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
              IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, lj)
              CALL BLR_BLOCK_APPLY(A_BLR%CACHE_A%SLOTS(IDX)%BLK, NRHS, CONE, XK_FULL, KW, CONE, &
                                    CONTRIB(RE+1,1), NROWTRAIL)
              RE = RE + MI
            END DO
          END IF
          ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
          END DO
          CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                             CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
          IF (NRHS_LOC .GT. 0) THEN
            RE = 0; RS = 1
            DO li = 1, A_BLR%NBLKR_LOC
              I = (li-1)*GRID%PR + GRID%MYROW + 1
              MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
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
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('4) Resolution ooc (pivot)', STEP, NSTEPS, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
        END IF
      END DO

      ! === 2) U*X = Y (substitution arriere) -- AUCUN pivot necessaire ===
      DO K = NBLK, 1, -1
        KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
        PROW_K = MOD(K-1, GRID%PR)
        PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)

        ALLOCATE(XK_FULL(KW,NRHS))
        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(DKK(KW,KW))
          IF (GRID%MYCOL .EQ. PCOL_K) THEN
            IDX_D = OOC_DIAG_CACHE_GET(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK, K)
            DKK = A_BLR%CACHE_D%SLOTS(IDX_D)%BLK%DENSE
          END IF
          CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)

          RS = LOCAL_ROW_OFFSET_OOC(A_BLR, GRID, K)
          IF (NRHS_LOC .GT. 0) THEN
            CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS_LOC, &
                        CONE, DKK, KW, B_LOCAL(RS,1), LDB)
          END IF
          DEALLOCATE(DKK)
          ALLOCATE(XK_LOC(KW,MAX(NRHS_LOC,1)))
          IF (NRHS_LOC .GT. 0) XK_LOC = B_LOCAL(RS:RS+KW-1,1:NRHS_LOC)
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
          END DO
          CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                               XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
          DEALLOCATE(XK_LOC)
        END IF
        CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

        lj = BLR_COL_LOCAL(A_BLR%CMAP, K, GRID)
        NROWTRAIL = 0
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LT. K) NROWTRAIL = NROWTRAIL + A_BLR%BEND(I)-A_BLR%BSTART(I)+1
        END DO
        IF (NROWTRAIL .GT. 0) THEN
          ALLOCATE(CONTRIB(NROWTRAIL,NRHS))
          IF (GRID%MYCOL .EQ. PCOL_K) THEN
            CONTRIB = CZERO
            RE = 0
            DO li = 1, A_BLR%NBLKR_LOC
              I = (li-1)*GRID%PR + GRID%MYROW + 1
              IF (I .GE. K) CYCLE
              MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
              IDX = OOC_CACHE_GET(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK, li, lj)
              CALL BLR_BLOCK_APPLY(A_BLR%CACHE_A%SLOTS(IDX)%BLK, NRHS, CONE, XK_FULL, KW, CONE, &
                                    CONTRIB(RE+1,1), NROWTRAIL)
              RE = RE + MI
            END DO
          END IF
          ALLOCATE(CONTRIB_LOC(NROWTRAIL,MAX(NRHS_LOC,1)))
          DO C = 0, GRID%PC-1
            RCOUNTS(C) = NROWTRAIL*RLOCS(C); RDISPLS(C) = NROWTRAIL*(RSTARTS(C)-1)
          END DO
          CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                             CONTRIB_LOC, NROWTRAIL*NRHS_LOC, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
          IF (NRHS_LOC .GT. 0) THEN
            RE = 0; RS = 1
            DO li = 1, A_BLR%NBLKR_LOC
              I = (li-1)*GRID%PR + GRID%MYROW + 1
              MI = A_BLR%BEND(I)-A_BLR%BSTART(I)+1
              IF (I .LT. K) THEN
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
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('4) Resolution ooc (pivot)', STEP, NSTEPS, PBAR_STATE)
        IF (PRESENT(PROGRESS_CALLBACK)) THEN
          IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
        END IF
      END DO

      DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
      CALL CONSUME_X_BATCH(RSTART, REND, GSTART, GEND, B_LOCAL, LDB, NRHS_LOC)
      DEALLOCATE(B_LOCAL)
    END DO

    ! Libere immediatement ce que le cache de resolution a accumule (voir
    ! commentaire de tete) -- jamais DIRTY (lecture seule), un flush ne
    ! reecrit donc rien sur disque, il se contente de liberer la RAM sans
    ! attendre l'eventuel MPI_BLR_MATRIX_OOC_FREE (qui peut arriver bien
    ! plus tard, voire jamais si l'appelant reutilise A_BLR autrement).
    CALL OOC_CACHE_FLUSH_ALL(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK)
  END SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D

END MODULE MPI_BLR_MATRIX_MOD
