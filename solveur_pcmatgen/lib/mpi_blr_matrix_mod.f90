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
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK, BCAST_BLR_PANEL, BCAST_NCALL, BCAST_ELEMS
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE, DENSE_PANEL_T
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_MILESTONE, PRINT_PROGRESS_BAR
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB, BLR_TRIM_HEAP
  IMPLICIT NONE

  ! Profil memoire de la factorisation (diagnostic, FACT_MEM_PROF) : a la
  ! fin de chaque etape K, MAX sur la grille (Mo) du pic de RSS de l'etape
  ! (VmHWM remis a zero au debut de l'etape), de la RSS, de la memoire
  ! allouee (mallinfo), de la matrice/facteurs, des panneaux recus (WRECV,
  ! LTRECV ; LRECV, URECV en LU), de la copie des L(I,K) chez le proprietaire de la colonne K
  ! et des espaces de travail (WS, DWORK). Affiche par le rang 0.
  LOGICAL :: FACT_MEM_PROF = .FALSE.

  ! Trace des attentes MPI de la factorisation (diagnostic, FACT_STEP_TRACE) :
  ! par processus, temps de communication ventile (1 debut d'etape : pivots
  ! et bloc diagonal, 2 panneaux le long des lignes comme emetteur, 3
  ! panneaux le long des colonnes, 4 anticipation, 5 panneaux le long des
  ! lignes comme recepteur), et par etape le Schur et les communications
  ! de chaque processus. Affiche par le rang 0.
  LOGICAL :: FACT_STEP_TRACE = .FALSE.


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
    INCLUDE 'mpif.h'
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
  ! Libere la memoire locale de A_BLR (blocs diagonaux et hors-diagonale,
  ! puis decoupage et pivots).
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

  ! Profil memoire : remet a zero le pic de RSS (VmHWM) du processus.
  SUBROUTINE MEMPROF_RESET_HWM()
    INTEGER :: U, IOS
    OPEN(NEWUNIT=U, FILE='/proc/self/clear_refs', STATUS='OLD', ACTION='WRITE', IOSTAT=IOS)
    IF (IOS .NE. 0) RETURN
    WRITE(U, '(A)', IOSTAT=IOS) '5'
    CLOSE(U)
  END SUBROUTINE MEMPROF_RESET_HWM

  ! Profil memoire : pic de RSS (VmHWM, Ko ; -1 si illisible).
  FUNCTION MEMPROF_HWM_KB() RESULT(HWM)
    INTEGER(8) :: HWM
    INTEGER :: U, IOS
    CHARACTER(LEN=256) :: LINE
    HWM = -1_8
    OPEN(NEWUNIT=U, FILE='/proc/self/status', STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) RETURN
    DO
      READ(U, '(A)', IOSTAT=IOS) LINE
      IF (IOS .NE. 0) EXIT
      IF (LINE(1:6) .EQ. 'VmHWM:') THEN
        READ(LINE(7:), *, IOSTAT=IOS) HWM
        EXIT
      END IF
    END DO
    CLOSE(U)
  END FUNCTION MEMPROF_HWM_KB

  ! Profil memoire : octets alloues par malloc et vivants (glibc mallinfo,
  ! champs int lus modulo 2^32 ; uordblks + hblkhd).
  FUNCTION MEMPROF_ALLOC_BYTES() RESULT(NB)
    USE ISO_C_BINDING, ONLY: C_INT
    INTEGER(8) :: NB
    TYPE, BIND(C) :: MALLINFO_T
      INTEGER(C_INT) :: ARENA, ORDBLKS, SMBLKS, HBLKS, HBLKHD, USMBLKS, FSMBLKS, UORDBLKS, FORDBLKS, KEEPCOST
    END TYPE MALLINFO_T
    INTERFACE
      FUNCTION MALLINFO() BIND(C, NAME='mallinfo')
        IMPORT :: MALLINFO_T
        TYPE(MALLINFO_T) :: MALLINFO
      END FUNCTION MALLINFO
    END INTERFACE
    TYPE(MALLINFO_T) :: MI
    MI = MALLINFO()
    NB = MODULO(INT(MI%UORDBLKS,8), 4294967296_8) + MODULO(INT(MI%HBLKHD,8), 4294967296_8)
  END FUNCTION MEMPROF_ALLOC_BYTES

  !-----------------------------------------------------------------------
  ! Trace d'une factorisation (FACT_STEP_TRACE ou BLR_TRACE_FACT=1), commune
  ! au LU et au LDLt. STEPT(1:6,K) par processus et par etape : temps de
  ! Schur, temps de communication, duree, travail estime du Schur, appels
  ! et elements des diffusions de blocs. TW(1:5) : communications par
  ! categorie. Rassemble sur le rang 0, qui affiche les lignes [TRACE].
  !-----------------------------------------------------------------------
  SUBROUTINE FACT_TRACE_REPORT(GRID, N, NBLK, STEPT, TW)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)     :: GRID
    INTEGER, INTENT(IN)          :: N, NBLK
    DOUBLE PRECISION, INTENT(IN) :: STEPT(6,NBLK), TW(5)
    DOUBLE PRECISION, ALLOCATABLE :: STEPG(:,:,:), TWG(:,:)
    DOUBLE PRECISION :: WSYNC, WMAXT, WMOY, WDIAG
    INTEGER :: NPG, IPG, K, IERR
      NPG = GRID%PR*GRID%PC
      ALLOCATE(STEPG(6,NBLK,MERGE(NPG,1,GRID%MYID.EQ.0)), TWG(5,MERGE(NPG,1,GRID%MYID.EQ.0)))
      CALL MPI_GATHER(STEPT, 6*NBLK, MPI_DOUBLE_PRECISION, STEPG, 6*NBLK, MPI_DOUBLE_PRECISION, 0, GRID%COMM, IERR)
      CALL MPI_GATHER(TW, 5, MPI_DOUBLE_PRECISION, TWG, 5, MPI_DOUBLE_PRECISION, 0, GRID%COMM, IERR)
      IF (GRID%MYID .EQ. 0) THEN
        ! Resume independant des temps : travail estime du Schur (WORK_COUNT)
        ! et diffusions de blocs par processus.
        WSYNC = 0.0D0; WMAXT = 0.0D0; WMOY = 0.0D0
        DO K = 1, NBLK
          WSYNC = WSYNC + MAXVAL(STEPG(4,K,:))
        END DO
        DO IPG = 1, NPG
          WMAXT = MAX(WMAXT, SUM(STEPG(4,:,IPG)))
          WMOY = WMOY + SUM(STEPG(4,:,IPG))/DBLE(NPG)
        END DO
        WRITE(6,'(A,I0,A,I0,A,I0,A,I0)') ' [TRACE] grille ', GRID%PR, ' x ', GRID%PC, ', N = ', N, ', blocs = ', NBLK
        WRITE(6,'(A,ES10.3,A,ES10.3,A,ES10.3)') ' [TRACE] travail du Schur par processus : moyenne ', WMOY, &
             ', maximum ', WMAXT, ', somme des maximums par etape ', WSYNC
        WRITE(6,'(A,F6.1,A,F6.1,A,F6.1,A)') ' [TRACE] efficacite : desequilibre global ', 100.0D0*(WMAXT/MAX(WMOY,1.0D-30) &
             - 1.0D0), ' %, synchronisation par etape ', 100.0D0*(WSYNC/MAX(WMAXT,1.0D-30) - 1.0D0), &
             ' %, borne d''acceleration ', DBLE(NPG)*WMOY/MAX(WSYNC,1.0D-30), ' (ideal = nombre de processus)'
        WRITE(6,'(A,F10.0,A,F10.0,A,ES10.3,A)') ' [TRACE] diffusions de blocs par processus : moyenne ', &
             SUM(STEPG(5,:,:))/DBLE(NPG), ' appels, maximum ', MAXVAL(SUM(STEPG(5,:,:), DIM=1)), ' ; ', &
             SUM(STEPG(6,:,:))/DBLE(NPG), ' elements en moyenne'
        WDIAG = WMOY/MAX(SUM(STEPG(5,:,:))/DBLE(NPG), 1.0D0)
        WRITE(6,'(A,ES10.3)') ' [TRACE] travail du Schur par appel de diffusion (moyenne) : ', WDIAG
        IF (NPG .LE. 8) THEN
          WRITE(6,'(A)') ' [TRACE] communications par processus (s) : debut d''etape, lignes (emetteur), colonnes, '// &
               'anticipation, lignes (recepteur) ; total Schur ; total communications'
          DO IPG = 1, NPG
            WRITE(6,'(A,I3,A,I0,A,I0,A,5F8.2,A,F8.2,A,F8.2)') ' [TRACE] proc', IPG-1, ' (', (IPG-1)/GRID%PC, ',', &
                 MOD(IPG-1,GRID%PC), ')', TWG(:,IPG), ' ;', SUM(STEPG(1,:,IPG)), ' ;', SUM(STEPG(2,:,IPG))
          END DO
          WRITE(6,'(A)') ' [TRACE]   K  duree | Schur par processus | communications par processus'
          DO K = 1, NBLK
            WRITE(6,'(A,I4,F7.2,A,*(F7.2))') ' [TRACE] ', K, MAXVAL(STEPG(3,K,:)), ' |', STEPG(1,K,:), STEPG(2,K,:)
          END DO
        END IF
      END IF
      DEALLOCATE(STEPG, TWG)
  END SUBROUTINE FACT_TRACE_REPORT

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
                                             SHOW_PROGRESS, SHOW_RSS_TRACE, USE_MOVE_ALLOC, PROGRESS_CALLBACK, &
                                             T_BREAKDOWN, LOOKAHEAD)
    INCLUDE 'mpif.h'
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
    ! T_BREAKDOWN (optionnel, MAX sur la grille, s) : (1) factorisation des
    ! blocs diagonaux, (2) echanges de lignes/colonnes du pivot, (3) operations
    ! de panneau (TRSM des panneaux L et U), (4) communications MPI (attentes
    ! comprises), (5) Schur des blocs diagonaux, (6) Schur hors-diagonale,
    ! (7) boucle complete. Chaque categorie est un MAX independant : leur
    ! somme peut differer de (7).
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_BREAKDOWN(7)
    ! LOOKAHEAD (defaut .TRUE. ; BLR_LU_LOOKAHEAD=0 le coupe) : anticipation
    ! de profondeur 1, comme le LDLt. A l'etape K, des que la colonne K+1 a
    ! recu sa mise a jour de Schur, le proprietaire du bloc (K+1,K+1) le
    ! factorise et lance des diffusions NON BLOQUANTES (MPI_IBCAST) de
    ! INFO/IPIV a toute la grille et du bloc a sa colonne de grille, qui
    ! traite aussitot son panneau L(I,K+1). Le reste du Schur de l'etape K
    ! recouvre ces diffusions ; l'etape K+1 commence par un MPI_WAITALL.
    ! Resultat numerique identique (memes operations, seul l'ordre change).
    LOGICAL, INTENT(IN), OPTIONAL     :: LOOKAHEAD

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    ! Liste des blocs d'un panneau a diffuser (BCAST_BLR_PANEL).
    INTEGER, ALLOCATABLE :: PIX(:), PRT(:), PMS(:), PNS(:)
    INTEGER :: NPB
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, NS_LOCAL, li_K, PBAR_STATE
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0, T_FACTOR_L, RSS_MAX_MB, TB(7), TT0, TLOOP0
    LOGICAL :: PROGRESS, RSS_TRACE
    ! tampons des diffusions anticipees (en vol entre deux etapes)
    INTEGER, ALLOCATABLE, ASYNCHRONOUS :: IBUF_N(:)
    COMPLEX(WP), ALLOCATABLE, ASYNCHRONOUS :: DKN_N(:)
    INTEGER :: REQ(2), NREQ
    LOGICAL :: LA, PENDING, LADONE
    ! trace, compteurs, profil memoire, releve du travail (voir le LDLt)
    DOUBLE PRECISION :: TW(5), TDW, TB0(7), TSTEP0, W0, WUPD
    DOUBLE PRECISION, ALLOCATABLE :: STEPT(:,:)
    INTEGER(8) :: NC0, NE0
    CHARACTER(LEN=8) :: ENVT
    LOGICAL :: TRACE, DUMPW
    INTEGER :: UDW
    CHARACTER(LEN=64) :: FDW

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
    NPB = MAX(A_BLR%NBLKR_LOC, A_BLR%NBLKC_LOC, 1)
    ALLOCATE(PIX(NPB), PRT(NPB), PMS(NPB), PNS(NPB))
    ALLOCATE(A_BLR%IPIV(A_BLR%NB, NBLK))
    A_BLR%IPIV = 0
    CALL BLR_WORKSPACE_ALLOC(WS, A_BLR%NB, FACT_WORKSPACE_RANK(A_BLR%NB))
    LA = .TRUE.
    IF (PRESENT(LOOKAHEAD)) LA = LOOKAHEAD
    CALL GET_ENVIRONMENT_VARIABLE('BLR_LU_LOOKAHEAD', ENVT)
    IF (TRIM(ENVT) .EQ. '0') LA = .FALSE.
    PENDING = .FALSE.
    NREQ = 0
    IF (LA) ALLOCATE(IBUF_N(0:A_BLR%NB), DKN_N(A_BLR%NB*A_BLR%NB))
    TB = 0.0D0
    TW = 0.0D0
    CALL GET_ENVIRONMENT_VARIABLE('BLR_TRACE_FACT', ENVT)
    TRACE = FACT_STEP_TRACE .OR. (TRIM(ENVT) .EQ. '1')
    IF (TRACE) WORK_COUNT = .TRUE.
    ! BLR_DUMP_WORK=1 : releve du travail de chaque mise a jour de Schur
    ! (K, I, J, travail), entree de tools/blr_sim/simulate.py.
    CALL GET_ENVIRONMENT_VARIABLE('BLR_DUMP_WORK', ENVT)
    DUMPW = (TRIM(ENVT) .EQ. '1')
    IF (DUMPW) THEN
      WORK_COUNT = .TRUE.
      WRITE(FDW,'(A,I0,A)') 'blr_work_rank', GRID%MYID, '.bin'
      OPEN(NEWUNIT=UDW, FILE=TRIM(FDW), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
    END IF
    ALLOCATE(STEPT(6,NBLK))
    STEPT = 0.0D0
    IF (FACT_MEM_PROF .AND. GRID%MYID .EQ. 0) WRITE(6,'(A)') ' [MEMFACT]  K   pic_etape   RSS  alloue  matrice  '// &
         'L_recus  U_recus  WS  (Mo, max par processus)'
    TLOOP0 = MPI_WTIME()

    DO K = 1, NBLK
      IF (FACT_MEM_PROF) CALL MEMPROF_RESET_HWM()
      TB0 = TB; TSTEP0 = MPI_WTIME(); W0 = SCHUR_WORK; NC0 = BCAST_NCALL; NE0 = BCAST_ELEMS
      LADONE = LA .AND. K .GT. 1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = BLR_COL_OWNER(A_BLR%CMAP, K, GRID)
      INFO_STEP = 0
      ALLOCATE(IPIV_K(KW))
      IPIV_K = 0

      ! 1) facteur diagonal AVEC PIVOT (CGETRF), local au proprietaire -- ou,
      !    en anticipation, deja fait a l'etape K-1 : on attend alors
      !    seulement la fin des diffusions non bloquantes.
      IF (PENDING) THEN
        TT0 = MPI_WTIME()
        CALL MPI_WAITALL(NREQ, REQ, MPI_STATUSES_IGNORE, IERR)
        TDW = MPI_WTIME() - TT0; TB(4) = TB(4) + TDW; TW(1) = TW(1) + TDW
        PENDING = .FALSE.
        INFO_STEP = IBUF_N(0)
        IPIV_K = IBUF_N(1:KW)
      ELSE
        IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
          TT0 = MPI_WTIME()
          CALL CGETRF(KW, KW, A_BLR%DIAGBLK(K)%D, KW, IPIV_K, IINFO)
          IF (IINFO .NE. 0) INFO_STEP = A_BLR%BSTART(K) - 1 + IINFO
          TB(1) = TB(1) + (MPI_WTIME() - TT0)
        END IF
        TT0 = MPI_WTIME()
        CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
        ! Diffusion GLOBALE (voir commentaire de tete) -- petite (KW entiers)
        IF (INFO_STEP .EQ. 0) CALL MPI_BCAST(IPIV_K, KW, MPI_INTEGER, PROW_K*GRID%PC+PCOL_K, GRID%COMM, IERR)
        TDW = MPI_WTIME() - TT0; TB(4) = TB(4) + TDW; TW(1) = TW(1) + TDW
      END IF
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        DEALLOCATE(LRECV, URECV, IPIV_K, STEPT)
        IF (DUMPW) CLOSE(UDW)
        CALL BLR_WORKSPACE_FREE(WS)
        RETURN
      END IF
      A_BLR%IPIV(1:KW,K) = IPIV_K

      ! 2) propage le pivot aux AUTRES blocs de la MEME ligne-bloc K, deja
      !    stockes localement (J<K finalises, J>K encore bruts) -- SEUL le
      !    proprietaire de la ligne PROW_K en possede.
      TT0 = MPI_WTIME()
      IF (GRID%MYROW .EQ. PROW_K) THEN
        li_K = (K-1)/GRID%PR + 1
        DO lj = 1, A_BLR%NBLKC_LOC
          J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
          IF (J .EQ. K) CYCLE
          CALL BLR_BLOCK_ROW_SWAP(A_BLR%OFFBLK(li_K,lj), KW, IPIV_K)
        END DO
      END IF
      TB(2) = TB(2) + (MPI_WTIME() - TT0)
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
      TT0 = MPI_WTIME()
      IF (GRID%MYROW .EQ. PROW_K) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      ! (colonne PCOL_K : deja recu et utilise a l'etape K-1 en anticipation)
      IF (GRID%MYCOL .EQ. PCOL_K .AND. .NOT. LADONE) CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      TDW = MPI_WTIME() - TT0; TB(4) = TB(4) + TDW; TW(1) = TW(1) + TDW
      TT0 = MPI_WTIME()

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
      !    Deja fait a l'etape K-1 en anticipation (POST_NEXT_DIAG).
      IF (GRID%MYCOL .EQ. PCOL_K .AND. .NOT. LADONE) THEN
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
      TB(3) = TB(3) + (MPI_WTIME() - TT0)
      TT0 = MPI_WTIME()

      ! 6) diffusion des panneaux L/U DEJA COMPRESSES : un envoi groupe par
      !    panneau (BCAST_BLR_PANEL) ; la racine copie d'abord ses blocs
      !    dans LRECV/URECV et diffuse depuis cette copie.
      NPB = 0
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        IF (I .LE. K) CYCLE
        NPB = NPB + 1
        PIX(NPB) = li; PRT(NPB) = PCOL_K
        PMS(NPB) = A_BLR%BEND(I)-A_BLR%BSTART(I)+1; PNS(NPB) = KW
        IF (GRID%MYCOL .EQ. PCOL_K) LRECV(li) = A_BLR%OFFBLK(li,BLR_COL_LOCAL(A_BLR%CMAP, K, GRID))
      END DO
      CALL BCAST_BLR_PANEL(GRID%ROW_COMM, NPB, PIX, PRT, PMS, PNS, LRECV)
      TDW = MPI_WTIME() - TT0; TB(4) = TB(4) + TDW
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        TW(2) = TW(2) + TDW
      ELSE
        TW(5) = TW(5) + TDW
      END IF
      TT0 = MPI_WTIME()

      NPB = 0
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NPB = NPB + 1
        PIX(NPB) = lj; PRT(NPB) = PROW_K
        PMS(NPB) = KW; PNS(NPB) = A_BLR%BEND(J)-A_BLR%BSTART(J)+1
        IF (GRID%MYROW .EQ. PROW_K) URECV(lj) = A_BLR%OFFBLK((K-1)/GRID%PR+1,lj)
      END DO
      CALL BCAST_BLR_PANEL(GRID%COL_COMM, NPB, PIX, PRT, PMS, PNS, URECV)
      TDW = MPI_WTIME() - TT0; TB(4) = TB(4) + TDW; TW(3) = TW(3) + TDW

      ! 7) mise a jour de Schur locale. En anticipation, des que la colonne
      !    K+1 est a jour, le bloc (K+1,K+1) est factorise et diffuse sans
      !    bloquer ; les autres colonnes de grille postent leur reception
      !    avant leur Schur.
      IF (LA .AND. K .LT. NBLK) THEN
        IF (GRID%MYCOL .NE. BLR_COL_OWNER(A_BLR%CMAP, K+1, GRID)) CALL POST_NEXT_DIAG(K+1)
      END IF
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (J .LE. K) CYCLE
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        DO li = 1, A_BLR%NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          IF (I .LE. K) CYCLE
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          TT0 = MPI_WTIME()
          IF (I .EQ. J) THEN
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LRECV(li), URECV(lj), A_BLR%DIAGBLK(I)%D, MI)
            ! travail estime du bloc diagonal : MI x NJ x rang du panneau (dense : KW)
            IF (WORK_COUNT) THEN
              WUPD = DBLE(MI)*DBLE(NJ)*DBLE(MERGE(LRECV(li)%RANK, KW, LRECV(li)%IS_LOWRANK))
              SCHUR_WORK = SCHUR_WORK + WUPD
              IF (DUMPW) WRITE(UDW) K, I, J, WUPD
            END IF
            TB(5) = TB(5) + (MPI_WTIME() - TT0)
          ELSE
            WUPD = SCHUR_WORK
            CALL SCHUR_UPDATE_BLR_BLOCK(A_BLR%OFFBLK(li,lj), MI, NJ, KW, LRECV(li), URECV(lj), TOL, WS, &
                                         USE_MOVE_ALLOC=USE_MOVE_ALLOC)
            IF (DUMPW) WRITE(UDW) K, I, J, SCHUR_WORK - WUPD
            TB(6) = TB(6) + (MPI_WTIME() - TT0)
          END IF
        END DO
        IF (LA .AND. J .EQ. K+1) CALL POST_NEXT_DIAG(K+1)
      END DO

      STEPT(1,K) = (TB(5) + TB(6)) - (TB0(5) + TB0(6))
      STEPT(2,K) = TB(4) - TB0(4)
      STEPT(3,K) = MPI_WTIME() - TSTEP0
      STEPT(4,K) = SCHUR_WORK - W0
      STEPT(5,K) = DBLE(BCAST_NCALL - NC0)
      STEPT(6,K) = DBLE(BCAST_ELEMS - NE0)
      IF (FACT_MEM_PROF) CALL MEMPROF_STEP()

      ! Panneaux recus des blocs-lignes / blocs-colonnes deja traites (I, J
      ! <= K+1) : plus jamais lus, liberes a chaque etape (voir le LDLt).
      DO li = 1, A_BLR%NBLKR_LOC
        IF ((li-1)*GRID%PR + GRID%MYROW + 1 .LE. K+1) CALL BLR_BLOCK_RELEASE(LRECV(li))
      END DO
      DO lj = 1, A_BLR%NBLKC_LOC
        IF (BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID) .LE. K+1) CALL BLR_BLOCK_RELEASE(URECV(lj))
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
    TB(7) = MPI_WTIME() - TLOOP0
    IF (PRESENT(T_BREAKDOWN)) CALL MPI_ALLREDUCE(TB, T_BREAKDOWN, 7, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)

    IF (TRACE) CALL FACT_TRACE_REPORT(GRID, A_BLR%N, NBLK, STEPT, TW)
    IF (DUMPW) CLOSE(UDW)
    DEALLOCATE(LRECV, URECV, STEPT)
    IF (ALLOCATED(IBUF_N)) DEALLOCATE(IBUF_N, DKN_N)
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

  CONTAINS

    ! Profil memoire de l'etape K (FACT_MEM_PROF), voir l'en-tete du module.
    SUBROUTINE MEMPROF_STEP()
      INTEGER(8) :: V(8), VM(8), BL, BU
      INTEGER :: LIP, LJP, IERRP
      BL = 0_8; BU = 0_8
      DO LIP = 1, A_BLR%NBLKR_LOC
        IF ((LIP-1)*GRID%PR + GRID%MYROW + 1 .GT. K) BL = BL + BLR_BLOCK_BYTES(LRECV(LIP))
      END DO
      DO LJP = 1, A_BLR%NBLKC_LOC
        IF (BLR_COL_GLOBAL(A_BLR%CMAP, LJP, GRID) .GT. K) BU = BU + BLR_BLOCK_BYTES(URECV(LJP))
      END DO
      V(1) = INT(K,8)
      V(2) = MEMPROF_HWM_KB()*1024_8
      V(3) = GET_RSS_KB()*1024_8
      V(4) = MEMPROF_ALLOC_BYTES()
      V(5) = MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR)
      V(6) = BL; V(7) = BU
      V(8) = BLR_WORKSPACE_BYTES(WS)
      CALL MPI_REDUCE(V, VM, 8, MPI_INTEGER8, MPI_MAX, 0, GRID%COMM, IERRP)
      IF (GRID%MYID .EQ. 0) WRITE(6,'(A,I4,7F9.1)') ' [MEMFACT] ', INT(VM(1)), REAL(VM(2:8),8)/1.0D6
    END SUBROUTINE MEMPROF_STEP

    ! Anticipation : factorise (proprietaire) le bloc KN deja a jour et
    ! lance les diffusions non bloquantes de (INFO, IPIV) sur la grille et
    ! du bloc factorise sur la colonne de grille PCOL_KN. Cette colonne
    ! attend ensuite le bloc et traite tout de suite son panneau L(I,KN)
    ! (TRSM), qui sort ainsi du chemin critique. Appelee une fois par etape
    ! par CHAQUE processus, dans le meme ordre relatif aux autres
    ! collectives de chaque communicateur.
    SUBROUTINE POST_NEXT_DIAG(KN)
      INTEGER, INTENT(IN) :: KN
      INTEGER :: KWN, PRN, PCN, IERRN, LIN, IN, LJN, IINFON
      COMPLEX(WP), ALLOCATABLE :: DKN(:,:)
      DOUBLE PRECISION :: TF0
      KWN = A_BLR%BEND(KN) - A_BLR%BSTART(KN) + 1
      PRN = MOD(KN-1, GRID%PR)
      PCN = BLR_COL_OWNER(A_BLR%CMAP, KN, GRID)
      IF (GRID%MYROW.EQ.PRN .AND. GRID%MYCOL.EQ.PCN) THEN
        TF0 = MPI_WTIME()
        IBUF_N(0) = 0
        CALL CGETRF(KWN, KWN, A_BLR%DIAGBLK(KN)%D, KWN, IBUF_N(1:KWN), IINFON)
        IF (IINFON .NE. 0) IBUF_N(0) = A_BLR%BSTART(KN) - 1 + IINFON
        DKN_N(1:KWN*KWN) = RESHAPE(A_BLR%DIAGBLK(KN)%D, [KWN*KWN])
        TB(1) = TB(1) + (MPI_WTIME() - TF0)
      END IF
      CALL MPI_IBCAST(IBUF_N, KWN+1, MPI_INTEGER, PRN*GRID%PC + PCN, GRID%COMM, REQ(1), IERRN)
      NREQ = 1
      IF (GRID%MYCOL .EQ. PCN) THEN
        CALL MPI_IBCAST(DKN_N, KWN*KWN, MPI_COMPLEX, PRN, GRID%COL_COMM, REQ(2), IERRN)
        NREQ = 2
        TF0 = MPI_WTIME()
        CALL MPI_WAITALL(NREQ, REQ, MPI_STATUSES_IGNORE, IERRN)
        TDW = MPI_WTIME() - TF0; TB(4) = TB(4) + TDW; TW(4) = TW(4) + TDW
        IF (IBUF_N(0) .EQ. 0) THEN
          TF0 = MPI_WTIME()
          ALLOCATE(DKN(KWN,KWN))
          DKN = RESHAPE(DKN_N(1:KWN*KWN), [KWN,KWN])
          LJN = BLR_COL_LOCAL(A_BLR%CMAP, KN, GRID)
          DO LIN = 1, A_BLR%NBLKR_LOC
            IN = (LIN-1)*GRID%PR + GRID%MYROW + 1
            IF (IN .LE. KN) CYCLE
            CALL BLR_TRSM_RIGHT_UPPER_NONUNIT(A_BLR%OFFBLK(LIN,LJN), DKN, KWN)
            NO_LOCAL = NO_LOCAL + 1
            IF (A_BLR%OFFBLK(LIN,LJN)%IS_LOWRANK) THEN
              NC_LOCAL = NC_LOCAL + 1
              RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(LIN,LJN)%RANK
            ELSE
              RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(LIN,LJN)%M, A_BLR%OFFBLK(LIN,LJN)%N)
            END IF
          END DO
          DEALLOCATE(DKN)
          TB(3) = TB(3) + (MPI_WTIME() - TF0)
        END IF
      END IF
      PENDING = .TRUE.
    END SUBROUTINE POST_NEXT_DIAG
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT

  !-----------------------------------------------------------------------
  ! Largeur effective des facteurs LU : stockage des blocs hors diagonale
  ! (rang x (M+N) en basse-rang, M x N en dense) rapporte au nombre de
  ! lignes de ces blocs. C'est le volume diffuse par ligne de matrice et
  ! par etape en resolution distribuee (DESCENT_MODE = 2), a comparer a
  ! NRHS, volume du mode 1. Valeur identique sur tous les processus.
  !-----------------------------------------------------------------------
  REAL(WP) FUNCTION LU_EFFECTIVE_WIDTH(GRID, A_BLR)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)           :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER :: li, lj, I, J, IERR
    DOUBLE PRECISION :: LOC(2), GLOB(2)
    LOC = 0.0D0
    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (I .EQ. J) CYCLE
        IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
          LOC(1) = LOC(1) + DBLE(A_BLR%OFFBLK(li,lj)%RANK)*DBLE(A_BLR%OFFBLK(li,lj)%M + A_BLR%OFFBLK(li,lj)%N)
        ELSE
          LOC(1) = LOC(1) + DBLE(A_BLR%OFFBLK(li,lj)%M)*DBLE(A_BLR%OFFBLK(li,lj)%N)
        END IF
        LOC(2) = LOC(2) + DBLE(A_BLR%OFFBLK(li,lj)%M)
      END DO
    END DO
    CALL MPI_ALLREDUCE(LOC, GLOB, 2, MPI_DOUBLE_PRECISION, MPI_SUM, GRID%COMM, IERR)
    LU_EFFECTIVE_WIDTH = REAL(GLOB(1)/MAX(GLOB(2),1.0D0), WP)
  END FUNCTION LU_EFFECTIVE_WIDTH

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
                                               PROGRESS_CALLBACK, DESCENT_MODE)
    INCLUDE 'mpif.h'
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
    ! DESCENT_MODE (optionnel, defaut 0 = automatique) : schema de la
    ! descente L et de la remontee U (meme principe que le LDLt).
    !   1 = "colonne" : X_K rassemble en pleine largeur sur la colonne de
    !       grille PCOL_K, qui calcule seule L(I,K)*X_K (ou U(I,K)*X_K) pour
    !       tous les seconds membres, puis disperse le long des lignes.
    !   2 = "distribuee" : X_K (tranche locale de seconds membres) diffuse le
    !       long de la colonne de grille ; la colonne K de L (ou de U) de mes
    !       lignes diffusee compressee le long de la ligne de grille
    !       (BCAST_BLR_PANEL) ; chaque processus l'applique a SES seconds
    !       membres. Toutes les colonnes de grille calculent.
    !   0 = automatique : mode 2 si NRHS > 1.5 x largeur effective des
    !       facteurs (LU_EFFECTIVE_WIDTH), seuil mesure pour le LDLt.
    ! Sans l'argument, la variable d'environnement BLR_LU_SOLVE_MODE (1 ou
    ! 2) force le schema (comparaisons).
    INTEGER, INTENT(IN), OPTIONAL :: DESCENT_MODE

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, KW, MI, RS, RE, NROWTRAIL, C, IERR
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, PBAR_STATE
    INTEGER :: DMODE, NPB
    INTEGER, ALLOCATABLE :: PIX(:), PRT(:), PMS(:), PNS(:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: CB(:)
    COMPLEX(WP), ALLOCATABLE :: XK(:,:)
    CHARACTER(LEN=8) :: ENVM
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
    DMODE = 0
    IF (PRESENT(DESCENT_MODE)) THEN
      DMODE = DESCENT_MODE
    ELSE
      CALL GET_ENVIRONMENT_VARIABLE('BLR_LU_SOLVE_MODE', ENVM)
      IF (TRIM(ENVM) .EQ. '1') DMODE = 1
      IF (TRIM(ENVM) .EQ. '2') DMODE = 2
    END IF
    IF (DMODE .NE. 1 .AND. DMODE .NE. 2) THEN
      DMODE = 1
      IF (REAL(NRHS,WP) .GT. 1.5_WP*LU_EFFECTIVE_WIDTH(GRID, A_BLR)) DMODE = 2
    END IF
    IF (DMODE .EQ. 2) THEN
      NPB = MAX(A_BLR%NBLKR_LOC, 1)
      ALLOCATE(PIX(NPB), PRT(NPB), PMS(NPB), PNS(NPB), CB(NPB))
    END IF
    ! tailles des tampons de seconds membres (pleine largeur) en entier 32 bits
    CALL BLR_CHECK_COUNT32('lignes locales x seconds membres', MAX(LDB,1), NRHS, GRID%COMM)
    CALL BLR_CHECK_COUNT32('taille de bloc x seconds membres', A_BLR%NB, NRHS, GRID%COMM)
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
      END IF
      IF (DMODE .EQ. 2) THEN
        CALL DIST_STEP(K, .TRUE.)
      ELSE
      IF (GRID%MYROW .EQ. PROW_K) THEN
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
      END IF
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
      END IF
      IF (DMODE .EQ. 2) THEN
        CALL DIST_STEP(K, .FALSE.)
      ELSE
      IF (GRID%MYROW .EQ. PROW_K) THEN
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
      END IF
      STEP = STEP + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Resolution (pivot)', STEP, NSTEPS, PBAR_STATE)
      IF (PRESENT(PROGRESS_CALLBACK)) THEN
        IF (STAGE_PROGRESS_SHOULD_FIRE(STEP, NSTEPS)) CALL PROGRESS_CALLBACK(STEP, NSTEPS)
      END IF
    END DO

    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  CONTAINS

    ! Etape K de la descente (LOWER : blocs I > K, colonne K de L) ou de la
    ! remontee (blocs I < K, colonne K de U) en mode distribue, apres la
    ! resolution du bloc diagonal sur la ligne de grille PROW_K.
    SUBROUTINE DIST_STEP(KS, LOWER)
      INTEGER, INTENT(IN) :: KS
      LOGICAL, INTENT(IN) :: LOWER
      INTEGER :: LIS, IS, MIS, RSS, LJS, B, IERRS
      ! X_K : ma tranche de seconds membres, le long de ma colonne de grille.
      IF (NRHS_LOC .GT. 0) THEN
        ALLOCATE(XK(KW,NRHS_LOC))
        IF (GRID%MYROW .EQ. PROW_K) THEN
          RSS = LOCAL_ROW_OFFSET(A_BLR, GRID, KS)
          XK = B_LOCAL(RSS:RSS+KW-1,1:NRHS_LOC)
        END IF
        CALL MPI_BCAST(XK, KW*NRHS_LOC, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERRS)
      END IF
      ! Colonne K des facteurs, pour mes lignes : un envoi groupe le long de
      ! la ligne de grille, depuis une copie chez le proprietaire.
      LJS = BLR_COL_LOCAL(A_BLR%CMAP, KS, GRID)
      NPB = 0
      DO LIS = 1, A_BLR%NBLKR_LOC
        IS = (LIS-1)*GRID%PR + GRID%MYROW + 1
        IF (IS .EQ. KS .OR. ((IS .GT. KS) .NEQV. LOWER)) CYCLE
        NPB = NPB + 1
        PIX(NPB) = LIS; PRT(NPB) = PCOL_K
        PMS(NPB) = A_BLR%BEND(IS)-A_BLR%BSTART(IS)+1; PNS(NPB) = KW
        IF (GRID%MYCOL .EQ. PCOL_K) CB(LIS) = A_BLR%OFFBLK(LIS,LJS)
      END DO
      CALL BCAST_BLR_PANEL(GRID%ROW_COMM, NPB, PIX, PRT, PMS, PNS, CB)
      ! B_I <- B_I - bloc(I,K) * X_K sur mes seconds membres.
      RSS = 1; B = 0
      DO LIS = 1, A_BLR%NBLKR_LOC
        IS = (LIS-1)*GRID%PR + GRID%MYROW + 1
        MIS = A_BLR%BEND(IS)-A_BLR%BSTART(IS)+1
        IF (IS .NE. KS .AND. ((IS .GT. KS) .EQV. LOWER)) THEN
          IF (NRHS_LOC .GT. 0) CALL BLR_BLOCK_APPLY(CB(LIS), NRHS_LOC, -CONE, XK, KW, CONE, B_LOCAL(RSS,1), LDB)
          CALL BLR_BLOCK_RELEASE(CB(LIS))
        END IF
        RSS = RSS + MIS
      END DO
      IF (ALLOCATED(XK)) DEALLOCATE(XK)
    END SUBROUTINE DIST_STEP

  END SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D

  !-----------------------------------------------------------------------
  ! Decoupage simple par blocs de NRHS sur PC (identique a RHS_LOCAL_RANGE
  ! de mpi_lu_solve_2d_mod.f90 -- duplique ici pour eviter de dependre du
  ! module de resolution DENSE dans ce module BLR ; meme formule).
  !-----------------------------------------------------------------------
  ! Garde-fou des compteurs MPI et des tailles en entier 32 bits : arrete
  ! tous les processus avec un message si NA*NB_ depasse HUGE(1) (au lieu
  ! d'un debordement silencieux). QUOI : grandeur concernee, pour le message.
  SUBROUTINE BLR_CHECK_COUNT32(QUOI, NA, NB_, COMM)
    INCLUDE 'mpif.h'
    CHARACTER(LEN=*), INTENT(IN) :: QUOI
    INTEGER, INTENT(IN) :: NA, NB_, COMM
    INTEGER :: IERR
    IF (INT(NA,8)*INT(NB_,8) .LE. INT(HUGE(1),8)) RETURN
    WRITE(6,'(A,A,A,I0,A,I0,A)') ' ERREUR BLR : ', QUOI, ' = ', NA, ' x ', NB_, &
         ' depasse 2^31-1 elements (compteur MPI 32 bits) : reduire le nombre de seconds membres par lot.'
    CALL MPI_ABORT(COMM, 3, IERR)
  END SUBROUTINE BLR_CHECK_COUNT32

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


END MODULE MPI_BLR_MATRIX_MOD
