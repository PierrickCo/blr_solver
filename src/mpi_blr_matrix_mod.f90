!===============================================================================
! SOLVEUR LU MPI, BLR, SANS PIVOT (etape 1 d'un plan en 2 etapes explicitement
! valide avec l'utilisateur -- le pivot viendra ensuite, sur cette base).
!
! DIFFERENCE avec TOUT ce qui existait jusqu'ici cote MPI (6c/6d/6e/P3) :
! ces modules stockent TOUJOURS A_LOCAL(MLOC,NLOC) comme un tableau DENSE
! (la compression n'accelere que la mise a jour de Schur, ou la diffusion
! des panneaux -- jamais le STOCKAGE local lui-meme). ICI, le stockage local
! de CHAQUE processus est un tableau de BLR_BLOCK_T (compresse des
! l'assemblage, EXACTEMENT comme blr_matrix_mod.f90/"A" en sequentiel) --
! "evite au maximum le passage au bloc plein", demande explicite.
!
! STRUCTURE DE DONNEES : au lieu de l'indexation ELEMENT-level habituelle
! (NUMROC/INDXL2G sur des LIGNES/COLONNES individuelles), on indexe ici au
! niveau BLOC : chaque processus possede NBLKR_LOC blocs-lignes locaux et
! NBLKC_LOC blocs-colonnes locaux (BLOCK_NUMROC, meme principe que NUMROC
! mais un "element" = un bloc NB x NB entier, cycle=1 bloc par processus).
! Bloc-ligne LOCAL li <-> bloc-ligne GLOBAL I = (li-1)*PR + MYROW + 1 (et
! symetriquement pour les colonnes) -- CHAQUE processus d'une meme ligne de
! grille (meme MYROW) partage EXACTEMENT le meme ensemble de lignes-bloc
! locales (la formule ne depend que de MYROW/PR), ce qui est la cle de
! l'algorithme de diffusion ci-dessous.
!
! MESSAGES MPI MINIMISES (demande explicite) : les panneaux L21/U12 sont
! diffuses DEJA COMPRESSES (technique de la proposition 3, BCAST_BLR_BLOCK)
! -- jamais dense-puis-compression (6e) ni meme dense du tout ici. Comme
! CHAQUE bloc du panneau garde sa propre representation (dense ou UFAC/VFAC
! de rang tres inferieur a NB), le volume envoye est proportionnel au rang
! reel de la matrice, pas a NB.
!
! MISE A JOUR DE SCHUR : reste compressee aussi souvent que possible via
! SCHUR_UPDATE_BLR_BLOCK (deja construit pour blr_matrix_mod.f90/"A") --
! meme mecanisme, applique ici bloc local par bloc local.
!
! RESOLUTION : reprend l'algorithme 2D complet de mpi_lu_solve_2d_mod.f90
! (ALLGATHERV+diffusion+SCATTERV sur B, AUCUNE replication residuelle) mais
! remplace le CGEMM dense par BLR_BLOCK_APPLY, bloc local par bloc local --
! exploite directement la compression des facteurs stockes.
!
! HORS SCOPE ICI (valide avec l'utilisateur, etape suivante) : PAS de pivot
! (matrices a diagonale dominante uniquement, comme TOUS les modules "sans
! pivot" du projet) ; PAS de memoire partagee MPI (MPI_Win/RMA) -- cette
! architecture ne l'empeche pas d'etre ajoutee plus tard (chaque processus
! garde une frontiere claire entre ses donnees locales et les messages
! recus), mais ce n'est pas construit maintenant.
!===============================================================================
MODULE MPI_BLR_MATRIX_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE, DENSE_PANEL_T
  USE PROGRESS_BAR_MOD, ONLY: PRINT_PROGRESS_MILESTONE
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE

  TYPE :: MPI_BLR_MATRIX_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
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
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(DENSE_PANEL_T), ALLOCATABLE      :: DIAGBLK(:)
    TYPE(BLR_BLOCK_ARENA_T), ALLOCATABLE  :: OFFBLK(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:,:)
    TYPE(BLR_ARENA_T) :: ARENA
  END TYPE MPI_BLR_MATRIX_ARENA_T

  !-----------------------------------------------------------------------
  ! Variante "CAPACITE FIXE" (demande explicite : "les blocs LR stoquent
  ! le maximum de rang que COMPRESS_WORTHWHILE_RATIO le permet ; les
  ! echanges MPI et les operations associees ne se font que sur le rang
  ! REEL") -- voir BLR_BLOCK_FIXEDCAP_T (blr_compress_mod.f90) pour le
  ! principe complet. Contrairement a MPI_BLR_MATRIX_ARENA_T, AUCUN champ
  ! ARENA ici -- chaque bloc gere sa propre capacite RMAX via de simples
  ! composantes ALLOCATABLE, aucun allocateur/pool partage, aucun attribut
  ! TARGET necessaire nulle part dans cette variante.
  !-----------------------------------------------------------------------
  TYPE :: MPI_BLR_MATRIX_FIXEDCAP_T
    INTEGER :: N = 0, NB = 0, NBLK = 0, NBLKR_LOC = 0, NBLKC_LOC = 0
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(DENSE_PANEL_T), ALLOCATABLE         :: DIAGBLK(:)
    TYPE(BLR_BLOCK_FIXEDCAP_T), ALLOCATABLE  :: OFFBLK(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:,:)
  END TYPE MPI_BLR_MATRIX_FIXEDCAP_T

CONTAINS

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
  ! ETAPE 1 : assemblage bloc par bloc + compression IMMEDIATE des blocs
  ! hors-diagonale locaux (meme esprit que ASSEMBLE_BLR_MATRIX/"A", mais
  ! chaque processus n'assemble que SA part -- GET_BLOCK recoit des
  ! indices GLOBAUX, comme en sequentiel).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, &
                                      T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                      SHOW_PROGRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_T), INTENT(OUT) :: A_BLR
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS

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
    LAST_PCT = -5

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%NBLKC_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PC, GRID%MYCOL)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(li,lj))
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
  ! ETAPE 3 : factorisation LU EN PLACE, SANS PIVOT, sur A_BLR -- reste
  ! compressee tout du long (voir commentaire de tete de fichier pour
  ! l'algorithme complet). N_COMPRESSED/N_OFFDIAG/RANK_SUM caracterisent
  ! les facteurs L/U FINAUX (sommes sur tous les processus).
  !-----------------------------------------------------------------------
  ! T_DIAG/T_BCAST_DKK/T_TRSM/T_BCAST_PANEL/T_ALLOC/T_MPICALL/T_SCHUR_PRODUCT/
  ! T_RECOMPRESS/T_SCHUR_DIAG (tous optionnels, INTENT(OUT), max/processus) :
  ! decoupage fin du temps de factorisation demande par l'utilisateur --
  ! T_BCAST_PANEL est le temps TOTAL de l'etape 5 (diffusion des panneaux
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
    INCLUDE 'mpif.h'
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
      PCOL_K = MOD(K-1, GRID%PC)
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
        lj = (K-1)/GRID%PC + 1
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

      ! 5) diffusion des panneaux L/U DEJA COMPRESSES (reduit au maximum
      !    le volume envoye -- technique de la proposition 3) : le
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
                                A_BLR%OFFBLK(li,(K-1)/GRID%PC+1), T_ALLOC_L, T_MPICALL_L)
          LRECV(li) = A_BLR%OFFBLK(li,(K-1)/GRID%PC+1)
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, &
                                A_BLR%OFFBLK((K-1)/GRID%PR+1,lj), T_ALLOC_L, T_MPICALL_L)
          URECV(lj) = A_BLR%OFFBLK((K-1)/GRID%PR+1,lj)
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj), &
                                T_ALLOC_L, T_MPICALL_L)
        END DO
      END IF
      T_BCASTPANEL_L = T_BCASTPANEL_L + (MPI_WTIME() - TT0)

      ! 6) mise a jour de Schur locale, bloc par bloc, en gardant chaque
      !    bloc hors-diagonale compresse aussi souvent que possible.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
  ! ETAPE 3bis : factorisation LU EN PLACE, AVEC PIVOT PARTIEL CONFINE AU
  ! BLOC DIAGONAL (meme philosophie que block_lu_mod.f90/etape 1 et
  ! mpi_lu_pivot_mod.f90/6d -- pivot recherche UNIQUEMENT dans les KW
  ! lignes du bloc diagonal courant, jamais globalement). Choisi plutot
  ! qu'une methode "ne decompresse jamais" (option initialement envisagee,
  ! haut risque, cf. discussion avec l'utilisateur) car permuter des LIGNES
  ! d'un bloc bas-rang ne touche QUE UFAC (cout O(rang), jamais de
  ! decompression) -- voir BLR_BLOCK_ROW_SWAP.
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
                                             SHOW_PROGRESS, SHOW_RSS_TRACE, USE_MOVE_ALLOC)
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
    ! DIAGNOSTIC TEMPORAIRE (a retirer une fois l'investigation terminee,
    ! voir la memoire du projet -- ecart RSS-vs-theorique apres factori-
    ! sation, causes MKL/malloc_trim ecartees) : affiche la RSS du
    ! processus 0 A CHAQUE etape K, pour voir la FORME de la croissance
    ! (reguliere = fragmentation continue plausible ; par sauts ponctuels
    ! = plutot un redimensionnement de workspace qui grossit sans jamais
    ! retrecir). N'affecte AUCUN resultat -- uniquement des WRITE.
    LOGICAL, INTENT(IN), OPTIONAL     :: SHOW_RSS_TRACE
    ! USE_MOVE_ALLOC (optionnel, defaut .FALSE.) : transmis tel quel a
    ! SCHUR_UPDATE_BLR_BLOCK (voir blr_compress_mod.f90) -- evite la
    ! reallocation implicite de "ABLK = NEWBLK" au profit d'un transfert
    ! de propriete (MOVE_ALLOC). Resultat numerique IDENTIQUE.
    LOGICAL, INTENT(IN), OPTIONAL     :: USE_MOVE_ALLOC

    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LRECV(:), URECV(:)
    TYPE(BLR_WORKSPACE_T) :: WS
    INTEGER, ALLOCATABLE :: IPIV_K(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, NJ, IINFO, INFO_STEP
    INTEGER :: PROW_K, PCOL_K, IERR, NC_LOCAL, NO_LOCAL, NS_LOCAL, li_K, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS, RSS_TRACE

    ! Boucle K=1..NBLK GLOBALEMENT synchronisee (MPI_BCAST a chaque etape,
    ! voir plus bas) -- contrairement a l'assemblage, K/NBLK est ici un
    ! pourcentage d'avancement EXACT (tous les processus traitent le meme K
    ! au meme moment), pas une approximation. Affiche uniquement par le
    ! processus 0 pour eviter des lignes entrelacees.
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
      PCOL_K = MOD(K-1, GRID%PC)
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
        lj = (K-1)/GRID%PC + 1
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
                                A_BLR%OFFBLK(li,(K-1)/GRID%PC+1))
          LRECV(li) = A_BLR%OFFBLK(li,(K-1)/GRID%PC+1)
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, &
                                A_BLR%OFFBLK((K-1)/GRID%PR+1,lj))
          URECV(lj) = A_BLR%OFFBLK((K-1)/GRID%PR+1,lj)
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      END IF

      ! 7) mise a jour de Schur locale -- inchange.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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

      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('3) Factorisation LU (pivot)', K, NBLK, LAST_PCT)
      IF (RSS_TRACE) WRITE(6,'(A,I4,A,F8.1,A)') '     [RSS trace] K=', K, ' : RSS = ', &
           REAL(GET_RSS_KB(),8)/1024.0D0, ' Mo'
    END DO

    DEALLOCATE(LRECV, URECV)
    CALL BLR_WORKSPACE_FREE(WS)
    IF (PRESENT(T_FACTOR)) T_FACTOR = MPI_WTIME() - T0

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NO_LOCAL, N_OFFDIAG, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RS_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    ! NS_LOCAL est deja IDENTIQUE sur tous les processus (chacun le calcule
    ! a partir de SA PROPRE copie de IPIV_K, diffusee identique a tous via
    ! GRID%COMM ci-dessus) -- aucune reduction necessaire, juste le lire.
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS_LOCAL
  END SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT

  !-----------------------------------------------------------------------
  ! ETAPE 4 : resolution A*X=B, DISTRIBUTION 2D COMPLETE de B (comme
  ! mpi_lu_solve_2d_mod.f90 -- ALLGATHERV+diffusion+SCATTERV, RSTARTS/RLOCS
  ! calcules localement sans communication), mais la contribution
  ! MI x NRHS de chaque bloc trainant local est calculee via
  ! BLR_BLOCK_APPLY (exploite directement la compression), PAS un CGEMM
  ! dense sur tout le panneau -- boucle bloc local par bloc local.
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)           :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)         :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, RS, RE, NROWTRAIL, C, IERR
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
      PCOL_K = MOD(K-1, GRID%PC)

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
      lj = (K-1)/GRID%PC + 1
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
      PCOL_K = MOD(K-1, GRID%PC)

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
      lj = (K-1)/GRID%PC + 1
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
  SUBROUTINE MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, LDB, SHOW_PROGRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)           :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)         :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL      :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, RS, RE, NROWTRAIL, C, IERR
    INTEGER :: PROW_K, PCOL_K, STEP, NSTEPS, LAST_PCT
    LOGICAL :: PROGRESS

    ! Une seule barre pour toute la resolution (substitution avant PUIS
    ! arriere, 2*NBLK etapes au total) -- comme la factorisation, chaque K
    ! est GLOBALEMENT synchronise (MPI_BCAST a chaque etape), K/NBLK (ou
    ! ici STEP/NSTEPS) est donc un pourcentage exact, pas une approximation.
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

    ! === 0) P*B : applique le swap de CHAQUE bloc a B_LOCAL EN UNE FOIS,
    !     AVANT toute substitution -- IMPORTANT : ceci doit precede TOUTE
    !     soustraction croisee entre blocs (etape 1 ci-dessous), sinon le
    !     swap de l'etape K permute aussi les contributions DEJA soustraites
    !     par les etapes m<K, ce qui est faux (P_K ne doit s'appliquer qu'au
    !     B ORIGINAL, jamais a une valeur partiellement mise a jour -- bug
    !     reel trouve et corrige : la version precedente appliquait le swap
    !     "en ligne" au debut de l'iteration K, trop tard). Chaque processus
    !     traite uniquement les blocs dont il possede des lignes
    !     (MYROW==PROW_K) ; IPIV deja diffuse a toute la grille, aucune
    !     communication necessaire ici. ===
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
      PCOL_K = MOD(K-1, GRID%PC)

      ALLOCATE(XK_FULL(KW,NRHS))
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
        DO C = 0, GRID%PC-1
          RCOUNTS(C) = KW*RLOCS(C); RDISPLS(C) = KW*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XK_LOC, KW*NRHS_LOC, MPI_COMPLEX, &
                             XK_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, GRID%ROW_COMM, IERR)
        DEALLOCATE(XK_LOC)
      END IF
      CALL MPI_BCAST(XK_FULL, KW*NRHS, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)

      lj = (K-1)/GRID%PC + 1
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('4) Resolution (pivot)', STEP, NSTEPS, LAST_PCT)
    END DO

    ! === 2) U*X = Y (substitution arriere) -- AUCUN pivot necessaire ===
    DO K = NBLK, 1, -1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      PROW_K = MOD(K-1, GRID%PR)
      PCOL_K = MOD(K-1, GRID%PC)

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

      lj = (K-1)/GRID%PC + 1
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
      IF (PROGRESS) CALL PRINT_PROGRESS_MILESTONE('4) Resolution (pivot)', STEP, NSTEPS, LAST_PCT)
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
                                            SHOW_PROGRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(OUT), TARGET :: A_BLR
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS

    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: li, lj, I, J, MI, NJ, NC_LOCAL, NO_LOCAL, IERR, NDONE, NTOTAL, LAST_PCT
    INTEGER(8) :: RS_LOCAL
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS .AND. (GRID%MYID.EQ.0)
    NDONE = 0
    LAST_PCT = -5

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%NBLKC_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PC, GRID%MYCOL)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          CALL COMPRESS_BLOCK_ARENA(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(li,lj), A_BLR%ARENA)
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
  ! DIAGNOSTIC TEMPORAIRE (a retirer une fois l'investigation terminee,
  ! voir la memoire du projet -- reste ~6x de gaspillage arena apres le
  ! correctif de reutilisation en place, ~x28 avant) : decompose
  ! l'empreinte OFFBLK locale (ARENA_BYTES(A_BLR%ARENA)) en 3 categories,
  ! pour verifier l'hypothese avant de coder quoi que ce soit de plus :
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
    INCLUDE 'mpif.h'
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
      PCOL_K = MOD(K-1, GRID%PC)
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
        lj = (K-1)/GRID%PC + 1
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
          CALL COPY_ARENA_TO_REGULAR(A_BLR%OFFBLK(li,(K-1)/GRID%PC+1), LRECV(li))
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL COPY_ARENA_TO_REGULAR(A_BLR%OFFBLK((K-1)/GRID%PR+1,lj), URECV(lj))
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      END IF

      ! 7) mise a jour de Schur -- SCHUR_UPDATE_BLR_BLOCK_ARENA.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)                 :: GRID
    TYPE(MPI_BLR_MATRIX_ARENA_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                      :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)               :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL            :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, RS, RE, NROWTRAIL, C, IERR
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
      PCOL_K = MOD(K-1, GRID%PC)

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

      lj = (K-1)/GRID%PC + 1
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
      PCOL_K = MOD(K-1, GRID%PC)

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

      lj = (K-1)/GRID%PC + 1
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
                                               SHOW_PROGRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB
    REAL(WP), INTENT(IN)            :: TOL
    PROCEDURE(GET_BLOCK_IFACE)      :: GET_BLOCK
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(OUT) :: A_BLR
    DOUBLE PRECISION, INTENT(OUT)   :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL   :: SHOW_PROGRESS

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

    T0 = MPI_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%NBLKC_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PC, GRID%MYCOL)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))
    NTOTAL = A_BLR%NBLKR_LOC * A_BLR%NBLKC_LOC

    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8

    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, TMPBLK)
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
  !     INCHANGEE -- le volume MPI ne depend QUE du rang REEL, jamais de
  !     RMAX (demande explicite de l'utilisateur).
  !   - etape 7 : SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP (ecrase en place, ZERO
  !     reallocation apres la premiere compression du bloc) -- pas de
  !     USE_MOVE_ALLOC ici, non pertinent (aucune reallocation a eviter).
  !-----------------------------------------------------------------------
  SUBROUTINE MPI_BLR_MATRIX_LU_FACTOR_PIVOT_FIXEDCAP(GRID, A_BLR, TOL, INFO, &
                                                      N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                                      SHOW_PROGRESS, SHOW_RSS_TRACE)
    INCLUDE 'mpif.h'
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
      PCOL_K = MOD(K-1, GRID%PC)
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
        lj = (K-1)/GRID%PC + 1
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
          CALL COPY_FIXEDCAP_TO_REGULAR(A_BLR%OFFBLK(li,(K-1)/GRID%PC+1), LRECV(li))
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
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL COPY_FIXEDCAP_TO_REGULAR(A_BLR%OFFBLK((K-1)/GRID%PR+1,lj), URECV(lj))
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      ELSE
        DO lj = 1, A_BLR%NBLKC_LOC
          J = (lj-1)*GRID%PC + GRID%MYCOL + 1
          IF (J .LE. K) CYCLE
          CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, A_BLR%BEND(J)-A_BLR%BSTART(J)+1, URECV(lj))
        END DO
      END IF

      ! 7) mise a jour de Schur -- SCHUR_UPDATE_BLR_BLOCK_FIXEDCAP.
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)                    :: GRID
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)                         :: NRHS, NRHS_LOC, LDB
    COMPLEX(WP), INTENT(INOUT)                  :: B_LOCAL(MAX(LDB,1),MAX(NRHS_LOC,1))
    LOGICAL, INTENT(IN), OPTIONAL               :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:), XK_LOC(:,:), XK_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: NBLK, K, li, lj, I, J, KW, MI, RS, RE, NROWTRAIL, C, IERR
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
      PCOL_K = MOD(K-1, GRID%PC)

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

      lj = (K-1)/GRID%PC + 1
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
      PCOL_K = MOD(K-1, GRID%PC)

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

      lj = (K-1)/GRID%PC + 1
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

END MODULE MPI_BLR_MATRIX_MOD
