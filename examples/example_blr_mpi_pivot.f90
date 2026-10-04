!===============================================================================
! EXEMPLE CANONIQUE -- BLR+MPI avec pivot partiel confine au bloc diagonal
! (MPI_BLR_MATRIX_LU_FACTOR_PIVOT), comme MUMPS/STRUMPACK : jamais de pivot
! global -- permuter les lignes d'un bloc compresse ne coute que O(rang).
! Tolere les matrices sans diagonale dominante (contrairement a la variante
! sans pivot) ; la matrice de demo est volontairement mal conditionnee pour
! forcer un vrai pivotage.
!
! But : montrer l'enchainement minimal (assembler, factoriser, resoudre,
! rassembler) sans l'instrumentation de test/validation croisee (voir
! test/test_mpi_blr_matrix_pivot.f90 pour cela).
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]
!         (--help/-h affiche la liste complete). PR PC obligatoires ;
!         N/NB/NRHS/TOL positionnels optionnels dans cet ordre (defauts
!         4000/200/4/1e-4) ; toutes les autres options sont NOMMEES,
!         --cle=valeur ou --cle seul (= --cle=1) pour les booleennes, dans
!         n'importe quel ordre : --use-move-alloc --eps --ratio --use-rsvd
!         --check-error --check-residual --profile-usage --dump-structure
!         --save-matrix --save-rhs --save-factored --save-solution
!         --load-matrix --load-factored --load-rhs --checkpoint-dir
!         --nrhs-batch --ntailbloc (N/NB/NRHS/TOL aussi dispo en --n/--nb/
!         --nrhs/--tol). Exemple :
!           mpirun -np 4 ./example_blr_mpi_pivot 2 2 4000 200 8 1e-4 \
!               --check-error --save-factored --checkpoint-dir=ckpt
!
! USE_MOVE_ALLOC (defaut 0) : SCHUR_UPDATE_BLR_BLOCK transfere la propriete
! via MOVE_ALLOC au lieu de "ABLK=NEWBLK" (copie implicite en moins) --
! resultat identique, gain RSS mesure negligeable a N=20000, garde par
! defaut desactive pour ne rien changer au comportement existant.
!
! EPS (bruit plein-rang, defaut 0) : ajoute EPS*HASH_NOISE(GI,GJ) a chaque
! entree hors anti-diagonale du noyau de Cauchy. EPS croissant degrade
! progressivement la compressibilite (plus de blocs en repli dense a TOL
! fixe) sans toucher l'anti-diagonale (pivot). Essayer 0/0.01/0.1/1.0 pour
! voir la transition compressible -> incompressible (voir le taux de
! compression REEL affiche apres assemblage/factorisation).
!
! RATIO (COMPRESS_WORTHWHILE_RATIO, defaut 1/3) : un bloc M x N reste
! bas-rang seulement si son rang R < RATIO*MIN(M,N), sinon repli DENSE
! meme si la SVD trouve un rang valide sous TOL. RATIO plus grand accepte
! plus de blocs "bas-rang" (UFAC/VFAC plus larges) ; plus petit est plus
! strict. Effet surtout visible pres du plafond (voir
! test/bench_blr_mpi_pivot_ratio.f90).
!
! USE_RSVD (defaut 1) : 0=SVD complete (CGESVD, rang optimal garanti) ;
! 1=SVD randomisee adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE) -- mesuree
! x34-54 plus rapide, erreur du meme ordre, mais heuristique (pas de
! preuve formelle) avec une zone morte pres du plafond RATIO (replie en
! dense la ou la SVD complete aurait compresse -- jamais de corruption,
! juste un gain manque).
!
! CHECK_ERROR (defaut 0) : si 1, B=A*Xexact (Xexact connu par formule) au
! lieu d'un B arbitraire -- permet d'afficher ||X-Xexact||/||Xexact|| en
! fin de resolution. Calcule par produit BRUT (jamais une 2e matrice en
! RAM) : cout calcul O(N^2*NRHS) supplementaire (peut prendre plusieurs
! dizaines de secondes a N=20000) mais aucun cout memoire au-dela de
! X_GLOBAL deja alloue pour l'affichage -- l'avantage memoire MPI/OOC
! reste entier.
!
! CHECK_RESIDUAL (defaut 0, independant de CHECK_ERROR) : ||A*X-B||/||B||
! apres resolution, avec le B REELLEMENT utilise -- cas d'usage general
! (verifier sans connaitre la solution exacte). Meme principe que
! CHECK_ERROR : produit brut, zero cout memoire supplementaire, cout
! calcul O(N^2*NRHS) cumulatif si les deux sont actifs.
!
! PROFILE_USAGE (defaut 0) : si 1, chaque processus enregistre a chaque
! etape (+ chaque palier de 5% dans assemblage/factorisation/resolution,
! via PROGRESS_CALLBACK) le temps ecoule, sa RSS et son %CPU (CPU_TIME/
! temps ecoule, comme `top`) dans usage_profile_rank<N>.csv. Tracer :
! python3 examples/plot_usage_profile.py [dossier] (matplotlib).
!
! DUMP_STRUCTURE (defaut 0) : si 1, chaque processus ecrit la nature de
! ses blocs locaux (diagonal/dense/bas-rang + rang) apres l'assemblage ET
! apres la factorisation (le fill-in de Schur peut changer certains blocs
! de bas-rang a dense) -- blr_structure_<etape>_rank<N>.csv. Tracer :
! python3 examples/plot_blr_structure.py [dossier] (matplotlib).
!
! SAUVEGARDE/REPRISE (SAVE_MATRIX/RHS/FACTORED/SOLUTION, LOAD_MATRIX/
! FACTORED/RHS, CHECKPOINT_DIR) : chaque SAVE_* ecrit un fichier PAR
! PROCESSUS dans CHECKPOINT_DIR (defaut '.', doit deja exister), STREAM
! non formate, sans communication MPI -- meme convention que les caches
! OOC.
!   SAVE_MATRIX/FACTORED : A_BLR apres assemblage / apres factorisation.
!   SAVE_RHS/SOLUTION     : B_LOCAL apres FILL_B / apres resolution.
!   LOAD_MATRIX/FACTORED  : relit et saute assemblage / assemblage+
!                           factorisation (LOAD_FACTORED prioritaire).
!   LOAD_RHS              : relit et saute FILL_B.
! Chaque fichier commence par un en-tete (N/NB/PR/PC/grille/rang) verifie
! au chargement -- une incompatibilite sur UN SEUL rang arrete TOUS les
! rangs (MPI_ALLREDUCE), jamais un sous-ensemble qui continue seul.
!
! NRHS_BATCH (defaut 0 = pas de decoupage) : borne le nombre de colonnes
! GLOBALES de second membre traitees a la fois -- assemblage/resolution/
! sauvegarde se font PAR LOT, jamais plus d'un lot de large NRHS_BATCH en
! RAM. Reutilise MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D tel quel (deja generique
! sur la largeur). Resultat mathematiquement identique decoupe ou non ;
! en pratique, de legers ecarts d'arrondi (simple precision, matrice mal
! conditionnee) sont attendus entre deux decoupages differents -- pas un
! signe d'erreur. LOAD_RHS ignore NRHS_BATCH et relit le decoupage stocke
! a la sauvegarde. LIMITE : CHECK_ERROR/CHECK_RESIDUAL rassemblent
! toujours X_GLOBAL(N,NRHS) en entier (cout inchange) -- a eviter avec un
! NRHS massif.
!
! NTAILBLOC (defaut 0 = desactive) : quand SAVE_RHS/SAVE_SOLUTION/LOAD_RHS
! sont actifs, ecrit/relit le contenu numerique de chaque lot par BLOCS DE
! NTAILBLOC LIGNES GLOBALES, un fichier par bloc (saved_rhs_batch<B>_
! rank<N>_bloc<K>.dat, idem solution) -- un seul WRITE/READ par bloc,
! jamais tout le lot en memoire (voir SAVE_ROWBLOCKS/LOAD_ROWBLOCKS). Ne
! s'applique pas a SAVE_MATRIX/FACTORED (facteur LU compresse, pas de
! decoupage naturel par lignes). LOAD_RHS doit utiliser EXACTEMENT le
! meme NTAILBLOC qu'a la sauvegarde (verifie bloc par bloc).
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_PIVOT
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL, EPS, RATIO
  CHARACTER(LEN=256) :: ARG
  CHARACTER(LEN=64)  :: KEYSTR
  CHARACTER(LEN=256) :: VALSTR
  INTEGER :: IARG, POS_IDX
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  LOGICAL :: WANT_FULL_X
  ! Apercu de X(1:5,1) sans assembler X_GLOBAL en entier -- voir son
  ! utilisation plus bas (CHECK_ERROR=0 et CHECK_RESIDUAL=0).
  INTEGER :: NPREVIEW, PREVIEW_RS, PREVIEW_MI, PREVIEW_GI
  COMPLEX(WP) :: PREVIEW_LOCAL(5), PREVIEW_GLOBAL(5)

  ! Second membre/resolution/sauvegarde PAR LOT (NRHS_BATCH, voir le
  ! commentaire de tete) -- bornent la RAM a un lot au lieu de NRHS entier.
  INTEGER :: NRHS_BATCH, NRHS_BATCH_EFF, NBATCHES, BATCH_IDX
  ! E/S rapide par blocs de lignes (voir SAVE_ROWBLOCKS/LOAD_ROWBLOCKS) --
  ! NTAILBLOC=0 (defaut) garde l'ancien format 1-fichier-par-lot inchange.
  INTEGER :: NTAILBLOC
  CHARACTER(LEN=8) :: BATCH_STR
  INTEGER :: RSTART_B, REND_B, WIDTH_B, LOCAL_OFFSET_B, GSTART_B, GEND_B, NRHS_LOC_B
  INTEGER :: SAVE_RHS_UNIT, SAVE_SOLUTION_UNIT, LOAD_RHS_UNIT
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  INTEGER :: USE_MOVE_ALLOC_INT, USE_RSVD_INT, CHECK_ERROR_INT, CHECK_RESIDUAL_INT, PROFILE_USAGE_INT, DUMP_STRUCTURE_INT
  LOGICAL :: USE_MOVE_ALLOC_FLAG, USE_RSVD_FLAG, CHECK_ERROR_FLAG, CHECK_RESIDUAL_FLAG, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG
  DOUBLE PRECISION :: RSS_MAX_B_MB, RSS_MAX_RESID_MB
  REAL(WP) :: ERR, RESID

  ! Sauvegarde/reprise (voir le commentaire de tete) -- un fichier par
  ! processus dans CHECKPOINT_DIR par objet sauvegarde/relu.
  INTEGER :: SAVE_MATRIX_INT, SAVE_RHS_INT, SAVE_FACTORED_INT, SAVE_SOLUTION_INT
  INTEGER :: LOAD_MATRIX_INT, LOAD_FACTORED_INT, LOAD_RHS_INT
  LOGICAL :: SAVE_MATRIX_FLAG, SAVE_RHS_FLAG, SAVE_FACTORED_FLAG, SAVE_SOLUTION_FLAG
  LOGICAL :: LOAD_MATRIX_FLAG, LOAD_FACTORED_FLAG, LOAD_RHS_FLAG, OK
  CHARACTER(LEN=256) :: CHECKPOINT_DIR

  ! Suivi memoire/CPU par etape (PROFILE_USAGE) -- tableaux de taille fixe,
  ! largement au-dela des ~8 etapes reellement utilisees ci-dessous.
  ! 32 suffisait pour un point par etape ; avec un echantillon tous les
  ! 5% DANS l'assemblage/la factorisation/la resolution (~20 chacune), il
  ! en faut nettement plus -- 128 laisse une bonne marge.
  INTEGER, PARAMETER :: MAX_PROFILE_SAMPLES = 128
  CHARACTER(LEN=40) :: PROF_STAGE(MAX_PROFILE_SAMPLES)
  DOUBLE PRECISION :: PROF_T(MAX_PROFILE_SAMPLES), PROF_RSS_MB(MAX_PROFILE_SAMPLES), PROF_CPU_PCT(MAX_PROFILE_SAMPLES)
  INTEGER :: PROF_N
  DOUBLE PRECISION :: PROF_T_START, PROF_T_PREV
  REAL :: PROF_CPU_PREV

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP; USE_MOVE_ALLOC_INT = 0; EPS = 0.0_WP
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 1; CHECK_ERROR_INT = 0; CHECK_RESIDUAL_INT = 0; PROFILE_USAGE_INT = 0
  DUMP_STRUCTURE_INT = 0
  SAVE_MATRIX_INT = 0; SAVE_RHS_INT = 0; SAVE_FACTORED_INT = 0; SAVE_SOLUTION_INT = 0
  LOAD_MATRIX_INT = 0; LOAD_FACTORED_INT = 0; LOAD_RHS_INT = 0; CHECKPOINT_DIR = '.'
  NRHS_BATCH = 0
  NTAILBLOC = 0
  PROF_N = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG)
    IF (TRIM(ARG).EQ.'--help' .OR. TRIM(ARG).EQ.'-h') THEN
      IF (MYID0.EQ.0) CALL PRINT_USAGE()
      CALL MPI_FINALIZE(IERR)
      STOP 0
    END IF
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) CALL PRINT_USAGE()
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC

  ! --- positionnels OPTIONNELS N, NB, NRHS, TOL, dans cet ordre fixe --
  !     s'arretent au premier "--" (bascule vers les options nommees) ou
  !     a la fin de la ligne de commande.
  POS_IDX = 0
  IARG = 3
  DO WHILE (IARG .LE. COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(IARG, ARG)
    IF (ARG(1:2) .EQ. '--') EXIT
    POS_IDX = POS_IDX + 1
    SELECT CASE (POS_IDX)
    CASE (1); READ(ARG,*) N
    CASE (2); READ(ARG,*) NB
    CASE (3); READ(ARG,*) NRHS
    CASE (4); READ(ARG,*) TOL
    CASE DEFAULT
      IF (MYID0.EQ.0) WRITE(6,'(A)') &
           ' ERREUR : trop d''arguments positionnels -- apres PR PC N NB NRHS TOL, '// &
           'toute option supplementaire doit utiliser --cle=valeur (voir --help).'
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END SELECT
    IARG = IARG + 1
  END DO

  ! --- options NOMMEES --cle=valeur (ou --cle seul = --cle=1), dans
  !     n'importe quel ordre (voir --help). Une option nommee l'emporte
  !     sur la meme valeur positionnelle.
  DO WHILE (IARG .LE. COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(IARG, ARG)
    IF (TRIM(ARG).EQ.'--help' .OR. TRIM(ARG).EQ.'-h') THEN
      IF (MYID0.EQ.0) CALL PRINT_USAGE()
      CALL MPI_FINALIZE(IERR)
      STOP 0
    END IF
    IF (ARG(1:2) .NE. '--') THEN
      IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : argument inattendu "', TRIM(ARG), &
           '" -- les options apres TOL doivent commencer par -- (voir --help).'
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    CALL SPLIT_KV(ARG, KEYSTR, VALSTR)
    SELECT CASE (TRIM(KEYSTR))
    CASE ('n');               READ(VALSTR,*) N
    CASE ('nb');              READ(VALSTR,*) NB
    CASE ('nrhs');            READ(VALSTR,*) NRHS
    CASE ('tol');             READ(VALSTR,*) TOL
    CASE ('use-move-alloc');  READ(VALSTR,*) USE_MOVE_ALLOC_INT
    CASE ('eps');             READ(VALSTR,*) EPS
    CASE ('ratio');           READ(VALSTR,*) RATIO
    CASE ('use-rsvd');        READ(VALSTR,*) USE_RSVD_INT
    CASE ('check-error');     READ(VALSTR,*) CHECK_ERROR_INT
    CASE ('check-residual');  READ(VALSTR,*) CHECK_RESIDUAL_INT
    CASE ('profile-usage');   READ(VALSTR,*) PROFILE_USAGE_INT
    CASE ('dump-structure');  READ(VALSTR,*) DUMP_STRUCTURE_INT
    CASE ('save-matrix');     READ(VALSTR,*) SAVE_MATRIX_INT
    CASE ('save-rhs');        READ(VALSTR,*) SAVE_RHS_INT
    CASE ('save-factored');   READ(VALSTR,*) SAVE_FACTORED_INT
    CASE ('save-solution');   READ(VALSTR,*) SAVE_SOLUTION_INT
    CASE ('load-matrix');     READ(VALSTR,*) LOAD_MATRIX_INT
    CASE ('load-factored');   READ(VALSTR,*) LOAD_FACTORED_INT
    CASE ('load-rhs');        READ(VALSTR,*) LOAD_RHS_INT
    CASE ('checkpoint-dir');  CHECKPOINT_DIR = VALSTR
    CASE ('nrhs-batch');      READ(VALSTR,*) NRHS_BATCH
    CASE ('ntailbloc');       READ(VALSTR,*) NTAILBLOC
    CASE DEFAULT
      IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : option inconnue "--', TRIM(KEYSTR), '" (voir --help).'
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END SELECT
    IARG = IARG + 1
  END DO

  USE_MOVE_ALLOC_FLAG = (USE_MOVE_ALLOC_INT .NE. 0)
  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
  CHECK_ERROR_FLAG = (CHECK_ERROR_INT .NE. 0)
  CHECK_RESIDUAL_FLAG = (CHECK_RESIDUAL_INT .NE. 0)
  PROFILE_USAGE_FLAG = (PROFILE_USAGE_INT .NE. 0)
  DUMP_STRUCTURE_FLAG = (DUMP_STRUCTURE_INT .NE. 0)
  SAVE_MATRIX_FLAG = (SAVE_MATRIX_INT .NE. 0)
  SAVE_RHS_FLAG = (SAVE_RHS_INT .NE. 0)
  SAVE_FACTORED_FLAG = (SAVE_FACTORED_INT .NE. 0)
  SAVE_SOLUTION_FLAG = (SAVE_SOLUTION_INT .NE. 0)
  LOAD_MATRIX_FLAG = (LOAD_MATRIX_INT .NE. 0)
  LOAD_FACTORED_FLAG = (LOAD_FACTORED_INT .NE. 0)
  LOAD_RHS_FLAG = (LOAD_RHS_INT .NE. 0)

  ! Surcharge COMPRESS_WORTHWHILE_RATIO avant tout assemblage (voir RATIO
  ! dans le commentaire de tete) ; le defaut (1/3) ne change rien.
  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,L1,A,ES9.2,A,F6.3)') &
         ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
         ' , USE_MOVE_ALLOC=', USE_MOVE_ALLOC_FLAG, ' , EPS(bruit)=', EPS, ' , RATIO=', RATIO
    IF (USE_RSVD_FLAG) THEN
      WRITE(6,'(A)') ' Compression : SVD randomisee adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE)'
    ELSE
      WRITE(6,'(A)') ' Compression : SVD complete (COMPRESS_BLOCK, CGESVD)'
    END IF
    IF (CHECK_ERROR_FLAG) THEN
      WRITE(6,'(A)') ' Verification : B=A*Xexact, erreur affichee en fin de resolution'
    END IF
    IF (CHECK_RESIDUAL_FLAG) THEN
      WRITE(6,'(A)') ' Verification : residu ||A*X-B||/||B|| affiche en fin de resolution'
    END IF
    IF (PROFILE_USAGE_FLAG) THEN
      WRITE(6,'(A)') ' Profilage RSS/CPU par etape active : usage_profile_rank<N>.csv par processus'
    END IF
    IF (DUMP_STRUCTURE_FLAG) THEN
      WRITE(6,'(A)') ' Structure BLR (dense/bas-rang) exportee : blr_structure_<etape>_rank<N>.csv par processus'
    END IF
    IF (LOAD_FACTORED_FLAG) THEN
      WRITE(6,'(A)') ' Reprise : matrice factorisee relue depuis disque (LOAD_FACTORED=1) -- '// &
           'assemblage et factorisation ignores'
      IF (LOAD_MATRIX_FLAG) WRITE(6,'(A)') ' (LOAD_MATRIX=1 ignore : LOAD_FACTORED est prioritaire)'
    ELSE IF (LOAD_MATRIX_FLAG) THEN
      WRITE(6,'(A)') ' Reprise : matrice relue depuis disque (LOAD_MATRIX=1) -- assemblage ignore'
    END IF
    IF (LOAD_RHS_FLAG) WRITE(6,'(A)') ' Reprise : second membre relu depuis disque (LOAD_RHS=1)'
    IF (SAVE_MATRIX_FLAG .OR. SAVE_RHS_FLAG .OR. SAVE_FACTORED_FLAG .OR. SAVE_SOLUTION_FLAG) &
         WRITE(6,'(A,A)') ' Sauvegarde(s) activee(s) dans : ', TRIM(CHECKPOINT_DIR)
    IF (NRHS_BATCH.GT.0 .AND. NRHS_BATCH.LT.NRHS .AND. .NOT.LOAD_RHS_FLAG) &
         WRITE(6,'(A,I0,A)') ' Second membre/resolution/sauvegarde par lots de ', NRHS_BATCH, ' colonnes (NRHS_BATCH)'
  END IF

  ! RSS de reference, capturee avant toute allocation -- isole le cout
  ! fixe du runtime deja charge (MPI/MKL) de la croissance due au code.
  PROF_T_START = MPI_WTIME()
  PROF_T_PREV = 0.0D0
  CALL CPU_TIME(PROF_CPU_PREV)
  CALL PROFILE_SAMPLE('debut')

  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  IF (LOAD_FACTORED_FLAG) THEN
    ! --- reprise : matrice DEJA FACTORISEE relue depuis disque -- assemblage
    !     ET factorisation entierement ignores (voir le commentaire de tete).
    CALL LOAD_BLR_MATRIX(TRIM(CHECKPOINT_DIR)//'/saved_factored', .TRUE., OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    T_FACTOR = 0.0D0
    N_SWAPS = -1
    CALL RECOMPUTE_COMPRESSION_STATS(N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (factorisation, charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
           ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
      WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (factorisation, charge) : ', &
           REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
  ELSE
    IF (LOAD_MATRIX_FLAG) THEN
      ! --- reprise : matrice ASSEMBLEE (pas encore factorisee) relue depuis
      !     disque -- assemblage ignore, factorisation normale ci-dessous.
      CALL LOAD_BLR_MATRIX(TRIM(CHECKPOINT_DIR)//'/saved_matrix', .FALSE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
      T_ASM = 0.0D0
      CALL RECOMPUTE_COMPRESSION_STATS(N_COMPRESSED, N_OFFDIAG, RANK_SUM)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage, charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage, charge) : ', &
             REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
    ELSE
      ! --- 1) assemblage BLR distribue ---
      ! SHOW_PROGRESS : barre affichee par le rang 0, approximative ici car
      ! l'assemblage n'est pas synchronise entre processus.
      CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                    N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., USE_RSVD=USE_RSVD_FLAG, &
                                    PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY)
      ! N_COMPRESSED/N_OFFDIAG sont deja des sommes globales (MPI_ALLREDUCE
      ! interne) -- un seul affichage (rang 0) suffit.
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        ! Taux REEL base sur le rang MOYEN (pas juste le compte de blocs,
        ! qui reste a 100% meme si chaque bloc n'economise presque rien).
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage) : ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
      IF (SAVE_MATRIX_FLAG) THEN
        CALL SAVE_BLR_MATRIX(TRIM(CHECKPOINT_DIR)//'/saved_matrix', .FALSE., OK)
        IF (.NOT. OK) THEN
          CALL MPI_FINALIZE(IERR)
          STOP 1
        END IF
      END IF
    END IF
    ! Empreinte theorique (blocs reellement stockes), pas mesuree via l'OS.
    CALL REPORT_MEMORY_MB(GRID, 'Matrice A_BLR (compressee)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
    ! RSS reelle (/proc/self/status) -- toujours >= la theorique (inclut
    ! workspaces SVD/QR, runtime MPI/MKL, allocateur...). Le ratio n'est
    ! pas constant : overhead fixe dominant a petite echelle (>60x mesure
    ! a N~2000), se stabilise vers 2.5-3x a N~20000 -- ce second regime
    ! sert de reference pour detecter une vraie anomalie (fuite, copie).
    CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
    CALL PROFILE_SAMPLE('assemblage')
    CALL DUMP_BLR_STRUCTURE('assembly')

    ! --- 2) factorisation LU, avec pivot confine au bloc diagonal ---
    ! SHOW_PROGRESS : ici K/NBLK est exact (boucle synchronisee par
    ! MPI_BCAST a chaque etape K).
    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                         N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                         SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE., &
                                         USE_MOVE_ALLOC=USE_MOVE_ALLOC_FLAG, &
                                         PROGRESS_CALLBACK=PROFILE_CALLBACK_FACTOR)
    IF (INFO .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
      WRITE(6,'(A,I0)') ' Nombre de permutations de lignes (pivot) : ', N_SWAPS
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
           ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
      WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (fin factorisation) : ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
    IF (SAVE_FACTORED_FLAG) THEN
      CALL SAVE_BLR_MATRIX(TRIM(CHECKPOINT_DIR)//'/saved_factored', .TRUE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
  END IF
  ! A_BLR est factorisee EN PLACE -- cette valeur peut differer de celle
  ! de l'assemblage (fill-in de la mise a jour de Schur).
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU (compresses)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('factorisation')
  CALL DUMP_BLR_STRUCTURE('factorization')

  ! --- 3+4) second membre B, resolution 2D et sauvegarde, PAR LOT d'au
  !     plus NRHS_BATCH colonnes (voir le commentaire de tete). NRHS_
  !     BATCH_EFF=NRHS (un seul lot) si NRHS_BATCH n'est pas positionne,
  !     LOAD_RHS=1, ou NRHS_BATCH>=NRHS.
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO

  IF (LOAD_RHS_FLAG) THEN
    ! LOAD_RHS ignore NRHS_BATCH : le decoupage stocke a la sauvegarde
    ! fait foi (comme LOAD_MATRIX/LOAD_FACTORED).
    CALL OPEN_RHS_CHECKPOINT_READ(TRIM(CHECKPOINT_DIR)//'/saved_rhs', LOAD_RHS_UNIT, NBATCHES, NRHS_BATCH_EFF, OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  ELSE
    NRHS_BATCH_EFF = NRHS_BATCH
    IF (NRHS_BATCH_EFF .LE. 0 .OR. NRHS_BATCH_EFF .GT. NRHS) NRHS_BATCH_EFF = NRHS
    NBATCHES = (NRHS + NRHS_BATCH_EFF - 1) / NRHS_BATCH_EFF
    IF (SAVE_RHS_FLAG) THEN
      CALL OPEN_RHS_CHECKPOINT_WRITE(TRIM(CHECKPOINT_DIR)//'/saved_rhs', SAVE_RHS_UNIT, &
                                      NBATCHES, NRHS_BATCH_EFF, MLOC, OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
  END IF
  IF (SAVE_SOLUTION_FLAG) THEN
    CALL OPEN_RHS_CHECKPOINT_WRITE(TRIM(CHECKPOINT_DIR)//'/saved_solution', SAVE_SOLUTION_UNIT, &
                                    NBATCHES, NRHS_BATCH_EFF, MLOC, OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  END IF

  WANT_FULL_X = (CHECK_ERROR_FLAG .OR. CHECK_RESIDUAL_FLAG)
  ! X_GLOBAL n'est necessaire que pour CHECK_ERROR/CHECK_RESIDUAL (compares
  ! sur toute sa largeur) -- rempli lot par lot (GATHER_2D_RHS directement
  ! dans X_GLOBAL(:,RSTART_B:REND_B)). Sans verification, l'allouer serait
  ! un cout memoire O(N*NRHS) pur gaspillage -- voir la branche ELSE plus
  ! bas (apercu 5 valeurs).
  IF (WANT_FULL_X) ALLOCATE(X_GLOBAL(N,NRHS))
  NPREVIEW = MIN(5, N)
  PREVIEW_LOCAL = (0.0_WP, 0.0_WP)

  RSS_MAX_B_MB = 0.0D0
  T_SOLVE = 0.0D0
  DO BATCH_IDX = 1, NBATCHES
    RSTART_B = (BATCH_IDX-1)*NRHS_BATCH_EFF + 1
    REND_B = MIN(BATCH_IDX*NRHS_BATCH_EFF, NRHS)
    WIDTH_B = REND_B - RSTART_B + 1
    CALL RHS_LOCAL_RANGE_MBM(WIDTH_B, GRID%PC, GRID%MYCOL, LOCAL_OFFSET_B, NRHS_LOC_B)
    GSTART_B = RSTART_B + LOCAL_OFFSET_B - 1
    GEND_B = GSTART_B + NRHS_LOC_B - 1

    WRITE(BATCH_STR,'(I4.4)') BATCH_IDX
    ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC_B,1)))
    IF (LOAD_RHS_FLAG) THEN
      IF (NTAILBLOC .GT. 0) THEN
        CALL LOAD_ROWBLOCKS(TRIM(CHECKPOINT_DIR)//'/saved_rhs_batch'//BATCH_STR, B_LOCAL, MLOC, NRHS_LOC_B, NTAILBLOC, OK)
      ELSE
        CALL LOAD_RHS_BATCH(LOAD_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
        OK = .TRUE.
      END IF
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    ELSE
      CALL FILL_B(A_BLR, GRID, NRHS_LOC_B, GSTART_B, B_LOCAL, MLOC)
      IF (SAVE_RHS_FLAG) THEN
        IF (NTAILBLOC .GT. 0) THEN
          CALL SAVE_ROWBLOCKS(TRIM(CHECKPOINT_DIR)//'/saved_rhs_batch'//BATCH_STR, B_LOCAL, MLOC, NRHS_LOC_B, NTAILBLOC, OK)
        ELSE
          CALL SAVE_RHS_BATCH(SAVE_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
          OK = .TRUE.
        END IF
        IF (.NOT. OK) THEN
          CALL MPI_FINALIZE(IERR)
          STOP 1
        END IF
      END IF
    END IF
    IF (NBATCHES.EQ.1) THEN
      CALL PROFILE_SAMPLE('second_membre')
    ELSE
      CALL PROFILE_SAMPLE(PCT_LABEL('second_membre', BATCH_IDX, NBATCHES))
    END IF

    ! Resolution 2D, AVEC pivot (applique la meme sequence de swaps que
    ! celle utilisee a la factorisation, deja diffusee a la grille) --
    ! MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D est deja generique sur la LARGEUR
    ! resolue (WIDTH_B/NRHS_LOC_B DU LOT, pas NRHS global) : aucune
    ! modification necessaire a la bibliotheque pour ce decoupage.
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0 = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, WIDTH_B, NRHS_LOC_B, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE., &
                                           PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE)
    T1 = MPI_WTIME()
    T_SOLVE = T_SOLVE + (T1 - T0)
    ! B_LOCAL contient maintenant X de ce lot (ecrase en place) -- meme
    ! routine de sauvegarde que le second membre.
    IF (SAVE_SOLUTION_FLAG) THEN
      IF (NTAILBLOC .GT. 0) THEN
        CALL SAVE_ROWBLOCKS(TRIM(CHECKPOINT_DIR)//'/saved_solution_batch'//BATCH_STR, B_LOCAL, MLOC, NRHS_LOC_B, NTAILBLOC, OK)
      ELSE
        CALL SAVE_RHS_BATCH(SAVE_SOLUTION_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
        OK = .TRUE.
      END IF
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
    IF (NBATCHES.EQ.1) THEN
      CALL PROFILE_SAMPLE('resolution')
    ELSE
      CALL PROFILE_SAMPLE(PCT_LABEL('resolution', BATCH_IDX, NBATCHES))
    END IF

    IF (WANT_FULL_X) THEN
      CALL GATHER_2D_RHS(GRID, N, NB, WIDTH_B, B_LOCAL, MAX(MLOC,1), NRHS_LOC_B, 0, X_GLOBAL(:,RSTART_B:REND_B))
    ELSE IF (GSTART_B .LE. 1 .AND. 1 .LE. GEND_B) THEN
      ! Ni CHECK_ERROR ni CHECK_RESIDUAL : inutile d'assembler X en entier.
      ! Recupere juste les 5 premieres valeurs de la colonne globale 1
      ! (forcement dans le premier lot) sans GATHER_2D_RHS -- chaque rang
      ! contribue sa valeur reelle s'il possede la ligne ET la colonne 1,
      ! zero sinon ; une somme reduite sur 5 complexes suffit.
      PREVIEW_RS = 1
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        PREVIEW_MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
        DO PREVIEW_GI = 1, NPREVIEW
          IF (PREVIEW_GI.GE.A_BLR%BSTART(I) .AND. PREVIEW_GI.LE.A_BLR%BEND(I)) THEN
            PREVIEW_LOCAL(PREVIEW_GI) = B_LOCAL(PREVIEW_RS + (PREVIEW_GI-A_BLR%BSTART(I)), 1-GSTART_B+1)
          END IF
        END DO
        PREVIEW_RS = PREVIEW_RS + PREVIEW_MI
      END DO
    END IF
    DEALLOCATE(B_LOCAL)
  END DO

  IF (LOAD_RHS_FLAG) THEN
    CLOSE(LOAD_RHS_UNIT)
  ELSE IF (SAVE_RHS_FLAG) THEN
    CLOSE(SAVE_RHS_UNIT)
  END IF
  IF (SAVE_SOLUTION_FLAG) CLOSE(SAVE_SOLUTION_UNIT)

  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Resolution : ', T_SOLVE, ' s'
  IF (CHECK_ERROR_FLAG .AND. MYID0.EQ.0) &
       WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant B=A*Xexact = ', RSS_MAX_B_MB, ' Mo'
  ! Empreinte reelle d'un lot (pas NRHS entier) -- pic de RAM reel,
  ! borne par NRHS_BATCH.
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B (par lot)', &
       8_8*INT(MLOC,8)*INT(RHS_NUMROC_MBM(NRHS_BATCH_EFF,GRID%PC,GRID%MYCOL),8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres resolution)', GET_RSS_KB()*1024_8)

  IF (WANT_FULL_X) THEN
    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
      DO I = 1, MIN(5, N)
        WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
      END DO
      IF (CHECK_ERROR_FLAG) THEN
        ERR = COMPUTE_ERR(X_GLOBAL)
        WRITE(6,'(A,ES12.4)') ' Erreur relative ||X-Xexact|| / ||Xexact|| = ', ERR
        CALL PROFILE_SAMPLE('verif_erreur')
      END IF
    END IF

    ! Appel COLLECTIF (MPI_BCAST/MPI_ALLREDUCE a l'interieur) -- tous les
    ! rangs doivent l'appeler, contrairement a COMPUTE_ERR ci-dessus.
    IF (CHECK_RESIDUAL_FLAG) THEN
      RSS_MAX_RESID_MB = 0.0D0
      RESID = COMPUTE_RESIDUAL(X_GLOBAL)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant le residu = ', RSS_MAX_RESID_MB, ' Mo'
        WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
      END IF
      CALL PROFILE_SAMPLE('verif_residu')
    END IF
  ELSE
    CALL MPI_REDUCE(PREVIEW_LOCAL, PREVIEW_GLOBAL, NPREVIEW, MPI_COMPLEX, MPI_SUM, 0, GRID%COMM, IERR)
    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
      DO I = 1, NPREVIEW
        WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(PREVIEW_GLOBAL(I)), ' + i*', AIMAG(PREVIEW_GLOBAL(I))
      END DO
      WRITE(6,'(A)') ' (X complet non assemble -- CHECK_ERROR=0 et CHECK_RESIDUAL=0, memoire economisee)'
    END IF
  END IF

  CALL PROFILE_SAMPLE('fin')
  CALL PROFILE_WRITE_CSV()

  IF (ALLOCATED(X_GLOBAL)) DEALLOCATE(X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Coupe "--cle=valeur" en KEYOUT="cle"/VALOUT="valeur" (le "--" initial
  ! est deja verifie par l'appelant). Sans "=" (option booleenne nue),
  ! VALOUT="1" -- equivalent a "--cle=1".
  SUBROUTINE SPLIT_KV(ARGSTR, KEYOUT, VALOUT)
    CHARACTER(LEN=*), INTENT(IN)  :: ARGSTR
    CHARACTER(LEN=*), INTENT(OUT) :: KEYOUT, VALOUT
    INTEGER :: EQPOS
    EQPOS = INDEX(ARGSTR, '=')
    IF (EQPOS .EQ. 0) THEN
      KEYOUT = ARGSTR(3:)
      VALOUT = '1'
    ELSE
      KEYOUT = ARGSTR(3:EQPOS-1)
      VALOUT = ARGSTR(EQPOS+1:)
    END IF
  END SUBROUTINE SPLIT_KV

  ! Message d'aide -- affiche sur --help/-h, ligne de commande invalide,
  ! ou option/argument non reconnu.
  SUBROUTINE PRINT_USAGE()
    WRITE(6,'(A)') 'Usage : mpirun -np PR*PC ./example_blr_mpi_pivot PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  PR PC          grille de processus MPI (obligatoire, PR*PC = nb de processus lances)'
    WRITE(6,'(A)') '  N NB NRHS TOL  positionnels optionnels, DANS CET ORDRE (defauts 4000 200 4 1e-4)'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Options nommees (apres les positionnels ci-dessus, dans N''IMPORTE QUEL ORDRE ;'
    WRITE(6,'(A)') '--cle=valeur, ou --cle seul pour une option booleenne, equivalent a --cle=1) :'
    WRITE(6,'(A)') '  --n=ENTIER              equivalent au positionnel N (utile pour ne changer que NRHS/TOL...)'
    WRITE(6,'(A)') '  --nb=ENTIER             equivalent au positionnel NB'
    WRITE(6,'(A)') '  --nrhs=ENTIER           equivalent au positionnel NRHS'
    WRITE(6,'(A)') '  --tol=REEL              equivalent au positionnel TOL'
    WRITE(6,'(A)') '  --use-move-alloc=0|1    (defaut 0) transfert MOVE_ALLOC au lieu de copie (Schur)'
    WRITE(6,'(A)') '  --eps=REEL              (defaut 0.0) bruit plein-rang (compressibilite)'
    WRITE(6,'(A)') '  --ratio=REEL            (defaut 0.333) COMPRESS_WORTHWHILE_RATIO'
    WRITE(6,'(A)') '  --use-rsvd=0|1          (defaut 1) SVD randomisee adaptative ; 0 = CGESVD complete'
    WRITE(6,'(A)') '  --check-error=0|1       (defaut 0) verification B=A*Xexact'
    WRITE(6,'(A)') '  --check-residual=0|1    (defaut 0) verification ||A*X-B||/||B||'
    WRITE(6,'(A)') '  --profile-usage=0|1     (defaut 0) trace RSS/CPU par etape (usage_profile_rank<N>.csv)'
    WRITE(6,'(A)') '  --dump-structure=0|1    (defaut 0) export structure BLR (blr_structure_<etape>_rank<N>.csv)'
    WRITE(6,'(A)') '  --save-matrix=0|1       (defaut 0) sauvegarde la matrice assemblee'
    WRITE(6,'(A)') '  --save-rhs=0|1          (defaut 0) sauvegarde le second membre'
    WRITE(6,'(A)') '  --save-factored=0|1     (defaut 0) sauvegarde la matrice factorisee'
    WRITE(6,'(A)') '  --save-solution=0|1     (defaut 0) sauvegarde la solution'
    WRITE(6,'(A)') '  --load-matrix=0|1       (defaut 0) relit la matrice sauvegardee (saute l''assemblage)'
    WRITE(6,'(A)') '  --load-factored=0|1     (defaut 0) relit la matrice factorisee (saute assemblage+factorisation,'
    WRITE(6,'(A)') '                          prioritaire sur --load-matrix)'
    WRITE(6,'(A)') '  --load-rhs=0|1          (defaut 0) relit le second membre sauvegarde'
    WRITE(6,'(A)') '  --checkpoint-dir=CHEMIN (defaut .) repertoire des fichiers de checkpoint (doit deja exister)'
    WRITE(6,'(A)') '  --nrhs-batch=ENTIER     (defaut 0 = pas de decoupage) colonnes RHS traitees par lot'
    WRITE(6,'(A)') '  --ntailbloc=ENTIER      (defaut 0 = desactive) SAVE_RHS/SAVE_SOLUTION/LOAD_RHS ecrivent/'
    WRITE(6,'(A)') '                          relisent par blocs de ce nombre de lignes, un fichier par bloc'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Exemple :'
    WRITE(6,'(A)') '  mpirun -np 4 ./example_blr_mpi_pivot 2 2 4000 200 8 1e-4 \'
    WRITE(6,'(A)') '      --check-error --save-factored --checkpoint-dir=ckpt'
  END SUBROUTINE PRINT_USAGE

  ! Recalcule N_COMPRESSED/N_OFFDIAG/RANK_SUM depuis A_BLR%OFFBLK (meme
  ! calcul qu'en interne dans MPI_ASSEMBLE_BLR_MATRIX/LU_FACTOR_PIVOT) --
  ! necessaire apres LOAD_MATRIX/LOAD_FACTORED, qui sautent justement
  ! l'appel qui produirait normalement ces stats.
  SUBROUTINE RECOMPUTE_COMPRESSION_STATS(NC_OUT, NO_OUT, RS_OUT)
    INTEGER, INTENT(OUT) :: NC_OUT, NO_OUT
    INTEGER(8), INTENT(OUT) :: RS_OUT
    INTEGER :: li2, lj2, I2, J2, NC_LOCAL, NO_LOCAL, IERR2
    INTEGER(8) :: RS_LOCAL
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8
    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        IF (I2 .EQ. J2) CYCLE
        NO_LOCAL = NO_LOCAL + 1
        IF (A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK) THEN
          NC_LOCAL = NC_LOCAL + 1
          RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li2,lj2)%RANK
        ELSE
          RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li2,lj2)%M, A_BLR%OFFBLK(li2,lj2)%N)
        END IF
      END DO
    END DO
    CALL MPI_ALLREDUCE(NC_LOCAL, NC_OUT, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(NO_LOCAL, NO_OUT, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(RS_LOCAL, RS_OUT, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR2)
  END SUBROUTINE RECOMPUTE_COMPRESSION_STATS

  ! Ecrit A_BLR (SAVE_MATRIX/SAVE_FACTORED) -- un fichier par processus,
  ! STREAM non formate, sans communication MPI. En-tete d'abord (verifie
  ! par LOAD_BLR_MATRIX), puis les blocs diagonaux localement possedes,
  ! puis tous les blocs hors-diagonale locaux (toujours "possedes" par
  ! construction du decoupage 2D block-cyclique). WITH_IPIV n'ecrit IPIV
  ! que si elle est allouee (apres factorisation).
  SUBROUTINE SAVE_BLR_MATRIX(PREFIX, WITH_IPIV, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: U, li2, lj2, I2, J2, ND, MI, IOS, IERR2
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    OK_LOCAL = (IOS .EQ. 0)
    IF (.NOT. OK_LOCAL) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible de creer ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- le repertoire CHECKPOINT_DIR doit deja exister.'
    END IF
    ! Un seul rang qui ne peut pas ecrire doit arreter tous les rangs --
    ! sinon les autres continueraient seuls vers un appel MPI collectif
    ! que ce rang n'atteindra jamais.
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U, STATUS='DELETE')
      RETURN
    END IF

    ! Marqueur de format (N > 0 dans l'ancien format, sans marqueur) +
    ! distribution des blocs-colonnes (CMAP, voir BLR_COL_OWNER).
    WRITE(U) -7411, A_BLR%CMAP
    WRITE(U) N, NB, A_BLR%NBLK, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, &
         A_BLR%NBLKR_LOC, A_BLR%NBLKC_LOC

    ND = 0
    DO I2 = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(I2)%D)) ND = ND + 1
    END DO
    WRITE(U) ND
    DO I2 = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(I2)%D)) THEN
        MI = SIZE(A_BLR%DIAGBLK(I2)%D, 1)
        WRITE(U) I2, MI
        WRITE(U) A_BLR%DIAGBLK(I2)%D
      END IF
    END DO

    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        ! Position diagonale d'OFFBLK : jamais compressee (I.EQ.J), rien
        ! a ecrire.
        IF (I2 .EQ. J2) CYCLE
        WRITE(U) I2, J2, A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK, A_BLR%OFFBLK(li2,lj2)%M, &
             A_BLR%OFFBLK(li2,lj2)%N, A_BLR%OFFBLK(li2,lj2)%RANK
        IF (A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK) THEN
          WRITE(U) A_BLR%OFFBLK(li2,lj2)%UFAC
          WRITE(U) A_BLR%OFFBLK(li2,lj2)%VFAC
        ELSE
          WRITE(U) A_BLR%OFFBLK(li2,lj2)%DENSE
        END IF
      END DO
    END DO

    IF (WITH_IPIV) WRITE(U) A_BLR%IPIV
    CLOSE(U)
  END SUBROUTINE SAVE_BLR_MATRIX

  ! Relit A_BLR ecrite par SAVE_BLR_MATRIX -- reconstruit BSTART/BEND/
  ! NBLKR_LOC/NBLKC_LOC/DIAGBLK/OFFBLK comme le ferait MPI_ASSEMBLE_BLR_
  ! MATRIX, donc la factorisation qui suit ne voit aucune difference.
  ! OK=.FALSE. (tous les rangs s'arretent) si le fichier est absent ou si
  ! son en-tete ne correspond pas au run courant.
  SUBROUTINE LOAD_BLR_MATRIX(PREFIX, WITH_IPIV, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: U, li2, lj2, I2, J2, ND, K2, MI, IOS, IERR2
    INTEGER :: F_CMAP, F_N, F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
    INTEGER :: R_I, R_J, R_M, R_N, R_RANK
    LOGICAL :: R_LOWRANK, OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    OK_LOCAL = .TRUE.
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- avez-vous sauvegarde avec SAVE_MATRIX/SAVE_FACTORED=1 auparavant ?'
      OK_LOCAL = .FALSE.
    ELSE
      ! Nouveau format : marqueur -7411 puis CMAP ; ancien format (premier
      ! entier = N > 0) : distribution cyclique pure (CMAP=0).
      READ(U) F_N
      F_CMAP = 0
      IF (F_N .EQ. -7411) READ(U) F_CMAP, F_N
      READ(U) F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
      IF (F_N.NE.N .OR. F_NB.NE.NB .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NB/PR/PC ou grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    ! Un seul rang incompatible doit arreter tous les rangs -- sinon les
    ! rangs corrects continueraient seuls vers un appel MPI collectif :
    ! blocage garanti plutot qu'un arret net.
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U)
      RETURN
    END IF

    A_BLR%N = N; A_BLR%NB = NB
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%CMAP = F_CMAP
    A_BLR%NBLKC_LOC = BLR_NBLKC_LOC(A_BLR%CMAP, A_BLR%NBLK, GRID)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))

    READ(U) ND
    DO K2 = 1, ND
      READ(U) I2, MI
      ALLOCATE(A_BLR%DIAGBLK(I2)%D(MI,MI))
      READ(U) A_BLR%DIAGBLK(I2)%D
    END DO

    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        ! Meme position diagonale sautee a l'ecriture -- rien a relire,
        ! l'element reste dans son etat par defaut.
        IF (I2 .EQ. J2) CYCLE
        READ(U) R_I, R_J, R_LOWRANK, R_M, R_N, R_RANK
        A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK = R_LOWRANK
        A_BLR%OFFBLK(li2,lj2)%M = R_M
        A_BLR%OFFBLK(li2,lj2)%N = R_N
        A_BLR%OFFBLK(li2,lj2)%RANK = R_RANK
        IF (R_LOWRANK) THEN
          ALLOCATE(A_BLR%OFFBLK(li2,lj2)%UFAC(R_M,R_RANK), A_BLR%OFFBLK(li2,lj2)%VFAC(R_RANK,R_N))
          READ(U) A_BLR%OFFBLK(li2,lj2)%UFAC
          READ(U) A_BLR%OFFBLK(li2,lj2)%VFAC
        ELSE
          ALLOCATE(A_BLR%OFFBLK(li2,lj2)%DENSE(R_M,R_N))
          READ(U) A_BLR%OFFBLK(li2,lj2)%DENSE
        END IF
      END DO
    END DO

    IF (WITH_IPIV) THEN
      ALLOCATE(A_BLR%IPIV(NB, A_BLR%NBLK))
      READ(U) A_BLR%IPIV
    END IF
    CLOSE(U)
  END SUBROUTINE LOAD_BLR_MATRIX

  ! Ouvre un checkpoint RHS/solution en ECRITURE et ecrit son en-tete une
  ! seule fois (N, NRHS, grille, MLOC, NRHS_BATCH_HDR, NBATCHES_HDR) --
  ! reste ouvert pendant toute la boucle de lots, ferme par l'appelant.
  ! NRHS_BATCH_HDR/NBATCHES_HDR permettent a LOAD_RHS de reproduire
  ! exactement le meme decoupage, sans dependre du NRHS_BATCH du run qui
  ! relit.
  SUBROUTINE OPEN_RHS_CHECKPOINT_WRITE(PREFIX, UNIT, NBATCHES_HDR, NRHS_BATCH_HDR, MLOC_HDR, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER, INTENT(OUT) :: UNIT
    INTEGER, INTENT(IN) :: NBATCHES_HDR, NRHS_BATCH_HDR, MLOC_HDR
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: IOS, IERR2
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=UNIT, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    OK_LOCAL = (IOS .EQ. 0)
    IF (.NOT. OK_LOCAL) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible de creer ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- le repertoire CHECKPOINT_DIR doit deja exister.'
    END IF
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(UNIT, STATUS='DELETE')
      RETURN
    END IF
    WRITE(UNIT) N, NRHS, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, MLOC_HDR, NRHS_BATCH_HDR, NBATCHES_HDR
  END SUBROUTINE OPEN_RHS_CHECKPOINT_WRITE

  ! Ecrit un lot (NLOC_B colonnes locales de BDATA) -- appelee pour
  ! SAVE_RHS (apres FILL_B) et SAVE_SOLUTION (apres la resolution) --
  ! meme routine, deux fichiers/unites distincts.
  SUBROUTINE SAVE_RHS_BATCH(UNIT, NLOC_B, BDATA, ML)
    INTEGER, INTENT(IN) :: UNIT, NLOC_B, ML
    COMPLEX(WP), INTENT(IN) :: BDATA(MAX(ML,1),MAX(NLOC_B,1))
    WRITE(UNIT) NLOC_B
    IF (ML.GT.0 .AND. NLOC_B.GT.0) WRITE(UNIT) BDATA(1:ML,1:NLOC_B)
  END SUBROUTINE SAVE_RHS_BATCH

  ! Ouvre un checkpoint RHS en LECTURE, verifie son en-tete contre le run
  ! courant, et retourne NBATCHES_OUT/NRHS_BATCH_OUT lus depuis le fichier
  ! -- c'est ce decoupage, pas le NRHS_BATCH de ce run, qui pilote la
  ! boucle de lots. Meme logique OK/MPI_ALLREDUCE que LOAD_BLR_MATRIX.
  SUBROUTINE OPEN_RHS_CHECKPOINT_READ(PREFIX, UNIT, NBATCHES_OUT, NRHS_BATCH_OUT, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER, INTENT(OUT) :: UNIT, NBATCHES_OUT, NRHS_BATCH_OUT
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: IOS, IERR2
    INTEGER :: F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    OK_LOCAL = .TRUE.
    NBATCHES_OUT = 0; NRHS_BATCH_OUT = 0
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=UNIT, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- avez-vous sauvegarde avec SAVE_RHS=1 auparavant ?'
      OK_LOCAL = .FALSE.
    ELSE
      READ(UNIT) F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, NRHS_BATCH_OUT, NBATCHES_OUT
      IF (F_N.NE.N .OR. F_NRHS.NE.NRHS .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0 .OR. F_MLOC.NE.MLOC) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NRHS/PR/PC/grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(UNIT)
    END IF
  END SUBROUTINE OPEN_RHS_CHECKPOINT_READ

  ! Relit un lot (NLOC_B colonnes locales attendues, connues de l'appelant
  ! via le meme decoupage RHS_LOCAL_RANGE_MBM que la sauvegarde) dans
  ! BDATA.
  SUBROUTINE LOAD_RHS_BATCH(UNIT, NLOC_B, BDATA, ML)
    INTEGER, INTENT(IN) :: UNIT, NLOC_B, ML
    COMPLEX(WP), INTENT(OUT) :: BDATA(MAX(ML,1),MAX(NLOC_B,1))
    INTEGER :: NLOC_B_FILE
    READ(UNIT) NLOC_B_FILE
    IF (ML.GT.0 .AND. NLOC_B_FILE.GT.0) READ(UNIT) BDATA(1:ML,1:NLOC_B_FILE)
  END SUBROUTINE LOAD_RHS_BATCH

  ! E/S rapide par blocs de lignes (NTAILBLOC, voir le commentaire de tete)
  ! -- ecrit DATA(1:ML,1:NCOLS) (second membre ou solution, jamais le
  ! facteur LU) en tranches de NTAILBLOC lignes globales, chaque tranche
  ! dans son propre fichier : un seul WRITE par bloc, chacun lisible
  ! independamment des autres. Le dernier bloc peut etre plus petit si ML
  ! n'est pas un multiple exact. Meme convention OK/MPI_ALLREDUCE que les
  ! autres routines de checkpoint.
  SUBROUTINE SAVE_ROWBLOCKS(PREFIX, DATA, ML, NCOLS, NTAILBLOC_ARG, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER, INTENT(IN)          :: ML, NCOLS, NTAILBLOC_ARG
    COMPLEX(WP), INTENT(IN)      :: DATA(MAX(ML,1),MAX(NCOLS,1))
    LOGICAL, INTENT(OUT)         :: OK
    INTEGER :: U, RSTART, REND, NROWS, IBLOC, IOS, IERR2
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR, BSTR
    CHARACTER(LEN=220) :: FNAME
    OK_LOCAL = .TRUE.
    WRITE(RSTR,'(I6.6)') MYID0
    IBLOC = 0
    RSTART = 1
    DO WHILE (RSTART .LE. ML)
      REND = MIN(RSTART+NTAILBLOC_ARG-1, ML)
      NROWS = REND - RSTART + 1
      IBLOC = IBLOC + 1
      WRITE(BSTR,'(I6.6)') IBLOC
      FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'_bloc'//TRIM(BSTR)//'.dat'
      OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
           STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible de creer ', TRIM(FNAME), &
             ' (rang ', MYID0, ') -- le repertoire CHECKPOINT_DIR doit deja exister.'
        OK_LOCAL = .FALSE.
        EXIT
      END IF
      WRITE(U) NROWS, NCOLS
      WRITE(U) DATA(RSTART:REND, 1:NCOLS)
      CLOSE(U)
      RSTART = RSTART + NTAILBLOC_ARG
    END DO
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
  END SUBROUTINE SAVE_ROWBLOCKS

  ! Relit ce que SAVE_ROWBLOCKS a ecrit -- NTAILBLOC_ARG doit correspondre
  ! exactement a celui utilise a la sauvegarde ; chaque bloc verifie en
  ! plus sa propre taille avant de lire, donc un NTAILBLOC different est
  ! detecte proprement plutot que de charger des donnees decalees.
  SUBROUTINE LOAD_ROWBLOCKS(PREFIX, DATA, ML, NCOLS, NTAILBLOC_ARG, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER, INTENT(IN)          :: ML, NCOLS, NTAILBLOC_ARG
    COMPLEX(WP), INTENT(OUT)     :: DATA(MAX(ML,1),MAX(NCOLS,1))
    LOGICAL, INTENT(OUT)         :: OK
    INTEGER :: U, RSTART, REND, NROWS, IBLOC, IOS, IERR2, F_NROWS, F_NCOLS
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR, BSTR
    CHARACTER(LEN=220) :: FNAME
    OK_LOCAL = .TRUE.
    WRITE(RSTR,'(I6.6)') MYID0
    IBLOC = 0
    RSTART = 1
    DO WHILE (RSTART .LE. ML)
      REND = MIN(RSTART+NTAILBLOC_ARG-1, ML)
      NROWS = REND - RSTART + 1
      IBLOC = IBLOC + 1
      WRITE(BSTR,'(I6.6)') IBLOC
      FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'_bloc'//TRIM(BSTR)//'.dat'
      OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
           STATUS='OLD', ACTION='READ', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
             ' (rang ', MYID0, ') -- meme NTAILBLOC qu''a la sauvegarde ?'
        OK_LOCAL = .FALSE.
        EXIT
      END IF
      READ(U) F_NROWS, F_NCOLS
      IF (F_NROWS.NE.NROWS .OR. F_NCOLS.NE.NCOLS) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : bloc ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, ') -- NTAILBLOC different de la sauvegarde.'
        OK_LOCAL = .FALSE.
        CLOSE(U)
        EXIT
      END IF
      READ(U) DATA(RSTART:REND, 1:NCOLS)
      CLOSE(U)
      RSTART = RSTART + NTAILBLOC_ARG
    END DO
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
  END SUBROUTINE LOAD_ROWBLOCKS

  ! Taux de compression REEL (stockage bas-rang effectif / dense
  ! equivalent) -- pas juste 2*NB*RANK_SUM/(NB^2*N_OFFDIAG) : RANK_SUM
  ! melange le rang reel des blocs compresses et un "rang virtuel"=NB pour
  ! chaque bloc dense (convention de comptage interne, pas un vrai
  ! stockage UV), qui occupe NB^2 elements et non 2*NB^2 (une formule
  ! naive donnerait a tort 200% quand tout est dense). La contribution des
  ! blocs denses est donc retiree de RANK_SUM avant de les recompter a
  ! leur cout reel NB^2.
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  ! Echantillon RSS/CPU pour l'etape STAGE_NAME (PROFILE_USAGE) -- no-op
  ! si desactive. CPU_TIME mesure le temps processeur (tous threads MKL/
  ! OpenMP confondus), MPI_WTIME le temps horloge ecoule -- leur rapport
  ! sur l'intervalle donne le nombre moyen de coeurs actifs, en %, comme
  ! `top`. Chaque processus echantillonne son propre etat, exploitable
  ! point par point dans plot_usage_profile.py.
  SUBROUTINE PROFILE_SAMPLE(STAGE_NAME)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_NAME
    REAL :: CPU_NOW
    DOUBLE PRECISION :: T_NOW, DT_WALL, DT_CPU
    IF (.NOT. PROFILE_USAGE_FLAG) RETURN
    IF (PROF_N .GE. MAX_PROFILE_SAMPLES) RETURN
    CALL CPU_TIME(CPU_NOW)
    T_NOW = MPI_WTIME() - PROF_T_START
    DT_WALL = T_NOW - PROF_T_PREV
    DT_CPU = REAL(CPU_NOW - PROF_CPU_PREV, 8)
    PROF_N = PROF_N + 1
    PROF_STAGE(PROF_N) = STAGE_NAME
    PROF_T(PROF_N) = T_NOW
    PROF_RSS_MB(PROF_N) = REAL(GET_RSS_KB(),8) / 1024.0D0
    IF (DT_WALL .GT. 1.0D-6) THEN
      PROF_CPU_PCT(PROF_N) = 100.0D0 * DT_CPU / DT_WALL
    ELSE
      PROF_CPU_PCT(PROF_N) = 0.0D0
    END IF
    PROF_T_PREV = T_NOW
    PROF_CPU_PREV = CPU_NOW
  END SUBROUTINE PROFILE_SAMPLE

  ! Etiquette "STAGE@NN%" pour un echantillon intra-etape. Le separateur
  ! "@" (absent des etiquettes de bornes comme 'assemblage') permet a
  ! plot_usage_profile.py de ne tracer un repere que sur les bornes.
  FUNCTION PCT_LABEL(STAGE, NDONE, NTOTAL) RESULT(LBL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    CHARACTER(LEN=40) :: LBL
    WRITE(LBL,'(A,A,I0,A)') TRIM(STAGE), '@', (100*NDONE)/MAX(NTOTAL,1), '%'
  END FUNCTION PCT_LABEL

  ! Les 3 callbacks suivants ont l'interface STAGE_PROGRESS_CALLBACK_IFACE
  ! -- appeles par la bibliotheque a chaque palier de 5% dans l'assemblage/
  ! la factorisation/la resolution (boucles internes a mpi_blr_matrix_
  ! mod.f90, invisibles d'ici sans ce mecanisme). Procedures internes
  ! passees comme argument (Fortran 2008) : acces par host-association a
  ! PROFILE_SAMPLE sans variable de module supplementaire.
  SUBROUTINE PROFILE_CALLBACK_ASSEMBLY(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    IF (.NOT. PROFILE_USAGE_FLAG) RETURN
    CALL PROFILE_SAMPLE(PCT_LABEL('assemblage', NDONE, NTOTAL))
  END SUBROUTINE PROFILE_CALLBACK_ASSEMBLY

  SUBROUTINE PROFILE_CALLBACK_FACTOR(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    IF (.NOT. PROFILE_USAGE_FLAG) RETURN
    CALL PROFILE_SAMPLE(PCT_LABEL('factorisation', NDONE, NTOTAL))
  END SUBROUTINE PROFILE_CALLBACK_FACTOR

  SUBROUTINE PROFILE_CALLBACK_SOLVE(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    IF (.NOT. PROFILE_USAGE_FLAG) RETURN
    CALL PROFILE_SAMPLE(PCT_LABEL('resolution', NDONE, NTOTAL))
  END SUBROUTINE PROFILE_CALLBACK_SOLVE

  ! Ecrit les echantillons de ce processus dans son propre CSV -- un
  ! fichier par rang, pas de rassemblement MPI.
  SUBROUTINE PROFILE_WRITE_CSV()
    INTEGER :: U, K
    CHARACTER(LEN=64) :: FNAME
    IF (.NOT. PROFILE_USAGE_FLAG) RETURN
    WRITE(FNAME,'(A,I0,A)') 'usage_profile_rank', MYID0, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'rank,stage,t_wall_s,rss_mb,cpu_pct'
    DO K = 1, PROF_N
      WRITE(U,'(I0,A,A,A,F12.4,A,F12.3,A,F8.2)') MYID0, ',', TRIM(PROF_STAGE(K)), ',', &
           PROF_T(K), ',', PROF_RSS_MB(K), ',', PROF_CPU_PCT(K)
    END DO
    CLOSE(U)
  END SUBROUTINE PROFILE_WRITE_CSV

  ! Ecrit la nature de chaque bloc local de A_BLR (DUMP_STRUCTURE) -- un
  ! fichier par processus, jamais rassemble ici : plot_blr_structure.py
  ! reconstitue la grille globale en lisant tous les fichiers d'un meme
  ! STAGE_LABEL.
  SUBROUTINE DUMP_BLR_STRUCTURE(STAGE_LABEL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_LABEL
    INTEGER :: U, li, lj, I, J
    CHARACTER(LEN=96) :: FNAME
    IF (.NOT. DUMP_STRUCTURE_FLAG) RETURN
    WRITE(FNAME,'(A,A,A,I0,A)') 'blr_structure_', TRIM(STAGE_LABEL), '_rank', MYID0, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'block_row,block_col,kind,rank'
    DO I = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(I)%D)) THEN
        WRITE(U,'(I0,A,I0,A,I0)') I, ',', I, ',diagonal,', A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END IF
    END DO
    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = BLR_COL_GLOBAL(A_BLR%CMAP, lj, GRID)
        IF (I .EQ. J) CYCLE
        IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
          WRITE(U,'(I0,A,I0,A,I0)') I, ',', J, ',lowrank,', A_BLR%OFFBLK(li,lj)%RANK
        ELSE
          WRITE(U,'(I0,A,I0,A,I0)') I, ',', J, ',dense,', A_BLR%BEND(J)-A_BLR%BSTART(J)+1
        END IF
      END DO
    END DO
    CLOSE(U)
  END SUBROUTINE DUMP_BLR_STRUCTURE

  ! Hash deterministe (I,J) -> [-1,1) -- pas de RNG a etat, recalculable
  ! independamment par n'importe quel processus pour n'importe quel bloc.
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  ! A(GI,GJ) = 50 sur une anti-diagonale confinee au bloc diagonal de GI
  ! (GJ_STRONG, inchangee quel que soit EPS), sinon un noyau de Cauchy
  ! amorti 1/(1+|GJ_STRONG-GJ|) + un bruit plein-rang EPS*HASH_NOISE.
  !
  ! L'anti-diagonale est confinee au bloc LOCAL de GI (plutot qu'une
  ! anti-diagonale globale GJ=N+1-GI) : une anti-diagonale globale serait
  ! injuste envers le pivotage confine au bloc diagonal (l'element fort
  ! pourrait tomber hors de portee du pivot local -> systeme
  ! catastrophiquement mal conditionne). Extrait en fonction independante
  ! pour etre reutilisable par FILL_B en mode CHECK_ERROR.
  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: DF
    INTEGER :: BSTART0, BEND0, GJ_STRONG
    BSTART0 = ((GI-1)/NB)*NB + 1
    BEND0 = MIN(BSTART0+NB-1, N)
    GJ_STRONG = BEND0 - (GI-BSTART0)
    IF (GJ .EQ. GJ_STRONG) THEN
      LOWRANK_A = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
    ELSE
      DF = 1.0_WP / (1.0_WP + ABS(REAL(GJ_STRONG-GJ,WP)))
      IF (EPS .NE. 0.0_WP) DF = DF + EPS*HASH_NOISE(GI,GJ)
      LOWRANK_A = CMPLX(DF, 0.3_WP*DF, KIND=WP)
    END IF
  END FUNCTION LOWRANK_A

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        BLK(II-ISTART+1, JJ-JSTART+1) = LOWRANK_A(II, JJ)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  ! Solution exacte deterministe, connue par formule -- utilisee seulement
  ! si CHECK_ERROR=1, pour construire B=A*Xexact (FILL_B) et mesurer
  ! l'erreur finale (COMPUTE_ERR).
  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  ! Erreur relative ||X-Xexact||/||Xexact|| -- rang 0 uniquement (X_GLOBAL
  ! y est deja rassemble, aucune allocation globale supplementaire).
  REAL(WP) FUNCTION COMPUTE_ERR(XG)
    COMPLEX(WP), INTENT(IN) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: B_REF(:,:)
    INTEGER :: GI, R
    ALLOCATE(B_REF(N,NRHS))
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = XEX_FORMULA(GI,R)
      END DO
    END DO
    COMPUTE_ERR = SQRT(SUM(ABS(XG-B_REF)**2)) / SQRT(SUM(ABS(B_REF)**2))
    DEALLOCATE(B_REF)
  END FUNCTION COMPUTE_ERR

  ! Valeur de B a l'entree (GI,R) pour le cas CHECK_ERROR=0 -- source
  ! unique partagee par FILL_B et COMPUTE_RESIDUAL, garantit que le residu
  ! mesure exactement le B utilise. Le cas CHECK_ERROR=1 (B=A*Xexact) est
  ! traite separement via CGEMM direct, jamais via cette fonction. Pour un
  ! probleme reel, cette fonction serait remplacee par la formule/lecture
  ! de B de l'utilisateur.
  COMPLEX(WP) FUNCTION B_ENTRY(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    B_ENTRY = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP), KIND=WP)
  END FUNCTION B_ENTRY

  ! Residu ||A*X-B||/||B|| -- cas general : ne suppose aucune connaissance
  ! de la solution exacte, utilisable meme si A et B viennent d'un vrai
  ! probleme utilisateur.
  !
  ! Deux optimisations vs une triple boucle scalaire naive (~12x plus
  ! lente, mesure a N=8000) :
  !   1) BLAS (CGEMM) : A materialisee par blocs de NB lignes (via
  !      GET_BLOCK, jamais la matrice entiere) puis multipliee par X (et
  !      Xexact si CHECK_ERROR=1, en un seul CGEMM via XCOMBINED pour
  !      eviter de reconstruire deux fois le meme bloc de A).
  !   2) Distribution MPI : les N lignes sont partagees entre tous les
  !      rangs (pas seulement le rang 0), les normes partielles reduites
  !      par MPI_ALLREDUCE. X doit d'abord etre diffuse (MPI_BCAST) car
  !      A*X(GI,:) necessite tout X.
  ! Devient un appel COLLECTIF (voir le site d'appel, hors du bloc
  ! "IF (MYID0.EQ.0)"). Cout memoire : un tampon transitoire NB x N par
  ! rang au lieu d'un tableau (N,NRHS) complet -- l'avantage memoire
  ! reste entier.
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG)
    COMPLEX(WP), INTENT(INOUT) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: BLKROWS(:,:), RESBLK(:,:), XCOMBINED(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, NROWS_MINE, NTOTAL_LOCAL
    INTEGER :: RSTART_CHUNK, REND_CHUNK, NCHUNK, GI, R, NCOLS_MULT, IERR2, NDONE, PBAR_STATE
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2
    COMPLEX(WP) :: BVAL

    ! X n'est rassemble que sur le rang 0 (GATHER_2D_RHS) -- les autres
    ! rangs n'ont rien de valide avant cette diffusion.
    CALL MPI_BCAST(XG, N*NRHS, MPI_COMPLEX, 0, GRID%COMM, IERR2)

    ! Repartition simple (pas bloc-cyclique, inutile pour ce calcul
    ! ponctuel de verification) : chaque rang recoit un intervalle
    ! contigu de lignes.
    ROWS_PER_PROC = (N + GRID%NPROCS - 1) / GRID%NPROCS
    RSTART_MINE = GRID%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GRID%MYID+1)*ROWS_PER_PROC, N)
    NROWS_MINE = MAX(REND_MINE - RSTART_MINE + 1, 0)
    NTOTAL_LOCAL = MAX(NROWS_MINE, 1)

    ! Si CHECK_ERROR=1, Xexact est concatene a X : un seul CGEMM par bloc
    ! calcule A*X (residu) ET A*Xexact (=B), evitant de reconstruire le
    ! meme bloc de A deux fois. Recalcule redondamment par chaque rang
    ! (formule pure O(N*NRHS), moins cher que de le diffuser).
    NCOLS_MULT = NRHS
    IF (CHECK_ERROR_FLAG) THEN
      NCOLS_MULT = 2*NRHS
      ALLOCATE(XCOMBINED(N,NCOLS_MULT))
      XCOMBINED(:,1:NRHS) = XG
      DO R = 1, NRHS
        DO GI = 1, N
          XCOMBINED(GI,NRHS+R) = XEX_FORMULA(GI,R)
        END DO
      END DO
    END IF

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP
    RSS_MAX_RESID_MB = REAL(GET_RSS_KB(),8)/1024.0D0
    NDONE = 0
    PBAR_STATE = -1
    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, NB
      REND_CHUNK = MIN(RSTART_CHUNK+NB-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(BLKROWS(NCHUNK,N), RESBLK(NCHUNK,NCOLS_MULT))
      CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, N, BLKROWS)
      IF (CHECK_ERROR_FLAG) THEN
        CALL CGEMM('N','N', NCHUNK, NCOLS_MULT, N, CONE, BLKROWS, NCHUNK, XCOMBINED, N, CZERO, RESBLK, NCHUNK)
      ELSE
        CALL CGEMM('N','N', NCHUNK, NCOLS_MULT, N, CONE, BLKROWS, NCHUNK, XG, N, CZERO, RESBLK, NCHUNK)
      END IF
      DO R = 1, NRHS
        DO GI = 1, NCHUNK
          IF (CHECK_ERROR_FLAG) THEN
            BVAL = RESBLK(GI,NRHS+R)
          ELSE
            BVAL = B_ENTRY(RSTART_CHUNK+GI-1, R)
          END IF
          NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(RESBLK(GI,R)-BVAL)**2,WP)
          DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BVAL)**2,WP)
        END DO
      END DO
      DEALLOCATE(BLKROWS, RESBLK)
      NDONE = NDONE + NCHUNK
      RSS_MAX_RESID_MB = MAX(RSS_MAX_RESID_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
      ! Approximation par le rang 0 : sa progression est representative
      ! sans communication supplementaire.
      IF (GRID%MYID.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul du residu ||A*X-B|| (verification)', NDONE, NTOTAL_LOCAL, PBAR_STATE)
    END DO
    IF (ALLOCATED(XCOMBINED)) DEALLOCATE(XCOMBINED)

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, RSS_MAX_RESID_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

  ! Meme optimisation BLAS que COMPUTE_RESIDUAL pour CHECK_ERROR=1
  ! (B=A*Xexact, seul cas ou une somme O(N) couteuse est necessaire --
  ! sinon B_ENTRY est O(1), boucle scalaire directe conservee). Pas de
  ! plafonnage memoire ici : ABLOCK(BSZ,N) utilise BSZ = taille de bloc
  ! BLR locale (deja <= NB, deja borne par construction).
  SUBROUTINE FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: BL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), BBLOCK(:,:), XEXACT_VEC(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: lli, II, GII, RR, GRR, RS, NDONE, BSZ, PBAR_STATE
    NDONE = 0
    PBAR_STATE = -1
    IF (CHECK_ERROR_FLAG .AND. NRL.GT.0) THEN
      ALLOCATE(XEXACT_VEC(N,NRL))
      DO RR = 1, NRL
        GRR = RSTART + RR - 1
        DO GII = 1, N
          XEXACT_VEC(GII,RR) = XEX_FORMULA(GII,GRR)
        END DO
      END DO
    END IF
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      BSZ = AB%BEND(II) - AB%BSTART(II) + 1
      IF (CHECK_ERROR_FLAG .AND. NRL.GT.0) THEN
        ALLOCATE(ABLOCK(BSZ,N), BBLOCK(BSZ,NRL))
        CALL GET_BLOCK(AB%BSTART(II), AB%BEND(II), 1, N, ABLOCK)
        CALL CGEMM('N','N', BSZ, NRL, N, CONE, ABLOCK, BSZ, XEXACT_VEC, N, CZERO, BBLOCK, BSZ)
        BL(RS:RS+BSZ-1, 1:NRL) = BBLOCK
        DEALLOCATE(ABLOCK, BBLOCK)
        NDONE = NDONE + BSZ
        RSS_MAX_B_MB = MAX(RSS_MAX_B_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
        IF (GD%MYID.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul de B=A*Xexact (verification)', NDONE, ML, PBAR_STATE)
      ELSE
        DO GII = AB%BSTART(II), AB%BEND(II)
          DO RR = 1, NRL
            GRR = RSTART + RR - 1
            BL(RS + (GII-AB%BSTART(II)), RR) = B_ENTRY(GII, GRR)
          END DO
        END DO
      END IF
      RS = RS + BSZ
    END DO
    IF (ALLOCATED(XEXACT_VEC)) DEALLOCATE(XEXACT_VEC)
  END SUBROUTINE FILL_B

END PROGRAM EXAMPLE_BLR_MPI_PIVOT
