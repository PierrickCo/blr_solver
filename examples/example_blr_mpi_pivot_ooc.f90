!===============================================================================
! VARIANTE OUT-OF-CORE de example_blr_mpi_pivot.f90 -- MEME algorithme
! (pivot partiel confine au bloc diagonal, MEME matrice de demonstration),
! MEME resultat numerique, pour permettre une comparaison A/B DIRECTE avec
! l'executable existant (example_blr_mpi_pivot, JAMAIS modifie) -- voir
! mpi_blr_matrix_mod.f90, type MPI_BLR_MATRIX_OOC_T, et blr_ooc_mod.f90.
!
! PRINCIPE : ne plus etre limite par la RAM disponible, ni pour la
! matrice, ni pour les facteurs LU apres remplissage de Schur, ni pour un
! second membre massif -- 4 etapes separees, chacune bornee en memoire
! par un parametre utilisateur :
!   1) assemblage de A       -> ecrit sur un fichier scratch PAR PROCESSUS
!      (jamais partage), borne par MAX_MB_PER_PROCESS (cache LRU).
!   2) assemblage du second membre B -> ECRIT SUR DISQUE LUI AUSSI, par
!      paquets d'au plus NRHS_BATCH colonnes -- jamais tout B en RAM a la
!      fois, meme etape separee de la resolution.
!   3) factorisation LU -> les blocs "finis" (panneau/diagonale de l'etape
!      K) repartent sur disque des qu'ils ne servent plus a la mise a jour
!      de Schur ; la sous-matrice traineante reste bornee par
!      MAX_MB_PER_PROCESS via le meme cache LRU qu'a l'assemblage.
!   4) resolution -> relit L/U depuis le disque, par paquets d'au plus
!      NRHS_BATCH seconds membres (CONSUME_X_BATCH_CB rassemble chaque
!      paquet resolu dans X_GLOBAL ICI, a des fins de verification -- un
!      vrai cas massivement second membre ecrirait plutot chaque paquet
!      directement sur disque au lieu de tout accumuler).
!
! Base sur l'algorithme CLASSIQUE (jamais l'arena -- le cycle allocate/
! deallocate individuel de BLR_BLOCK_T se prete naturellement a l'eviction
! disque, contrairement a l'arena qui ne libere jamais un bloc
! individuellement).
!
! CHOISIR MAX_MB_PER_PROCESS (voir le commentaire detaille sur MPI_BLR_
! MATRIX_OOC_T dans mpi_blr_matrix_mod.f90, mesures a l'appui) : NE PAS le
! prendre "le plus petit possible" pour economiser de la RAM -- un budget
! sous l'ensemble de travail reel du processus ne degrade pas gracieusement,
! il fait THRASHER le cache (mesure : 68% du temps de factorisation perdu
! en acces cache/disque a un budget trop juste, /13.5 en portant le budget
! au-dessus de l'ensemble de travail reel, pour une RAM quasi inchangee).
! Choisir large par defaut, et ne descendre que si la RAM disponible
! l'impose vraiment.
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot_ooc PR PC [N] [NB] [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [EPS] [RATIO] [USE_RSVD] [CHECK_ERROR] [CHECK_RESIDUAL] [PROFILE_USAGE] [DUMP_STRUCTURE] [SAVE_MATRIX] [SAVE_RHS] [SAVE_FACTORED] [SAVE_SOLUTION] [LOAD_MATRIX] [LOAD_FACTORED] [LOAD_RHS] [CHECKPOINT_DIR]
!         (defauts : N=4000 NB=200 NRHS=4 TOL=1e-4 MAX_MB_PER_PROCESS=64 NRHS_BATCH=2 SCRATCH_DIR=./blr_ooc_scratch EPS=0.0 RATIO=1/3 USE_RSVD=1 CHECK_ERROR=0 CHECK_RESIDUAL=0 PROFILE_USAGE=0 DUMP_STRUCTURE=0 SAVE_MATRIX=0 SAVE_RHS=0 SAVE_FACTORED=0 SAVE_SOLUTION=0 LOAD_MATRIX=0 LOAD_FACTORED=0 LOAD_RHS=0 CHECKPOINT_DIR=.)
!
! EPS (BRUIT PLEIN-RANG, POUR EXPERIMENTER FACILEMENT AVEC LA COMPRESSIBI-
! LITE -- IDENTIQUE a example_blr_mpi_pivot.f90, voir son commentaire de
! tete pour le detail) : EPS=0 (defaut) = noyau tres compressible ; EPS
! croissant degrade progressivement la compressibilite, de plus en plus de
! blocs tombant en repli dense a TOL fixe. Passer la MEME valeur d'EPS aux
! deux executables pour comparer classique/OOC sur EXACTEMENT la meme
! matrice.
!
! RATIO (COMPRESS_WORTHWHILE_RATIO -- IDENTIQUE a example_blr_mpi_pivot.f90,
! voir son commentaire de tete pour le detail) : plafond de rang au-dela
! duquel un bloc est repli en DENSE meme si sa SVD trouve un rang valide
! sous TOL. Passer la MEME valeur aux deux executables pour comparer.
!
! USE_RSVD (methode de compression a l'assemblage -- IDENTIQUE a example_
! blr_mpi_pivot.f90, voir son commentaire de tete pour le detail) : 0
! (defaut) = SVD complete (CGESVD) ; 1 = SVD randomisee adaptative.
! Passer la MEME valeur aux deux executables pour comparer.
!
! CHECK_ERROR (verification numerique -- IDENTIQUE a example_blr_mpi_
! pivot.f90, voir son commentaire de tete pour le detail complet) : 0
! (defaut) = comportement INCHANGE (B arbitraire). 1 = B=A*Xexact par
! produit BRUT (JAMAIS via une deuxieme matrice assemblee ou chargee en
! RAM -- l'avantage memoire de l'OOC reste ENTIER, voir plus bas), erreur
! affichee en fin de resolution. Barre de progression + RSS max dedies
! pendant le calcul de B (etape potentiellement longue a grand N, cout
! CALCUL O(N^2*NRHS) mais ZERO cout memoire supplementaire).
!
! CHECK_RESIDUAL (verification SANS solution exacte connue -- INDEPENDANT
! de CHECK_ERROR) : calcule ||A*X-B||/||B|| APRES resolution avec le B
! REELLEMENT utilise -- le cas d'usage GENERAL (A et B fournis par
! l'utilisateur, solution inconnue).
!
! Pour eviter de rappeler GET_BLOCK/FILL_RHS une deuxieme fois (couteux si
! ces callbacks le sont vraiment, cas d'usage reel), A (dense, telle que
! renvoyee par GET_BLOCK) et B (telle qu'ecrite par FILL_RHS, donc deja
! coherente avec CHECK_ERROR) sont mises en cache sur disque UNE SEULE
! FOIS, DES QUE POSSIBLE (juste apres assemblage / apres assemblage du
! RHS, donc AVANT que la factorisation n'ecrase A par L/U en place) --
! acheck_rank_NNNNNN.dat et bcheck_rank_NNNNNN.dat, meme convention
! 1-fichier-par-processus, meme decoupage par blocs de lignes que la
! verification elle-meme (borne par MAX_MB_PER_PROCESS) -- voir
! example_user_matrix_ooc.f90 pour le detail complet de la technique
! (identique ici). Des qu'un bloc de B a servi au calcul du residu, il
! n'est plus utile : X (deja disponible dans X_GLOBAL) est ECRIT A SA
! PLACE, DANS LE MEME FICHIER (reecriture en place par POS=, comme
! OOC_WRITE_BLOCK_FIXED dans blr_ooc_mod.f90) -- a la fin,
! acheck_rank_*.dat est supprime, bcheck_rank_*.dat survit et contient
! desormais X (copie disque de la solution, en plus de X_GLOBAL en RAM).
!
! PROFILE_USAGE / DUMP_STRUCTURE : voir example_blr_mpi_pivot.f90 pour le
! detail complet (memes parametres, memes scripts de trace). Ici,
! l'echantillonnage tous les 5% DANS la resolution redemarre a chaque
! PAQUET de second membre (NRHS_BATCH) -- chaque paquet est sa propre
! "etape" du point de vue du profilage. Sans CHECK_ERROR ni
! CHECK_RESIDUAL, X_GLOBAL(N,NRHS) n'est PAS accumule (voir CONSUME_X_
! BATCH) -- seul un apercu de X(1:5,1) est recupere, sans le cout memoire
! ni communication du rassemblement complet a chaque paquet.
!
! SAVE_MATRIX/SAVE_RHS/SAVE_FACTORED/SAVE_SOLUTION et LOAD_MATRIX/
! LOAD_FACTORED/LOAD_RHS (+ CHECKPOINT_DIR, defaut '.', doit deja exister)
! -- meme fonctionnalite que example_blr_mpi_pivot.f90 (voir son
! commentaire de tete pour la description generale), adaptee au stockage
! OOC : UN FICHIER PAR PROCESSUS dans CHECKPOINT_DIR (STREAM non formate),
! DISTINCT des fichiers scratch internes a_rank_*.dat/b_rank_*.dat
! (SCRATCH_DIR) -- ceux-ci restent du scratch PUR, supprimes en fin
! d'execution (OOC_STORE_FREE), alors que les fichiers de checkpoint
! (prefixe saved_..._ooc_) sont PERSISTANTS, jamais supprimes
! automatiquement :
!   - saved_matrix_ooc_rank<N>.dat / saved_factored_ooc_rank<N>.dat :
!     en-tete (N/NB/grille) puis, pour chaque bloc LOCAL de A_BLR (relu
!     depuis le fichier scratch OOC via OOC_READ_BLOCK_RAW -- les caches
!     CACHE_A/CACHE_D sont d'abord vidangees via OOC_CACHE_FLUSH_ALL/
!     OOC_DIAG_CACHE_FLUSH_ALL pour garantir un OFFSET valide sur TOUS les
!     blocs, y compris ceux jamais encore evinces), les memes champs que
!     la version RAM (IS_LOWRANK/M/N/RANK + UFAC/VFAC ou DENSE) ; IPIV en
!     plus pour la version factorisee. LOAD_* reconstruit A_BLR (BSTART/
!     BEND/NBLKR_LOC/NBLKC_LOC via BLOCK_PARTITION/BLOCK_NUMROC, comme une
!     assemblage neuve), REINITIALISE le stockage OOC (OOC_STORE_INIT/
!     OOC_CACHE_INIT/OOC_DIAG_CACHE_INIT, comme MPI_ASSEMBLE_BLR_MATRIX_
!     OOC), puis reecrit chaque bloc lu sur le NOUVEAU fichier scratch via
!     OOC_WRITE_BLOCK_FIXED -- exactement comme si le bloc venait d'etre
!     compresse pour la premiere fois (RAM residente reste bornee : un
!     seul bloc a la fois, jamais toute la matrice).
!   - saved_rhs_ooc_rank<N>.dat : en-tete (N/NRHS/grille/MLOC/NRHS_BATCH/
!     NBATCHES) puis chaque paquet RHS (relu via OOC_READ_BLOCK_RAW, aucun
!     cache a vidanger ici -- RHS%STORE_B est toujours ecrit directement,
!     voir blr_ooc_mod.f90). LOAD_RHS_OOC reconstruit OOC_RHS_T a
!     l'identique (meme STORE_B/BATCH neufs, remplis via OOC_WRITE_BLOCK_
!     RAW dans le meme ordre).
!   - saved_solution_ooc_rank<N>.dat : PAS de contrepartie RAM a relire
!     (comme dans la version RAM, aucun LOAD_SOLUTION) -- X n'existe nulle
!     part en un seul bloc pour l'OOC (voir CONSUME_X_BATCH), donc chaque
!     paquet resolu est ecrit AU FIL DE L'EAU, directement depuis CONSUME_
!     X_BATCH, des qu'il est calcule (jamais accumule au-dela d'un seul
!     paquet a la fois -- l'avantage memoire de l'OOC reste entier).
! Comme la version RAM : en-tete verifie (N/NB ou NRHS/grille) avant tout
! chargement, arret COLLECTIF propre (MPI_ALLREDUCE) si un seul rang
! detecte une incoherence, et LOAD_FACTORED prioritaire sur LOAD_MATRIX si
! les deux sont positionnes.
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_PIVOT_OOC
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS, EPS, RATIO
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_BLR
  TYPE(OOC_RHS_T) :: RHS
  COMPLEX(WP), ALLOCATABLE :: X_GLOBAL(:,:)
  INTEGER :: INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR
  INTEGER :: USE_RSVD_INT, CHECK_ERROR_INT, CHECK_RESIDUAL_INT, PROFILE_USAGE_INT, DUMP_STRUCTURE_INT
  LOGICAL :: USE_RSVD_FLAG, CHECK_ERROR_FLAG, CHECK_RESIDUAL_FLAG, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG
  DOUBLE PRECISION :: RSS_MAX_RHS_MB, RSS_MAX_RESID_MB
  REAL(WP) :: ERR, RESID

  ! Sauvegarde/reprise (voir le commentaire de tete) -- un fichier par
  ! processus dans CHECKPOINT_DIR par objet sauvegarde/relu, DISTINCT des
  ! fichiers scratch internes (SCRATCH_DIR).
  INTEGER :: SAVE_MATRIX_INT, SAVE_RHS_INT, SAVE_FACTORED_INT, SAVE_SOLUTION_INT
  INTEGER :: LOAD_MATRIX_INT, LOAD_FACTORED_INT, LOAD_RHS_INT
  LOGICAL :: SAVE_MATRIX_FLAG, SAVE_RHS_FLAG, SAVE_FACTORED_FLAG, SAVE_SOLUTION_FLAG
  LOGICAL :: LOAD_MATRIX_FLAG, LOAD_FACTORED_FLAG, LOAD_RHS_FLAG, OK
  CHARACTER(LEN=256) :: CHECKPOINT_DIR
  ! Unite ouverte pendant TOUTE la resolution quand SAVE_SOLUTION_FLAG :
  ! chaque paquet resolu (CONSUME_X_BATCH) y ecrit sa part au fil de
  ! l'eau, jamais d'accumulation au-dela d'un seul paquet a la fois.
  INTEGER :: SAVE_SOLUTION_UNIT

  ! Apercu de X(1:5,1) sans accumuler X_GLOBAL en entier (CHECK_ERROR=0 et
  ! CHECK_RESIDUAL=0) -- accumule DANS CONSUME_X_BATCH (voir plus bas),
  ! reduit UNE FOIS apres la resolution complete (tous paquets confondus).
  LOGICAL :: WANT_FULL_X
  INTEGER :: NPREVIEW
  COMPLEX(WP) :: PREVIEW_LOCAL(5), PREVIEW_GLOBAL(5)

  ! Suivi memoire/CPU par etape (PROFILE_USAGE) -- voir example_blr_mpi_
  ! pivot.f90 pour le detail complet de la technique.
  INTEGER, PARAMETER :: MAX_PROFILE_SAMPLES = 128
  CHARACTER(LEN=40) :: PROF_STAGE(MAX_PROFILE_SAMPLES)
  DOUBLE PRECISION :: PROF_T(MAX_PROFILE_SAMPLES), PROF_RSS_MB(MAX_PROFILE_SAMPLES), PROF_CPU_PCT(MAX_PROFILE_SAMPLES)
  INTEGER :: PROF_N
  DOUBLE PRECISION :: PROF_T_START, PROF_T_PREV
  REAL :: PROF_CPU_PREV

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP
  MAX_MB_PER_PROCESS = 64.0_WP; NRHS_BATCH = 2; SCRATCH_DIR = './blr_ooc_scratch'; EPS = 0.0_WP
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 1; CHECK_ERROR_INT = 0; CHECK_RESIDUAL_INT = 0
  PROFILE_USAGE_INT = 0; DUMP_STRUCTURE_INT = 0
  SAVE_MATRIX_INT = 0; SAVE_RHS_INT = 0; SAVE_FACTORED_INT = 0; SAVE_SOLUTION_INT = 0
  LOAD_MATRIX_INT = 0; LOAD_FACTORED_INT = 0; LOAD_RHS_INT = 0; CHECKPOINT_DIR = '.'
  PROF_N = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./example_blr_mpi_pivot_ooc PR PC [N] [NB] [NRHS] [TOL] '// &
         '[MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [EPS] [RATIO] [USE_RSVD] '// &
         '[CHECK_ERROR] [CHECK_RESIDUAL] [PROFILE_USAGE] [DUMP_STRUCTURE] [SAVE_MATRIX] [SAVE_RHS] '// &
         '[SAVE_FACTORED] [SAVE_SOLUTION] [LOAD_MATRIX] [LOAD_FACTORED] [LOAD_RHS] [CHECKPOINT_DIR]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) MAX_MB_PER_PROCESS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN
    CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) NRHS_BATCH
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) THEN
    CALL GET_COMMAND_ARGUMENT(9, SCRATCH_DIR)
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN
    CALL GET_COMMAND_ARGUMENT(10, ARG); READ(ARG,*) EPS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN
    CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) RATIO
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 12) THEN
    CALL GET_COMMAND_ARGUMENT(12, ARG); READ(ARG,*) USE_RSVD_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 13) THEN
    CALL GET_COMMAND_ARGUMENT(13, ARG); READ(ARG,*) CHECK_ERROR_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 14) THEN
    CALL GET_COMMAND_ARGUMENT(14, ARG); READ(ARG,*) CHECK_RESIDUAL_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 15) THEN
    CALL GET_COMMAND_ARGUMENT(15, ARG); READ(ARG,*) PROFILE_USAGE_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 16) THEN
    CALL GET_COMMAND_ARGUMENT(16, ARG); READ(ARG,*) DUMP_STRUCTURE_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 17) THEN
    CALL GET_COMMAND_ARGUMENT(17, ARG); READ(ARG,*) SAVE_MATRIX_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 18) THEN
    CALL GET_COMMAND_ARGUMENT(18, ARG); READ(ARG,*) SAVE_RHS_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 19) THEN
    CALL GET_COMMAND_ARGUMENT(19, ARG); READ(ARG,*) SAVE_FACTORED_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 20) THEN
    CALL GET_COMMAND_ARGUMENT(20, ARG); READ(ARG,*) SAVE_SOLUTION_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 21) THEN
    CALL GET_COMMAND_ARGUMENT(21, ARG); READ(ARG,*) LOAD_MATRIX_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 22) THEN
    CALL GET_COMMAND_ARGUMENT(22, ARG); READ(ARG,*) LOAD_FACTORED_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 23) THEN
    CALL GET_COMMAND_ARGUMENT(23, ARG); READ(ARG,*) LOAD_RHS_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 24) THEN
    CALL GET_COMMAND_ARGUMENT(24, CHECKPOINT_DIR)
  END IF
  SAVE_MATRIX_FLAG = (SAVE_MATRIX_INT .NE. 0)
  SAVE_RHS_FLAG = (SAVE_RHS_INT .NE. 0)
  SAVE_FACTORED_FLAG = (SAVE_FACTORED_INT .NE. 0)
  SAVE_SOLUTION_FLAG = (SAVE_SOLUTION_INT .NE. 0)
  LOAD_MATRIX_FLAG = (LOAD_MATRIX_INT .NE. 0)
  LOAD_FACTORED_FLAG = (LOAD_FACTORED_INT .NE. 0)
  LOAD_RHS_FLAG = (LOAD_RHS_INT .NE. 0)
  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
  CHECK_ERROR_FLAG = (CHECK_ERROR_INT .NE. 0)
  CHECK_RESIDUAL_FLAG = (CHECK_RESIDUAL_INT .NE. 0)
  PROFILE_USAGE_FLAG = (PROFILE_USAGE_INT .NE. 0)
  DUMP_STRUCTURE_FLAG = (DUMP_STRUCTURE_INT .NE. 0)
  WANT_FULL_X = (CHECK_ERROR_FLAG .OR. CHECK_RESIDUAL_FLAG)

  ! Voir example_blr_mpi_pivot.f90 pour le detail (SET_COMPRESS_WORTHWHILE_
  ! RATIO, blr_compress_mod.f90).
  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1,A,I0,A,A,A,ES9.2,A,F6.3)') &
         ' [OOC] Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
         ' , MAX_MB_PER_PROCESS=', MAX_MB_PER_PROCESS, ' , NRHS_BATCH=', NRHS_BATCH, &
         ' , SCRATCH_DIR=', TRIM(SCRATCH_DIR), ' , EPS(bruit)=', EPS, ' , RATIO=', RATIO
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
  END IF

  PROF_T_START = MPI_WTIME()
  PROF_T_PREV = 0.0D0
  CALL CPU_TIME(PROF_CPU_PREV)
  CALL PROFILE_SAMPLE('debut')

  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  ! --- 1) assemblage BLR distribue, out-of-core -> disque (ou reprise) ---
  IF (LOAD_FACTORED_FLAG) THEN
    CALL LOAD_BLR_MATRIX_OOC(TRIM(CHECKPOINT_DIR)//'/saved_factored_ooc', .TRUE., OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    T_ASM = 0.0D0
    CALL RECOMPUTE_COMPRESSION_STATS_OOC(N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (factorisation, charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
           ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
      WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (factorisation, charge) : ', &
           REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
  ELSE
    IF (LOAD_MATRIX_FLAG) THEN
      CALL LOAD_BLR_MATRIX_OOC(TRIM(CHECKPOINT_DIR)//'/saved_matrix_ooc', .FALSE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
      T_ASM = 0.0D0
      CALL RECOMPUTE_COMPRESSION_STATS_OOC(N_COMPRESSED, N_OFFDIAG, RANK_SUM)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage, charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage, charge) : ', &
             REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
    ELSE
      CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                        T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., &
                                        USE_RSVD=USE_RSVD_FLAG, PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        ! Taux de compression REEL, base sur le RANG MOYEN -- voir
        ! REAL_COMPRESSION_PCT ci-dessous (meme correctif que example_blr_mpi_
        ! pivot.f90 : un bloc en repli dense compte NB dans RANK_SUM, PAS 2*NB,
        ! la formule naive donnait a tort 200% quand tout est dense).
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage) : ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
      IF (SAVE_MATRIX_FLAG) CALL SAVE_BLR_MATRIX_OOC(TRIM(CHECKPOINT_DIR)//'/saved_matrix_ooc', .FALSE.)
    END IF
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'RAM residente (cache, apres assemblage)', MPI_BLR_MATRIX_OOC_RESIDENT_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('assemblage')
  CALL DUMP_BLR_STRUCTURE('assembly')

  ! Cache de A pour CHECK_RESIDUAL : AVANT la factorisation, qui ecrase A
  ! par L/U en place -- GET_BLOCK n'est plus jamais rappelee ensuite.
  IF (CHECK_RESIDUAL_FLAG) CALL CACHE_A_FOR_CHECK(GRID, N, NB)

  ! --- 2) assemblage du second membre, SEPAREMENT, aussi sur disque (ou
  !     reprise) -- INDEPENDANT du choix ci-dessus (assemble/charge/
  !     factorisee chargee), toujours execute --------------------------
  RSS_MAX_RHS_MB = 0.0D0
  IF (LOAD_RHS_FLAG) THEN
    CALL LOAD_RHS_OOC(TRIM(CHECKPOINT_DIR)//'/saved_rhs_ooc', OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  ELSE
    CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
    IF (SAVE_RHS_FLAG) CALL SAVE_RHS_OOC(TRIM(CHECKPOINT_DIR)//'/saved_rhs_ooc')
  END IF
  IF (CHECK_ERROR_FLAG .AND. MYID0.EQ.0) &
       WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant B=A*Xexact = ', RSS_MAX_RHS_MB, ' Mo'
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage B)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('second_membre')

  ! Cache de B pour CHECK_RESIDUAL -- deja coherente avec CHECK_ERROR
  ! puisque FILL_RHS (rappelee ici) applique la meme logique qu'a
  ! l'assemblage reel du second membre.
  IF (CHECK_RESIDUAL_FLAG) CALL CACHE_B_FOR_CHECK(GRID, N, NB, NRHS)

  ! --- 3) factorisation LU, AVEC pivot, out-of-core (sautee si LOAD_
  !     FACTORED -- deja disponible, voir ci-dessus) -----------------------
  IF (.NOT. LOAD_FACTORED_FLAG) THEN
    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                             N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                             SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE., &
                                             PROGRESS_CALLBACK=PROFILE_CALLBACK_FACTOR)
    IF (INFO .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO
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
    IF (SAVE_FACTORED_FLAG) CALL SAVE_BLR_MATRIX_OOC(TRIM(CHECKPOINT_DIR)//'/saved_factored_ooc', .TRUE.)
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'RAM residente (cache, apres factorisation)', MPI_BLR_MATRIX_OOC_RESIDENT_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('factorisation')
  CALL DUMP_BLR_STRUCTURE('factorization')

  ! --- 4) resolution, out-of-core, par paquets de second membre -----------
  ! X_GLOBAL(N,NRHS) n'est accumule (par CONSUME_X_BATCH) que si
  ! CHECK_ERROR ou CHECK_RESIDUAL le necessite -- sinon CONSUME_X_BATCH ne
  ! recupere qu'un apercu de X(1:5,1) (PREVIEW_LOCAL/PREVIEW_GLOBAL, voir
  ! le commentaire de tete et la fonction elle-meme), sans jamais assembler
  ! X en entier ni meme le paquet contenant la colonne RHS 1.
  IF (WANT_FULL_X) ALLOCATE(X_GLOBAL(N,NRHS))
  NPREVIEW = MIN(5, N)
  PREVIEW_LOCAL = (0.0_WP, 0.0_WP)
  ! X n'existe nulle part en un seul bloc pour l'OOC (voir CONSUME_X_
  ! BATCH) -- SAVE_SOLUTION_FLAG ecrit donc chaque paquet resolu au fil de
  ! l'eau, DANS le callback, des qu'il est disponible (jamais accumule).
  IF (SAVE_SOLUTION_FLAG) THEN
    BLOCK
      CHARACTER(LEN=16) :: RSTR
      WRITE(RSTR,'(I6.6)') MYID0
      OPEN(NEWUNIT=SAVE_SOLUTION_UNIT, FILE=TRIM(CHECKPOINT_DIR)//'/saved_solution_ooc_rank'//TRIM(RSTR)//'.dat', &
           ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
      WRITE(SAVE_SOLUTION_UNIT) N, NRHS, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, &
           RHS%MLOC, RHS%NRHS_BATCH, RHS%NBATCHES
    END BLOCK
  END IF
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH, SHOW_PROGRESS=.TRUE., &
                                             PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE)
  IF (SAVE_SOLUTION_FLAG) CLOSE(SAVE_SOLUTION_UNIT)
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres resolution)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('resolution')

  ! PAS de verification independante ICI (voir test/test_mpi_blr_matrix_
  ! pivot_ooc.f90 pour la validation rigoureuse, X_exact connu, quel que
  ! soit NRHS_BATCH/MAX_MB_PER_PROCESS -- meme philosophie que les autres
  ! exemples canoniques : demonstration de l'enchainement d'appels, pas un
  ! test de non-regression). NOTE : sur CETTE matrice de demonstration
  ! (chaotique, mal conditionnee), la valeur numerique de X(:,1) peut
  ! varier legerement avec NRHS_BATCH (simple precision + formes d'appels
  ! BLAS differentes amplifiees par le mauvais conditionnement) -- ce n'est
  ! pas un signe d'erreur.
  IF (.NOT. WANT_FULL_X) THEN
    CALL MPI_REDUCE(PREVIEW_LOCAL, PREVIEW_GLOBAL, NPREVIEW, MPI_COMPLEX, MPI_SUM, 0, GRID%COMM, IERR)
  END IF
  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
    BLOCK
      INTEGER :: I
      IF (WANT_FULL_X) THEN
        DO I = 1, NPREVIEW
          WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
        END DO
      ELSE
        DO I = 1, NPREVIEW
          WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(PREVIEW_GLOBAL(I)), ' + i*', AIMAG(PREVIEW_GLOBAL(I))
        END DO
        WRITE(6,'(A)') ' (X complet non accumule -- CHECK_ERROR=0 et CHECK_RESIDUAL=0, memoire economisee)'
      END IF
    END BLOCK
    IF (CHECK_ERROR_FLAG) THEN
      ERR = COMPUTE_ERR(X_GLOBAL)
      WRITE(6,'(A,ES12.4)') ' Erreur relative ||X-Xexact|| / ||Xexact|| = ', ERR
      CALL PROFILE_SAMPLE('verif_erreur')
    END IF
  END IF

  ! APPEL COLLECTIF (MPI_BCAST/MPI_ALLREDUCE a l'interieur) -- TOUS les
  ! rangs doivent l'appeler, pas seulement le rang 0. Consomme et supprime
  ! acheck_rank_*.dat au passage ; bcheck_rank_*.dat survit seul,
  ! contenant desormais X (voir COMPUTE_RESIDUAL).
  IF (CHECK_RESIDUAL_FLAG) THEN
    RSS_MAX_RESID_MB = 0.0D0
    RESID = COMPUTE_RESIDUAL(X_GLOBAL)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant le residu = ', RSS_MAX_RESID_MB, ' Mo'
      WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
      WRITE(6,'(A)') ' acheck_rank_*.dat supprime ; bcheck_rank_*.dat contient desormais X '// &
           '(economie d''espace disque).'
    END IF
    CALL PROFILE_SAMPLE('verif_residu')
  END IF

  CALL PROFILE_SAMPLE('fin')
  CALL PROFILE_WRITE_CSV()

  IF (ALLOCATED(X_GLOBAL)) DEALLOCATE(X_GLOBAL)
  CALL OOC_RHS_FREE(RHS)
  CALL MPI_BLR_MATRIX_OOC_FREE(A_BLR)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Recalcule N_COMPRESSED/N_OFFDIAG/RANK_SUM depuis l'etat COURANT de
  ! A_BLR%OFFBLK -- necessaire apres LOAD_MATRIX_OOC/LOAD_FACTORED_OOC
  ! puisque ces stats ne sont normalement produites QUE par un appel reel
  ! a l'assemblage/la factorisation, que l'on vient justement de sauter.
  ! CONTRAIREMENT a la version RAM, AUCUNE lecture disque necessaire ici :
  ! IS_LOWRANK/M/N/RANK sont des CHAMPS DU DESCRIPTEUR (OOC_BLOCK_DESC_T),
  ! toujours resident en RAM et valides quel que soit l'etat du cache (voir
  ! DUMP_BLR_STRUCTURE ci-dessous, meme principe).
  SUBROUTINE RECOMPUTE_COMPRESSION_STATS_OOC(NC_OUT, NO_OUT, RS_OUT)
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
  END SUBROUTINE RECOMPUTE_COMPRESSION_STATS_OOC

  ! Ecrit A_BLR (SAVE_MATRIX ou SAVE_FACTORED, voir le commentaire de
  ! tete) -- vidange D'ABORD les caches (OOC_CACHE_FLUSH_ALL/OOC_DIAG_
  ! CACHE_FLUSH_ALL) pour garantir un OFFSET valide sur TOUS les blocs
  ! locaux, y compris ceux jamais encore evinces vers le fichier scratch
  ! (normalement deja fait en fin d'assemblage/factorisation, voir
  ! mpi_blr_matrix_mod.f90 -- appel defensif ici, sans effet si deja
  ! vidangees). Chaque bloc est ensuite RELU depuis le disque UN A LA FOIS
  ! (OOC_READ_BLOCK_RAW) et ecrit sur le fichier de checkpoint -- jamais
  ! plus d'un bloc resident a la fois, l'avantage memoire de l'OOC reste
  ! entier meme pendant la sauvegarde.
  SUBROUTINE SAVE_BLR_MATRIX_OOC(PREFIX, WITH_IPIV)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    TYPE(BLR_BLOCK_T) :: BLK
    INTEGER :: U, li2, lj2, I2, J2, ND
    CHARACTER(LEN=16) :: RSTR

    CALL OOC_CACHE_FLUSH_ALL(A_BLR%CACHE_A, A_BLR%STORE_A, A_BLR%OFFBLK)
    CALL OOC_DIAG_CACHE_FLUSH_ALL(A_BLR%CACHE_D, A_BLR%STORE_A, A_BLR%DIAGBLK)

    WRITE(RSTR,'(I6.6)') MYID0
    OPEN(NEWUNIT=U, FILE=TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
    ! Marqueur de format (N > 0 dans l'ancien format, sans marqueur) +
    ! distribution des blocs-colonnes (CMAP, voir BLR_COL_OWNER).
    WRITE(U) -7411, A_BLR%CMAP
    WRITE(U) N, NB, A_BLR%NBLK, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, &
         A_BLR%NBLKR_LOC, A_BLR%NBLKC_LOC

    ND = 0
    DO I2 = 1, A_BLR%NBLK
      IF (A_BLR%DIAGBLK(I2)%OFFSET .GE. 0_8) ND = ND + 1
    END DO
    WRITE(U) ND
    DO I2 = 1, A_BLR%NBLK
      IF (A_BLR%DIAGBLK(I2)%OFFSET .GE. 0_8) THEN
        CALL OOC_READ_BLOCK_RAW(A_BLR%STORE_A, A_BLR%DIAGBLK(I2), BLK)
        WRITE(U) I2, BLK%M
        WRITE(U) BLK%DENSE
      END IF
    END DO

    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        IF (I2 .EQ. J2) CYCLE
        CALL OOC_READ_BLOCK_RAW(A_BLR%STORE_A, A_BLR%OFFBLK(li2,lj2), BLK)
        WRITE(U) I2, J2, BLK%IS_LOWRANK, BLK%M, BLK%N, BLK%RANK
        IF (BLK%IS_LOWRANK) THEN
          WRITE(U) BLK%UFAC
          WRITE(U) BLK%VFAC
        ELSE
          WRITE(U) BLK%DENSE
        END IF
      END DO
    END DO

    IF (WITH_IPIV) WRITE(U) A_BLR%IPIV
    CLOSE(U)
  END SUBROUTINE SAVE_BLR_MATRIX_OOC

  ! Relit A_BLR ecrite par SAVE_BLR_MATRIX_OOC -- reconstruit BSTART/BEND/
  ! NBLKR_LOC/NBLKC_LOC (BLOCK_PARTITION/BLOCK_NUMROC, comme une
  ! assemblage neuve), REINITIALISE le stockage OOC (OOC_STORE_INIT/
  ! OOC_CACHE_INIT/OOC_DIAG_CACHE_INIT -- meme sequence que MPI_ASSEMBLE_
  ! BLR_MATRIX_OOC, sur un NOUVEAU fichier scratch dans SCRATCH_DIR), puis
  ! reecrit chaque bloc lu via OOC_WRITE_BLOCK_FIXED -- exactement comme
  ! s'il venait d'etre compresse pour la premiere fois (OFFSET=-1 par
  ! defaut sur les descripteurs neufs, meme convention de reservation
  ! pire-cas-dense qu'un bloc frais). Un seul bloc resident a la fois.
  ! OK=.FALSE. (et TOUS les rangs s'arretent -- reduction MPI_ALLREDUCE,
  ! voir le commentaire de tete) si le fichier est absent OU si son
  ! en-tete ne correspond pas au run COURANT (N/NB/PR/PC/grille).
  SUBROUTINE LOAD_BLR_MATRIX_OOC(PREFIX, WITH_IPIV, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    LOGICAL, INTENT(OUT) :: OK
    TYPE(BLR_BLOCK_T) :: BLK
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
    CALL OOC_STORE_INIT(A_BLR%STORE_A, TRIM(SCRATCH_DIR), 'a', MYID0)
    CALL OOC_CACHE_INIT(A_BLR%CACHE_A, MAX_MB_PER_PROCESS, A_BLR%NBLKR_LOC*A_BLR%NBLKC_LOC+4, &
                         MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1))
    CALL OOC_DIAG_CACHE_INIT(A_BLR%CACHE_D, MAX_MB_PER_PROCESS, A_BLR%NBLKR_LOC+2, A_BLR%NBLK)
    A_BLR%MAX_MB_PER_PROCESS = MAX_MB_PER_PROCESS

    READ(U) ND
    DO K2 = 1, ND
      READ(U) I2, MI
      BLK%IS_LOWRANK = .FALSE.; BLK%M = MI; BLK%N = MI; BLK%RANK = MI
      ALLOCATE(BLK%DENSE(MI,MI))
      READ(U) BLK%DENSE
      CALL OOC_WRITE_BLOCK_FIXED(A_BLR%STORE_A, A_BLR%DIAGBLK(I2), BLK)
      DEALLOCATE(BLK%DENSE)
    END DO

    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        IF (I2 .EQ. J2) CYCLE
        READ(U) R_I, R_J, R_LOWRANK, R_M, R_N, R_RANK
        BLK%IS_LOWRANK = R_LOWRANK; BLK%M = R_M; BLK%N = R_N; BLK%RANK = R_RANK
        IF (R_LOWRANK) THEN
          ALLOCATE(BLK%UFAC(R_M,R_RANK), BLK%VFAC(R_RANK,R_N))
          READ(U) BLK%UFAC
          READ(U) BLK%VFAC
        ELSE
          ALLOCATE(BLK%DENSE(R_M,R_N))
          READ(U) BLK%DENSE
        END IF
        CALL OOC_WRITE_BLOCK_FIXED(A_BLR%STORE_A, A_BLR%OFFBLK(li2,lj2), BLK)
        IF (ALLOCATED(BLK%UFAC)) DEALLOCATE(BLK%UFAC)
        IF (ALLOCATED(BLK%VFAC)) DEALLOCATE(BLK%VFAC)
        IF (ALLOCATED(BLK%DENSE)) DEALLOCATE(BLK%DENSE)
      END DO
    END DO

    IF (WITH_IPIV) THEN
      ALLOCATE(A_BLR%IPIV(NB, A_BLR%NBLK))
      READ(U) A_BLR%IPIV
    END IF
    CLOSE(U)
  END SUBROUTINE LOAD_BLR_MATRIX_OOC

  ! Ecrit RHS (SAVE_RHS) -- pas de cache a vidanger ici (RHS%STORE_B est
  ! toujours ecrit directement des l'assemblage, voir blr_ooc_mod.f90),
  ! chaque paquet est simplement RELU (OOC_READ_BLOCK_RAW) puis ecrit sur
  ! le fichier de checkpoint, un paquet a la fois.
  SUBROUTINE SAVE_RHS_OOC(PREFIX)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    TYPE(BLR_BLOCK_T) :: BLK
    INTEGER :: U, B
    CHARACTER(LEN=16) :: RSTR
    WRITE(RSTR,'(I6.6)') MYID0
    OPEN(NEWUNIT=U, FILE=TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
    WRITE(U) N, NRHS, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, RHS%MLOC, RHS%NRHS_BATCH, RHS%NBATCHES
    DO B = 1, RHS%NBATCHES
      CALL OOC_READ_BLOCK_RAW(RHS%STORE_B, RHS%BATCH(B), BLK)
      WRITE(U) BLK%M, BLK%N
      WRITE(U) BLK%DENSE
    END DO
    CLOSE(U)
  END SUBROUTINE SAVE_RHS_OOC

  ! Relit RHS ecrit par SAVE_RHS_OOC -- reconstruit OOC_RHS_T a
  ! l'identique : REINITIALISE STORE_B (OOC_STORE_INIT, nouveau fichier
  ! scratch dans SCRATCH_DIR) puis reecrit chaque paquet lu via OOC_WRITE_
  ! BLOCK_RAW, dans le MEME ordre (B=1..NBATCHES) -- la resolution
  ! retrouve exactement la meme structure qu'un OOC_RHS_T fraichement
  ! assemble. Meme logique OK/MPI_ALLREDUCE que LOAD_BLR_MATRIX_OOC.
  SUBROUTINE LOAD_RHS_OOC(PREFIX, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(OUT) :: OK
    TYPE(BLR_BLOCK_T) :: BLK
    INTEGER :: U, B, BM, BN, IOS, IERR2
    INTEGER :: F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, F_NRHS_BATCH, F_NBATCHES
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    OK_LOCAL = .TRUE.
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- avez-vous sauvegarde avec SAVE_RHS=1 auparavant ?'
      OK_LOCAL = .FALSE.
    ELSE
      READ(U) F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, F_NRHS_BATCH, F_NBATCHES
      IF (F_N.NE.N .OR. F_NRHS.NE.NRHS .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NRHS/PR/PC ou grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U)
      RETURN
    END IF

    RHS%NRHS = NRHS
    RHS%NRHS_BATCH = F_NRHS_BATCH
    RHS%NBATCHES = F_NBATCHES
    RHS%MLOC = F_MLOC
    CALL OOC_STORE_INIT(RHS%STORE_B, TRIM(SCRATCH_DIR), 'b', MYID0)
    ALLOCATE(RHS%BATCH(RHS%NBATCHES))

    DO B = 1, RHS%NBATCHES
      READ(U) BM, BN
      BLK%IS_LOWRANK = .FALSE.; BLK%M = BM; BLK%N = BN; BLK%RANK = MIN(BM,BN)
      ALLOCATE(BLK%DENSE(BM,BN))
      READ(U) BLK%DENSE
      CALL OOC_WRITE_BLOCK_RAW(RHS%STORE_B, RHS%BATCH(B), BLK)
      DEALLOCATE(BLK%DENSE)
    END DO
    CLOSE(U)
  END SUBROUTINE LOAD_RHS_OOC

  ! Echantillon RSS/CPU pour l'etape STAGE_NAME (PROFILE_USAGE) -- voir
  ! example_blr_mpi_pivot.f90 pour le detail complet de la technique.
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

  FUNCTION PCT_LABEL(STAGE, NDONE, NTOTAL) RESULT(LBL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    CHARACTER(LEN=40) :: LBL
    WRITE(LBL,'(A,A,I0,A)') TRIM(STAGE), '@', (100*NDONE)/MAX(NTOTAL,1), '%'
  END FUNCTION PCT_LABEL

  ! Callbacks PROGRESS_CALLBACK (STAGE_PROGRESS_CALLBACK_IFACE) -- pour la
  ! resolution, NDONE/NTOTAL redemarrent a chaque PAQUET de second membre
  ! (voir MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D), donc l'etiquette
  ! "resolution@NN%" se repete une fois par paquet -- attendu, voir le
  ! commentaire de tete.
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

  ! Ecrit les echantillons de CE processus (usage_profile_rank<N>.csv) --
  ! voir example_blr_mpi_pivot.f90 pour le detail complet.
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

  ! Ecrit la nature de chaque bloc LOCAL de A_BLR (DUMP_STRUCTURE) -- voir
  ! example_blr_mpi_pivot.f90 pour le detail complet de la technique.
  ! DIFFERENCE avec la version RAM : A_BLR%DIAGBLK(I) est ici un
  ! OOC_BLOCK_DESC_T (descripteur LEGER, pas un tableau dense) -- son
  ! proprietaire se detecte via OFFSET>=0 (ecrit sur disque par
  ! OOC_WRITE_BLOCK_FIXED, voir blr_ooc_mod.f90), au lieu de
  ! ALLOCATED(...%D). A_BLR%OFFBLK(li,lj)%IS_LOWRANK/RANK restent, eux,
  ! toujours valides des l'assemblage (poses directement sur le
  ! descripteur en RAM), qu'ils aient ete evinces sur disque ou non.
  SUBROUTINE DUMP_BLR_STRUCTURE(STAGE_LABEL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_LABEL
    INTEGER :: U, li, lj, I, J
    CHARACTER(LEN=96) :: FNAME
    IF (.NOT. DUMP_STRUCTURE_FLAG) RETURN
    WRITE(FNAME,'(A,A,A,I0,A)') 'blr_structure_', TRIM(STAGE_LABEL), '_rank', MYID0, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'block_row,block_col,kind,rank'
    DO I = 1, A_BLR%NBLK
      IF (A_BLR%DIAGBLK(I)%OFFSET .GE. 0_8) THEN
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

  ! Voir example_blr_mpi_pivot.f90 pour le detail du correctif.
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  ! Hash deterministe (I,J) -> [-1,1) -- meme recette que example_blr_mpi_
  ! pivot.f90/test/bench_blr_mpi_pivot_compressibility.f90.
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  ! Matrice IDENTIQUE a example_blr_mpi_pivot.f90 (comparaison A/B directe,
  ! y compris le bruit EPS -- passer la MEME valeur d'EPS aux deux
  ! executables). Anti-diagonale CONFINEE au bloc diagonal de GI (voir le
  ! commentaire detaille dans example_blr_mpi_pivot.f90/LOWRANK_A pour la
  ! raison du correctif -- meme correction appliquee ICI a l'identique).
  ! Extrait en fonction independante (au lieu d'etre inline dans GET_BLOCK)
  ! pour etre reutilisable par FILL_RHS en mode CHECK_ERROR.
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

  ! Solution EXACTE deterministe -- IDENTIQUE a example_blr_mpi_pivot.f90,
  ! utilisee UNIQUEMENT si CHECK_ERROR=1.
  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  ! Erreur relative ||X-Xexact||/||Xexact|| -- RANG 0 UNIQUEMENT (X_GLOBAL
  ! y est deja rassemble par CONSUME_X_BATCH, aucune allocation supple-
  ! mentaire de portee globale : l'avantage memoire de l'OOC reste entier).
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

  ! Valeur de B a l'entree (GI,R) -- SOURCE UNIQUE partagee par FILL_RHS
  ! ET COMPUTE_RESIDUAL, voir example_blr_mpi_pivot.f90/B_ENTRY pour le
  ! detail (meme raisonnement, identique ici).
  COMPLEX(WP) FUNCTION B_ENTRY(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    INTEGER :: JJ
    COMPLEX(WP) :: ACC
    COMPLEX(WP), PARAMETER :: CZERO = (0.0_WP,0.0_WP)
    IF (CHECK_ERROR_FLAG) THEN
      ACC = CZERO
      DO JJ = 1, N
        ACC = ACC + LOWRANK_A(GI,JJ) * XEX_FORMULA(JJ,R)
      END DO
      B_ENTRY = ACC
    ELSE
      B_ENTRY = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP), KIND=WP)
    END IF
  END FUNCTION B_ENTRY

  ! Ecrit UNE FOIS le cache dense de A (acheck_rank_NNNNNN.dat), par blocs
  ! de lignes bornes par MAX_MB_PER_PROCESS -- meme decoupage que
  ! COMPUTE_RESIDUAL, pour que la relecture s'y cale exactement. Appelee
  ! AVANT la factorisation (qui ecrase A par L/U en place). Voir
  ! example_user_matrix_ooc.f90/CACHE_A_FOR_CHECK pour le detail complet
  ! (technique identique).
  SUBROUTINE CACHE_A_FOR_CHECK(GD, NN, NNB)
    TYPE(GRID_T), INTENT(IN) :: GD
    INTEGER, INTENT(IN)      :: NN, NNB
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK, CR
    INTEGER :: UNIT_A
    CHARACTER(LEN=16) :: RSTR

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)
    CR = MAX(1, MIN(NNB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(NN,WP)))))

    WRITE(RSTR,'(I6.6)') GD%MYID
    OPEN(NEWUNIT=UNIT_A, FILE=TRIM(SCRATCH_DIR)//'/acheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')

    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CR
      REND_CHUNK = MIN(RSTART_CHUNK+CR-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(ABLOCK(NCHUNK,NN))
      CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, NN, ABLOCK)
      WRITE(UNIT_A) RSTART_CHUNK, REND_CHUNK, NCHUNK
      WRITE(UNIT_A) ABLOCK(1:NCHUNK,1:NN)
      DEALLOCATE(ABLOCK)
    END DO
    CLOSE(UNIT_A)
  END SUBROUTINE CACHE_A_FOR_CHECK

  ! Ecrit UNE FOIS le cache de B (bcheck_rank_NNNNNN.dat), par bloc de
  ! lignes -- meme decoupage que COMPUTE_RESIDUAL. B est ici deja LARGE
  ! de NRHS colonnes en une seule passe (pas de sous-paquets : X_GLOBAL
  ! lui-meme n'est jamais borne dans ce fichier de demonstration, voir le
  ! commentaire de COMPUTE_RESIDUAL). Rappelle FILL_RHS -- deja coherente
  ! avec CHECK_ERROR puisque c'est la MEME fonction qui a servi a
  ! l'assemblage reel du second membre.
  SUBROUTINE CACHE_B_FOR_CHECK(GD, NN, NNB, NR)
    TYPE(GRID_T), INTENT(IN) :: GD
    INTEGER, INTENT(IN)      :: NN, NNB, NR
    COMPLEX(WP), ALLOCATABLE :: BBLOCK(:,:)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK, CR
    INTEGER :: UNIT_B
    CHARACTER(LEN=16) :: RSTR

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)
    CR = MAX(1, MIN(NNB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(NN,WP)))))

    WRITE(RSTR,'(I6.6)') GD%MYID
    OPEN(NEWUNIT=UNIT_B, FILE=TRIM(SCRATCH_DIR)//'/bcheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')

    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CR
      REND_CHUNK = MIN(RSTART_CHUNK+CR-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(BBLOCK(NCHUNK,NR))
      CALL FILL_RHS(RSTART_CHUNK, REND_CHUNK, 1, NR, BBLOCK)
      WRITE(UNIT_B) RSTART_CHUNK, REND_CHUNK, NCHUNK
      WRITE(UNIT_B) BBLOCK(1:NCHUNK,1:NR)
      DEALLOCATE(BBLOCK)
    END DO
    CLOSE(UNIT_B)
  END SUBROUTINE CACHE_B_FOR_CHECK

  ! Residu ||A*X-B||/||B|| -- CAS GENERAL, voir example_blr_mpi_pivot.f90/
  ! COMPUTE_RESIDUAL pour le detail COMPLET du correctif de plafonnage
  ! memoire (identique ici) : BLAS (CGEMM par blocs de lignes) +
  ! DISTRIBUTION MPI sur TOUS les rangs (X diffuse via MPI_BCAST, normes
  ! partielles reduites via MPI_ALLREDUCE) -- APPEL COLLECTIF (voir le
  ! site d'appel).
  !
  ! A et B sont RELUES depuis le cache disque (acheck_rank_*.dat,
  ! bcheck_rank_*.dat -- ecrit par CACHE_A_FOR_CHECK/CACHE_B_FOR_CHECK,
  ! AVANT la factorisation) au lieu de rappeler GET_BLOCK/FILL_RHS/
  ! B_ENTRY -- B etant deja coherente avec CHECK_ERROR (meme FILL_RHS),
  ! XCOMBINED/NCOLS_MULT ne sont plus necessaires : la comparaison se fait
  ! directement contre le B cache, que CHECK_ERROR soit actif ou non.
  ! Des qu'un bloc de B a servi, X (deja disponible dans XG) l'ecrase EN
  ! PLACE dans bcheck_rank_*.dat (POS=, meme forme) -- economie d'espace
  ! disque, voir le commentaire de tete. acheck_rank_*.dat est supprime en
  ! fin de fonction.
  !
  ! LIMITE RESTANTE, HONNETE (prexistante, DOCUMENTEE) : X_GLOBAL(N,NRHS)
  ! reste un tableau COMPLET non borne par MAX_MB_PER_PROCESS -- MEME
  ! limite deja assumee par le commentaire de tete de ce fichier pour
  ! l'affichage des 5 premieres valeurs de X ("un vrai cas massivement
  ! second membre ecrirait plutot chaque paquet directement sur disque").
  ! Cette verification (comme l'affichage) est un outil de DEMONSTRATION,
  ! pas le coeur de l'algorithme OOC (assemblage/factorisation/resolution,
  ! ceux-la restent strictement bornes) -- rester honnete sur cette
  ! distinction plutot que de le cacher.
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG)
    COMPLEX(WP), INTENT(INOUT) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), AXBLOCK(:,:), BBLOCK(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, NROWS_MINE, NTOTAL_LOCAL
    INTEGER :: RSTART_CHUNK, REND_CHUNK, NCHUNK, GI, R, IERR2, NDONE, PBAR_STATE
    INTEGER :: CHUNK_ROWS, UNIT_A, UNIT_B
    INTEGER :: HA_RSTART, HA_REND, HA_NCHUNK, HB_RSTART, HB_REND, HB_NCHUNK
    INTEGER(8) :: PAYLOAD_POS
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2
    CHARACTER(LEN=16) :: RSTR

    CALL MPI_BCAST(XG, N*NRHS, MPI_COMPLEX, 0, GRID%COMM, IERR2)

    ROWS_PER_PROC = (N + GRID%NPROCS - 1) / GRID%NPROCS
    RSTART_MINE = GRID%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GRID%MYID+1)*ROWS_PER_PROC, N)
    NROWS_MINE = MAX(REND_MINE - RSTART_MINE + 1, 0)
    NTOTAL_LOCAL = MAX(NROWS_MINE, 1)

    ! Taille de bloc de lignes bornee par MAX_MB_PER_PROCESS -- MEME
    ! decoupage que CACHE_A_FOR_CHECK/CACHE_B_FOR_CHECK (relecture calee
    ! exactement dessus).
    CHUNK_ROWS = MAX(1, MIN(NB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(N,WP)))))

    WRITE(RSTR,'(I6.6)') GRID%MYID
    OPEN(NEWUNIT=UNIT_A, FILE=TRIM(SCRATCH_DIR)//'/acheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')
    OPEN(NEWUNIT=UNIT_B, FILE=TRIM(SCRATCH_DIR)//'/bcheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READWRITE')

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP
    RSS_MAX_RESID_MB = REAL(GET_RSS_KB(),8)/1024.0D0
    NDONE = 0
    PBAR_STATE = -1
    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CHUNK_ROWS
      REND_CHUNK = MIN(RSTART_CHUNK+CHUNK_ROWS-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1

      READ(UNIT_A) HA_RSTART, HA_REND, HA_NCHUNK
      ALLOCATE(ABLOCK(NCHUNK,N))
      READ(UNIT_A) ABLOCK(1:NCHUNK,1:N)

      READ(UNIT_B) HB_RSTART, HB_REND, HB_NCHUNK
      INQUIRE(UNIT=UNIT_B, POS=PAYLOAD_POS)
      ALLOCATE(BBLOCK(NCHUNK,NRHS))
      READ(UNIT_B) BBLOCK(1:NCHUNK,1:NRHS)

      ALLOCATE(AXBLOCK(NCHUNK,NRHS))
      CALL CGEMM('N','N', NCHUNK, NRHS, N, CONE, ABLOCK, NCHUNK, XG, N, CZERO, AXBLOCK, NCHUNK)
      DO R = 1, NRHS
        DO GI = 1, NCHUNK
          NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(AXBLOCK(GI,R)-BBLOCK(GI,R))**2,WP)
          DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BBLOCK(GI,R))**2,WP)
        END DO
      END DO

      ! B vient de servir : plus utile -- X ecrase sa place sur disque.
      WRITE(UNIT_B, POS=PAYLOAD_POS) XG(RSTART_CHUNK:REND_CHUNK, 1:NRHS)

      DEALLOCATE(ABLOCK, AXBLOCK, BBLOCK)
      NDONE = NDONE + NCHUNK
      RSS_MAX_RESID_MB = MAX(RSS_MAX_RESID_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
      IF (GRID%MYID.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul du residu ||A*X-B|| (verification)', NDONE, NTOTAL_LOCAL, PBAR_STATE)
    END DO
    CLOSE(UNIT_B)
    CLOSE(UNIT_A, STATUS='DELETE')

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, RSS_MAX_RESID_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

  ! B(GI,GR) : formule arbitraire (defaut, CHECK_ERROR=0, INCHANGEE, O(1)
  ! par entree -- boucle scalaire directe conservee) ou produit BRUT
  ! B=A*Xexact (CHECK_ERROR=1, voir commentaire de tete).
  !
  ! Meme optimisation BLAS que COMPUTE_RESIDUAL pour le cas CHECK_ERROR=1,
  ! et meme plafonnage par MAX_MB_PER_PROCESS -- [ISTART,IEND] (fourni par
  ! MPI_ASSEMBLE_RHS_OOC) n'est pas garanti respecter le budget a lui seul
  ! (sa largeur depend de la structure BLR, pas du budget memoire
  ! utilisateur), donc le tampon ABLOCK(CHUNK,N) est encore sous-decoupe
  ! si besoin -- meme raisonnement que dans COMPUTE_RESIDUAL (voir son
  ! commentaire) : ne jamais laisser un tampon BLAS depasser le budget
  ! annonce par l'utilisateur.
  SUBROUTINE FILL_RHS(ISTART, IEND, RSTART, REND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RSTART, REND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, REND-RSTART+1)
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), BBLOCK(:,:), XEXACT_VEC(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: II, RR, NDONE, NCOLS, CHUNK_ROWS, RSTART_CHUNK, REND_CHUNK, NCHUNK, PBAR_STATE
    NDONE = 0
    PBAR_STATE = -1
    IF (CHECK_ERROR_FLAG) THEN
      NCOLS = REND - RSTART + 1
      ALLOCATE(XEXACT_VEC(N,NCOLS))
      DO RR = RSTART, REND
        DO II = 1, N
          XEXACT_VEC(II,RR-RSTART+1) = XEX_FORMULA(II,RR)
        END DO
      END DO
      CHUNK_ROWS = MAX(1, MIN(IEND-ISTART+1, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(N,WP)))))
      DO RSTART_CHUNK = ISTART, IEND, CHUNK_ROWS
        REND_CHUNK = MIN(RSTART_CHUNK+CHUNK_ROWS-1, IEND)
        NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
        ALLOCATE(ABLOCK(NCHUNK,N), BBLOCK(NCHUNK,NCOLS))
        CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, N, ABLOCK)
        CALL CGEMM('N','N', NCHUNK, NCOLS, N, CONE, ABLOCK, NCHUNK, XEXACT_VEC, N, CZERO, BBLOCK, NCHUNK)
        BLK(RSTART_CHUNK-ISTART+1:REND_CHUNK-ISTART+1, 1:NCOLS) = BBLOCK
        DEALLOCATE(ABLOCK, BBLOCK)
        NDONE = NDONE + NCHUNK
        RSS_MAX_RHS_MB = MAX(RSS_MAX_RHS_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
        IF (MYID0.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul de B=A*Xexact (verification)', NDONE, IEND-ISTART+1, PBAR_STATE)
      END DO
      DEALLOCATE(XEXACT_VEC)
    ELSE
      DO RR = RSTART, REND
        DO II = ISTART, IEND
          BLK(II-ISTART+1, RR-RSTART+1) = B_ENTRY(II, RR)
        END DO
      END DO
    END IF
  END SUBROUTINE FILL_RHS

  ! Rassemble CE paquet resolu dans X_GLOBAL SI CHECK_ERROR/CHECK_RESIDUAL
  ! le necessite (verification -- un vrai cas massivement second membre
  ! ecrirait plutot CE paquet directement sur disque ici, voir commentaire
  ! de tete). SINON, extrait juste l'apercu de X(1:5,1) SANS assembler quoi
  ! que ce soit -- SEUL le paquet qui possede effectivement la colonne RHS
  ! globale 1 (GSTART<=1<=GEND) contribue, les autres ne font rien ; meme
  ! technique que example_blr_mpi_pivot.f90 (X_LOCAL joue ici le role de
  ! B_LOCAL la-bas, PREVIEW_LOCAL est reduit UNE FOIS apres la resolution
  ! complete, voir le site d'appel).
  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    COMPLEX(WP), ALLOCATABLE :: X_BATCH_GLOBAL(:,:)
    INTEGER :: WIDTH, li, I, PREVIEW_RS, PREVIEW_MI, PREVIEW_GI, LOCAL_COL
    IF (SAVE_SOLUTION_FLAG) THEN
      WRITE(SAVE_SOLUTION_UNIT) RHS%MLOC, NRHS_LOC_B
      IF (RHS%MLOC.GT.0 .AND. NRHS_LOC_B.GT.0) WRITE(SAVE_SOLUTION_UNIT) X_LOCAL(1:RHS%MLOC,1:NRHS_LOC_B)
    END IF
    IF (WANT_FULL_X) THEN
      WIDTH = BATCH_REND - BATCH_RSTART + 1
      ALLOCATE(X_BATCH_GLOBAL(N,WIDTH))
      CALL GATHER_2D_RHS(GRID, N, NB, WIDTH, X_LOCAL, LDX, NRHS_LOC_B, 0, X_BATCH_GLOBAL)
      IF (MYID0 .EQ. 0) X_GLOBAL(:,BATCH_RSTART:BATCH_REND) = X_BATCH_GLOBAL
      DEALLOCATE(X_BATCH_GLOBAL)
    ELSE IF (GSTART .LE. 1 .AND. 1 .LE. GEND) THEN
      LOCAL_COL = 1 - GSTART + 1
      PREVIEW_RS = 1
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        PREVIEW_MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
        DO PREVIEW_GI = 1, NPREVIEW
          IF (PREVIEW_GI.GE.A_BLR%BSTART(I) .AND. PREVIEW_GI.LE.A_BLR%BEND(I)) THEN
            PREVIEW_LOCAL(PREVIEW_GI) = X_LOCAL(PREVIEW_RS + (PREVIEW_GI-A_BLR%BSTART(I)), LOCAL_COL)
          END IF
        END DO
        PREVIEW_RS = PREVIEW_RS + PREVIEW_MI
      END DO
    END IF
  END SUBROUTINE CONSUME_X_BATCH

END PROGRAM EXAMPLE_BLR_MPI_PIVOT_OOC
