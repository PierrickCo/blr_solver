!===============================================================================
! EXEMPLE CANONIQUE -- reference BLR+MPI, AVEC PIVOT PARTIEL CONFINE AU
! BLOC DIAGONAL (mpi_blr_matrix_mod.f90, MPI_BLR_MATRIX_LU_FACTOR_PIVOT) :
! meme philosophie que MUMPS/STRUMPACK (pivotage local au front, jamais
! global) -- permuter les lignes d'un bloc compresse ne coute que
! O(rang) (seul le facteur UFAC est touche), jamais de decompression.
!
! Contrairement a la variante sans pivot (example_blr_mpi_sans_pivot.f90),
! celle-ci tolere des matrices SANS diagonale dominante -- la matrice de
! demonstration ci-dessous est volontairement "anti-diagonale forte /
! diagonale naturelle faible" pour forcer un vrai pivotage.
!
! But de ce fichier : montrer l'enchainement minimal (assembler,
! factoriser AVEC pivot, resoudre AVEC pivot, rassembler) sans toute
! l'instrumentation de mesure/validation croisee des drivers de test (voir
! test/test_mpi_blr_matrix_pivot.f90 pour la version qui verifie aussi la
! correction contre une reference independante).
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot PR PC [N] [NB] [NRHS] [TOL] [USE_MOVE_ALLOC] [EPS] [RATIO] [USE_RSVD] [CHECK_ERROR] [CHECK_RESIDUAL] [PROFILE_USAGE] [DUMP_STRUCTURE] [SAVE_MATRIX] [SAVE_RHS] [SAVE_FACTORED] [SAVE_SOLUTION] [LOAD_MATRIX] [LOAD_FACTORED] [LOAD_RHS] [CHECKPOINT_DIR]
!         (defauts : N=4000 NB=200 NRHS=4 TOL=1e-4 USE_MOVE_ALLOC=0 EPS=0.0 RATIO=1/3 USE_RSVD=0 CHECK_ERROR=0 CHECK_RESIDUAL=0 PROFILE_USAGE=0 DUMP_STRUCTURE=0 SAVE_MATRIX=0 SAVE_RHS=0 SAVE_FACTORED=0 SAVE_SOLUTION=0 LOAD_MATRIX=0 LOAD_FACTORED=0 LOAD_RHS=0 CHECKPOINT_DIR=.)
!         USE_MOVE_ALLOC=1 : SCHUR_UPDATE_BLR_BLOCK transfere la propriete
!         des tableaux via MOVE_ALLOC au lieu de "ABLK=NEWBLK" (qui
!         declenche une reallocation+copie implicite supplementaire) --
!         resultat numerique IDENTIQUE, bonne pratique defensive
!         independante du gain memoire (mesure a N=20000 : effet sur la
!         RSS negligeable, mais reste correct et sans risque -- garde par
!         defaut desactive pour ne rien changer au comportement existant).
!
! EPS (BRUIT PLEIN-RANG, POUR EXPERIMENTER FACILEMENT AVEC LA COMPRESSIBI-
! LITE) : ajoute EPS*HASH_NOISE(GI,GJ) (hash deterministe dans [-1,1),
! meme recette que test/bench_blr_mpi_pivot_compressibility.f90) a chaque
! entree hors anti-diagonale du noyau de Cauchy amorti ci-dessous.
!   EPS=0 (defaut)    : noyau de Cauchy amorti pur -- TRES compressible.
!   EPS croissant     : degrade PROGRESSIVEMENT la compressibilite -- de
!                       plus en plus de blocs tombent en repli dense a TOL
!                       fixe (voir "Taux de compression" affiche a deux
!                       endroits ci-dessous : apres l'assemblage et apres
!                       la factorisation). L'anti-diagonale (pivot) reste
!                       INCHANGEE quel que soit EPS -- seule la
!                       compressibilite du reste du noyau varie, sans
!                       confondre les deux effets.
!   Essayer : EPS=0, 0.01, 0.1, 1.0 pour voir la transition compressible
!   -> incompressible (le taux de compression REEL, base sur le rang
!   moyen et non le simple compte de blocs, est le plus parlant).
!
! RATIO (COMPRESS_WORTHWHILE_RATIO, blr_compress_mod.f90) : PLAFOND DE
! RANG au-dela duquel un bloc est juge "pas assez compressible pour valoir
! le coup" et repli en DENSE, meme si sa SVD trouve un rang techniquement
! valide sous TOL. Un bloc M x N reste bas-rang SEULEMENT si son rang
! R < RATIO*MIN(M,N) -- ex. RATIO=1/3 (defaut) : un bloc de taille NB=200
! bascule en dense des que son rang depasse ~66. RATIO plus GRAND accepte
! des rangs plus eleves comme "bas-rang" (plus de blocs compresses, mais
! chacun eventuellement moins avantageux, UFAC/VFAC plus larges) ; RATIO
! plus PETIT est plus strict (rejette plus tot en dense, sacrifie des
! opportunites de compression memoire-rentables entre RATIO et 0.5 contre
! une marge de securite sur le temps de calcul -- voir la mesure detaillee
! dans le commentaire de COMPRESS_WORTHWHILE_RATIO, blr_compress_mod.f90).
! Effet surtout visible dans le regime MARGINAL (rang proche du plafond) --
! voir test/bench_blr_mpi_pivot_ratio.f90 pour une etude dediee, combinee
! a EPS proche du point de transition compressible/incompressible.
!
! USE_RSVD (methode de compression a l'assemblage, blr_compress_mod.f90) :
!   0 (defaut) : SVD complete (COMPRESS_BLOCK, CGESVD) -- exacte, rang
!                optimal garanti par la theorie (Eckart-Young).
!   1          : SVD RANDOMISEE adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE,
!                rang decouvert par doublement, TOL seul) -- MESURE (voir
!                test/bench_rsvd_compress.f90) x34-54 plus rapide sur un
!                bloc bien compressible, erreur du meme ordre de grandeur.
!                GARANTIE HEURISTIQUE (pas une preuve formelle comme la
!                SVD complete) et une ZONE MORTE CARACTERISEE pres du
!                plafond RATIO ci-dessus (replie en dense la ou la SVD
!                complete aurait compresse -- comportement sur, jamais de
!                corruption, juste une occasion de gain manquee).
!   Le taux de compression REEL affiche ci-dessous (assemblage et fin de
!   factorisation) permet de comparer directement les deux methodes a
!   TOL/RATIO/EPS fixes.
!
! CHECK_ERROR (defaut 0 = comportement INCHANGE, B = formule arbitraire
! comme avant) : si 1, remplace B par un second membre EXACT, B=A*Xexact
! avec Xexact connu par formule deterministe -- permet d'afficher a la
! fin ||X_calcule-Xexact||/||Xexact||, une vraie mesure de correction (pas
! seulement "ca tourne sans planter"). CALCULE PAR PRODUIT BRUT (somme sur
! les N colonnes via LOWRANK_A, JAMAIS via une deuxieme matrice assemblee
! en RAM) -- cout CALCUL supplementaire O(N^2*NRHS) (voir la mesure faite
! sur bench_blr_mpi_pivot_ooc_geometric.f90 : a N=20000 ceci peut prendre
! plusieurs dizaines de secondes), mais AUCUN cout memoire supplementaire
! au-dela de X_GLOBAL(N,NRHS) deja alloue sur le rang 0 pour l'affichage
! des 5 premieres valeurs -- l'avantage memoire de MPI (et de la variante
! OOC, voir example_blr_mpi_pivot_ooc.f90) reste donc entier. Barre de
! progression + RSS max dedies pendant ce calcul (potentiellement long a
! grand N), voir plus bas.
!
! CHECK_RESIDUAL (defaut 0, INDEPENDANT de CHECK_ERROR -- voir ci-dessus) :
! calcule ||A*X-B||/||B|| APRES resolution, avec le B REELLEMENT utilise
! (arbitraire si CHECK_ERROR=0, ou B=A*Xexact si CHECK_ERROR=1) -- c'est
! le cas d'usage GENERAL : verifier une resolution SANS connaitre la
! solution exacte (A et B fournis par l'utilisateur, contrairement a
! CHECK_ERROR qui exige un Xexact connu par construction). Comme
! CHECK_ERROR : produit BRUT (A*X via LOWRANK_A, B reconstruit via la
! MEME formule que FILL_B -- fonction B_ENTRY partagee, jamais stocke a
! part), ZERO cout memoire au-dela de X_GLOBAL deja gathere sur le rang 0.
! Cout CALCUL O(N^2*NRHS) -- cumulatif avec CHECK_ERROR si les deux sont
! actifs (chacun refait sa propre passe O(N^2*NRHS), voir plus bas).
!
! PROFILE_USAGE (defaut 0) : si 1, chaque processus enregistre a chaque
! etape du solveur (debut, assemblage, factorisation, second membre,
! resolution, verifications eventuelles, fin) le temps ecoule, sa RSS
! (Mo) et son taux d'utilisation CPU depuis l'echantillon precedent, en
! % d'un coeur (CPU_TIME/temps ecoule -- peut depasser 100% si plusieurs
! threads MKL/OpenMP sont actifs simultanement, comme le %CPU de `top`).
! DANS l'assemblage, la factorisation ET la resolution, un echantillon
! supplementaire est pris a chaque palier de 5% de progression (via
! PROGRESS_CALLBACK, voir mpi_blr_matrix_mod.f90) -- pas seulement au
! debut/fin de chacune -- pour voir l'evolution FINE a l'interieur de
! ces 3 etapes, pas seulement leurs bornes. Ecrit un fichier
! usage_profile_rank<N>.csv par processus (jamais lu par ce programme,
! cout d'echantillonnage negligeable devant le reste). Tracer
! l'evolution : python3 examples/plot_usage_profile.py [dossier]
! (necessite matplotlib ; lit tous les usage_profile_rank*.csv trouves).
!
! DUMP_STRUCTURE (defaut 0) : si 1, chaque processus ecrit la nature de
! CHACUN de ses blocs LOCAUX (diagonal/dense/bas-rang, et le rang associe)
! a la fin de l'assemblage ET a la fin de la factorisation -- deux
! instantanes distincts, puisque le fill-in de la mise a jour de Schur
! peut faire passer certains blocs de bas-rang a dense entre les deux
! (voir SHOW_RSS_TRACE dans MPI_BLR_MATRIX_LU_FACTOR_PIVOT). Ecrit
! blr_structure_assembly_rank<N>.csv et blr_structure_factorization_
! rank<N>.csv par processus (jamais relu par ce programme). Tracer :
! python3 examples/plot_blr_structure.py [dossier] (necessite matplotlib) --
! une figure a deux panneaux (assemblage / factorisation), un pixel par
! bloc GLOBAL, colore par rang (bas-rang) ou par categorie fixe
! (diagonal, dense).
!
! SAUVEGARDE/REPRISE (SAVE_MATRIX, SAVE_RHS, SAVE_FACTORED, SAVE_SOLUTION,
! LOAD_MATRIX, LOAD_FACTORED, LOAD_RHS, CHECKPOINT_DIR) : chaque SAVE_*
! ecrit un fichier PAR PROCESSUS dans CHECKPOINT_DIR (defaut '.', doit
! deja exister), meme convention 1-fichier-par-processus/STREAM non
! formate que les caches OOC (voir example_blr_mpi_pivot_ooc.f90) --
! aucune communication MPI necessaire a l'ecriture NI a la relecture.
!   SAVE_MATRIX   : A_BLR juste apres l'assemblage (avant factorisation)
!                   -> saved_matrix_rank<N>.dat
!   SAVE_FACTORED : A_BLR (facteurs L/U + IPIV) juste apres la
!                   factorisation -> saved_factored_rank<N>.dat
!   SAVE_RHS      : B_LOCAL juste apres FILL_B (avant resolution)
!                   -> saved_rhs_rank<N>.dat
!   SAVE_SOLUTION : B_LOCAL juste apres la resolution (= X, meme memoire
!                   que B -- voir plus haut) -> saved_solution_rank<N>.dat
!   LOAD_MATRIX   : ignore l'assemblage, RELIT saved_matrix_rank<N>.dat --
!                   la factorisation s'execute normalement ensuite.
!   LOAD_FACTORED : ignore assemblage ET factorisation, RELIT
!                   saved_factored_rank<N>.dat (prioritaire sur
!                   LOAD_MATRIX si les deux sont actifs).
!   LOAD_RHS      : ignore FILL_B, RELIT saved_rhs_rank<N>.dat.
! CHAQUE fichier commence par un en-tete (N, NB, PR, PC, MYROW, MYCOL,
! rang) VERIFIE au chargement contre le run COURANT -- une grille/N/NB
! different produirait un chargement incoherent (chaque rang ne
! possederait plus les memes blocs GLOBAUX), donc REFUSE avec un message
! clair plutot que de charger silencieusement des donnees incoherentes.
! Une INCOMPATIBILITE sur UN SEUL rang arrete TOUS les rangs (reduction
! MPI_ALLREDUCE interne) -- jamais un sous-ensemble de rangs qui continue
! pendant que les autres se sont deja arretes (source de blocage MPI).
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
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  ! Apercu de X(1:5,1) sans assembler X_GLOBAL en entier -- voir son
  ! utilisation plus bas (CHECK_ERROR=0 et CHECK_RESIDUAL=0).
  INTEGER :: NPREVIEW, PREVIEW_RS, PREVIEW_MI, PREVIEW_GI
  COMPLEX(WP) :: PREVIEW_LOCAL(5), PREVIEW_GLOBAL(5)
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
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 0; CHECK_ERROR_INT = 0; CHECK_RESIDUAL_INT = 0; PROFILE_USAGE_INT = 0
  DUMP_STRUCTURE_INT = 0
  SAVE_MATRIX_INT = 0; SAVE_RHS_INT = 0; SAVE_FACTORED_INT = 0; SAVE_SOLUTION_INT = 0
  LOAD_MATRIX_INT = 0; LOAD_FACTORED_INT = 0; LOAD_RHS_INT = 0; CHECKPOINT_DIR = '.'
  PROF_N = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./example_blr_mpi_pivot PR PC [N] [NB] [NRHS] [TOL] '// &
         '[USE_MOVE_ALLOC] [EPS] [RATIO] [USE_RSVD] [CHECK_ERROR] [CHECK_RESIDUAL] [PROFILE_USAGE] [DUMP_STRUCTURE] '// &
         '[SAVE_MATRIX] [SAVE_RHS] [SAVE_FACTORED] [SAVE_SOLUTION] [LOAD_MATRIX] [LOAD_FACTORED] [LOAD_RHS] [CHECKPOINT_DIR]'
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
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) USE_MOVE_ALLOC_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN
    CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) EPS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) THEN
    CALL GET_COMMAND_ARGUMENT(9, ARG); READ(ARG,*) RATIO
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN
    CALL GET_COMMAND_ARGUMENT(10, ARG); READ(ARG,*) USE_RSVD_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN
    CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) CHECK_ERROR_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 12) THEN
    CALL GET_COMMAND_ARGUMENT(12, ARG); READ(ARG,*) CHECK_RESIDUAL_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 13) THEN
    CALL GET_COMMAND_ARGUMENT(13, ARG); READ(ARG,*) PROFILE_USAGE_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 14) THEN
    CALL GET_COMMAND_ARGUMENT(14, ARG); READ(ARG,*) DUMP_STRUCTURE_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 15) THEN
    CALL GET_COMMAND_ARGUMENT(15, ARG); READ(ARG,*) SAVE_MATRIX_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 16) THEN
    CALL GET_COMMAND_ARGUMENT(16, ARG); READ(ARG,*) SAVE_RHS_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 17) THEN
    CALL GET_COMMAND_ARGUMENT(17, ARG); READ(ARG,*) SAVE_FACTORED_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 18) THEN
    CALL GET_COMMAND_ARGUMENT(18, ARG); READ(ARG,*) SAVE_SOLUTION_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 19) THEN
    CALL GET_COMMAND_ARGUMENT(19, ARG); READ(ARG,*) LOAD_MATRIX_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 20) THEN
    CALL GET_COMMAND_ARGUMENT(20, ARG); READ(ARG,*) LOAD_FACTORED_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 21) THEN
    CALL GET_COMMAND_ARGUMENT(21, ARG); READ(ARG,*) LOAD_RHS_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 22) THEN
    CALL GET_COMMAND_ARGUMENT(22, CHECKPOINT_DIR)
  END IF
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

  ! Surcharge COMPRESS_WORTHWHILE_RATIO (blr_compress_mod.f90) AVANT tout
  ! assemblage/factorisation -- voir SET_COMPRESS_WORTHWHILE_RATIO et le
  ! commentaire RATIO de tete. Un run qui garde le defaut (1/3) reproduit
  ! EXACTEMENT le comportement d'avant l'ajout de ce parametre.
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
  END IF

  ! RSS DE REFERENCE : capturee AVANT toute allocation de donnees (A_BLR,
  ! B_LOCAL...) -- isole le cout FIXE du runtime deja charge en memoire a
  ! ce stade (Intel MPI, MKL, le Fortran runtime lui-meme) de la
  ! croissance imputable a NOTRE code. Comparer (RSS_apres - RSS_reference)
  ! a la valeur theorique est un test plus honnete que RSS_apres seule,
  ! qui melange les deux.
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
      ! --- 1) assemblage BLR distribue (matrice "anti-diagonale forte / diagonale
      !     naturelle faible" -- force un vrai pivotage, chaotique/peu compressible :
      !     la plupart des blocs hors-diagonale tombent en repli dense) ---
      ! SHOW_PROGRESS=.TRUE. : barre de progression affichee par le processus 0
      ! (voir mpi_blr_matrix_mod.f90 -- approximative ici, l'assemblage n'est
      ! pas synchronise entre processus, contrairement aux 2 etapes suivantes).
      CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                    N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., USE_RSVD=USE_RSVD_FLAG, &
                                    PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY)
      ! N_COMPRESSED/N_OFFDIAG sont deja des sommes GLOBALES (MPI_ALLREDUCE
      ! interne a MPI_ASSEMBLE_BLR_MATRIX) -- identiques sur tous les
      ! processus, un seul affichage (rang 0) suffit, pas une approximation.
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        ! Taux de compression REEL (stockage bas-rang effectif / stockage
        ! dense equivalent, base sur le RANG MOYEN -- pas juste le compte de
        ! blocs, qui peut rester a 100% meme quand chaque bloc n'economise
        ! presque rien) : 2*NB*rang_moyen / NB^2. Plus parlant que le %
        ! ci-dessus pour juger l'AMPLEUR reelle du gain memoire.
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage) : ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
      IF (SAVE_MATRIX_FLAG) CALL SAVE_BLR_MATRIX(TRIM(CHECKPOINT_DIR)//'/saved_matrix', .FALSE.)
    END IF
    ! Empreinte THEORIQUE (blocs BLR_BLOCK_T reellement stockes), PAS
    ! mesuree via l'OS.
    CALL REPORT_MEMORY_MB(GRID, 'Matrice A_BLR (compressee)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
    ! RSS REELLE (mesuree, /proc/self/status) -- forcement >= la valeur
    ! theorique ci-dessus : elle inclut aussi tout le reste du processus
    ! (workspaces SVD/QR transitoires par thread, runtime MPI/MKL charge en
    ! memoire, l'executable lui-meme, l'allocateur qui ne rend jamais toute
    ! la memoire liberee a l'OS...). Le RATIO RSS/theorique n'est PAS
    ! constant : a petite echelle (N~2000), l'overhead FIXE (runtime
    ! MPI/MKL, ~50-70 Mo/processus rien qu'au demarrage) domine largement
    ! et donne un ratio enorme (>60x, mesure) ; a plus grande echelle
    ! (N~20000), les donnees elles-memes dominent et le ratio se stabilise
    ! autour de 2.5-3x (mesure) -- c'est ce second regime qui sert de
    ! reference pour detecter une VRAIE anomalie (fuite, copie oubliee).
    CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
    CALL PROFILE_SAMPLE('assemblage')
    CALL DUMP_BLR_STRUCTURE('assembly')

    ! --- 2) factorisation LU, AVEC pivot confine au bloc diagonal ----------
    ! SHOW_PROGRESS=.TRUE. : ici K/NBLK est un pourcentage EXACT (boucle
    ! globalement synchronisee par les MPI_BCAST a chaque etape K).
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
    IF (SAVE_FACTORED_FLAG) CALL SAVE_BLR_MATRIX(TRIM(CHECKPOINT_DIR)//'/saved_factored', .TRUE.)
  END IF
  ! A_BLR est factorisee EN PLACE -- cette valeur peut differer de celle
  ! de l'assemblage (fill-in de la mise a jour de Schur).
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU (compresses)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('factorisation')
  CALL DUMP_BLR_STRUCTURE('factorization')

  ! --- 3) second membre B, distribue par blocs-lignes ET par colonnes sur
  !     NRHS -- formule deterministe, sans communication ---
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  IF (LOAD_RHS_FLAG) THEN
    CALL LOAD_RHS_VECTOR(TRIM(CHECKPOINT_DIR)//'/saved_rhs', OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  ELSE
    RSS_MAX_B_MB = 0.0D0
    CALL FILL_B(A_BLR, GRID, NRHS_LOC, RSTART_MINE, B_LOCAL, MLOC)
    IF (CHECK_ERROR_FLAG .AND. MYID0.EQ.0) &
         WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant B=A*Xexact = ', RSS_MAX_B_MB, ' Mo'
    IF (SAVE_RHS_FLAG) CALL SAVE_RHS_VECTOR(TRIM(CHECKPOINT_DIR)//'/saved_rhs')
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B', 8_8*INT(MLOC,8)*INT(NRHS_LOC,8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres allocation B)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('second_membre')

  ! --- 4) resolution 2D, AVEC pivot (applique la meme sequence de swaps
  !     que celle utilisee a la factorisation, deja diffusee a la grille) --
  ! SHOW_PROGRESS=.TRUE. : une seule barre pour les 2 sous-etapes
  ! (substitution avant + arriere), pourcentage exact (meme raison que la
  ! factorisation -- boucle globalement synchronisee).
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE., &
                                         PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE)
  T1 = MPI_WTIME()
  T_SOLVE = T1 - T0
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Resolution : ', T_SOLVE, ' s'
  CALL PROFILE_SAMPLE('resolution')
  ! B_LOCAL contient maintenant X (ecrase EN PLACE par la resolution) --
  ! meme routine de sauvegarde que le second membre, forme identique.
  IF (SAVE_SOLUTION_FLAG) CALL SAVE_RHS_VECTOR(TRIM(CHECKPOINT_DIR)//'/saved_solution')

  IF (CHECK_ERROR_FLAG .OR. CHECK_RESIDUAL_FLAG) THEN
    ! X_GLOBAL(N,NRHS) n'est reellement necessaire qu'ici : COMPUTE_ERR
    ! et COMPUTE_RESIDUAL comparent X a une reference sur toute sa
    ! largeur. Sans verification demandee, l'assembler serait un cout
    ! memoire O(N*NRHS) pur gaspillage -- voir la branche ELSE plus bas.
    ALLOCATE(X_GLOBAL(N,NRHS))
    CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

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

    ! APPEL COLLECTIF (MPI_BCAST/MPI_ALLREDUCE a l'interieur, voir son
    ! commentaire) -- TOUS les rangs doivent l'appeler, pas seulement le
    ! rang 0 (contrairement a COMPUTE_ERR ci-dessus, qui n'a aucun appel MPI).
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
    ! Ni CHECK_ERROR ni CHECK_RESIDUAL : X n'est jamais compare a une
    ! reference, inutile de l'assembler en entier. On recupere quand meme
    ! les 5 premieres valeurs de la colonne 1 pour l'affichage, sans
    ! GATHER_2D_RHS -- chaque rang contribue la valeur reelle s'il
    ! possede a la fois cette ligne (son A_BLR%BSTART/BEND local) et la
    ! colonne RHS 1 (son propre RSTART_MINE/NRHS_LOC), zero sinon ; une
    ! somme reduite sur 5 complexes suffit (cout independant de N).
    NPREVIEW = MIN(5, N)
    PREVIEW_LOCAL = (0.0_WP, 0.0_WP)
    IF (RSTART_MINE .LE. 1 .AND. 1 .LE. RSTART_MINE+NRHS_LOC-1) THEN
      PREVIEW_RS = 1
      DO li = 1, A_BLR%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        PREVIEW_MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
        DO PREVIEW_GI = 1, NPREVIEW
          IF (PREVIEW_GI.GE.A_BLR%BSTART(I) .AND. PREVIEW_GI.LE.A_BLR%BEND(I)) THEN
            PREVIEW_LOCAL(PREVIEW_GI) = B_LOCAL(PREVIEW_RS + (PREVIEW_GI-A_BLR%BSTART(I)), 1)
          END IF
        END DO
        PREVIEW_RS = PREVIEW_RS + PREVIEW_MI
      END DO
    END IF
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

  DEALLOCATE(B_LOCAL)
  IF (ALLOCATED(X_GLOBAL)) DEALLOCATE(X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Recalcule N_COMPRESSED/N_OFFDIAG/RANK_SUM depuis l'etat COURANT de
  ! A_BLR%OFFBLK (meme calcul que celui fait EN INTERNE par MPI_ASSEMBLE_
  ! BLR_MATRIX/MPI_BLR_MATRIX_LU_FACTOR_PIVOT) -- necessaire apres
  ! LOAD_MATRIX/LOAD_FACTORED puisque ces deux stats ne sont normalement
  ! produites QUE par un appel reel a l'assemblage/la factorisation, que
  ! l'on vient justement de sauter.
  SUBROUTINE RECOMPUTE_COMPRESSION_STATS(NC_OUT, NO_OUT, RS_OUT)
    INTEGER, INTENT(OUT) :: NC_OUT, NO_OUT
    INTEGER(8), INTENT(OUT) :: RS_OUT
    INTEGER :: li2, lj2, I2, J2, NC_LOCAL, NO_LOCAL, IERR2
    INTEGER(8) :: RS_LOCAL
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8
    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = (lj2-1)*GRID%PC + GRID%MYCOL + 1
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

  ! Ecrit A_BLR (SAVE_MATRIX ou SAVE_FACTORED, voir le commentaire de
  ! tete) -- un fichier PAR PROCESSUS, STREAM non formate, jamais de
  ! communication MPI. En-tete d'abord (permet a LOAD_BLR_MATRIX de
  ! verifier la coherence avec le run qui relit), puis les blocs
  ! diagonaux LOCALEMENT possedes (DIAGBLK(I)%D alloue seulement chez son
  ! proprietaire), puis TOUS les blocs hors-diagonale locaux (li,lj) --
  ! toujours "possedes" par construction du decoupage 2D block-cyclique,
  ! pas besoin de tester une allocation. WITH_IPIV n'ecrit A_BLR%IPIV que
  ! si elle est reellement allouee (apres factorisation seulement).
  SUBROUTINE SAVE_BLR_MATRIX(PREFIX, WITH_IPIV)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    INTEGER :: U, li2, lj2, I2, J2, ND, MI
    CHARACTER(LEN=16) :: RSTR
    WRITE(RSTR,'(I6.6)') MYID0
    OPEN(NEWUNIT=U, FILE=TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
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
        J2 = (lj2-1)*GRID%PC + GRID%MYCOL + 1
        ! Position diagonale du maillage OFFBLK : jamais compressee par
        ! MPI_ASSEMBLE_BLR_MATRIX (voir I.EQ.J dans mpi_blr_matrix_mod.f90)
        ! -- DENSE/UFAC/VFAC y restent non alloues, rien a ecrire.
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
  ! NBLKR_LOC/NBLKC_LOC/DIAGBLK/OFFBLK exactement comme le ferait
  ! MPI_ASSEMBLE_BLR_MATRIX (memes BLOCK_PARTITION/BLOCK_NUMROC), donc la
  ! factorisation qui suit ne voit aucune difference. OK=.FALSE. (et TOUS
  ! les rangs s'arretent -- reduction MPI_ALLREDUCE ci-dessous, voir le
  ! commentaire de tete) si le fichier est absent OU si son en-tete ne
  ! correspond pas au run COURANT (N/NB/PR/PC/grille).
  SUBROUTINE LOAD_BLR_MATRIX(PREFIX, WITH_IPIV, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: U, li2, lj2, I2, J2, ND, K2, MI, IOS, IERR2
    INTEGER :: F_N, F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
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
      READ(U) F_N, F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
      IF (F_N.NE.N .OR. F_NB.NE.NB .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NB/PR/PC ou grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    ! Un SEUL rang incompatible doit arreter TOUS les rangs -- sinon les
    ! rangs corrects continueraient seuls vers un appel MPI collectif
    ! (factorisation/resolution) que le(s) rang(s) arrete(s) n'atteindront
    ! jamais : blocage garanti plutot qu'un arret net.
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U)
      RETURN
    END IF

    A_BLR%N = N; A_BLR%NB = NB
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%NBLKC_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PC, GRID%MYCOL)
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
        J2 = (lj2-1)*GRID%PC + GRID%MYCOL + 1
        ! Meme position diagonale sautee a l'ecriture (SAVE_BLR_MATRIX) --
        ! rien a relire, A_BLR%OFFBLK(li2,lj2) reste dans son etat par
        ! defaut (IS_LOWRANK=.FALSE., M=N=RANK=0, DENSE/UFAC/VFAC non
        ! alloues), identique a une sortie fraiche de l'assemblage reel.
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

  ! Ecrit B_LOCAL (SAVE_RHS avant resolution, ou SAVE_SOLUTION apres --
  ! meme forme dans les deux cas, B_LOCAL devenant X en place). En-tete
  ! (N, NRHS, grille, MLOC/NRHS_LOC/RSTART_MINE) verifie par LOAD_RHS_
  ! VECTOR pour la meme raison que LOAD_BLR_MATRIX.
  SUBROUTINE SAVE_RHS_VECTOR(PREFIX)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER :: U
    CHARACTER(LEN=16) :: RSTR
    WRITE(RSTR,'(I6.6)') MYID0
    OPEN(NEWUNIT=U, FILE=TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
    WRITE(U) N, NRHS, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, MLOC, NRHS_LOC, RSTART_MINE
    IF (MLOC.GT.0 .AND. NRHS_LOC.GT.0) WRITE(U) B_LOCAL(1:MLOC,1:NRHS_LOC)
    CLOSE(U)
  END SUBROUTINE SAVE_RHS_VECTOR

  ! Relit B_LOCAL ecrit par SAVE_RHS_VECTOR (utilise pour LOAD_RHS) --
  ! MLOC/NRHS_LOC/RSTART_MINE doivent deja avoir ete calcules (via
  ! A_BLR/GRID, avant cet appel, comme pour FILL_B) puisqu'ils
  ! determinent la forme de B_LOCAL deja alloue par l'appelant. Meme
  ! logique OK/MPI_ALLREDUCE que LOAD_BLR_MATRIX.
  SUBROUTINE LOAD_RHS_VECTOR(PREFIX, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: U, IOS, IERR2
    INTEGER :: F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, F_NRHS_LOC, F_RSTART
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
      READ(U) F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, F_NRHS_LOC, F_RSTART
      IF (F_N.NE.N .OR. F_NRHS.NE.NRHS .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0 .OR. &
          F_MLOC.NE.MLOC .OR. F_NRHS_LOC.NE.NRHS_LOC .OR. F_RSTART.NE.RSTART_MINE) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NRHS/PR/PC/grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U)
      RETURN
    END IF

    IF (MLOC.GT.0 .AND. NRHS_LOC.GT.0) READ(U) B_LOCAL(1:MLOC,1:NRHS_LOC)
    CLOSE(U)
  END SUBROUTINE LOAD_RHS_VECTOR

  ! Taux de compression REEL (stockage bas-rang effectif / stockage dense
  ! equivalent) -- pas juste 2*NB*RANK_SUM/(NB^2*N_OFFDIAG) : RANK_SUM
  ! (MPI_ASSEMBLE_BLR_MATRIX) melange le rang reel des blocs compresses et
  ! un "rang virtuel" = NB pour chaque bloc en repli dense (convention
  ! interne de comptage, pas un vrai stockage UV) -- un bloc dense compte
  ! donc NB dans RANK_SUM mais occupe NB^2 elements (pas 2*NB^2 comme le
  ! sous-entendrait une formule naive, qui donnerait a tort 200% quand
  ! tout est dense, ex. EPS=1.0). La contribution des N_DENSE blocs denses
  ! (NB chacun) est donc d'abord retiree de RANK_SUM pour isoler le rang
  ! reel des seuls blocs compresses, puis les blocs denses sont recomptes
  ! a leur cout reel (NB^2, pas 2*NB^2).
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  ! Echantillon RSS/CPU pour l'etape STAGE_NAME (PROFILE_USAGE, voir le
  ! commentaire de tete) -- no-op si desactive, pour ne rien couter par
  ! defaut. CPU_TIME (intrinseque Fortran standard) mesure le temps
  ! PROCESSEUR consomme par CE processus (tous threads MKL/OpenMP
  ! confondus), contrairement a MPI_WTIME qui mesure le temps HORLOGE
  ! ecoule -- leur rapport sur l'intervalle depuis l'echantillon precedent
  ! donne directement le nombre moyen de coeurs actifs (ex: 3.5 coeurs en
  ! moyenne s'affiche comme 350%, meme convention que %CPU de `top`).
  ! Chaque processus echantillonne SON PROPRE etat -- c'est la variation
  ! ENTRE processus (deja visible ci-dessus via REPORT_MEMORY_MB, mais
  ! seulement au travers de son max/total reduits) que ce suivi rend
  ! exploitable point par point dans le graphe genere par
  ! plot_usage_profile.py.
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

  ! Etiquette "STAGE@NN%" pour un echantillon intra-etape (voir les 3
  ! callbacks ci-dessous). Le separateur "@" (absent des etiquettes de
  ! bornes d'etape comme 'assemblage'/'factorisation') permet a
  ! plot_usage_profile.py de ne tracer un trait de repere que sur les
  ! bornes, sans encombrer le graphe d'un trait tous les 5%.
  FUNCTION PCT_LABEL(STAGE, NDONE, NTOTAL) RESULT(LBL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    CHARACTER(LEN=40) :: LBL
    WRITE(LBL,'(A,A,I0,A)') TRIM(STAGE), '@', (100*NDONE)/MAX(NTOTAL,1), '%'
  END FUNCTION PCT_LABEL

  ! Les 3 callbacks ci-dessous ont l'interface STAGE_PROGRESS_CALLBACK_IFACE
  ! (mpi_blr_matrix_mod.f90) : appeles PAR LA BIBLIOTHEQUE au franchissement
  ! de chaque palier de 5% DANS l'assemblage/la factorisation/la resolution
  ! (voir PROGRESS_CALLBACK dans les 3 appels ci-dessus) -- c'est ce qui
  ! rend possible un echantillon RSS/CPU a une fraction donnee d'une etape,
  ! impossible depuis ce fichier seul puisque les boucles concernees sont
  ! internes a mpi_blr_matrix_mod.f90. Procedures INTERNES (CONTAINS)
  ! passees comme argument reel -- autorise depuis Fortran 2008, garde
  ! l'acces par host-association a PROFILE_SAMPLE/PROFILE_USAGE_FLAG sans
  ! variable de module supplementaire.
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

  ! Ecrit les echantillons de CE processus dans son propre CSV -- un
  ! fichier par rang (pas de rassemblement MPI) : simple, jamais de
  ! collision, chaque rang reste responsable du sien.
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

  ! Ecrit la nature de chaque bloc LOCAL de A_BLR (DUMP_STRUCTURE, voir le
  ! commentaire de tete) -- un fichier par processus, jamais rassemble ici
  ! (meme principe que PROFILE_WRITE_CSV) : plot_blr_structure.py
  ! reconstitue la grille GLOBALE NBLK x NBLK en lisant tous les fichiers
  ! d'un meme STAGE_LABEL. Un bloc diagonal (I,I) n'est possede que par UN
  ! SEUL processus (DIAGBLK(I)%D alloue seulement chez son proprietaire) ;
  ! les blocs hors-diagonale locaux sont parcourus via (li,lj) -> (I,J)
  ! comme partout ailleurs dans ce fichier (voir MPI_ASSEMBLE_BLR_MATRIX).
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
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
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
  ! independamment pour n'importe quel bloc par n'importe quel processus
  ! (meme recette que test/bench_blr_mpi_pivot_compressibility.f90).
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  ! A(GI,GJ) = 50 sur une anti-diagonale CONFINEE AU BLOC DIAGONAL de GI
  ! (GJ_STRONG, INCHANGEE quel que soit EPS), sinon un noyau de Cauchy
  ! amorti 1/(1+|GJ_STRONG-GJ|) + un bruit plein-rang EPS*HASH_NOISE (voir
  ! EPS dans le commentaire de tete : seul parametre qui degrade la
  ! compressibilite du reste du noyau).
  !
  ! L'anti-diagonale est deliberement confinee au bloc local de GI (plutot
  ! qu'une anti-diagonale globale GJ=N+1-GI) : une anti-diagonale globale
  ! est injuste envers le pivotage confine au bloc diagonal, puisque
  ! l'element fort peut tomber hors du bloc que le pivot local peut voir,
  ! ce qui rend le systeme catastrophiquement mal conditionne. Confiner
  ! l'anti-diagonale au bloc local (meme recette que test_mpi_blr_matrix_
  ! pivot.f90/bench_blr_mpi_pivot_compressibility.f90) force un pivotage
  ! tout aussi reel tout en restant bien conditionne.
  ! Extrait en fonction independante (au lieu d'etre inline dans GET_BLOCK)
  ! pour etre reutilisable par FILL_B en mode CHECK_ERROR (produit brut
  ! B=A*Xexact, voir commentaire de tete).
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

  ! Solution EXACTE deterministe, connue par formule (meme convention que
  ! test/bench_blr_mpi_pivot_ooc_compare.f90) -- utilisee UNIQUEMENT si
  ! CHECK_ERROR=1, pour construire B=A*Xexact (FILL_B) et pour mesurer
  ! l'erreur finale (COMPUTE_ERR).
  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  ! Erreur relative ||X-Xexact||/||Xexact|| -- RANG 0 UNIQUEMENT (X_GLOBAL
  ! y est deja rassemble pour l'affichage des 5 premieres valeurs, aucune
  ! allocation supplementaire de portee globale : l'avantage memoire MPI/
  ! OOC reste entier, voir le commentaire CHECK_ERROR de tete).
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

  ! Valeur de B a l'entree (GI,R) pour le cas CHECK_ERROR=0 -- SOURCE
  ! UNIQUE partagee par FILL_B (construction du second membre reellement
  ! resolu) ET COMPUTE_RESIDUAL (reconstruction de B pour le residu, voir
  ! plus bas) : garantit que le residu mesure EXACTEMENT le B utilise. Le
  ! cas CHECK_ERROR=1 (B=A*Xexact) est traite separement par ses deux
  ! appelants via un produit CGEMM direct (LOWRANK_A/XEX_FORMULA), jamais
  ! via cette fonction. Pour un probleme REEL fourni par l'utilisateur,
  ! cette fonction serait remplacee par SA propre formule/lecture de B.
  COMPLEX(WP) FUNCTION B_ENTRY(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    B_ENTRY = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP), KIND=WP)
  END FUNCTION B_ENTRY

  ! Residu ||A*X-B||/||B|| -- CAS GENERAL : ne suppose AUCUNE connaissance
  ! de la solution exacte (contrairement a COMPUTE_ERR), utilisable meme
  ! si A et B viennent d'un vrai probleme utilisateur.
  !
  ! DEUX optimisations combinees par rapport a une triple boucle scalaire
  ! naive (mesuree ~12x plus lente a N=8000) --
  !   1) BLAS (CGEMM) au lieu d'une triple boucle scalaire : A est
  !      materialisee PAR BLOCS DE LIGNES (NB a la fois, via GET_BLOCK deja
  !      existant -- PAS la matrice entiere, juste un tampon transitoire
  !      NB x N, libere apres chaque bloc) puis multipliee par X (et par
  !      Xexact si CHECK_ERROR=1, EN UN SEUL appel CGEMM en les concatenant
  !      dans XCOMBINED -- evite de reconstruire le meme bloc de A deux
  !      fois) via l'appel BLAS optimise (vectorise, MKL).
  !   2) DISTRIBUTION MPI : les N lignes sont partagees a parts egales
  !      entre TOUS les rangs de la grille (pas seulement le rang 0 comme
  !      avant) -- chaque rang calcule sa propre tranche, les normes
  !      partielles sont ensuite reduites (MPI_ALLREDUCE/SUM) pour obtenir
  !      le residu global. X doit d'abord etre DIFFUSE a tous les rangs
  !      (MPI_BCAST) : le calcul de A*X(GI,:) necessite TOUT X, pas
  !      seulement les lignes locales de GI.
  ! DEVIENT UN APPEL COLLECTIF (tous les rangs doivent l'appeler, plus
  ! seulement le rang 0) -- voir le site d'appel, deplace hors du bloc
  ! "IF (MYID0.EQ.0)". Cout memoire : un tampon BLAS transitoire NB x N
  ! par rang (quelques dizaines de Mo, libere a chaque bloc) au lieu d'un
  ! tableau (N,NRHS) complet -- PAS une deuxieme matrice N x N complete,
  ! l'avantage memoire reste entier.
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG)
    COMPLEX(WP), INTENT(INOUT) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: BLKROWS(:,:), RESBLK(:,:), XCOMBINED(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, NROWS_MINE, NTOTAL_LOCAL
    INTEGER :: RSTART_CHUNK, REND_CHUNK, NCHUNK, GI, R, NCOLS_MULT, IERR2, NDONE, PBAR_STATE
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2
    COMPLEX(WP) :: BVAL

    ! X est deja RASSEMBLE sur le rang 0 uniquement (GATHER_2D_RHS) -- les
    ! autres rangs n'ont rien de valide dedans avant cette diffusion.
    CALL MPI_BCAST(XG, N*NRHS, MPI_COMPLEX, 0, GRID%COMM, IERR2)

    ! Repartition SIMPLE (pas bloc-cyclique -- inutile ici, calcul ponctuel
    ! de verification, pas la structure de stockage du solveur) : chaque
    ! rang recoit un intervalle contigu de lignes.
    ROWS_PER_PROC = (N + GRID%NPROCS - 1) / GRID%NPROCS
    RSTART_MINE = GRID%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GRID%MYID+1)*ROWS_PER_PROC, N)
    NROWS_MINE = MAX(REND_MINE - RSTART_MINE + 1, 0)
    NTOTAL_LOCAL = MAX(NROWS_MINE, 1)

    ! Si CHECK_ERROR=1, Xexact est concatene a X : UN SEUL CGEMM par bloc
    ! calcule a la fois A*X (residu) ET A*Xexact (=B, evite de reconstruire
    ! le meme bloc de A une seconde fois). Recalcule INDEPENDAMMENT et
    ! REDONDANMENT par CHAQUE rang (formule pure, O(N*NRHS), sans
    ! communication -- moins cher que de le diffuser).
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
      ! Approximation par le rang 0 (comme l'assemblage) : chaque rang a
      ! une part comparable, sa progression est representative sans
      ! communication supplementaire.
      IF (GRID%MYID.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul du residu ||A*X-B|| (verification)', NDONE, NTOTAL_LOCAL, PBAR_STATE)
    END DO
    IF (ALLOCATED(XCOMBINED)) DEALLOCATE(XCOMBINED)

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, RSS_MAX_RESID_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

  ! Meme optimisation BLAS que COMPUTE_RESIDUAL pour le cas CHECK_ERROR=1
  ! (B=A*Xexact, seul cas ou B_ENTRY
  ! fait une somme O(N) couteuse -- sinon B_ENTRY est O(1), boucle scalaire
  ! directe conservee). Pas de plafonnage par un budget memoire ICI (le
  ! solveur classique n'a pas de notion de MAX_MB_PER_PROCESS) : le tampon
  ! ABLOCK(BSZ,N) utilise directement BSZ = taille de bloc BLR LOCALE
  ! (AB%BSTART/BEND, deja <= NB, une constante fixe independante de N) --
  ! deja borne par construction, rien a plafonner de plus.
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
