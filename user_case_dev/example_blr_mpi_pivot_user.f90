!===============================================================================
! Solveur BLR+MPI avec pivot -- version utilisateur, ou USER_GET_BLOCK/
! USER_FILL_B (matrix/user_matrix_mod.f90) sont les deux routines a
! remplacer par vos propres matrice/second membre ; rien ICI n'a besoin
! d'etre modifie. Ce fichier ne contient que l'orchestration (parsing
! CLI, assemblage/factorisation/resolution/verification) ; les routines
! reutilisables vivent dans leurs propres sous-dossiers :
!   cli/cli_options_mod.f90          : parsing CLI
!   checkpoint/checkpoint_io_mod.f90 : sauvegarde/reprise, stats de compression
!   profiling/profiling_mod.f90      : PROFILE_USAGE, DUMP_STRUCTURE
!   residual/residual_mod.f90        : CHECK_RESIDUAL (||A*X-B||/||B||, par lot)
!   matrix/user_matrix_mod.f90       : matrice/second membre (a personnaliser)
!   blr_lib/                         : bibliotheque BLR+MPI (vendorisee)
!   scalapack/                       : PCMATGEN + generateur "a saut"
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot_user PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]
!         (--help/-h affiche la liste complete). PR PC obligatoires ;
!         N/NB/NRHS/TOL positionnels optionnels dans cet ordre (defauts
!         4000/200/4/1e-4) ; toutes les autres options sont NOMMEES,
!         --cle=valeur ou --cle seul (= --cle=1) pour les booleennes,
!         dans n'importe quel ordre : --ratio --use-rsvd --check-residual
!         --profile-usage --dump-structure --save-matrix --save-rhs
!         --save-factored --save-solution --load-matrix --load-factored
!         --load-rhs --checkpoint-dir --nrhs-batch --reload-solution
!         (N/NB/NRHS/TOL aussi dispo en --n/--nb/--nrhs/--tol).
!
! CHECK_RESIDUAL (defaut 0) : ||A*X-B||/||B|| avec le B REELLEMENT resolu
! -- ne suppose pas de solution exacte connue. Calcule PAR LOT (voir
! residual/residual_mod.f90) : gathere seulement les colonnes du lot
! courant, jamais NRHS entier -- l'avantage memoire de NRHS_BATCH reste
! entier meme avec la verification active.
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_PIVOT_USER
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  USE USER_MATRIX_MOD, ONLY: USER_MATRIX_INIT, USER_GET_BLOCK, USER_FILL_B
  USE CLI_OPTIONS_MOD, ONLY: SPLIT_KV, PRINT_USAGE
  USE CHECKPOINT_IO_MOD
  USE PROFILING_MOD
  USE RESIDUAL_MOD, ONLY: ACCUMULATE_RESIDUAL
  USE MPI
  IMPLICIT NONE

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL, RATIO
  CHARACTER(LEN=256) :: ARG
  CHARACTER(LEN=64)  :: KEYSTR
  CHARACTER(LEN=256) :: VALSTR
  INTEGER :: IARG, POS_IDX
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), B_ORIG_LOCAL(:,:), XG_BATCH(:,:)
  INTEGER :: MLOC, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  ! Apercu de X(1:5,1) sans jamais assembler X en entier -- cout independant de N.
  INTEGER :: NPREVIEW, PREVIEW_RS, PREVIEW_MI, PREVIEW_GI
  COMPLEX(WP) :: PREVIEW_LOCAL(5), PREVIEW_GLOBAL(5)

  ! Second membre/resolution/sauvegarde PAR LOT (NRHS_BATCH) -- borne la
  ! RAM a un lot au lieu de NRHS entier.
  INTEGER :: NRHS_BATCH, NRHS_BATCH_EFF, NBATCHES, BATCH_IDX
  INTEGER :: RSTART_B, REND_B, WIDTH_B, LOCAL_OFFSET_B, GSTART_B, GEND_B, NRHS_LOC_B
  INTEGER :: SAVE_RHS_UNIT, SAVE_SOLUTION_UNIT, LOAD_RHS_UNIT
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  INTEGER :: USE_RSVD_INT, CHECK_RESIDUAL_INT, PROFILE_USAGE_INT, DUMP_STRUCTURE_INT
  LOGICAL :: USE_RSVD_FLAG, CHECK_RESIDUAL_FLAG, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG
  REAL(WP) :: RESID
  REAL(WP) :: NUM2_ACC, DEN2_ACC

  ! Sauvegarde/reprise -- un fichier par processus dans CHECKPOINT_DIR par
  ! objet sauvegarde/relu.
  INTEGER :: SAVE_MATRIX_INT, SAVE_RHS_INT, SAVE_FACTORED_INT, SAVE_SOLUTION_INT
  INTEGER :: LOAD_MATRIX_INT, LOAD_FACTORED_INT, LOAD_RHS_INT
  LOGICAL :: SAVE_MATRIX_FLAG, SAVE_RHS_FLAG, SAVE_FACTORED_FLAG, SAVE_SOLUTION_FLAG
  LOGICAL :: LOAD_MATRIX_FLAG, LOAD_FACTORED_FLAG, LOAD_RHS_FLAG, OK
  CHARACTER(LEN=256) :: CHECKPOINT_DIR

  ! --reload-solution : en toute fin de run, libere A_BLR (l'essentiel de
  ! la RAM restante) puis relit la solution depuis disque, SEUL le rang 0
  ! -- demontre/verifie qu'un post-traitement leger (un seul processus,
  ! sans la grille MPI ni A_BLR en memoire) peut recuperer X apres coup.
  ! Implique --save-solution (force a .TRUE. si l'utilisateur ne l'a pas
  ! deja mis -- sans quoi il n'y a rien a relire).
  INTEGER :: RELOAD_SOLUTION_INT
  LOGICAL :: RELOAD_SOLUTION_FLAG
  COMPLEX(WP), ALLOCATABLE :: X_RELOADED(:,:)

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 0; CHECK_RESIDUAL_INT = 0; PROFILE_USAGE_INT = 0
  DUMP_STRUCTURE_INT = 0
  SAVE_MATRIX_INT = 0; SAVE_RHS_INT = 0; SAVE_FACTORED_INT = 0; SAVE_SOLUTION_INT = 0
  LOAD_MATRIX_INT = 0; LOAD_FACTORED_INT = 0; LOAD_RHS_INT = 0; CHECKPOINT_DIR = '.'
  NRHS_BATCH = 0
  RELOAD_SOLUTION_INT = 0

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
  !     n'importe quel ordre (voir --help).
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
    CASE ('ratio');           READ(VALSTR,*) RATIO
    CASE ('use-rsvd');        READ(VALSTR,*) USE_RSVD_INT
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
    CASE ('reload-solution'); READ(VALSTR,*) RELOAD_SOLUTION_INT
    CASE DEFAULT
      IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : option inconnue "--', TRIM(KEYSTR), '" (voir --help).'
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END SELECT
    IARG = IARG + 1
  END DO

  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
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
  RELOAD_SOLUTION_FLAG = (RELOAD_SOLUTION_INT .NE. 0)
  ! --reload-solution n'a rien a relire sans une sauvegarde -- force
  ! SAVE_SOLUTION_FLAG si l'utilisateur ne l'a pas deja demande (le
  ! message correspondant est affiche plus bas, avec les autres options).
  IF (RELOAD_SOLUTION_FLAG) SAVE_SOLUTION_FLAG = .TRUE.

  ! Transmet N/NB/NRHS a USER_GET_BLOCK/USER_FILL_B (interface fixe
  ! imposee par la bibliotheque).
  CALL USER_MATRIX_INIT(N, NB, NRHS)

  ! Surcharge COMPRESS_WORTHWHILE_RATIO avant tout assemblage ; le defaut
  ! (1/3) ne change rien.
  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F6.3)') &
         ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, ' , RATIO=', RATIO
    IF (USE_RSVD_FLAG) THEN
      WRITE(6,'(A)') ' Compression : SVD randomisee adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE)'
    ELSE
      WRITE(6,'(A)') ' Compression : SVD complete (COMPRESS_BLOCK, CGESVD)'
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
    IF (RELOAD_SOLUTION_FLAG) THEN
      WRITE(6,'(A)') ' Fin de run : A_BLR sera libere puis la solution relue depuis disque (rang 0 seul)'
      IF (SAVE_SOLUTION_INT.EQ.0) WRITE(6,'(A)') ' (--reload-solution force --save-solution=1 : rien a relire sinon)'
    END IF
    IF (NRHS_BATCH.GT.0 .AND. NRHS_BATCH.LT.NRHS .AND. .NOT.LOAD_RHS_FLAG) &
         WRITE(6,'(A,I0,A)') ' Second membre/resolution/sauvegarde par lots de ', NRHS_BATCH, ' colonnes (NRHS_BATCH)'
  END IF

  CALL PROFILING_INIT(MYID0, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG)

  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  IF (LOAD_FACTORED_FLAG) THEN
    CALL LOAD_BLR_MATRIX(GRID, A_BLR, N, NB, PR, PC, MYID0, TRIM(CHECKPOINT_DIR)//'/saved_factored', .TRUE., OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    T_FACTOR = 0.0D0
    N_SWAPS = -1
    CALL RECOMPUTE_COMPRESSION_STATS(GRID, A_BLR, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (factorisation, charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
           ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
      WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (factorisation, charge) : ', &
           REAL_COMPRESSION_PCT(NB, N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
  ELSE
    IF (LOAD_MATRIX_FLAG) THEN
      CALL LOAD_BLR_MATRIX(GRID, A_BLR, N, NB, PR, PC, MYID0, TRIM(CHECKPOINT_DIR)//'/saved_matrix', .FALSE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
      T_ASM = 0.0D0
      CALL RECOMPUTE_COMPRESSION_STATS(GRID, A_BLR, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage, charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage, charge) : ', &
             REAL_COMPRESSION_PCT(NB, N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
    ELSE
      ! --- 1) assemblage BLR distribue, via VOTRE USER_GET_BLOCK ---
      CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, USER_GET_BLOCK, A_BLR, T_ASM, &
                                    N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., USE_RSVD=USE_RSVD_FLAG, &
                                    PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
             ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage) : ', REAL_COMPRESSION_PCT(NB, N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
      IF (SAVE_MATRIX_FLAG) THEN
        CALL SAVE_BLR_MATRIX(GRID, A_BLR, N, NB, PR, PC, MYID0, TRIM(CHECKPOINT_DIR)//'/saved_matrix', .FALSE., OK)
        IF (.NOT. OK) THEN
          CALL MPI_FINALIZE(IERR)
          STOP 1
        END IF
      END IF
    END IF
    CALL REPORT_MEMORY_MB(GRID, 'Matrice A_BLR (compressee)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
    CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
    CALL PROFILE_SAMPLE('assemblage')
    CALL DUMP_BLR_STRUCTURE(GRID, A_BLR, 'assembly')

    ! --- 2) factorisation LU, avec pivot confine au bloc diagonal ---
    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                         N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                         SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE., &
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
      WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (fin factorisation) : ', REAL_COMPRESSION_PCT(NB, N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
    IF (SAVE_FACTORED_FLAG) THEN
      CALL SAVE_BLR_MATRIX(GRID, A_BLR, N, NB, PR, PC, MYID0, TRIM(CHECKPOINT_DIR)//'/saved_factored', .TRUE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU (compresses)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('factorisation')
  CALL DUMP_BLR_STRUCTURE(GRID, A_BLR, 'factorization')

  ! --- 3) second membre B (via VOTRE USER_FILL_B) et resolution 2D, avec
  !     sauvegarde, PAR LOT d'au plus NRHS_BATCH colonnes GLOBALES.
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO

  IF (LOAD_RHS_FLAG) THEN
    CALL OPEN_RHS_CHECKPOINT_READ(GRID, N, NRHS, PR, PC, MYID0, MLOC, TRIM(CHECKPOINT_DIR)//'/saved_rhs', &
                                   LOAD_RHS_UNIT, NBATCHES, NRHS_BATCH_EFF, OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  ELSE
    NRHS_BATCH_EFF = NRHS_BATCH
    IF (NRHS_BATCH_EFF .LE. 0 .OR. NRHS_BATCH_EFF .GT. NRHS) NRHS_BATCH_EFF = NRHS
    NBATCHES = (NRHS + NRHS_BATCH_EFF - 1) / NRHS_BATCH_EFF
    IF (SAVE_RHS_FLAG) THEN
      CALL OPEN_RHS_CHECKPOINT_WRITE(GRID, N, NRHS, PR, PC, MYID0, TRIM(CHECKPOINT_DIR)//'/saved_rhs', SAVE_RHS_UNIT, &
                                      NBATCHES, NRHS_BATCH_EFF, MLOC, OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
  END IF
  IF (SAVE_SOLUTION_FLAG) THEN
    CALL OPEN_RHS_CHECKPOINT_WRITE(GRID, N, NRHS, PR, PC, MYID0, TRIM(CHECKPOINT_DIR)//'/saved_solution', SAVE_SOLUTION_UNIT, &
                                    NBATCHES, NRHS_BATCH_EFF, MLOC, OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  END IF

  NPREVIEW = MIN(5, N)
  PREVIEW_LOCAL = (0.0_WP, 0.0_WP)
  NUM2_ACC = 0.0_WP; DEN2_ACC = 0.0_WP
  T_SOLVE = 0.0D0
  DO BATCH_IDX = 1, NBATCHES
    RSTART_B = (BATCH_IDX-1)*NRHS_BATCH_EFF + 1
    REND_B = MIN(BATCH_IDX*NRHS_BATCH_EFF, NRHS)
    WIDTH_B = REND_B - RSTART_B + 1
    CALL RHS_LOCAL_RANGE_MBM(WIDTH_B, GRID%PC, GRID%MYCOL, LOCAL_OFFSET_B, NRHS_LOC_B)
    GSTART_B = RSTART_B + LOCAL_OFFSET_B - 1
    GEND_B = GSTART_B + NRHS_LOC_B - 1

    ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC_B,1)))
    IF (LOAD_RHS_FLAG) THEN
      CALL LOAD_RHS_BATCH(LOAD_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
    ELSE
      CALL USER_FILL_B(A_BLR, GRID, NRHS_LOC_B, GSTART_B, B_LOCAL, MLOC)
      IF (SAVE_RHS_FLAG) CALL SAVE_RHS_BATCH(SAVE_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
    END IF
    IF (CHECK_RESIDUAL_FLAG) THEN
      ! B_LOCAL sera ecrase EN PLACE par la resolution -- garder une copie
      ! du B REELLEMENT resolu (pas recalcule par formule) pour le residu.
      ALLOCATE(B_ORIG_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC_B,1)))
      B_ORIG_LOCAL = B_LOCAL
    END IF
    IF (NBATCHES.EQ.1) THEN
      CALL PROFILE_SAMPLE('second_membre')
    ELSE
      CALL PROFILE_SAMPLE(PCT_LABEL('second_membre', BATCH_IDX, NBATCHES))
    END IF

    ! Resolution 2D, avec pivot -- deja generique sur la LARGEUR resolue
    ! (WIDTH_B/NRHS_LOC_B DU LOT, pas NRHS global).
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0 = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, WIDTH_B, NRHS_LOC_B, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE., &
                                           PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE)
    T1 = MPI_WTIME()
    T_SOLVE = T_SOLVE + (T1 - T0)
    ! B_LOCAL contient maintenant X de ce lot (ecrase en place).
    IF (SAVE_SOLUTION_FLAG) CALL SAVE_RHS_BATCH(SAVE_SOLUTION_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)

    IF (CHECK_RESIDUAL_FLAG) THEN
      ! Gathere SEULEMENT les WIDTH_B colonnes de ce lot (jamais NRHS
      ! entier) -- l'avantage memoire de NRHS_BATCH reste entier.
      ALLOCATE(XG_BATCH(N,WIDTH_B))
      CALL GATHER_2D_RHS(GRID, N, NB, WIDTH_B, B_LOCAL, MAX(MLOC,1), NRHS_LOC_B, 0, XG_BATCH)
      IF (BATCH_IDX.EQ.1 .AND. MYID0.EQ.0) THEN
        WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
        DO I = 1, MIN(5, N)
          WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(XG_BATCH(I,1)), ' + i*', AIMAG(XG_BATCH(I,1))
        END DO
      END IF
      CALL ACCUMULATE_RESIDUAL(GRID, USER_GET_BLOCK, N, NB, MLOC, XG_BATCH, B_ORIG_LOCAL, WIDTH_B, NRHS_LOC_B, NUM2_ACC, DEN2_ACC)
      DEALLOCATE(XG_BATCH, B_ORIG_LOCAL)
    ELSE IF (GSTART_B .LE. 1 .AND. 1 .LE. GEND_B) THEN
      ! Apercu de X(1:5,1) sans assembler X en entier -- chaque rang
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
  ! Empreinte reelle d'un lot (pas NRHS entier) -- pic de RAM reel, borne
  ! par NRHS_BATCH.
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B (par lot)', &
       8_8*INT(MLOC,8)*INT(RHS_NUMROC_MBM(NRHS_BATCH_EFF,GRID%PC,GRID%MYCOL),8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres resolution)', GET_RSS_KB()*1024_8)

  ! --- 4) verification : residu relatif ||A*X-B||/||B|| (si --check-residual) ---
  IF (CHECK_RESIDUAL_FLAG) THEN
    RESID = SQRT(NUM2_ACC) / SQRT(DEN2_ACC)
    IF (MYID0.EQ.0) WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
    CALL PROFILE_SAMPLE('verif_residu')
  ELSE
    CALL MPI_REDUCE(PREVIEW_LOCAL, PREVIEW_GLOBAL, NPREVIEW, MPI_COMPLEX, MPI_SUM, 0, GRID%COMM, IERR)
    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
      DO I = 1, NPREVIEW
        WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(PREVIEW_GLOBAL(I)), ' + i*', AIMAG(PREVIEW_GLOBAL(I))
      END DO
      WRITE(6,'(A)') ' (X complet non assemble -- CHECK_RESIDUAL=0, memoire economisee)'
    END IF
  END IF

  CALL PROFILE_SAMPLE('fin')
  CALL PROFILE_WRITE_CSV()

  ! --- 5) --reload-solution : liberation d'un maximum de RAM (A_BLR --
  !     B_LOCAL est deja libere a chaque lot, ligne ~426) puis relecture
  !     de la solution depuis disque par le SEUL rang 0, pour verifier
  !     qu'un post-traitement leger (sans A_BLR ni la grille en memoire)
  !     peut recuperer X apres coup. LOAD_SOLUTION est sequentielle : elle
  !     parcourt elle-meme tous les PR*PC fichiers saved_solution_rank*.dat
  !     -- aucun autre rang n'a besoin d'y participer.
  IF (RELOAD_SOLUTION_FLAG) THEN
    CALL REPORT_MEMORY_MB(GRID, 'RSS avant liberation de A_BLR', GET_RSS_KB()*1024_8)
    CALL MPI_BLR_MATRIX_FREE(A_BLR)
    CALL REPORT_MEMORY_MB(GRID, 'RSS apres liberation de A_BLR', GET_RSS_KB()*1024_8)

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MYID0.EQ.0) THEN
      ALLOCATE(X_RELOADED(N,NRHS))
      CALL LOAD_SOLUTION(CHECKPOINT_DIR, N, NRHS, PR, PC, NB, 1, N, 1, NRHS, X_RELOADED, OK)
      IF (.NOT. OK) THEN
        WRITE(6,'(A)') ' ERREUR : relecture de la solution (--reload-solution) a echoue.'
      ELSE
        WRITE(6,'(A)') ' Solution relue depuis disque (rang 0 seul, A_BLR libere) :'
        DO I = 1, MIN(5, N)
          WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_RELOADED(I,1)), ' + i*', AIMAG(X_RELOADED(I,1))
        END DO
      END IF
      DEALLOCATE(X_RELOADED)
    END IF
    CALL REPORT_MEMORY_MB(GRID, 'RSS apres relecture solution (rang 0)', GET_RSS_KB()*1024_8)
  END IF

  CALL MPI_FINALIZE(IERR)

END PROGRAM EXAMPLE_BLR_MPI_PIVOT_USER
