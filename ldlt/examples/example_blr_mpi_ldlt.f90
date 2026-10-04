!===============================================================================
! EXEMPLE -- LDL^T BLR+MPI pour matrice COMPLEXE SYMETRIQUE (A^T = A), avec
! pivot de Bunch-Kaufman borne confine au bloc diagonal (CSYTRF_RK, pivots
! 1x1/2x2). Pendant symetrique de examples/example_blr_mpi_pivot.f90 (LU) :
! memes options (hors sauvegarde/reprise), stockage du seul triangle
! inferieur (~2x moins de memoire que le LU).
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_ldlt PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]
!         (--help/-h pour la liste). PR PC obligatoires ; N/NB/NRHS/TOL
!         positionnels optionnels dans cet ordre (defauts 4000/200/4/1e-4) ;
!         les autres options sont NOMMEES (--cle=valeur, ou --cle seul =
!         --cle=1), dans n'importe quel ordre. Exemple :
!           mpirun -np 4 ./example_blr_mpi_ldlt 2 2 4000 200 8 1e-4 --check-error --use-rsvd
!
! Options (meme sens que pour le LU, voir son commentaire de tete) :
!   --use-move-alloc  MOVE_ALLOC dans la mise a jour de Schur (resultat identique)
!   --eps             bruit plein-rang symetrique (degrade la compressibilite)
!   --ratio           COMPRESS_WORTHWHILE_RATIO (defaut 1/3)
!   --use-rsvd        SVD randomisee adaptative a l'assemblage au lieu de CGESVD
!   --check-error     B = A*Xexact (produit brut), affiche ||X-Xexact||/||Xexact||
!   --check-residual  affiche ||A*X-B||/||B|| avec le B reellement utilise
!   --profile-usage   RSS/CPU par etape -> usage_profile_rank<N>.csv
!                     (trace : python3 ../examples/plot_usage_profile.py)
!   --dump-structure  nature des blocs locaux (triangle inferieur) apres
!                     assemblage/factorisation -> blr_structure_<etape>_rank<N>.csv
!   --nrhs-batch      colonnes de second membre traitees par lot (borne aussi
!                     la copie pleine largeur XFULL de la remontee L^T)
! Sauvegarde / reprise (un fichier par processus dans --checkpoint-dir, qui
! doit exister ; en-tete verifie a la relecture, un rang incompatible arrete
! tous les rangs) :
!   --save-matrix / --load-matrix      matrice assemblee (saute l'assemblage)
!   --save-factored / --load-factored  facteur LDLt (saute assemblage +
!                                      factorisation, prioritaire)
!   --save-rhs / --load-rhs            second membre (relu avec le decoupage
!                                      en lots de la sauvegarde)
!   --save-solution                    solution X
!   --ntailbloc       ecrit/relit second membre et solution par blocs de ce
!                     nombre de lignes, un fichier par bloc
! Parametres de la matrice de demo (sym_testmat_mod.f90) :
!   --wavenumber (defaut 8), --strong (defaut 50, 0 = sans anti-diagonale),
!   --diag (defaut 0.01) -- choisissent le regime de pivots 1x1/2x2.
! La reprise suppose la MEME matrice : --load-* ne relit pas les parametres
! de la matrice de demo, ils doivent etre identiques a la sauvegarde pour
! que --check-error/--check-residual aient un sens.
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_LDLT
  USE MPI_BLR_LDLT_MOD
  USE SYM_TESTMAT_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  TYPE(MPI_BLR_SYM_T) :: A
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0, IARG, POS_IDX
  REAL(WP) :: TOL, EPS, RATIO, WAVENUMBER, STRONG, DIAG
  CHARACTER(LEN=256) :: ARG, VALSTR
  CHARACTER(LEN=64)  :: KEYSTR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, INFO, N_SWAPS, N_2X2, N_COMPRESSED, N_OFFDIAG, li, I
  INTEGER(8) :: RANK_SUM
  INTEGER :: NRHS_BATCH, NRHS_BATCH_EFF, NBATCHES, BATCH_IDX
  INTEGER :: RSTART_B, REND_B, WIDTH_B, LOCAL_OFFSET_B, GSTART_B, GEND_B, NRHS_LOC_B
  INTEGER :: NPREVIEW, PREVIEW_RS, PREVIEW_MI, PREVIEW_GI
  COMPLEX(WP) :: PREVIEW_LOCAL(5), PREVIEW_GLOBAL(5)
  LOGICAL :: WANT_FULL_X
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0
  INTEGER :: USE_MOVE_ALLOC_INT, USE_RSVD_INT, CHECK_ERROR_INT, CHECK_RESIDUAL_INT, PROFILE_USAGE_INT, DUMP_STRUCTURE_INT
  LOGICAL :: USE_MOVE_ALLOC_FLAG, USE_RSVD_FLAG, CHECK_ERROR_FLAG, CHECK_RESIDUAL_FLAG, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG
  REAL(WP) :: ERR, RESID
  INTEGER :: SAVE_MATRIX_INT, SAVE_RHS_INT, SAVE_FACTORED_INT, SAVE_SOLUTION_INT
  INTEGER :: LOAD_MATRIX_INT, LOAD_FACTORED_INT, LOAD_RHS_INT, NTAILBLOC
  LOGICAL :: SAVE_MATRIX_FLAG, SAVE_RHS_FLAG, SAVE_FACTORED_FLAG, SAVE_SOLUTION_FLAG
  LOGICAL :: LOAD_MATRIX_FLAG, LOAD_FACTORED_FLAG, LOAD_RHS_FLAG, OK
  CHARACTER(LEN=256) :: CHECKPOINT_DIR
  CHARACTER(LEN=8) :: BATCH_STR
  INTEGER :: SAVE_RHS_UNIT, SAVE_SOLUTION_UNIT, LOAD_RHS_UNIT

  INTEGER, PARAMETER :: MAX_PROFILE_SAMPLES = 128
  CHARACTER(LEN=40) :: PROF_STAGE(MAX_PROFILE_SAMPLES)
  DOUBLE PRECISION :: PROF_T(MAX_PROFILE_SAMPLES), PROF_RSS_MB(MAX_PROFILE_SAMPLES), PROF_CPU_PCT(MAX_PROFILE_SAMPLES)
  INTEGER :: PROF_N
  DOUBLE PRECISION :: PROF_T_START, PROF_T_PREV
  REAL :: PROF_CPU_PREV

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP; EPS = 0.0_WP; RATIO = 1.0_WP/3.0_WP
  WAVENUMBER = 8.0_WP; STRONG = 50.0_WP; DIAG = 0.01_WP
  USE_MOVE_ALLOC_INT = 0; USE_RSVD_INT = 1; CHECK_ERROR_INT = 0; CHECK_RESIDUAL_INT = 0
  PROFILE_USAGE_INT = 0; DUMP_STRUCTURE_INT = 0; NRHS_BATCH = 0; PROF_N = 0
  SAVE_MATRIX_INT = 0; SAVE_RHS_INT = 0; SAVE_FACTORED_INT = 0; SAVE_SOLUTION_INT = 0
  LOAD_MATRIX_INT = 0; LOAD_FACTORED_INT = 0; LOAD_RHS_INT = 0; NTAILBLOC = 0; CHECKPOINT_DIR = '.'

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
      IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : trop d''arguments positionnels (voir --help).'
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END SELECT
    IARG = IARG + 1
  END DO
  DO WHILE (IARG .LE. COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(IARG, ARG)
    IF (TRIM(ARG).EQ.'--help' .OR. TRIM(ARG).EQ.'-h') THEN
      IF (MYID0.EQ.0) CALL PRINT_USAGE()
      CALL MPI_FINALIZE(IERR)
      STOP 0
    END IF
    IF (ARG(1:2) .NE. '--') THEN
      IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : argument inattendu "', TRIM(ARG), '" (voir --help).'
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
    CASE ('nrhs-batch');      READ(VALSTR,*) NRHS_BATCH
    CASE ('wavenumber');      READ(VALSTR,*) WAVENUMBER
    CASE ('strong');          READ(VALSTR,*) STRONG
    CASE ('diag');            READ(VALSTR,*) DIAG
    CASE ('save-matrix');     READ(VALSTR,*) SAVE_MATRIX_INT
    CASE ('save-rhs');        READ(VALSTR,*) SAVE_RHS_INT
    CASE ('save-factored');   READ(VALSTR,*) SAVE_FACTORED_INT
    CASE ('save-solution');   READ(VALSTR,*) SAVE_SOLUTION_INT
    CASE ('load-matrix');     READ(VALSTR,*) LOAD_MATRIX_INT
    CASE ('load-factored');   READ(VALSTR,*) LOAD_FACTORED_INT
    CASE ('load-rhs');        READ(VALSTR,*) LOAD_RHS_INT
    CASE ('checkpoint-dir');  CHECKPOINT_DIR = VALSTR
    CASE ('ntailbloc');       READ(VALSTR,*) NTAILBLOC
    CASE DEFAULT
      IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : option inconnue "--', TRIM(KEYSTR), '" (voir --help).'
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END SELECT
    IARG = IARG + 1
  END DO
  USE_MOVE_ALLOC_FLAG = (USE_MOVE_ALLOC_INT .NE. 0)
  USE_RSVD_FLAG       = (USE_RSVD_INT .NE. 0)
  CHECK_ERROR_FLAG    = (CHECK_ERROR_INT .NE. 0)
  CHECK_RESIDUAL_FLAG = (CHECK_RESIDUAL_INT .NE. 0)
  PROFILE_USAGE_FLAG  = (PROFILE_USAGE_INT .NE. 0)
  DUMP_STRUCTURE_FLAG = (DUMP_STRUCTURE_INT .NE. 0)
  SAVE_MATRIX_FLAG    = (SAVE_MATRIX_INT .NE. 0)
  SAVE_RHS_FLAG       = (SAVE_RHS_INT .NE. 0)
  SAVE_FACTORED_FLAG  = (SAVE_FACTORED_INT .NE. 0)
  SAVE_SOLUTION_FLAG  = (SAVE_SOLUTION_INT .NE. 0)
  LOAD_MATRIX_FLAG    = (LOAD_MATRIX_INT .NE. 0)
  LOAD_FACTORED_FLAG  = (LOAD_FACTORED_INT .NE. 0)
  LOAD_RHS_FLAG       = (LOAD_RHS_INT .NE. 0)

  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)
  CALL SYM_TESTMAT_SETUP(N, NB, WAVENUMBER=WAVENUMBER, EPS=EPS, STRONG=STRONG, DIAG=DIAG)

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
    WRITE(6,'(A,F7.2,A,F7.2,A,F7.3)') ' Matrice complexe symetrique : WAVENUMBER=', WAVENUMBER, &
         ' , STRONG=', STRONG, ' , DIAG=', DIAG
    IF (USE_RSVD_FLAG) THEN
      WRITE(6,'(A)') ' Compression : SVD randomisee adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE)'
    ELSE
      WRITE(6,'(A)') ' Compression : SVD complete (COMPRESS_BLOCK, CGESVD)'
    END IF
    IF (CHECK_ERROR_FLAG) WRITE(6,'(A)') ' Verification : B=A*Xexact, erreur affichee en fin de resolution'
    IF (CHECK_RESIDUAL_FLAG) WRITE(6,'(A)') ' Verification : residu ||A*X-B||/||B|| affiche en fin de resolution'
    IF (PROFILE_USAGE_FLAG) WRITE(6,'(A)') ' Profilage RSS/CPU par etape active : usage_profile_rank<N>.csv par processus'
    IF (DUMP_STRUCTURE_FLAG) WRITE(6,'(A)') ' Structure BLR exportee : blr_structure_<etape>_rank<N>.csv par processus'
    IF (NRHS_BATCH.GT.0 .AND. NRHS_BATCH.LT.NRHS .AND. .NOT.LOAD_RHS_FLAG) &
         WRITE(6,'(A,I0,A)') ' Second membre/resolution par lots de ', NRHS_BATCH, ' colonnes (NRHS_BATCH)'
    IF (LOAD_FACTORED_FLAG) THEN
      WRITE(6,'(A)') ' Reprise : facteur LDLt relu (LOAD_FACTORED=1) -- assemblage et factorisation ignores'
      IF (LOAD_MATRIX_FLAG) WRITE(6,'(A)') ' (LOAD_MATRIX=1 ignore : LOAD_FACTORED est prioritaire)'
    ELSE IF (LOAD_MATRIX_FLAG) THEN
      WRITE(6,'(A)') ' Reprise : matrice relue (LOAD_MATRIX=1) -- assemblage ignore'
    END IF
    IF (LOAD_RHS_FLAG) WRITE(6,'(A)') ' Reprise : second membre relu (LOAD_RHS=1)'
    IF (SAVE_MATRIX_FLAG .OR. SAVE_RHS_FLAG .OR. SAVE_FACTORED_FLAG .OR. SAVE_SOLUTION_FLAG) &
         WRITE(6,'(A,A)') ' Sauvegarde(s) activee(s) dans : ', TRIM(CHECKPOINT_DIR)
  END IF

  PROF_T_START = MPI_WTIME()
  PROF_T_PREV = 0.0D0
  CALL CPU_TIME(PROF_CPU_PREV)
  CALL PROFILE_SAMPLE('debut')
  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  IF (LOAD_FACTORED_FLAG) THEN
    CALL MPI_BLR_SYM_LOAD(GRID, A, N, NB, TRIM(CHECKPOINT_DIR)//'/saved_ldlt_factored', .TRUE., OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    T_FACTOR = 0.0D0
    CALL MPI_BLR_SYM_COMPRESSION_STATS(GRID, A, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
    CALL COUNT_ALL_PIVOTS(N_SWAPS, N_2X2)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,I0,A,I0)') ' Echanges de pivot (charges) : ', N_SWAPS, ' , pivots 2x2 : ', N_2X2
      WRITE(6,'(A,I0,A,I0,A,F6.2,A)') ' Blocs L compresses (charges) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
           ' , taux REEL ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
  ELSE
    IF (LOAD_MATRIX_FLAG) THEN
      CALL MPI_BLR_SYM_LOAD(GRID, A, N, NB, TRIM(CHECKPOINT_DIR)//'/saved_ldlt_matrix', .FALSE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
      T_ASM = 0.0D0
      CALL MPI_BLR_SYM_COMPRESSION_STATS(GRID, A, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
      IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,F6.2,A)') ' Blocs compresses (matrice chargee) : ', N_COMPRESSED, &
           ' / ', N_OFFDIAG, ' , taux REEL ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    ELSE
      ! --- 1) assemblage BLR distribue, triangle inferieur seulement ---
      CALL MPI_ASSEMBLE_BLR_SYM(GRID, N, NB, TOL, SYM_GET_BLOCK, A, T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                SHOW_PROGRESS=.TRUE., USE_RSVD=USE_RSVD_FLAG, PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY)
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
        WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage, triangle inferieur) : ', N_COMPRESSED, ' / ', &
             N_OFFDIAG, ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
        WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage) : ', &
             REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
      END IF
      IF (SAVE_MATRIX_FLAG) THEN
        CALL MPI_BLR_SYM_SAVE(GRID, A, TRIM(CHECKPOINT_DIR)//'/saved_ldlt_matrix', .FALSE., OK)
        IF (.NOT. OK) THEN
          CALL MPI_FINALIZE(IERR)
          STOP 1
        END IF
      END IF
    END IF
    CALL REPORT_MEMORY_MB(GRID, 'Matrice A (triangle inferieur, compressee)', MPI_BLR_SYM_MEMORY_BYTES(A))
    CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
    CALL PROFILE_SAMPLE('assemblage')
    CALL DUMP_BLR_STRUCTURE('assembly')

    ! --- 2) factorisation LDL^T, Bunch-Kaufman confine au bloc diagonal ---
    CALL MPI_BLR_SYM_LDLT_FACTOR(GRID, A, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, N_2X2, T_FACTOR, &
                                 SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE., USE_MOVE_ALLOC=USE_MOVE_ALLOC_FLAG, &
                                 PROGRESS_CALLBACK=PROFILE_CALLBACK_FACTOR)
    IF (INFO .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_SYM_LDLT_FACTOR (D singuliere), INFO=', INFO
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,F8.3,A)') ' Factorisation LDLt : ', T_FACTOR, ' s'
      WRITE(6,'(A,I0,A,I0)') ' Echanges de pivot : ', N_SWAPS, ' , pivots 2x2 : ', N_2X2
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs L compresses (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
           ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
      WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (fin factorisation) : ', &
           REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
    END IF
    IF (SAVE_FACTORED_FLAG) THEN
      CALL MPI_BLR_SYM_SAVE(GRID, A, TRIM(CHECKPOINT_DIR)//'/saved_ldlt_factored', .TRUE., OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs L + D (compresses)', MPI_BLR_SYM_MEMORY_BYTES(A))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('factorisation')
  CALL DUMP_BLR_STRUCTURE('factorization')

  ! --- 3) second membre et resolution, PAR LOT de NRHS_BATCH colonnes ---
  MLOC = MPI_BLR_SYM_MLOC(A, GRID)
  IF (LOAD_RHS_FLAG) THEN
    ! le decoupage stocke a la sauvegarde fait foi
    CALL OPEN_RHS_CHECKPOINT_READ(TRIM(CHECKPOINT_DIR)//'/saved_ldlt_rhs', LOAD_RHS_UNIT, NBATCHES, NRHS_BATCH_EFF, OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  ELSE
    NRHS_BATCH_EFF = NRHS_BATCH
    IF (NRHS_BATCH_EFF .LE. 0 .OR. NRHS_BATCH_EFF .GT. NRHS) NRHS_BATCH_EFF = NRHS
    NBATCHES = (NRHS + NRHS_BATCH_EFF - 1) / NRHS_BATCH_EFF
    IF (SAVE_RHS_FLAG) THEN
      CALL OPEN_RHS_CHECKPOINT_WRITE(TRIM(CHECKPOINT_DIR)//'/saved_ldlt_rhs', SAVE_RHS_UNIT, NBATCHES, &
                                     NRHS_BATCH_EFF, MLOC, OK)
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    END IF
  END IF
  IF (SAVE_SOLUTION_FLAG) THEN
    CALL OPEN_RHS_CHECKPOINT_WRITE(TRIM(CHECKPOINT_DIR)//'/saved_ldlt_solution', SAVE_SOLUTION_UNIT, NBATCHES, &
                                   NRHS_BATCH_EFF, MLOC, OK)
    IF (.NOT. OK) THEN
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
  END IF
  WANT_FULL_X = (CHECK_ERROR_FLAG .OR. CHECK_RESIDUAL_FLAG)
  IF (WANT_FULL_X) ALLOCATE(X_GLOBAL(N,NRHS))
  NPREVIEW = MIN(5, N)
  PREVIEW_LOCAL = (0.0_WP, 0.0_WP)
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
        CALL LOAD_ROWBLOCKS(TRIM(CHECKPOINT_DIR)//'/saved_ldlt_rhs_batch'//BATCH_STR, B_LOCAL, MLOC, NRHS_LOC_B, &
                            NTAILBLOC, OK)
      ELSE
        CALL LOAD_RHS_BATCH(LOAD_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
        OK = .TRUE.
      END IF
      IF (.NOT. OK) THEN
        CALL MPI_FINALIZE(IERR)
        STOP 1
      END IF
    ELSE
      CALL FILL_B(NRHS_LOC_B, GSTART_B, B_LOCAL)
      IF (SAVE_RHS_FLAG) THEN
        IF (NTAILBLOC .GT. 0) THEN
          CALL SAVE_ROWBLOCKS(TRIM(CHECKPOINT_DIR)//'/saved_ldlt_rhs_batch'//BATCH_STR, B_LOCAL, MLOC, NRHS_LOC_B, &
                              NTAILBLOC, OK)
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
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0 = MPI_WTIME()
    CALL MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, A, WIDTH_B, NRHS_LOC_B, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE., &
                                   PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE)
    T_SOLVE = T_SOLVE + (MPI_WTIME() - T0)
    IF (SAVE_SOLUTION_FLAG) THEN
      IF (NTAILBLOC .GT. 0) THEN
        CALL SAVE_ROWBLOCKS(TRIM(CHECKPOINT_DIR)//'/saved_ldlt_solution_batch'//BATCH_STR, B_LOCAL, MLOC, &
                            NRHS_LOC_B, NTAILBLOC, OK)
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
      PREVIEW_RS = 1
      DO li = 1, A%NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        PREVIEW_MI = A%BEND(I) - A%BSTART(I) + 1
        DO PREVIEW_GI = 1, NPREVIEW
          IF (PREVIEW_GI.GE.A%BSTART(I) .AND. PREVIEW_GI.LE.A%BEND(I)) &
               PREVIEW_LOCAL(PREVIEW_GI) = B_LOCAL(PREVIEW_RS + (PREVIEW_GI-A%BSTART(I)), 1-GSTART_B+1)
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
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B (par lot)', &
       8_8*INT(MLOC,8)*INT(RHS_NUMROC_MBM(NRHS_BATCH_EFF,GRID%PC,GRID%MYCOL),8))
  CALL REPORT_MEMORY_MB(GRID, 'Copie pleine largeur XFULL (remontee L^T, par lot)', &
       8_8*INT(MLOC,8)*INT(NRHS_BATCH_EFF,8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres resolution)', GET_RSS_KB()*1024_8)

  IF (WANT_FULL_X) THEN
    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
      DO I = 1, NPREVIEW
        WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
      END DO
      IF (CHECK_ERROR_FLAG) THEN
        ERR = COMPUTE_ERR(X_GLOBAL)
        WRITE(6,'(A,ES12.4)') ' Erreur relative ||X-Xexact|| / ||Xexact|| = ', ERR
        CALL PROFILE_SAMPLE('verif_erreur')
      END IF
    END IF
    IF (CHECK_RESIDUAL_FLAG) THEN
      RESID = COMPUTE_RESIDUAL(X_GLOBAL)     ! collectif
      IF (MYID0.EQ.0) WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
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
  CALL MPI_BLR_SYM_FREE(A)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Echanges / pivots 2x2 recomptes depuis A%IPIV (apres LOAD_FACTORED).
  SUBROUTINE COUNT_ALL_PIVOTS(NSW, N2)
    INTEGER, INTENT(OUT) :: NSW, N2
    INTEGER :: KK, KW2, NS_K, N2_K
    NSW = 0; N2 = 0
    DO KK = 1, A%NBLK
      KW2 = A%BEND(KK) - A%BSTART(KK) + 1
      CALL COUNT_PIVOTS(KW2, A%IPIV(1:KW2,KK), A%DSUB(1:KW2,KK), NS_K, N2_K)
      NSW = NSW + NS_K; N2 = N2 + N2_K
    END DO
  END SUBROUTINE COUNT_ALL_PIVOTS

  ! --- Checkpoints du second membre / de la solution : reprises a
  !     l'identique de examples/example_blr_mpi_pivot.f90 (memes formats).
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

  SUBROUTINE PRINT_USAGE()
    WRITE(6,'(A)') 'Usage : mpirun -np PR*PC ./example_blr_mpi_ldlt PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  PR PC          grille de processus MPI (obligatoire, PR*PC = nb de processus lances)'
    WRITE(6,'(A)') '  N NB NRHS TOL  positionnels optionnels, DANS CET ORDRE (defauts 4000 200 4 1e-4)'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Options nommees (--cle=valeur, ou --cle seul = --cle=1, dans n''importe quel ordre) :'
    WRITE(6,'(A)') '  --n= --nb= --nrhs= --tol=   equivalents aux positionnels'
    WRITE(6,'(A)') '  --use-move-alloc=0|1    (defaut 0) transfert MOVE_ALLOC au lieu de copie (Schur)'
    WRITE(6,'(A)') '  --eps=REEL              (defaut 0.0) bruit plein-rang symetrique (compressibilite)'
    WRITE(6,'(A)') '  --ratio=REEL            (defaut 0.333) COMPRESS_WORTHWHILE_RATIO'
    WRITE(6,'(A)') '  --use-rsvd=0|1          (defaut 1) SVD randomisee adaptative ; 0 = CGESVD complete'
    WRITE(6,'(A)') '  --check-error=0|1       (defaut 0) verification B=A*Xexact'
    WRITE(6,'(A)') '  --check-residual=0|1    (defaut 0) verification ||A*X-B||/||B||'
    WRITE(6,'(A)') '  --profile-usage=0|1     (defaut 0) trace RSS/CPU par etape (usage_profile_rank<N>.csv)'
    WRITE(6,'(A)') '  --dump-structure=0|1    (defaut 0) export structure BLR (blr_structure_<etape>_rank<N>.csv)'
    WRITE(6,'(A)') '  --nrhs-batch=ENTIER     (defaut 0 = pas de decoupage) colonnes RHS traitees par lot'
    WRITE(6,'(A)') '  --wavenumber=REEL       (defaut 8) nombre d''onde du noyau de Helmholtz'
    WRITE(6,'(A)') '  --strong=REEL           (defaut 50) anti-diagonale forte par bloc (0 = aucune)'
    WRITE(6,'(A)') '  --diag=REEL             (defaut 0.01) diagonale (grande -> pivots 1x1, faible -> 2x2)'
    WRITE(6,'(A)') '  --save-matrix=0|1 --save-factored=0|1 --save-rhs=0|1 --save-solution=0|1'
    WRITE(6,'(A)') '  --load-matrix=0|1 --load-factored=0|1 --load-rhs=0|1   sauvegarde / reprise'
    WRITE(6,'(A)') '  --checkpoint-dir=CHEMIN (defaut .) repertoire des checkpoints (doit deja exister)'
    WRITE(6,'(A)') '  --ntailbloc=ENTIER      (defaut 0) second membre/solution ecrits par blocs de lignes'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Exemple :'
    WRITE(6,'(A)') '  mpirun -np 4 ./example_blr_mpi_ldlt 2 2 4000 200 8 1e-4 --check-error --use-rsvd'
  END SUBROUTINE PRINT_USAGE

  ! Taux de compression REEL (meme formule que l'exemple LU) : blocs
  ! denses retires de RANK_SUM (rang virtuel NB) puis comptes a NB^2.
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

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

  ! Meme format CSV que l'exemple LU (plot_blr_structure.py) ; seuls les
  ! blocs du triangle inferieur existent.
  SUBROUTINE DUMP_BLR_STRUCTURE(STAGE_LABEL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_LABEL
    INTEGER :: U, lli, llj, II, JJ
    CHARACTER(LEN=96) :: FNAME
    IF (.NOT. DUMP_STRUCTURE_FLAG) RETURN
    WRITE(FNAME,'(A,A,A,I0,A)') 'blr_structure_', TRIM(STAGE_LABEL), '_rank', MYID0, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'block_row,block_col,kind,rank'
    DO II = 1, A%NBLK
      IF (ALLOCATED(A%DIAGP(II)%AP)) WRITE(U,'(I0,A,I0,A,I0)') II, ',', II, ',diagonal,', A%BEND(II)-A%BSTART(II)+1
    END DO
    DO lli = 1, A%NBLKR_LOC
      II = (lli-1)*GRID%PR + GRID%MYROW + 1
      DO llj = 1, A%NBLKC_LOC
        JJ = BLR_COL_GLOBAL(A%CMAP, llj, GRID)
        IF (II .LE. JJ) CYCLE
        IF (A%OFFBLK(lli,llj)%IS_LOWRANK) THEN
          WRITE(U,'(I0,A,I0,A,I0)') II, ',', JJ, ',lowrank,', A%OFFBLK(lli,llj)%RANK
        ELSE
          WRITE(U,'(I0,A,I0,A,I0)') II, ',', JJ, ',dense,', A%BEND(JJ)-A%BSTART(JJ)+1
        END IF
      END DO
    END DO
    CLOSE(U)
  END SUBROUTINE DUMP_BLR_STRUCTURE

  ! B par defaut (CHECK_ERROR=0) -- source unique pour FILL_B et le residu.
  COMPLEX(WP) FUNCTION B_ENTRY(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    B_ENTRY = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP), KIND=WP)
  END FUNCTION B_ENTRY

  ! B_LOCAL du lot : B_ENTRY, ou A*Xexact (produit brut sur la formule de A,
  ! par blocs de lignes -- jamais les blocs compresses) si CHECK_ERROR=1.
  SUBROUTINE FILL_B(NRL, RSTART, BL)
    INTEGER, INTENT(IN)      :: NRL, RSTART
    COMPLEX(WP), INTENT(OUT) :: BL(MAX(MLOC,1),MAX(NRL,1))
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: AROWS(:,:), XEX(:,:)
    INTEGER :: lli, II, BSZ, RS, GI, RR
    IF (NRL .LE. 0) RETURN
    IF (CHECK_ERROR_FLAG) THEN
      ALLOCATE(XEX(N,NRL))
      DO RR = 1, NRL
        DO GI = 1, N
          XEX(GI,RR) = SYM_XEX(GI, RSTART+RR-1)
        END DO
      END DO
    END IF
    RS = 1
    DO lli = 1, A%NBLKR_LOC
      II = (lli-1)*GRID%PR + GRID%MYROW + 1
      BSZ = A%BEND(II) - A%BSTART(II) + 1
      IF (CHECK_ERROR_FLAG) THEN
        ALLOCATE(AROWS(BSZ,N))
        CALL SYM_GET_BLOCK(A%BSTART(II), A%BEND(II), 1, N, AROWS)
        CALL CGEMM('N', 'N', BSZ, NRL, N, CONE, AROWS, BSZ, XEX, N, CZERO, BL(RS,1), MAX(MLOC,1))
        DEALLOCATE(AROWS)
      ELSE
        DO GI = A%BSTART(II), A%BEND(II)
          DO RR = 1, NRL
            BL(RS + (GI-A%BSTART(II)), RR) = B_ENTRY(GI, RSTART+RR-1)
          END DO
        END DO
      END IF
      RS = RS + BSZ
    END DO
    IF (ALLOCATED(XEX)) DEALLOCATE(XEX)
  END SUBROUTINE FILL_B

  REAL(WP) FUNCTION COMPUTE_ERR(XG)
    COMPLEX(WP), INTENT(IN) :: XG(N,NRHS)
    REAL(WP) :: NUM, DEN
    INTEGER :: GI, R
    NUM = 0.0_WP; DEN = 0.0_WP
    DO R = 1, NRHS
      DO GI = 1, N
        NUM = NUM + ABS(XG(GI,R) - SYM_XEX(GI,R))**2
        DEN = DEN + ABS(SYM_XEX(GI,R))**2
      END DO
    END DO
    COMPUTE_ERR = SQRT(NUM/DEN)
  END FUNCTION COMPUTE_ERR

  ! ||A*X-B||/||B|| (collectif) : X diffuse, lignes reparties entre tous
  ! les rangs, A materialisee par blocs de NB lignes (formule exacte).
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG)
    COMPLEX(WP), INTENT(INOUT) :: XG(N,NRHS)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: AROWS(:,:), AX(:,:), BEX(:,:), XEX(:,:)
    INTEGER :: ROWS_PER_PROC, R0, R1, RC0, RC1, NCH, GI, R, IERR2
    REAL(WP) :: NUM2_L, DEN2_L, NUM2, DEN2
    COMPLEX(WP) :: BVAL
    CALL MPI_BCAST(XG, N*NRHS, MPI_COMPLEX, 0, GRID%COMM, IERR2)
    ROWS_PER_PROC = (N + GRID%NPROCS - 1) / GRID%NPROCS
    R0 = GRID%MYID*ROWS_PER_PROC + 1
    R1 = MIN((GRID%MYID+1)*ROWS_PER_PROC, N)
    IF (CHECK_ERROR_FLAG) THEN
      ALLOCATE(XEX(N,NRHS))
      DO R = 1, NRHS
        DO GI = 1, N
          XEX(GI,R) = SYM_XEX(GI,R)
        END DO
      END DO
    END IF
    NUM2_L = 0.0_WP; DEN2_L = 0.0_WP
    DO RC0 = R0, R1, NB
      RC1 = MIN(RC0+NB-1, R1)
      NCH = RC1 - RC0 + 1
      ALLOCATE(AROWS(NCH,N), AX(NCH,NRHS))
      CALL SYM_GET_BLOCK(RC0, RC1, 1, N, AROWS)
      CALL CGEMM('N', 'N', NCH, NRHS, N, CONE, AROWS, NCH, XG, N, CZERO, AX, NCH)
      IF (CHECK_ERROR_FLAG) THEN
        ALLOCATE(BEX(NCH,NRHS))
        CALL CGEMM('N', 'N', NCH, NRHS, N, CONE, AROWS, NCH, XEX, N, CZERO, BEX, NCH)
      END IF
      DO R = 1, NRHS
        DO GI = 1, NCH
          IF (CHECK_ERROR_FLAG) THEN
            BVAL = BEX(GI,R)
          ELSE
            BVAL = B_ENTRY(RC0+GI-1, R)
          END IF
          NUM2_L = NUM2_L + ABS(AX(GI,R)-BVAL)**2
          DEN2_L = DEN2_L + ABS(BVAL)**2
        END DO
      END DO
      DEALLOCATE(AROWS, AX)
      IF (ALLOCATED(BEX)) DEALLOCATE(BEX)
    END DO
    IF (ALLOCATED(XEX)) DEALLOCATE(XEX)
    CALL MPI_ALLREDUCE(NUM2_L, NUM2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_L, DEN2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2/DEN2)
  END FUNCTION COMPUTE_RESIDUAL

END PROGRAM EXAMPLE_BLR_MPI_LDLT
