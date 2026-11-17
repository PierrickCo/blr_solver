!===============================================================================
! Solveur BLR + MPI, LU avec pivot ou LDLt (choix a l'execution, --solver),
! sur une matrice A et un second membre B generes par PCMATGEN (generateur
! de test ScaLAPACK, voir matrix_gen_mod.f90).
!
! Etapes : 1) assemblage compresse de A (ou relecture), 2) factorisation
! (ou relecture), 3) generation de B et resolution par lots de seconds
! membres, 4) verification optionnelle ||A*X-B||/||B||.
!
! Usage : mpirun -np PR*PC ./solveur_blr PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]
!         (--help pour la liste des options).
!===============================================================================
PROGRAM SOLVEUR_BLR
  USE MPI_BLR_MATRIX_MOD
  USE MPI_BLR_LDLT_MOD, ONLY: MPI_BLR_SYM_T, MPI_ASSEMBLE_BLR_SYM, MPI_BLR_SYM_LDLT_FACTOR, &
                              MPI_BLR_SYM_LDLT_SOLVE_2D, MPI_BLR_SYM_SAVE, MPI_BLR_SYM_LOAD, &
                              MPI_BLR_SYM_COMPRESSION_STATS, MPI_BLR_SYM_MEMORY_BYTES, &
                              MPI_BLR_SYM_MLOC, MPI_BLR_SYM_FREE, COUNT_PIVOTS
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: PANEL_PACK_MODE
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  USE MATRIX_GEN_MOD, ONLY: GEN_INIT, GEN_GET_A_BLOCK, GEN_FILL_B
  USE OPTIONS_MOD, ONLY: OPTIONS_T, PARSE_OPTIONS, SET_ENV
  USE CHECKPOINT_IO_MOD
  USE PROFILING_MOD
  USE RESIDUAL_MOD, ONLY: ACCUMULATE_RESIDUAL
  USE MPI
  IMPLICIT NONE

  TYPE(OPTIONS_T) :: OPT
  TYPE(GRID_T) :: GRID
  TYPE(MPI_BLR_MATRIX_T) :: A_LU
  TYPE(MPI_BLR_SYM_T)    :: A_SYM
  INTEGER :: IERR, MYID0, STAT, INFO, N_SWAPS, N_2X2, N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0
  LOGICAL :: OK
  CHARACTER(LEN=8) :: SNAME
  CHARACTER(LEN=300) :: PFX

  ! Second membre / resolution par lots
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), B_ORIG_LOCAL(:,:), XG_BATCH(:,:), X_RELOADED(:,:)
  INTEGER :: MLOC, NBLKR_LOC, NRHS_BATCH_EFF, NBATCHES, BATCH_IDX
  INTEGER :: RSTART_B, REND_B, WIDTH_B, LOCAL_OFFSET_B, GSTART_B, GEND_B, NRHS_LOC_B
  INTEGER :: SAVE_RHS_UNIT, SAVE_SOLUTION_UNIT, LOAD_RHS_UNIT, I, li, RS, MI
  INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
  COMPLEX(WP) :: PREVIEW_LOCAL(5), PREVIEW_GLOBAL(5)
  INTEGER :: NPREVIEW
  REAL(WP) :: NUM2_ACC, DEN2_ACC

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  CALL PARSE_OPTIONS(MYID0, OPT, STAT)
  IF (STAT .NE. 0) THEN
    CALL MPI_FINALIZE(IERR)
    IF (STAT .EQ. 1) STOP
    ERROR STOP 1
  END IF

  ! --- options transmises a la bibliotheque ---
  ORTHO_RECOMP = OPT%ORTHO_RECOMP
  PANEL_PACK_MODE = MERGE(1, 0, OPT%PANEL_PACK)
  FACT_STEP_TRACE = OPT%TRACE_FACT
  IF (OPT%DUMP_WORK) CALL SET_ENV('BLR_DUMP_WORK', '1')
  IF (.NOT. OPT%MALLOC_TRIM) CALL SET_ENV('BLR_MALLOC_TRIM', '0')

  CALL GEN_INIT(OPT%N, OPT%NB, OPT%NRHS, OPT%LDLT)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, OPT%PR, OPT%PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : PR*PC doit egaler le nombre de processus lances.'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  SNAME = MERGE('LDLt    ', 'LU      ', OPT%LDLT)
  IF (MYID0.EQ.0) CALL PRINT_CONFIG()
  CALL PROFILING_INIT(MYID0, OPT%PROFILE_USAGE, OPT%DUMP_STRUCTURE)
  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  !=============================================================================
  ! 1) et 2) : matrice A -- assemblage (ou relecture) puis factorisation (ou
  ! relecture du facteur).
  !=============================================================================
  IF (OPT%LDLT) THEN
    CALL SETUP_LDLT()
  ELSE
    CALL SETUP_LU()
  END IF

  !=============================================================================
  ! 3) second membre B (PCMATGEN) et resolution, par lots de NRHS_BATCH
  ! colonnes globales.
  !=============================================================================
  IF (OPT%LDLT) THEN
    NBLKR_LOC = A_SYM%NBLKR_LOC
    ALLOCATE(BSTART(SIZE(A_SYM%BSTART)), BEND(SIZE(A_SYM%BEND)))
    BSTART = A_SYM%BSTART; BEND = A_SYM%BEND
    MLOC = MPI_BLR_SYM_MLOC(A_SYM, GRID)
  ELSE
    NBLKR_LOC = A_LU%NBLKR_LOC
    ALLOCATE(BSTART(SIZE(A_LU%BSTART)), BEND(SIZE(A_LU%BEND)))
    BSTART = A_LU%BSTART; BEND = A_LU%BEND
    MLOC = 0
    DO li = 1, NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + (BEND(I)-BSTART(I)+1)
    END DO
  END IF

  IF (OPT%LOAD_RHS) THEN
    ! Le decoupage en lots de la sauvegarde fait foi.
    CALL OPEN_RHS_CHECKPOINT_READ(GRID, OPT%N, OPT%NRHS, OPT%PR, OPT%PC, MYID0, MLOC, &
                                  TRIM(OPT%CHECKPOINT_DIR)//'/saved_rhs', LOAD_RHS_UNIT, NBATCHES, NRHS_BATCH_EFF, OK)
    IF (.NOT. OK) CALL ABORT_RUN()
  ELSE
    NRHS_BATCH_EFF = OPT%NRHS_BATCH
    IF (NRHS_BATCH_EFF .LE. 0 .OR. NRHS_BATCH_EFF .GT. OPT%NRHS) NRHS_BATCH_EFF = OPT%NRHS
    NBATCHES = (OPT%NRHS + NRHS_BATCH_EFF - 1) / NRHS_BATCH_EFF
    IF (OPT%SAVE_RHS) THEN
      CALL OPEN_RHS_CHECKPOINT_WRITE(GRID, OPT%N, OPT%NRHS, OPT%PR, OPT%PC, MYID0, &
                                     TRIM(OPT%CHECKPOINT_DIR)//'/saved_rhs', SAVE_RHS_UNIT, &
                                     NBATCHES, NRHS_BATCH_EFF, MLOC, OK)
      IF (.NOT. OK) CALL ABORT_RUN()
    END IF
  END IF
  IF (OPT%SAVE_SOLUTION) THEN
    CALL OPEN_RHS_CHECKPOINT_WRITE(GRID, OPT%N, OPT%NRHS, OPT%PR, OPT%PC, MYID0, &
                                   TRIM(OPT%CHECKPOINT_DIR)//'/saved_solution', SAVE_SOLUTION_UNIT, &
                                   NBATCHES, NRHS_BATCH_EFF, MLOC, OK)
    IF (.NOT. OK) CALL ABORT_RUN()
  END IF

  NPREVIEW = MIN(5, OPT%N)
  PREVIEW_LOCAL = (0.0_WP, 0.0_WP)
  NUM2_ACC = 0.0_WP; DEN2_ACC = 0.0_WP
  T_SOLVE = 0.0D0
  DO BATCH_IDX = 1, NBATCHES
    RSTART_B = (BATCH_IDX-1)*NRHS_BATCH_EFF + 1
    REND_B = MIN(BATCH_IDX*NRHS_BATCH_EFF, OPT%NRHS)
    WIDTH_B = REND_B - RSTART_B + 1
    CALL RHS_LOCAL_RANGE_MBM(WIDTH_B, GRID%PC, GRID%MYCOL, LOCAL_OFFSET_B, NRHS_LOC_B)
    GSTART_B = RSTART_B + LOCAL_OFFSET_B - 1
    GEND_B = GSTART_B + NRHS_LOC_B - 1

    ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC_B,1)))
    B_LOCAL = (0.0_WP, 0.0_WP)
    IF (OPT%LOAD_RHS) THEN
      CALL LOAD_RHS_BATCH(LOAD_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
    ELSE
      CALL GEN_FILL_B(BSTART, BEND, NBLKR_LOC, GRID%PR, GRID%MYROW, NRHS_LOC_B, GSTART_B, B_LOCAL, MLOC)
      IF (OPT%SAVE_RHS) CALL SAVE_RHS_BATCH(SAVE_RHS_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)
    END IF
    IF (OPT%CHECK_RESIDUAL) THEN
      ! B_LOCAL est ecrase par X : garder le B reellement resolu.
      ALLOCATE(B_ORIG_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC_B,1)))
      B_ORIG_LOCAL = B_LOCAL
    END IF
    IF (NBATCHES.EQ.1) THEN
      CALL PROFILE_SAMPLE('second_membre')
    ELSE
      CALL PROFILE_SAMPLE(PCT_LABEL('second_membre', BATCH_IDX, NBATCHES))
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0 = MPI_WTIME()
    IF (OPT%LDLT) THEN
      CALL MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, A_SYM, WIDTH_B, NRHS_LOC_B, B_LOCAL, MAX(MLOC,1), &
                                     SHOW_PROGRESS=.TRUE., PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE, &
                                     DESCENT_MODE=OPT%DESCENT_MODE)
    ELSE
      CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_LU, WIDTH_B, NRHS_LOC_B, B_LOCAL, MAX(MLOC,1), &
                                            SHOW_PROGRESS=.TRUE., PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE, &
                                            DESCENT_MODE=OPT%DESCENT_MODE)
    END IF
    T_SOLVE = T_SOLVE + (MPI_WTIME() - T0)
    ! B_LOCAL contient maintenant X pour ce lot.
    IF (OPT%SAVE_SOLUTION) CALL SAVE_RHS_BATCH(SAVE_SOLUTION_UNIT, NRHS_LOC_B, B_LOCAL, MLOC)

    IF (OPT%CHECK_RESIDUAL) THEN
      ALLOCATE(XG_BATCH(OPT%N, WIDTH_B))
      CALL GATHER_2D_RHS(GRID, OPT%N, OPT%NB, WIDTH_B, B_LOCAL, MAX(MLOC,1), NRHS_LOC_B, 0, XG_BATCH)
      IF (BATCH_IDX.EQ.1 .AND. MYID0.EQ.0) PREVIEW_GLOBAL(1:NPREVIEW) = XG_BATCH(1:NPREVIEW,1)
      CALL ACCUMULATE_RESIDUAL(GRID, GEN_GET_A_BLOCK, OPT%N, OPT%NB, MLOC, XG_BATCH, B_ORIG_LOCAL, &
                               WIDTH_B, NRHS_LOC_B, NUM2_ACC, DEN2_ACC)
      DEALLOCATE(XG_BATCH, B_ORIG_LOCAL)
    ELSE IF (GSTART_B .LE. 1 .AND. 1 .LE. GEND_B) THEN
      ! Apercu de X(1:5,1) sans rassembler X : chaque rang fournit les
      ! valeurs qu'il possede, zero sinon, puis somme reduite.
      RS = 1
      DO li = 1, NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        MI = BEND(I) - BSTART(I) + 1
        IF (BSTART(I) .LE. NPREVIEW) &
             PREVIEW_LOCAL(BSTART(I):MIN(BEND(I),NPREVIEW)) = &
             B_LOCAL(RS:RS+MIN(BEND(I),NPREVIEW)-BSTART(I), 1-GSTART_B+1)
        RS = RS + MI
      END DO
    END IF
    DEALLOCATE(B_LOCAL)
  END DO

  IF (OPT%LOAD_RHS) THEN
    CLOSE(LOAD_RHS_UNIT)
  ELSE IF (OPT%SAVE_RHS) THEN
    CLOSE(SAVE_RHS_UNIT)
  END IF
  IF (OPT%SAVE_SOLUTION) CLOSE(SAVE_SOLUTION_UNIT)

  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Resolution : ', T_SOLVE, ' s'
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B (par lot)', &
       8_8*INT(MLOC,8)*INT(RHS_NUMROC_MBM(NRHS_BATCH_EFF,GRID%PC,GRID%MYCOL),8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres resolution)', GET_RSS_KB()*1024_8)
  CALL PROFILE_SAMPLE('resolution')

  !=============================================================================
  ! 4) apercu de la solution et residu
  !=============================================================================
  IF (.NOT. OPT%CHECK_RESIDUAL) &
       CALL MPI_REDUCE(PREVIEW_LOCAL, PREVIEW_GLOBAL, NPREVIEW, MPI_COMPLEX, MPI_SUM, 0, GRID%COMM, IERR)
  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Premieres valeurs de X(:,1) :'
    DO I = 1, NPREVIEW
      WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(PREVIEW_GLOBAL(I)), ' + i*', AIMAG(PREVIEW_GLOBAL(I))
    END DO
  END IF
  IF (OPT%CHECK_RESIDUAL) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', SQRT(NUM2_ACC)/SQRT(DEN2_ACC)
    CALL PROFILE_SAMPLE('verif_residu')
  END IF

  CALL PROFILE_SAMPLE('fin')
  CALL PROFILE_WRITE_CSV()

  ! --- 5) --reload-solution : libere A puis relit X depuis disque (rang 0) ---
  IF (OPT%LDLT) THEN
    CALL MPI_BLR_SYM_FREE(A_SYM)
  ELSE
    CALL MPI_BLR_MATRIX_FREE(A_LU)
  END IF
  IF (OPT%RELOAD_SOLUTION) THEN
    CALL REPORT_MEMORY_MB(GRID, 'RSS apres liberation de A', GET_RSS_KB()*1024_8)
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MYID0.EQ.0) THEN
      ALLOCATE(X_RELOADED(OPT%N, OPT%NRHS))
      CALL LOAD_SOLUTION(OPT%CHECKPOINT_DIR, OPT%N, OPT%NRHS, OPT%PR, OPT%PC, OPT%NB, &
                         1, OPT%N, 1, OPT%NRHS, X_RELOADED, OK)
      IF (.NOT. OK) THEN
        WRITE(6,'(A)') ' ERREUR : relecture de la solution (--reload-solution) impossible.'
      ELSE
        WRITE(6,'(A)') ' Solution relue depuis disque (rang 0 seul) :'
        DO I = 1, NPREVIEW
          WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_RELOADED(I,1)), ' + i*', AIMAG(X_RELOADED(I,1))
        END DO
      END IF
      DEALLOCATE(X_RELOADED)
    END IF
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  SUBROUTINE ABORT_RUN()
    INTEGER :: IE
    CALL MPI_FINALIZE(IE)
    STOP 1
  END SUBROUTINE ABORT_RUN

  SUBROUTINE PRINT_CONFIG()
    WRITE(6,'(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)') ' Solveur ', TRIM(SNAME), ' , grille ', OPT%PR, 'x', OPT%PC, &
         ' , N=', OPT%N, ' , NB=', OPT%NB, ' , NRHS=', OPT%NRHS, ' , TOL=', OPT%TOL
    IF (OPT%LDLT) THEN
      WRITE(6,'(A)') ' Matrice A : PCMATGEN symetrique, diagonale dominante ; B : PCMATGEN'
    ELSE
      WRITE(6,'(A)') ' Matrice A : PCMATGEN generale, diagonale dominante ; B : PCMATGEN'
    END IF
    WRITE(6,'(A,L1,A,L1,A,L1,A,L1,A,L1,A,I0)') ' use-rsvd=', OPT%USE_RSVD, ' ortho-recomp=', OPT%ORTHO_RECOMP, &
         ' diag-spread=', OPT%DIAG_SPREAD, ' lookahead=', OPT%LOOKAHEAD, ' panel-pack=', OPT%PANEL_PACK, &
         ' descent-mode=', OPT%DESCENT_MODE
    IF (OPT%NRHS_BATCH.GT.0 .AND. OPT%NRHS_BATCH.LT.OPT%NRHS .AND. .NOT.OPT%LOAD_RHS) &
         WRITE(6,'(A,I0,A)') ' Seconds membres traites par lots de ', OPT%NRHS_BATCH, ' colonnes'
    IF (OPT%CHECK_RESIDUAL) WRITE(6,'(A)') ' Verification : residu ||A*X-B||/||B|| en fin de resolution'
    IF (OPT%PROFILE_USAGE)  WRITE(6,'(A)') ' Profilage RSS/CPU : usage_profile_rank<N>.csv'
    IF (OPT%DUMP_STRUCTURE) WRITE(6,'(A)') ' Structure BLR : blr_structure_<etape>_rank<N>.csv'
    IF (OPT%DUMP_WORK)      WRITE(6,'(A)') ' Releve du travail de Schur : blr_work_rank<N>.bin'
    IF (OPT%LOAD_FACTORED) THEN
      WRITE(6,'(A)') ' Reprise : facteur relu -- assemblage et factorisation ignores'
    ELSE IF (OPT%LOAD_MATRIX) THEN
      WRITE(6,'(A)') ' Reprise : matrice relue -- assemblage ignore'
    END IF
    IF (OPT%LOAD_RHS) WRITE(6,'(A)') ' Reprise : second membre relu'
    IF (OPT%SAVE_MATRIX .OR. OPT%SAVE_RHS .OR. OPT%SAVE_FACTORED .OR. OPT%SAVE_SOLUTION) &
         WRITE(6,'(A,A)') ' Sauvegardes dans : ', TRIM(OPT%CHECKPOINT_DIR)
  END SUBROUTINE PRINT_CONFIG

  SUBROUTINE PRINT_STATS(LABEL)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    IF (MYID0.NE.0) RETURN
    WRITE(6,'(A,A,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' Blocs compresses (', LABEL, ') : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %), taux de compression reel ', &
         REAL_COMPRESSION_PCT(OPT%NB, N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
  END SUBROUTINE PRINT_STATS

  !-----------------------------------------------------------------------------
  SUBROUTINE SETUP_LU()
    IF (OPT%LOAD_FACTORED) THEN
      PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_factored'
      CALL LOAD_BLR_MATRIX(GRID, A_LU, OPT%N, OPT%NB, OPT%PR, OPT%PC, MYID0, TRIM(PFX), .TRUE., OK)
      IF (.NOT. OK) CALL ABORT_RUN()
      CALL RECOMPUTE_COMPRESSION_STATS(GRID, A_LU, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
      CALL PRINT_STATS('facteurs relus')
    ELSE
      IF (OPT%LOAD_MATRIX) THEN
        PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_matrix'
        CALL LOAD_BLR_MATRIX(GRID, A_LU, OPT%N, OPT%NB, OPT%PR, OPT%PC, MYID0, TRIM(PFX), .FALSE., OK)
        IF (.NOT. OK) CALL ABORT_RUN()
        CALL RECOMPUTE_COMPRESSION_STATS(GRID, A_LU, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
        CALL PRINT_STATS('matrice relue')
      ELSE
        CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, OPT%N, OPT%NB, OPT%TOL, GEN_GET_A_BLOCK, A_LU, T_ASM, &
                                     N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., &
                                     USE_RSVD=OPT%USE_RSVD, PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY, &
                                     DIAG_SPREAD=OPT%DIAG_SPREAD)
        IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Assemblage + compression : ', T_ASM, ' s'
        CALL PRINT_STATS('assemblage')
        IF (OPT%SAVE_MATRIX) THEN
          PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_matrix'
          CALL SAVE_BLR_MATRIX(GRID, A_LU, OPT%N, OPT%NB, OPT%PR, OPT%PC, MYID0, TRIM(PFX), .FALSE., OK)
          IF (.NOT. OK) CALL ABORT_RUN()
        END IF
      END IF
      CALL REPORT_MEMORY_MB(GRID, 'Matrice A (compressee)', MPI_BLR_MATRIX_MEMORY_BYTES(A_LU))
      CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
      CALL PROFILE_SAMPLE('assemblage')
      CALL DUMP_BLR_STRUCTURE(GRID, A_LU, 'assembly')

      CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_LU, OPT%TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                          N_SWAPS, T_FACTOR, SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE., &
                                          PROGRESS_CALLBACK=PROFILE_CALLBACK_FACTOR, LOOKAHEAD=OPT%LOOKAHEAD)
      IF (INFO .NE. 0) THEN
        IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR : factorisation LU, INFO=', INFO
        CALL ABORT_RUN()
      END IF
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
        WRITE(6,'(A,I0)') ' Permutations de lignes (pivot) : ', N_SWAPS
      END IF
      CALL PRINT_STATS('facteurs')
      IF (OPT%SAVE_FACTORED) THEN
        PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_factored'
        CALL SAVE_BLR_MATRIX(GRID, A_LU, OPT%N, OPT%NB, OPT%PR, OPT%PC, MYID0, TRIM(PFX), .TRUE., OK)
        IF (.NOT. OK) CALL ABORT_RUN()
      END IF
    END IF
    CALL REPORT_MEMORY_MB(GRID, 'Facteurs L et U (compresses)', MPI_BLR_MATRIX_MEMORY_BYTES(A_LU))
    CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
    CALL PROFILE_SAMPLE('factorisation')
    CALL DUMP_BLR_STRUCTURE(GRID, A_LU, 'factorization')
  END SUBROUTINE SETUP_LU

  !-----------------------------------------------------------------------------
  SUBROUTINE SETUP_LDLT()
    INTEGER :: KK, KW, NS_K, N2_K
    IF (OPT%LOAD_FACTORED) THEN
      PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_ldlt_factored'
      CALL MPI_BLR_SYM_LOAD(GRID, A_SYM, OPT%N, OPT%NB, TRIM(PFX), .TRUE., OK)
      IF (.NOT. OK) CALL ABORT_RUN()
      CALL MPI_BLR_SYM_COMPRESSION_STATS(GRID, A_SYM, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
      N_SWAPS = 0; N_2X2 = 0
      DO KK = 1, A_SYM%NBLK
        KW = A_SYM%BEND(KK) - A_SYM%BSTART(KK) + 1
        CALL COUNT_PIVOTS(KW, A_SYM%IPIV(1:KW,KK), A_SYM%DSUB(1:KW,KK), NS_K, N2_K)
        N_SWAPS = N_SWAPS + NS_K; N_2X2 = N_2X2 + N2_K
      END DO
      IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0)') ' Echanges de pivot (relus) : ', N_SWAPS, ' , pivots 2x2 : ', N_2X2
      CALL PRINT_STATS('facteur relu')
    ELSE
      IF (OPT%LOAD_MATRIX) THEN
        PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_ldlt_matrix'
        CALL MPI_BLR_SYM_LOAD(GRID, A_SYM, OPT%N, OPT%NB, TRIM(PFX), .FALSE., OK)
        IF (.NOT. OK) CALL ABORT_RUN()
        CALL MPI_BLR_SYM_COMPRESSION_STATS(GRID, A_SYM, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
        CALL PRINT_STATS('matrice relue')
      ELSE
        CALL MPI_ASSEMBLE_BLR_SYM(GRID, OPT%N, OPT%NB, OPT%TOL, GEN_GET_A_BLOCK, A_SYM, T_ASM, &
                                  N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., &
                                  USE_RSVD=OPT%USE_RSVD, PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY, &
                                  DIAG_SPREAD=OPT%DIAG_SPREAD)
        IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Assemblage + compression (triangle inferieur) : ', T_ASM, ' s'
        CALL PRINT_STATS('assemblage')
        IF (OPT%SAVE_MATRIX) THEN
          PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_ldlt_matrix'
          CALL MPI_BLR_SYM_SAVE(GRID, A_SYM, TRIM(PFX), .FALSE., OK)
          IF (.NOT. OK) CALL ABORT_RUN()
        END IF
      END IF
      CALL REPORT_MEMORY_MB(GRID, 'Matrice A (triangle inferieur, compressee)', MPI_BLR_SYM_MEMORY_BYTES(A_SYM))
      CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)
      CALL PROFILE_SAMPLE('assemblage')
      CALL DUMP_BLR_STRUCTURE(GRID, A_SYM, 'assembly')

      CALL MPI_BLR_SYM_LDLT_FACTOR(GRID, A_SYM, OPT%TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                   N_SWAPS, N_2X2, T_FACTOR, SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE., &
                                   PROGRESS_CALLBACK=PROFILE_CALLBACK_FACTOR, LOOKAHEAD=OPT%LOOKAHEAD)
      IF (INFO .NE. 0) THEN
        IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR : factorisation LDLt (D singuliere), INFO=', INFO
        CALL ABORT_RUN()
      END IF
      IF (MYID0.EQ.0) THEN
        WRITE(6,'(A,F8.3,A)') ' Factorisation LDLt : ', T_FACTOR, ' s'
        WRITE(6,'(A,I0,A,I0)') ' Echanges de pivot : ', N_SWAPS, ' , pivots 2x2 : ', N_2X2
      END IF
      CALL PRINT_STATS('facteur L')
      IF (OPT%SAVE_FACTORED) THEN
        PFX = TRIM(OPT%CHECKPOINT_DIR)//'/saved_ldlt_factored'
        CALL MPI_BLR_SYM_SAVE(GRID, A_SYM, TRIM(PFX), .TRUE., OK)
        IF (.NOT. OK) CALL ABORT_RUN()
      END IF
    END IF
    CALL REPORT_MEMORY_MB(GRID, 'Facteurs L et D (compresses)', MPI_BLR_SYM_MEMORY_BYTES(A_SYM))
    CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)
    CALL PROFILE_SAMPLE('factorisation')
    CALL DUMP_BLR_STRUCTURE(GRID, A_SYM, 'factorization')
  END SUBROUTINE SETUP_LDLT

END PROGRAM SOLVEUR_BLR
