!===============================================================================
! Suivi RSS/CPU par etape (PROFILE_USAGE) et export de la structure BLR
! (DUMP_STRUCTURE).
!
! ETAT DE MODULE (contrairement a checkpoint_io_mod.f90) : PROFILE_
! CALLBACK_ASSEMBLY/FACTOR/SOLVE sont passees a la bibliotheque comme
! PROCEDURE(STAGE_PROGRESS_CALLBACK_IFACE) -- une interface FIXE (NDONE,
! NTOTAL uniquement, imposee par mpi_blr_matrix_mod.f90), qui ne laisse
! aucune place pour un argument supplementaire. Le seul moyen de leur
! donner acces a PROFILE_USAGE_FLAG/aux tableaux d'echantillons est donc
! un etat de module (SAVE), initialise UNE FOIS par PROFILING_INIT avant
! l'assemblage.
!===============================================================================
MODULE PROFILING_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: MPI_BLR_MATRIX_T, GRID_T, WP, BLR_COL_GLOBAL
  USE MPI_BLR_LDLT_MOD, ONLY: MPI_BLR_SYM_T
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  USE MPI
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: PROFILING_INIT, PROFILE_SAMPLE, PCT_LABEL
  PUBLIC :: PROFILE_CALLBACK_ASSEMBLY, PROFILE_CALLBACK_FACTOR, PROFILE_CALLBACK_SOLVE
  PUBLIC :: PROFILE_WRITE_CSV, DUMP_BLR_STRUCTURE

  ! Meme nom pour les deux solveurs (LU : MPI_BLR_MATRIX_T, LDLt : MPI_BLR_SYM_T).
  INTERFACE DUMP_BLR_STRUCTURE
    MODULE PROCEDURE DUMP_BLR_STRUCTURE_LU, DUMP_BLR_STRUCTURE_SYM
  END INTERFACE DUMP_BLR_STRUCTURE

  INTEGER, PARAMETER :: MAX_PROFILE_SAMPLES = 128
  INTEGER, SAVE :: MYID0_SAVED = 0
  LOGICAL, SAVE :: PROFILE_USAGE_FLAG_SAVED = .FALSE.
  LOGICAL, SAVE :: DUMP_STRUCTURE_FLAG_SAVED = .FALSE.
  CHARACTER(LEN=40), SAVE :: PROF_STAGE(MAX_PROFILE_SAMPLES)
  DOUBLE PRECISION, SAVE :: PROF_T(MAX_PROFILE_SAMPLES), PROF_RSS_MB(MAX_PROFILE_SAMPLES), PROF_CPU_PCT(MAX_PROFILE_SAMPLES)
  INTEGER, SAVE :: PROF_N = 0
  DOUBLE PRECISION, SAVE :: PROF_T_START = 0.0D0, PROF_T_PREV = 0.0D0
  REAL, SAVE :: PROF_CPU_PREV = 0.0

CONTAINS

  ! A appeler UNE FOIS, avant l'assemblage -- fixe MYID0/les indicateurs
  ! PROFILE_USAGE/DUMP_STRUCTURE pour tout le reste du run, demarre
  ! l'horloge de reference et enregistre le premier echantillon ('debut').
  SUBROUTINE PROFILING_INIT(MYID0, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG)
    INTEGER, INTENT(IN) :: MYID0
    LOGICAL, INTENT(IN) :: PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG
    MYID0_SAVED = MYID0
    PROFILE_USAGE_FLAG_SAVED = PROFILE_USAGE_FLAG
    DUMP_STRUCTURE_FLAG_SAVED = DUMP_STRUCTURE_FLAG
    PROF_N = 0
    PROF_T_START = MPI_WTIME()
    PROF_T_PREV = 0.0D0
    CALL CPU_TIME(PROF_CPU_PREV)
    CALL PROFILE_SAMPLE('debut')
  END SUBROUTINE PROFILING_INIT

  ! Echantillon RSS/CPU pour l'etape STAGE_NAME -- no-op si desactive.
  ! CPU_TIME mesure le temps processeur (tous threads MKL/OpenMP
  ! confondus), MPI_WTIME le temps horloge ecoule -- leur rapport sur
  ! l'intervalle donne le nombre moyen de coeurs actifs, en %, comme
  ! `top`. Chaque processus echantillonne son propre etat, exploitable
  ! point par point dans plot_usage_profile.py.
  SUBROUTINE PROFILE_SAMPLE(STAGE_NAME)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_NAME
    REAL :: CPU_NOW
    DOUBLE PRECISION :: T_NOW, DT_WALL, DT_CPU
    IF (.NOT. PROFILE_USAGE_FLAG_SAVED) RETURN
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
  ! la factorisation/la resolution.
  SUBROUTINE PROFILE_CALLBACK_ASSEMBLY(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    IF (.NOT. PROFILE_USAGE_FLAG_SAVED) RETURN
    CALL PROFILE_SAMPLE(PCT_LABEL('assemblage', NDONE, NTOTAL))
  END SUBROUTINE PROFILE_CALLBACK_ASSEMBLY

  SUBROUTINE PROFILE_CALLBACK_FACTOR(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    IF (.NOT. PROFILE_USAGE_FLAG_SAVED) RETURN
    CALL PROFILE_SAMPLE(PCT_LABEL('factorisation', NDONE, NTOTAL))
  END SUBROUTINE PROFILE_CALLBACK_FACTOR

  SUBROUTINE PROFILE_CALLBACK_SOLVE(NDONE, NTOTAL)
    INTEGER, INTENT(IN) :: NDONE, NTOTAL
    IF (.NOT. PROFILE_USAGE_FLAG_SAVED) RETURN
    CALL PROFILE_SAMPLE(PCT_LABEL('resolution', NDONE, NTOTAL))
  END SUBROUTINE PROFILE_CALLBACK_SOLVE

  ! Ecrit les echantillons de ce processus dans son propre CSV -- un
  ! fichier par rang, pas de rassemblement MPI.
  SUBROUTINE PROFILE_WRITE_CSV()
    INTEGER :: U, K
    CHARACTER(LEN=64) :: FNAME
    IF (.NOT. PROFILE_USAGE_FLAG_SAVED) RETURN
    WRITE(FNAME,'(A,I0,A)') 'usage_profile_rank', MYID0_SAVED, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'rank,stage,t_wall_s,rss_mb,cpu_pct'
    DO K = 1, PROF_N
      WRITE(U,'(I0,A,A,A,F12.4,A,F12.3,A,F8.2)') MYID0_SAVED, ',', TRIM(PROF_STAGE(K)), ',', &
           PROF_T(K), ',', PROF_RSS_MB(K), ',', PROF_CPU_PCT(K)
    END DO
    CLOSE(U)
  END SUBROUTINE PROFILE_WRITE_CSV

  ! Ecrit la nature de chaque bloc local de A_BLR -- un fichier par
  ! processus, jamais rassemble ici : plot_blr_structure.py reconstitue la
  ! grille globale en lisant tous les fichiers d'un meme STAGE_LABEL.
  SUBROUTINE DUMP_BLR_STRUCTURE_LU(GRID, A_BLR, STAGE_LABEL)
    TYPE(GRID_T), INTENT(IN) :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_LABEL
    INTEGER :: U, li, lj, I, J
    CHARACTER(LEN=96) :: FNAME
    IF (.NOT. DUMP_STRUCTURE_FLAG_SAVED) RETURN
    WRITE(FNAME,'(A,A,A,I0,A)') 'blr_structure_', TRIM(STAGE_LABEL), '_rank', MYID0_SAVED, '.csv'
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
  END SUBROUTINE DUMP_BLR_STRUCTURE_LU

  ! Variante LDLt : meme format CSV, seuls les blocs du triangle inferieur
  ! existent (blocs diagonaux stockes compacts dans DIAGP).
  SUBROUTINE DUMP_BLR_STRUCTURE_SYM(GRID, A, STAGE_LABEL)
    TYPE(GRID_T), INTENT(IN) :: GRID
    TYPE(MPI_BLR_SYM_T), INTENT(IN) :: A
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_LABEL
    INTEGER :: U, li, lj, I, J
    CHARACTER(LEN=96) :: FNAME
    IF (.NOT. DUMP_STRUCTURE_FLAG_SAVED) RETURN
    WRITE(FNAME,'(A,A,A,I0,A)') 'blr_structure_', TRIM(STAGE_LABEL), '_rank', MYID0_SAVED, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'block_row,block_col,kind,rank'
    DO I = 1, A%NBLK
      IF (ALLOCATED(A%DIAGP(I)%AP)) WRITE(U,'(I0,A,I0,A,I0)') I, ',', I, ',diagonal,', A%BEND(I)-A%BSTART(I)+1
    END DO
    DO li = 1, A%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A%NBLKC_LOC
        J = BLR_COL_GLOBAL(A%CMAP, lj, GRID)
        IF (I .LE. J) CYCLE
        IF (A%OFFBLK(li,lj)%IS_LOWRANK) THEN
          WRITE(U,'(I0,A,I0,A,I0)') I, ',', J, ',lowrank,', A%OFFBLK(li,lj)%RANK
        ELSE
          WRITE(U,'(I0,A,I0,A,I0)') I, ',', J, ',dense,', A%BEND(J)-A%BSTART(J)+1
        END IF
      END DO
    END DO
    CLOSE(U)
  END SUBROUTINE DUMP_BLR_STRUCTURE_SYM

END MODULE PROFILING_MOD
