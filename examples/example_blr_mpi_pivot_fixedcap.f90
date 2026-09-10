!===============================================================================
! VARIANTE FIXEDCAP de example_blr_mpi_pivot.f90 -- MEME algorithme (pivot
! partiel confine au bloc diagonal, MEME matrice de demonstration), MEME
! enchainement (assembler, factoriser, resoudre, rassembler), pour
! permettre une comparaison A/B DIRECTE avec l'executable existant
! (example_blr_mpi_pivot, JAMAIS modifie) ET avec la variante arena
! (example_blr_mpi_pivot_arena.f90) -- voir mpi_blr_matrix_mod.f90, type
! MPI_BLR_MATRIX_FIXEDCAP_T.
!
! PRINCIPE : les blocs bas-rang
! stockent le rang MAXIMAL que COMPRESS_WORTHWHILE_RATIO autorise (RMAX),
! alloue UNE SEULE FOIS -- JAMAIS realloue ensuite, quelle que soit la
! fluctuation du rang REEL au fil des recompressions (voir
! BLR_BLOCK_FIXEDCAP_T, blr_compress_mod.f90). Contrairement a l'arena
! (allocateur dedie, buffers "carve-ou-reutilise"), cette variante utilise
! l'ALLOCATE/DEALLOCATE Fortran ORDINAIRE : puisque la capacite est
! TOUJOURS la MEME pour tous les blocs d'un meme calibre (M fixe par la
! partition), glibc peut reutiliser efficacement les chunks liberes
! (repli dense) via ses free-lists internes SANS custom allocator.
!
! Les echanges MPI (diffusion des panneaux L/U) et les operations
! associees (TRSM, GEMM, APPLY) restent bornes au rang REEL a chaque
! instant, JAMAIS a RMAX -- voir COPY_FIXEDCAP_TO_REGULAR et les vues
! tranchees BLK%VFAC(1:RANK,:) dans blr_compress_mod.f90.
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot_fixedcap PR PC [N] [NB] [NRHS] [TOL]
!         (defauts : N=4000 NB=200 NRHS=4 TOL=1e-4)
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_PIVOT_FIXEDCAP
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_FIXEDCAP_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./example_blr_mpi_pivot_fixedcap PR PC [N] [NB] [NRHS] [TOL]'
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

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)') &
       ' [FIXEDCAP] Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL

  ! RSS de reference -- voir example_blr_mpi_pivot.f90 pour la justification.
  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  ! --- 1) assemblage BLR distribue, fixedcap -----------------------------
  CALL MPI_ASSEMBLE_BLR_MATRIX_FIXEDCAP(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                         N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE.)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Taux de compression (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF
  ! Empreinte THEORIQUE = capacite REELLEMENT allouee (M*RMAX+RMAX*N, PAS
  ! M*RANK+RANK*N) -- voir MPI_BLR_MATRIX_FIXEDCAP_MEMORY_BYTES.
  CALL REPORT_MEMORY_MB(GRID, 'Matrice A_BLR fixedcap (compressee)', MPI_BLR_MATRIX_FIXEDCAP_MEMORY_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, '  dont donnee utile (rang reel)', MPI_BLR_MATRIX_FIXEDCAP_USED_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)

  ! --- 2) factorisation LU, AVEC pivot, fixedcap --------------------------
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_FIXEDCAP(GRID, A_BLR, TOL, INFO, &
                                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                                SHOW_PROGRESS=.FALSE., SHOW_RSS_TRACE=.TRUE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_FIXEDCAP, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
    WRITE(6,'(A,I0)') ' Nombre de permutations de lignes (pivot) : ', N_SWAPS
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Taux de compression (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU fixedcap (compresses)', MPI_BLR_MATRIX_FIXEDCAP_MEMORY_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, '  dont donnee utile (rang reel)', MPI_BLR_MATRIX_FIXEDCAP_USED_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)

  ! --- 3) second membre B ------------------------------------------------
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_B(A_BLR, GRID, NRHS_LOC, RSTART_MINE, B_LOCAL, MLOC)
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B', 8_8*INT(MLOC,8)*INT(NRHS_LOC,8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres allocation B)', GET_RSS_KB()*1024_8)

  ! --- 4) resolution 2D, fixedcap ------------------------------------------
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_FIXEDCAP_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1), &
                                                  SHOW_PROGRESS=.TRUE.)

  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
    DO I = 1, MIN(5, N)
      WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
    END DO
  END IF

  DEALLOCATE(B_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Matrice IDENTIQUE a example_blr_mpi_pivot.f90/example_blr_mpi_pivot_arena.f90
  ! (comparaison A/B directe).
  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    REAL(WP) :: DF
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        IF (JJ .EQ. N+1-II) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
        ELSE
          DF = 1.0_WP / (1.0_WP + ABS(REAL(N+1-II-JJ,WP)))
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(DF, 0.3_WP*DF, KIND=WP)
        END IF
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  SUBROUTINE FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_FIXEDCAP_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)                    :: GD
    INTEGER, INTENT(IN)                         :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)                    :: BL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, GRR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          GRR = RSTART + RR - 1
          BL(RS + (GII-AB%BSTART(II)), RR) = CMPLX(REAL(GII,WP)*0.01_WP, -0.3_WP*REAL(GRR,WP), KIND=WP)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_B

END PROGRAM EXAMPLE_BLR_MPI_PIVOT_FIXEDCAP
