!===============================================================================
! BENCHMARK -- verifie que la RSS de la variante OUT-OF-CORE reste BORNEE
! independamment de NRHS, meme pour un second membre MASSIF (des centaines
! a des milliers de colonnes) -- c'est l'un des buts explicites de la
! conception OOC (voir mpi_blr_matrix_mod.f90, MPI_ASSEMBLE_RHS_OOC/
! MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D).
!
! DEUX CHOIX DE CONCEPTION, distincts de bench_blr_mpi_pivot_ooc_compare.f90 :
!   1) FILL_RHS est une FORMULE DIRECTE, O(1) par entree (PAS le produit
!      brut A*Xexact utilise ailleurs pour la validation rigoureuse sous
!      batching arbitraire -- ce dernier coute O(N^2*NRHS), imprati-
!      cable ici : a NRHS=1000, N=20000, ce serait ~800 MILLIARDS
!      d'operations). Consequence : PAS de verification ||X-Xexact||
!      independante ici (voir test_mpi_blr_matrix_pivot_ooc.f90 pour la
!      correction, deja etablie a petite echelle) -- ce banc mesure
!      UNIQUEMENT la RSS/le temps, pas la precision numerique.
!   2) CONSUME_X_BATCH n'ACCUMULE JAMAIS les paquets resolus dans un
!      tableau X_GLOBAL(N,NRHS) (ce que fait l'exemple canonique, a des
!      fins de demonstration/verification a echelle modeste, PAS
!      representatif d'un vrai cas massivement second membre) -- il
!      calcule juste une somme locale ||X_LOCAL||^2, cumulee dans une
!      variable SCALAIRE partagee par association d'hote. Memoire
!      supplementaire de ce chemin : O(1), independante de NRHS.
!
! Lance le MEME systeme (meme N/NB/MAX_MB_PER_PROCESS/NRHS_BATCH) pour
! PLUSIEURS valeurs de NRHS successives (meme processus, memes fichiers
! scratch a chaque fois) et affiche la RSS a chaque etape -- la RSS de la
! resolution ne doit PAS croitre avec NRHS (bornee par NRHS_BATCH), la
! RSS de l'assemblage/factorisation NE DEPEND PAS de NRHS du tout.
!
! Usage : mpirun -np <PR*PC> ./bench_blr_mpi_pivot_ooc_massive_rhs PR PC N NB \
!         TOL MAX_MB_PER_PROCESS NRHS_BATCH [SCRATCH_DIR] NRHS1 [NRHS2 ...]
!===============================================================================
PROGRAM BENCH_BLR_MPI_PIVOT_OOC_MASSIVE_RHS
  USE MPI_BLR_MATRIX_MOD
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_BLR
  TYPE(OOC_RHS_T) :: RHS
  INTEGER :: INFO, N_SWAPS, N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_ASM_B, T_FACTOR, T0, T1
  INTEGER :: NARGS, IARG, NRHS, NRUNS, IRUN
  INTEGER, ALLOCATABLE :: NRHS_LIST(:)
  DOUBLE PRECISION :: CHECKSUM_LOCAL, CHECKSUM_GLOBAL

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 8) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./bench_blr_mpi_pivot_ooc_massive_rhs PR PC N NB TOL '// &
         'MAX_MB_PER_PROCESS NRHS_BATCH [SCRATCH_DIR] NRHS1 [NRHS2 ...]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) TOL
  CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) MAX_MB_PER_PROCESS
  CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) NRHS_BATCH

  NARGS = COMMAND_ARGUMENT_COUNT()
  CALL GET_COMMAND_ARGUMENT(8, ARG)
  IF (VERIFY(TRIM(ARG), '0123456789') .EQ. 0) THEN
    ! 8e argument purement numerique -- pas de SCRATCH_DIR fourni, c'est le
    ! premier NRHS de la liste.
    SCRATCH_DIR = './blr_ooc_scratch_massive_rhs'
    IARG = 8
  ELSE
    SCRATCH_DIR = TRIM(ARG)
    IARG = 9
  END IF
  NRUNS = NARGS - IARG + 1
  IF (NRUNS .LT. 1) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : au moins une valeur de NRHS requise apres SCRATCH_DIR.'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  ALLOCATE(NRHS_LIST(NRUNS))
  DO IRUN = 1, NRUNS
    CALL GET_COMMAND_ARGUMENT(IARG+IRUN-1, ARG); READ(ARG,*) NRHS_LIST(IRUN)
  END DO

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,F8.1,A,I0)') &
       ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, &
       ' , MAX_MB_PER_PROCESS=', MAX_MB_PER_PROCESS, ' , NRHS_BATCH=', NRHS_BATCH

  ! === 1+3) assemblage de A et factorisation -- UNE SEULE FOIS, reutilise
  ! pour toutes les valeurs de NRHS (independant du second membre) ===
  CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                    T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM)
  CALL REPORT_MEMORY_MB(GRID, 'RSS apres assemblage A', GET_RSS_KB()*1024_8)
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Temps assemblage A : ', T_ASM, ' s'

  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                           N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'RSS apres factorisation', GET_RSS_KB()*1024_8)
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Temps factorisation : ', T_FACTOR, ' s'

  ! === 2+4) pour CHAQUE NRHS demande : assemblage de B PUIS resolution,
  ! avec la MEME matrice/facteurs deja en place -- montre que la RSS de
  ! ces deux etapes ne depend QUE de NRHS_BATCH, jamais de NRHS lui-meme.
  DO IRUN = 1, NRUNS
    NRHS = NRHS_LIST(IRUN)
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' '
    IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A)') ' --- NRHS=', NRHS, ' (', &
         (NRHS+NRHS_BATCH-1)/NRHS_BATCH, ' paquets) ---'

    T0 = MPI_WTIME()
    CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
    T1 = MPI_WTIME()
    T_ASM_B = T1 - T0
    CALL REPORT_MEMORY_MB(GRID, 'RSS apres assemblage B', GET_RSS_KB()*1024_8)
    IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Temps assemblage B : ', T_ASM_B, ' s'

    CHECKSUM_LOCAL = 0.0D0
    T0 = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH)
    T1 = MPI_WTIME()
    CALL REPORT_MEMORY_MB(GRID, 'RSS apres resolution   ', GET_RSS_KB()*1024_8)
    IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Temps resolution : ', T1-T0, ' s'

    CALL MPI_ALLREDUCE(CHECKSUM_LOCAL, CHECKSUM_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
    IF (MYID0.EQ.0) WRITE(6,'(A,ES14.6)') ' Checksum ||X||^2 (global, non accumule) = ', CHECKSUM_GLOBAL

    CALL OOC_RHS_FREE(RHS)
  END DO

  CALL MPI_BLR_MATRIX_OOC_FREE(A_BLR)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Matrice de demonstration, IDENTIQUE a example_blr_mpi_pivot_ooc.f90.
  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ, NN)
    INTEGER, INTENT(IN) :: GI, GJ, NN
    REAL(WP) :: DF
    IF (GJ .EQ. NN+1-GI) THEN
      LOWRANK_A = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
    ELSE
      DF = 1.0_WP / (1.0_WP + ABS(REAL(NN+1-GI-GJ,WP)))
      LOWRANK_A = CMPLX(DF, 0.3_WP*DF, KIND=WP)
    END IF
  END FUNCTION LOWRANK_A

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        BLK(II-ISTART+1, JJ-JSTART+1) = LOWRANK_A(II, JJ, N)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  ! Formule DIRECTE, O(1) par entree -- PAS de produit brut A*Xexact (voir
  ! commentaire de tete : impraticable a NRHS massif). Aucune verification
  ! de precision n'est faite ici, uniquement RSS/temps.
  SUBROUTINE FILL_RHS(ISTART, IEND, RSTART, REND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RSTART, REND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, REND-RSTART+1)
    INTEGER :: II, RR
    DO RR = RSTART, REND
      DO II = ISTART, IEND
        BLK(II-ISTART+1, RR-RSTART+1) = CMPLX(REAL(II,WP)*0.01_WP, -0.3_WP*REAL(RR,WP), KIND=WP)
      END DO
    END DO
  END SUBROUTINE FILL_RHS

  ! N'ACCUMULE JAMAIS le paquet resolu -- juste une somme locale ajoutee a
  ! CHECKSUM_LOCAL (partagee par association d'hote), memoire O(1)
  ! independante de NRHS. C'est le comportement d'un VRAI cas massivement
  ! second membre (voir commentaire de tete) -- l'exemple canonique
  ! (example_blr_mpi_pivot_ooc.f90) accumule au contraire tout dans
  ! X_GLOBAL(N,NRHS) a des fins de demonstration/verification, ce qui
  ! grandit AVEC NRHS et ne conviendrait pas ici.
  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    IF (NRHS_LOC_B .GT. 0) CHECKSUM_LOCAL = CHECKSUM_LOCAL + SUM(ABS(X_LOCAL(:,1:NRHS_LOC_B))**2)
  END SUBROUTINE CONSUME_X_BATCH

END PROGRAM BENCH_BLR_MPI_PIVOT_OOC_MASSIVE_RHS
