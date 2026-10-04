!===============================================================================
! BENCHMARK COMPARATIF -- solveur BLR+MPI+pivot, VERSION CLASSIQUE (tout en
! RAM, mpi_blr_matrix_mod.f90 / MPI_BLR_MATRIX_T) vs VERSION OUT-OF-CORE
! (MPI_BLR_MATRIX_OOC_T, voir blr_ooc_mod.f90) -- MEME matrice, MEME
! Xexact, MEME grille, pour permettre une comparaison directe temps CPU +
! RSS a CHAQUE etape (assemblage, factorisation, resolution).
!
! MODE (10e argument, OBLIGATOIRE, 'classical' ou 'ooc') choisit UN SEUL
! chemin par execution -- PAS les deux dans le meme processus : la RSS
! (/proc/self/status) ne redescend PAS apres un DEALLOCATE, l'allocateur
! (glibc/MKL) garde la memoire liberee pour la reutiliser plutot que de la
! rendre a l'OS (meme constat deja documente dans la memoire du projet pour
! l'ecart RSS-vs-theorique) -- enchainer les deux chemins dans UN processus
! contaminerait la RSS du second par les allocations du premier, JAMAIS
! liberees au sens du systeme. Lancer donc DEUX mpirun separes (memes
! parametres N/NB/NRHS/..., MODE different) et comparer leurs sorties.
!
! MATRICE : IDENTIQUE dans les deux chemins -- anti-diagonale forte CONFINEE
! au bloc diagonal (force un vrai pivotage QUE LE PIVOT CONFINE PEUT VOIR,
! voir test_mpi_blr_matrix_pivot.f90/_ooc.f90) + noyau de Cauchy amorti +
! bruit plein-rang EPS*HASH_NOISE (meme recette EXACTE que bench_blr_mpi_
! pivot_compressibility.f90 -- EPS controle la compressibilite : EPS=0 =
! noyau de Cauchy pur (tres compressible) ; EPS croissant degrade
! progressivement le rang effectif des blocs hors-diagonale, jusqu'a un
! repli dense pour beaucoup d'entre eux si TOL fixe. Xexact connu par
! formule deterministe (XEX_FORMULA) dans les deux cas.
!
! Cote OOC, le second membre B est calcule par PRODUIT BRUT (voir
! test_mpi_blr_matrix_pivot_ooc.f90 pour la justification : la distribution
! colonnes<->processus d'un paquet NRHS_BATCH peut differer de celle du
! NRHS complet utilise cote classique -- un calcul base uniquement sur les
! indices globaux demandes reste correct quelle que soit cette distribution).
! Cote classique, B = A_BLR*Xexact via APPLY_MPI_BLR_MATRIX (plus rapide,
! deja valide par test_mpi_blr_matrix_pivot.f90) -- les deux methodes
! calculent la MEME quantite mathematique, seule la voie de calcul differe.
!
! Usage : mpirun -np <PR*PC> ./bench_blr_mpi_pivot_ooc_compare PR PC MODE [N] [NB] \
!         [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [EPS]
!         MODE = 'classical' ou 'ooc' (voir note ci-dessus -- obligatoire).
!         (defauts : N=8000 NB=400 NRHS=6 TOL=1e-4 MAX_MB_PER_PROCESS=16
!          NRHS_BATCH=2 SCRATCH_DIR=./blr_ooc_scratch_bench EPS=0.0)
!===============================================================================
PROGRAM BENCH_BLR_MPI_PIVOT_OOC_COMPARE
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS, EPS
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  CHARACTER(LEN=16) :: MODE

  ! --- chemin CLASSIQUE (tout en RAM) ---
  TYPE(MPI_BLR_MATRIX_T) :: A_CL
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), XEX_LOCAL(:,:), X_GLOBAL_CL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO_CL, N_SWAPS_CL
  INTEGER :: N_COMPRESSED_CL, N_OFFDIAG_CL, li, I
  INTEGER(8) :: RANK_SUM_CL
  DOUBLE PRECISION :: T_ASM_CL, T_FACTOR_CL, T_SOLVE_CL
  DOUBLE PRECISION :: RSS_ASM_CL, RSS_FACT_CL, RSS_SOLVE_CL
  REAL(WP) :: ERR_CL

  ! --- chemin OUT-OF-CORE ---
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_OOC
  TYPE(OOC_RHS_T) :: RHS
  COMPLEX(WP), ALLOCATABLE :: X_GLOBAL_OOC(:,:)
  INTEGER :: INFO_OOC, N_SWAPS_OOC, N_COMPRESSED_OOC, N_OFFDIAG_OOC
  INTEGER(8) :: RANK_SUM_OOC
  DOUBLE PRECISION :: T_ASM_A_OOC, T_ASM_B_OOC, T_FACTOR_OOC, T_SOLVE_OOC
  DOUBLE PRECISION :: RSS_ASMA_OOC, RSS_ASMB_OOC, RSS_FACT_OOC, RSS_SOLVE_OOC
  DOUBLE PRECISION :: T_CACHE_GET, T_CACHE_EVICT, T_COMPUTE, T_MPICALL
  DOUBLE PRECISION :: T_GETBLK_ASM, T_COMPRESS_ASM, T_IO_ASM
  REAL(WP) :: ERR_OOC
  DOUBLE PRECISION :: T0B, T1B

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  N = 8000; NB = 400; NRHS = 6; TOL = 1.0E-4_WP
  MAX_MB_PER_PROCESS = 16.0_WP; NRHS_BATCH = 2; SCRATCH_DIR = './blr_ooc_scratch_bench'
  EPS = 0.0_WP

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./bench_blr_mpi_pivot_ooc_compare PR PC MODE [N] [NB] [NRHS] [TOL] '// &
         '[MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR]   (MODE = classical|ooc, obligatoire -- voir '// &
         'commentaire de tete : ne PAS lancer les deux chemins dans un seul processus)'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(3, MODE)
  IF (MODE .NE. 'classical' .AND. MODE .NE. 'ooc') THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : MODE doit etre "classical" ou "ooc", recu : ', TRIM(MODE)
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN
    CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) MAX_MB_PER_PROCESS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) THEN
    CALL GET_COMMAND_ARGUMENT(9, ARG); READ(ARG,*) NRHS_BATCH
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN
    CALL GET_COMMAND_ARGUMENT(10, SCRATCH_DIR)
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN
    CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) EPS
  END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1,A,I0)') &
       ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
       ' , MAX_MB_PER_PROCESS(OOC)=', MAX_MB_PER_PROCESS, ' , NRHS_BATCH(OOC)=', NRHS_BATCH

  IF (MODE .EQ. 'classical') THEN
    !===========================================================================
    ! CHEMIN CLASSIQUE (tout en RAM)
    !===========================================================================
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_CL, T_ASM_CL, &
                                  N_COMPRESSED_CL, N_OFFDIAG_CL, RANK_SUM_CL, &
                                  T_GET_BLOCK=T_GETBLK_ASM, T_COMPRESS=T_COMPRESS_ASM)
    RSS_ASM_CL = RSS_NOW_MB()

    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_CL, TOL, INFO_CL, &
                                         N_COMPRESSED_CL, N_OFFDIAG_CL, RANK_SUM_CL, N_SWAPS_CL, T_FACTOR_CL)
    IF (INFO_CL .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO_CL
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    RSS_FACT_CL = RSS_NOW_MB()

    MLOC = 0
    DO li = 1, A_CL%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      MLOC = MLOC + (A_CL%BEND(I)-A_CL%BSTART(I)+1)
    END DO
    NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
    CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
    ALLOCATE(XEX_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)), B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
    CALL FILL_XEX(A_CL, GRID, NRHS_LOC, RSTART_MINE, XEX_LOCAL, MLOC)

    ! B = A*Xexact -- construction hors mesure (meme convention que les autres
    ! bench_blr_mpi_pivot_*.f90 : cette matrice AUXILIAIRE sert uniquement a
    ! generer un second membre coherent pour la validation, ce n'est pas une
    ! etape de l'algorithme lui-meme).
    BLOCK
      TYPE(MPI_BLR_MATRIX_T) :: A_CL_FOR_B
      DOUBLE PRECISION :: T_ASM_DUMMY
      INTEGER :: NC_DUMMY, NO_DUMMY
      INTEGER(8) :: RS_DUMMY
      CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_CL_FOR_B, T_ASM_DUMMY, &
                                    NC_DUMMY, NO_DUMMY, RS_DUMMY)
      CALL APPLY_MPI_BLR_MATRIX(GRID, A_CL_FOR_B, NRHS, NRHS_LOC, XEX_LOCAL, B_LOCAL, MLOC)
    END BLOCK

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0B = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_CL, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1))
    T1B = MPI_WTIME()
    T_SOLVE_CL = T1B - T0B
    RSS_SOLVE_CL = RSS_NOW_MB()

    ALLOCATE(X_GLOBAL_CL(N,NRHS))
    CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL_CL)
    ERR_CL = COMPUTE_ERR(X_GLOBAL_CL)
    DEALLOCATE(B_LOCAL, XEX_LOCAL)

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' '
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A)') '  [CLASSIQUE] Temps (s) et RSS max/processus (Mo) par etape'
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Assemblage    : temps = ', T_ASM_CL,   ' s  , RSS = ', RSS_ASM_CL,   ' Mo'
      WRITE(6,'(A,F10.3,A)') '    dont GET_BLOCK (evaluation du noyau) = ', T_GETBLK_ASM, ' s'
      WRITE(6,'(A,F10.3,A)') '    dont COMPRESS_BLOCK (SVD complete, CGESVD) = ', T_COMPRESS_ASM, ' s'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Factorisation : temps = ', T_FACTOR_CL,' s  , RSS = ', RSS_FACT_CL,  ' Mo'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Resolution    : temps = ', T_SOLVE_CL, ' s  , RSS = ', RSS_SOLVE_CL, ' Mo'
      WRITE(6,'(A,F10.3,A)') '  TOTAL         : temps = ', T_ASM_CL+T_FACTOR_CL+T_SOLVE_CL, ' s'
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Compression : ', N_COMPRESSED_CL, ' / ', N_OFFDIAG_CL, ' (', &
           100.0_WP*REAL(N_COMPRESSED_CL,WP)/REAL(MAX(N_OFFDIAG_CL,1),WP), ' %)'
      WRITE(6,'(A,I0)') '  Swaps de pivot (doit etre > 0) = ', N_SWAPS_CL
      WRITE(6,'(A,ES12.4)') '  Erreur relative ||X-Xexact|| = ', ERR_CL
      IF (ERR_CL.LT.1.0E-3_WP .AND. N_SWAPS_CL.GT.0) THEN
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_COMPARE [CLASSIQUE] : OK ==='
      ELSE
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_COMPARE [CLASSIQUE] : ECHEC ==='
      END IF
    END IF
    DEALLOCATE(X_GLOBAL_CL)

  ELSE
    !===========================================================================
    ! CHEMIN OUT-OF-CORE
    !===========================================================================
    CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_OOC, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                      T_ASM_A_OOC, N_COMPRESSED_OOC, N_OFFDIAG_OOC, RANK_SUM_OOC, &
                                      T_GET_BLOCK=T_GETBLK_ASM, T_COMPRESS=T_COMPRESS_ASM, T_IO=T_IO_ASM)
    RSS_ASMA_OOC = RSS_NOW_MB()

    T0B = MPI_WTIME()
    CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_OOC, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
    T1B = MPI_WTIME()
    T_ASM_B_OOC = T1B - T0B
    RSS_ASMB_OOC = RSS_NOW_MB()

    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_OOC, TOL, INFO_OOC, &
                                             N_COMPRESSED_OOC, N_OFFDIAG_OOC, RANK_SUM_OOC, N_SWAPS_OOC, T_FACTOR_OOC, &
                                             T_CACHE_GET=T_CACHE_GET, T_CACHE_EVICT=T_CACHE_EVICT, &
                                             T_COMPUTE=T_COMPUTE, T_MPICALL=T_MPICALL)
    IF (INFO_OOC .NE. 0) THEN
      IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO_OOC
      CALL MPI_FINALIZE(IERR)
      STOP 1
    END IF
    RSS_FACT_OOC = RSS_NOW_MB()

    ALLOCATE(X_GLOBAL_OOC(N,NRHS))
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    T0B = MPI_WTIME()
    CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_OOC, RHS, CONSUME_X_BATCH)
    T1B = MPI_WTIME()
    T_SOLVE_OOC = T1B - T0B
    RSS_SOLVE_OOC = RSS_NOW_MB()

    ERR_OOC = COMPUTE_ERR(X_GLOBAL_OOC)

    IF (MYID0 .EQ. 0) THEN
      WRITE(6,'(A)') ' '
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A)') '  [OOC] Temps (s) et RSS max/processus (Mo) par etape'
      WRITE(6,'(A)') ' ================================================================'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Assemblage A  : temps = ', T_ASM_A_OOC, ' s  , RSS = ', RSS_ASMA_OOC, ' Mo'
      WRITE(6,'(A,F10.3,A)') '    dont GET_BLOCK (evaluation du noyau) = ', T_GETBLK_ASM, ' s'
      WRITE(6,'(A,F10.3,A)') '    dont COMPRESS_BLOCK (SVD complete, CGESVD) = ', T_COMPRESS_ASM, ' s'
      WRITE(6,'(A,F10.3,A)') '    dont E/S disque (ecriture diag + cache + flush) = ', T_IO_ASM, ' s'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Assemblage B  : temps = ', T_ASM_B_OOC, ' s  , RSS = ', RSS_ASMB_OOC, ' Mo'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Factorisation : temps = ', T_FACTOR_OOC,' s  , RSS = ', RSS_FACT_OOC, ' Mo'
      WRITE(6,'(A,F10.3,A)')         '    dont cache GET (lookup+lecture disque) = ', T_CACHE_GET, ' s'
      WRITE(6,'(A,F10.3,A)')         '    dont eviction panneau (ecriture disque, inconditionnelle) = ', T_CACHE_EVICT, ' s'
      WRITE(6,'(A,F10.3,A)')         '    dont calcul pur (CGETRF/TRSM/Schur/swap) = ', T_COMPUTE, ' s'
      WRITE(6,'(A,F10.3,A)')         '    dont communications MPI (BCAST) = ', T_MPICALL, ' s'
      WRITE(6,'(A,F10.3,A)')         '    somme des 4 categories = ', T_CACHE_GET+T_CACHE_EVICT+T_COMPUTE+T_MPICALL, ' s'
      WRITE(6,'(A,F10.3,A,F10.1,A)') '  Resolution    : temps = ', T_SOLVE_OOC, ' s  , RSS = ', RSS_SOLVE_OOC,' Mo'
      WRITE(6,'(A,F10.3,A)') '  TOTAL         : temps = ', T_ASM_A_OOC+T_ASM_B_OOC+T_FACTOR_OOC+T_SOLVE_OOC, ' s'
      WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '  Compression : ', N_COMPRESSED_OOC, ' / ', N_OFFDIAG_OOC, ' (', &
           100.0_WP*REAL(N_COMPRESSED_OOC,WP)/REAL(MAX(N_OFFDIAG_OOC,1),WP), ' %)'
      WRITE(6,'(A,I0)') '  Swaps de pivot (doit etre > 0) = ', N_SWAPS_OOC
      WRITE(6,'(A,ES12.4)') '  Erreur relative ||X-Xexact|| = ', ERR_OOC
      IF (ERR_OOC.LT.1.0E-3_WP .AND. N_SWAPS_OOC.GT.0) THEN
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_COMPARE [OOC] : OK ==='
      ELSE
        WRITE(6,*) '=== BENCH_BLR_MPI_PIVOT_OOC_COMPARE [OOC] : ECHEC ==='
      END IF
    END IF
    DEALLOCATE(X_GLOBAL_OOC)
    CALL OOC_RHS_FREE(RHS)
    CALL MPI_BLR_MATRIX_OOC_FREE(A_OOC)
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! RSS reelle du processus, max sur toute la grille -- releve juste apres
  ! chaque etape (voir la description de tete).
  DOUBLE PRECISION FUNCTION RSS_NOW_MB()
    DOUBLE PRECISION :: TMP_MB
    INTEGER :: IERR3
    TMP_MB = REAL(GET_RSS_KB(),8)/1024.0D0
    CALL MPI_ALLREDUCE(TMP_MB, RSS_NOW_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR3)
  END FUNCTION RSS_NOW_MB

  REAL(WP) FUNCTION COMPUTE_ERR(XG)
    COMPLEX(WP), INTENT(IN) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: B_REF(:,:)
    INTEGER :: GI, R
    COMPUTE_ERR = -1.0_WP
    IF (MYID0 .NE. 0) RETURN
    ALLOCATE(B_REF(N,NRHS))
    DO R = 1, NRHS
      DO GI = 1, N
        B_REF(GI,R) = XEX_FORMULA(GI,R)
      END DO
    END DO
    COMPUTE_ERR = SQRT(SUM(ABS(XG-B_REF)**2)) / SQRT(SUM(ABS(B_REF)**2))
    DEALLOCATE(B_REF)
  END FUNCTION COMPUTE_ERR

  ! A(GI,GJ) = 50 sur une anti-diagonale CONFINEE AU BLOC DIAGONAL de GI --
  ! IDENTIQUE a test_mpi_blr_matrix_pivot.f90/_ooc.f90 (force un vrai
  ! pivotage QUE LE PIVOT CONFINE PEUT VOIR, meme matrice sur les DEUX
  ! chemins -- comparaison directe).
  ! Hash deterministe (I,J) -> [-1,1) -- meme recette EXACTE que
  ! bench_blr_mpi_pivot_compressibility.f90/_arena.f90 (pas de RNG a etat
  ! partage, chaque bloc recalculable independamment).
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  COMPLEX(WP) FUNCTION LOWRANK_A(GI, GJ, NN)
    INTEGER, INTENT(IN) :: GI, GJ, NN
    REAL(WP) :: DF
    INTEGER :: BSTART0, BEND0, GJ_STRONG
    BSTART0 = ((GI-1)/NB)*NB + 1
    BEND0 = MIN(BSTART0+NB-1, NN)
    GJ_STRONG = BEND0 - (GI-BSTART0)
    IF (GJ .EQ. GJ_STRONG) THEN
      LOWRANK_A = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
    ELSE
      DF = 1.0_WP / (1.0_WP + ABS(REAL(GJ_STRONG-GJ,WP))) + EPS*HASH_NOISE(GI,GJ)
      LOWRANK_A = CMPLX(DF, 0.3_WP*DF, KIND=WP)
    END IF
  END FUNCTION LOWRANK_A

  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

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

  ! B = A*Xexact, PAR PRODUIT BRUT -- INDEPENDANT de toute distribution
  ! colonnes<->processus (voir commentaire de tete et test_mpi_blr_matrix_
  ! pivot_ooc.f90) : correct pour n'importe quel decoupage NRHS_BATCH.
  SUBROUTINE FILL_RHS(ISTART, IEND, RSTART, REND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RSTART, REND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, REND-RSTART+1)
    INTEGER :: II, JJ, RR
    COMPLEX(WP) :: ACC
    COMPLEX(WP), PARAMETER :: CZERO = (0.0_WP,0.0_WP)
    DO RR = RSTART, REND
      DO II = ISTART, IEND
        ACC = CZERO
        DO JJ = 1, N
          ACC = ACC + LOWRANK_A(II,JJ,N) * XEX_FORMULA(JJ,RR)
        END DO
        BLK(II-ISTART+1, RR-RSTART+1) = ACC
      END DO
    END DO
  END SUBROUTINE FILL_RHS

  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    COMPLEX(WP), ALLOCATABLE :: X_BATCH_GLOBAL(:,:)
    INTEGER :: WIDTH
    WIDTH = BATCH_REND - BATCH_RSTART + 1
    ALLOCATE(X_BATCH_GLOBAL(N,WIDTH))
    CALL GATHER_2D_RHS(GRID, N, NB, WIDTH, X_LOCAL, LDX, NRHS_LOC_B, 0, X_BATCH_GLOBAL)
    IF (MYID0 .EQ. 0) X_GLOBAL_OOC(:,BATCH_RSTART:BATCH_REND) = X_BATCH_GLOBAL
    DEALLOCATE(X_BATCH_GLOBAL)
  END SUBROUTINE CONSUME_X_BATCH

  SUBROUTINE FILL_XEX(AB, GD, NRL, RSTART, XL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: XL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, GRR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          GRR = RSTART + RR - 1
          XL(RS + (GII-AB%BSTART(II)), RR) = XEX_FORMULA(GII, GRR)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_XEX

  SUBROUTINE APPLY_MPI_BLR_MATRIX(GD, AB, NR, NRL, XL, YL, ML)
    TYPE(GRID_T), INTENT(IN)           :: GD
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    INTEGER, INTENT(IN)                :: NR, NRL, ML
    COMPLEX(WP), INTENT(IN)            :: XL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), INTENT(OUT)           :: YL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), PARAMETER :: CONE=(1.0_WP,0.0_WP), CZERO=(0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: XJ_LOC(:,:), XJ_FULL(:,:), CONTRIB(:,:), CONTRIB_LOC(:,:)
    INTEGER, ALLOCATABLE :: RSTARTS(:), RLOCS(:), RCOUNTS(:), RDISPLS(:)
    INTEGER :: J, lj, li, I, NJ, MI, RS, RE, C, IERR2, PROWJ, PCOLJ, MLOC_ALL

    IF (NRL .GT. 0) YL = CZERO

    ALLOCATE(RSTARTS(0:GD%PC-1), RLOCS(0:GD%PC-1), RCOUNTS(0:GD%PC-1), RDISPLS(0:GD%PC-1))
    DO C = 0, GD%PC-1
      CALL RHS_LOCAL_RANGE_MBM(NR, GD%PC, C, RSTARTS(C), RLOCS(C))
    END DO

    DO J = 1, AB%NBLK
      NJ = AB%BEND(J)-AB%BSTART(J)+1
      PROWJ = MOD(J-1, GD%PR); PCOLJ = BLR_COL_OWNER(AB%CMAP, J, GD)

      ALLOCATE(XJ_FULL(NJ,NR))
      IF (GD%MYROW .EQ. PROWJ) THEN
        RS = LOCAL_ROW_OFFSET(AB, GD, J)
        ALLOCATE(XJ_LOC(NJ,MAX(NRL,1)))
        IF (NRL .GT. 0) XJ_LOC = XL(RS:RS+NJ-1,1:NRL)
        DO C = 0, GD%PC-1
          RCOUNTS(C) = NJ*RLOCS(C); RDISPLS(C) = NJ*(RSTARTS(C)-1)
        END DO
        CALL MPI_ALLGATHERV(XJ_LOC, NJ*NRL, MPI_COMPLEX, XJ_FULL, RCOUNTS, RDISPLS, MPI_COMPLEX, &
                             GD%ROW_COMM, IERR2)
        DEALLOCATE(XJ_LOC)
      END IF
      CALL MPI_BCAST(XJ_FULL, NJ*NR, MPI_COMPLEX, PROWJ, GD%COL_COMM, IERR2)

      MLOC_ALL = 0
      DO li = 1, AB%NBLKR_LOC
        I = (li-1)*GD%PR + GD%MYROW + 1
        MLOC_ALL = MLOC_ALL + (AB%BEND(I)-AB%BSTART(I)+1)
      END DO
      ALLOCATE(CONTRIB(MLOC_ALL,NR))
      IF (GD%MYCOL .EQ. PCOLJ) THEN
        lj = BLR_COL_LOCAL(AB%CMAP, J, GD)
        CONTRIB = CZERO
        RE = 0
        DO li = 1, AB%NBLKR_LOC
          I = (li-1)*GD%PR + GD%MYROW + 1
          MI = AB%BEND(I)-AB%BSTART(I)+1
          IF (I .EQ. J) THEN
            CALL CGEMM('No transpose','No transpose', MI, NR, NJ, CONE, AB%DIAGBLK(I)%D, MI, &
                        XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
          ELSE
            CALL BLR_BLOCK_APPLY(AB%OFFBLK(li,lj), NR, CONE, XJ_FULL, NJ, CZERO, CONTRIB(RE+1,1), MLOC_ALL)
          END IF
          RE = RE + MI
        END DO
      END IF
      ALLOCATE(CONTRIB_LOC(MLOC_ALL,MAX(NRL,1)))
      DO C = 0, GD%PC-1
        RCOUNTS(C) = MLOC_ALL*RLOCS(C); RDISPLS(C) = MLOC_ALL*(RSTARTS(C)-1)
      END DO
      CALL MPI_SCATTERV(CONTRIB, RCOUNTS, RDISPLS, MPI_COMPLEX, CONTRIB_LOC, MLOC_ALL*NRL, &
                         MPI_COMPLEX, PCOLJ, GD%ROW_COMM, IERR2)
      IF (NRL .GT. 0) YL(1:MLOC_ALL,1:NRL) = YL(1:MLOC_ALL,1:NRL) + CONTRIB_LOC(1:MLOC_ALL,1:NRL)
      DEALLOCATE(CONTRIB, CONTRIB_LOC)
      DEALLOCATE(XJ_FULL)
    END DO
    DEALLOCATE(RSTARTS, RLOCS, RCOUNTS, RDISPLS)
  END SUBROUTINE APPLY_MPI_BLR_MATRIX

END PROGRAM BENCH_BLR_MPI_PIVOT_OOC_COMPARE
