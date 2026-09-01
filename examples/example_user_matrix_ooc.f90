!===============================================================================
! EXEMPLE MINIMAL -- point de depart pour brancher VOTRE PROPRE systeme
! lineaire A*X=B sur le solveur BLR+MPI OUT-OF-CORE (RAM bornee par
! MAX_MB_PER_PROCESS, quelle que soit la taille de A). Version tout-en-RAM :
! voir example_user_matrix.f90, meme structure.
!
! CE QU'IL FAUT REMPLACER (marque ci-dessous) :
!   - GET_BLOCK : renvoie un bloc de VOTRE matrice A -- appelee PAR BLOCS,
!     jamais pour A entiere (c'est ce qui permet a l'assemblage OOC de ne
!     JAMAIS avoir toute A en RAM, meme transitoirement).
!   - FILL_RHS : renvoie VOTRE second membre b, meme principe.
! RIEN D'AUTRE n'a besoin d'etre modifie.
!
! Le placeholder ci-dessous (identique a example_user_matrix.f90) n'est
! qu'un exemple RUNNABLE pour verifier que tout compile et tourne AVANT de
! brancher vos vraies donnees.
!
! CHOISIR MAX_MB_PER_PROCESS : NE PAS prendre le plus petit possible -- un
! budget sous l'ensemble de travail reel du processus fait THRASHER le
! cache (temps tres degrade pour une RAM presque inchangee). Voir
! examples/README.md pour la mesure detaillee et la methode pour choisir
! une valeur raisonnable (augmenter jusqu'a ce que la RSS reelle cesse de
! croitre).
!
! Usage : mpirun -np <PR*PC> ./example_user_matrix_ooc PR PC [N] [NB] [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR]
!         (defauts : N=4000 NB=200 NRHS=1 TOL=1e-4 MAX_MB_PER_PROCESS=64 NRHS_BATCH=2 SCRATCH_DIR=./blr_ooc_scratch)
!===============================================================================
PROGRAM EXAMPLE_USER_MATRIX_OOC
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_BLR
  TYPE(OOC_RHS_T) :: RHS
  COMPLEX(WP), ALLOCATABLE :: X_GLOBAL(:,:)
  INTEGER :: INFO, N_SWAPS, N_COMPRESSED, N_OFFDIAG, I
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR
  REAL(WP) :: RESID

  N = 4000; NB = 200; NRHS = 1; TOL = 1.0E-4_WP
  MAX_MB_PER_PROCESS = 64.0_WP; NRHS_BATCH = 2; SCRATCH_DIR = './blr_ooc_scratch'

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./example_user_matrix_ooc PR PC [N] [NB] [NRHS] [TOL] '// &
         '[MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR]'
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

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1)') &
       ' [OOC] Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
       ' , MAX_MB_PER_PROCESS=', MAX_MB_PER_PROCESS

  ! --- 1) assemblage BLR distribue, out-of-core -> disque -------------------
  CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                    T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE.)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF

  ! --- 2) second membre, SEPAREMENT, aussi sur disque ------------------------
  CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)

  ! --- 3) factorisation LU, avec pivot, out-of-core --------------------------
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                           N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                           SHOW_PROGRESS=.TRUE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'

  ! --- 4) resolution, out-of-core, par paquets de second membre --------------
  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH, SHOW_PROGRESS=.TRUE.)
  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
    DO I = 1, MIN(5, N)
      WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
    END DO
  END IF

  ! --- 5) verification ||A*X-B||/||B|| -- APPEL COLLECTIF (tous les rangs) --
  RESID = COMPUTE_RESIDUAL(X_GLOBAL, GRID, N, NB, NRHS)
  IF (MYID0.EQ.0) WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID

  DEALLOCATE(X_GLOBAL)
  CALL OOC_RHS_FREE(RHS)
  CALL MPI_BLR_MATRIX_OOC_FREE(A_BLR)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  !===========================================================================
  ! REMPLACEZ CETTE FONCTION par votre propre matrice A.
  ! BLK(i,j) doit contenir A(ISTART+i-1, JSTART+j-1) -- appelee par blocs,
  ! jamais pour la matrice entiere (c'est ce qui garde l'assemblage OOC
  ! borne en memoire). Placeholder RUNNABLE, sans signification physique.
  !===========================================================================
  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        IF (II .EQ. JJ) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        ELSE
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(1.0_WP/(1.0_WP+ABS(REAL(II-JJ,WP))), 0.0_WP, KIND=WP)
        END IF
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  !===========================================================================
  ! REMPLACEZ CETTE FONCTION par votre propre second membre b.
  ! BLK(i,r) doit contenir b(ISTART+i-1, RSTART+r-1) -- appelee par blocs
  ! ET par paquets de colonnes (jamais tout B en RAM d'un coup, meme
  ! principe que GET_BLOCK). Placeholder : b = vecteur de 1.
  !===========================================================================
  SUBROUTINE FILL_RHS(ISTART, IEND, RSTART, REND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RSTART, REND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, REND-RSTART+1)
    BLK = CMPLX(1.0_WP, 0.0_WP, KIND=WP)
  END SUBROUTINE FILL_RHS

  ! Rassemble CE paquet resolu dans X_GLOBAL (verification/demonstration --
  ! un vrai cas massivement second membre ecrirait plutot chaque paquet
  ! directement sur disque au lieu de tout accumuler, voir CONSUME_X_BATCH_
  ! IFACE dans mpi_blr_matrix_mod.f90).
  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    COMPLEX(WP), ALLOCATABLE :: X_BATCH_GLOBAL(:,:)
    INTEGER :: WIDTH
    WIDTH = BATCH_REND - BATCH_RSTART + 1
    ALLOCATE(X_BATCH_GLOBAL(N,WIDTH))
    CALL GATHER_2D_RHS(GRID, N, NB, WIDTH, X_LOCAL, LDX, NRHS_LOC_B, 0, X_BATCH_GLOBAL)
    IF (MYID0 .EQ. 0) X_GLOBAL(:,BATCH_RSTART:BATCH_REND) = X_BATCH_GLOBAL
    DEALLOCATE(X_BATCH_GLOBAL)
  END SUBROUTINE CONSUME_X_BATCH

  !===========================================================================
  ! ||A*X-B||/||B|| -- APPEL COLLECTIF (MPI_BCAST/MPI_ALLREDUCE a
  ! l'interieur, tous les rangs doivent l'appeler). A*X calcule par blocs de
  ! lignes via GET_BLOCK + CGEMM (BLAS), lignes reparties entre TOUS les
  ! rangs. B est recalcule via FILL_RHS (fonction PURE des indices, jamais
  ! besoin de sauvegarder une copie contrairement au chemin classique ou la
  ! resolution ecrase B en place -- ici B reste sur disque, jamais touche).
  ! Tampon BLAS BORNE PAR MAX_MB_PER_PROCESS (jamais NB seul) : voir
  ! example_blr_mpi_pivot_ooc.f90/COMPUTE_RESIDUAL pour le detail complet
  ! de ce correctif, applique ici a l'identique.
  !===========================================================================
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG, GD, NN, NNB, NR)
    TYPE(GRID_T), INTENT(IN)   :: GD
    INTEGER, INTENT(IN)        :: NN, NNB, NR
    COMPLEX(WP), INTENT(INOUT) :: XG(NN,NR)
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), AXBLOCK(:,:), BBLOCK(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK, CHUNK_ROWS
    INTEGER :: GI, R, IERR2
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2

    CALL MPI_BCAST(XG, NN*NR, MPI_COMPLEX, 0, GD%COMM, IERR2)

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)

    CHUNK_ROWS = MAX(1, MIN(NNB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(NN,WP)))))

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP
    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CHUNK_ROWS
      REND_CHUNK = MIN(RSTART_CHUNK+CHUNK_ROWS-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(ABLOCK(NCHUNK,NN), AXBLOCK(NCHUNK,NR), BBLOCK(NCHUNK,NR))
      CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, NN, ABLOCK)
      CALL FILL_RHS(RSTART_CHUNK, REND_CHUNK, 1, NR, BBLOCK)
      CALL CGEMM('N','N', NCHUNK, NR, NN, CONE, ABLOCK, NCHUNK, XG, NN, CZERO, AXBLOCK, NCHUNK)
      DO R = 1, NR
        DO GI = 1, NCHUNK
          NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(AXBLOCK(GI,R)-BBLOCK(GI,R))**2,WP)
          DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BBLOCK(GI,R))**2,WP)
        END DO
      END DO
      DEALLOCATE(ABLOCK, AXBLOCK, BBLOCK)
    END DO

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GD%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GD%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

END PROGRAM EXAMPLE_USER_MATRIX_OOC
