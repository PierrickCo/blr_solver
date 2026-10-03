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
! RATIO, USE_RSVD : memes parametres, memes valeurs par defaut, que
! example_user_matrix.f90 (version RAM) -- voir son commentaire de tete
! pour le detail complet, notamment sur la lecture des deux indicateurs
! de compression ci-dessous.
!
! X (SOLUTION) : ECRITE SUR DISQUE, jamais accumulee en RAM -- coherent
! avec A et B. CONSUME_X_BATCH
! ecrit CHAQUE paquet directement dans UN SEUL FICHIER PAR PROCESSUS
! (x_rank_NNNNNN.dat dans SCRATCH_DIR, meme convention que a_rank_*/
! b_rank_*), SANS AUCUNE communication MPI (chaque processus ecrit SA
! PROPRE part locale, jamais rassemblee a l'ecriture) -- le mode a
! utiliser pour un VRAI cas massivement multi-second-membre (NRHS
! enorme), ou un tableau (N,NRHS) complet ne tiendrait pas en RAM.
! Ecriture SEQUENTIELLE (fichier STREAM, aucun POS= explicite -- chaque
! WRITE reprend automatiquement ou le precedent s'est arrete), un
! enregistrement AUTO-DESCRIPTIF par paquet (6 entiers d'en-tete :
! BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B, puis les
! donnees LDX*NRHS_LOC_B).
!
! VERIFICATION SANS SOLUTION CONNUE (CHECK_RESIDUAL) : ||A*X-B||/||B||.
! A (dense, telle que renvoyee par GET_BLOCK, AVANT compression BLR) et B
! sont mis en cache sur disque UNE SEULE FOIS, DES QUE POSSIBLE (juste
! apres l'assemblage / apres l'assemblage du RHS, donc AVANT que la
! factorisation n'ecrase A par L/U en place) -- acheck_rank_NNNNNN.dat et
! bcheck_rank_NNNNNN.dat, meme convention 1-fichier-par-processus, meme
! decoupage par blocs de lignes que la verification elle-meme (borne par
! MAX_MB_PER_PROCESS). COMPUTE_RESIDUAL les RELIT ensuite au lieu de
! rappeler GET_BLOCK/FILL_RHS -- utile si ces callbacks sont couteux (mesh,
! quadrature, ressource externe...), pour ne jamais les evaluer deux fois.
! X_LOCAL est relu depuis x_rank_NNNNNN.dat (meme technique que
! B_ORIG_LOCAL gathere dans example_user_matrix.f90/COMPUTE_RESIDUAL --
! GATHER_2D_RHS puis MPI_BCAST). Des qu'un bloc de B a servi au calcul du
! residu, il n'est plus utile : X est ECRIT A SA PLACE, DANS LE MEME
! FICHIER (meme forme, reecriture en place par POS=, comme
! OOC_WRITE_BLOCK_FIXED dans blr_ooc_mod.f90) -- economie d'espace disque
! demandee explicitement : a la fin, x_rank_*.dat et acheck_rank_*.dat
! sont supprimes, seul bcheck_rank_*.dat survit et contient desormais X.
! Memoire bornee : jamais plus d'UN SEUL bloc de lignes de A, ni plus
! d'UN SEUL paquet de X pleinement rassemble, a la fois.
!
! Usage : mpirun -np <PR*PC> ./example_user_matrix_ooc PR PC [N] [NB] [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [RATIO] [USE_RSVD] [CHECK_RESIDUAL] [PROFILE_USAGE] [DUMP_STRUCTURE]
!         (defauts : N=4000 NB=200 NRHS=1 TOL=1e-4 MAX_MB_PER_PROCESS=64 NRHS_BATCH=2 SCRATCH_DIR=./blr_ooc_scratch RATIO=1/3 USE_RSVD=1 CHECK_RESIDUAL=0 PROFILE_USAGE=0 DUMP_STRUCTURE=0)
!
! PROFILE_USAGE / DUMP_STRUCTURE : voir example_blr_mpi_pivot.f90 pour le
! detail complet (memes parametres, memes scripts de trace). Comme pour
! example_blr_mpi_pivot_ooc.f90, l'echantillonnage tous les 5% DANS la
! resolution redemarre a chaque PAQUET de second membre. X est deja ecrite
! directement sur disque par paquet (voir plus haut) -- aucune economie
! supplementaire n'est necessaire ici (contrairement a example_blr_mpi_
! pivot_ooc.f90, qui accumulait X_GLOBAL en RAM).
!===============================================================================
PROGRAM EXAMPLE_USER_MATRIX_OOC
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS, RATIO
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_BLR
  TYPE(OOC_RHS_T) :: RHS
  INTEGER :: INFO, N_SWAPS, N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR
  REAL(WP) :: RESID
  INTEGER :: USE_RSVD_INT, CHECK_RESIDUAL_INT, PROFILE_USAGE_INT, DUMP_STRUCTURE_INT
  LOGICAL :: USE_RSVD_FLAG, CHECK_RESIDUAL_FLAG, PROFILE_USAGE_FLAG, DUMP_STRUCTURE_FLAG
  INTEGER :: UNIT_X, NBATCH_WRITTEN
  CHARACTER(LEN=16) :: RANKSTR

  ! Suivi memoire/CPU par etape (PROFILE_USAGE) -- voir example_blr_mpi_
  ! pivot.f90 pour le detail complet de la technique.
  INTEGER, PARAMETER :: MAX_PROFILE_SAMPLES = 128
  CHARACTER(LEN=40) :: PROF_STAGE(MAX_PROFILE_SAMPLES)
  DOUBLE PRECISION :: PROF_T(MAX_PROFILE_SAMPLES), PROF_RSS_MB(MAX_PROFILE_SAMPLES), PROF_CPU_PCT(MAX_PROFILE_SAMPLES)
  INTEGER :: PROF_N
  DOUBLE PRECISION :: PROF_T_START, PROF_T_PREV
  REAL :: PROF_CPU_PREV

  N = 4000; NB = 200; NRHS = 1; TOL = 1.0E-4_WP
  MAX_MB_PER_PROCESS = 64.0_WP; NRHS_BATCH = 2; SCRATCH_DIR = './blr_ooc_scratch'
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 1; CHECK_RESIDUAL_INT = 0
  PROFILE_USAGE_INT = 0; DUMP_STRUCTURE_INT = 0
  PROF_N = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./example_user_matrix_ooc PR PC [N] [NB] [NRHS] [TOL] '// &
         '[MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [RATIO] [USE_RSVD] [CHECK_RESIDUAL] '// &
         '[PROFILE_USAGE] [DUMP_STRUCTURE]'
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
    CALL GET_COMMAND_ARGUMENT(10, ARG); READ(ARG,*) RATIO
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN
    CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) USE_RSVD_INT
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
  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
  CHECK_RESIDUAL_FLAG = (CHECK_RESIDUAL_INT .NE. 0)
  PROFILE_USAGE_FLAG = (PROFILE_USAGE_INT .NE. 0)
  DUMP_STRUCTURE_FLAG = (DUMP_STRUCTURE_INT .NE. 0)

  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1,A,F6.3)') &
         ' [OOC] Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
         ' , MAX_MB_PER_PROCESS=', MAX_MB_PER_PROCESS, ' , RATIO=', RATIO
    IF (USE_RSVD_FLAG) THEN
      WRITE(6,'(A)') ' Compression : SVD randomisee adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE)'
    ELSE
      WRITE(6,'(A)') ' Compression : SVD complete (COMPRESS_BLOCK, CGESVD)'
    END IF
    IF (PROFILE_USAGE_FLAG) THEN
      WRITE(6,'(A)') ' Profilage RSS/CPU par etape active : usage_profile_rank<N>.csv par processus'
    END IF
    IF (DUMP_STRUCTURE_FLAG) THEN
      WRITE(6,'(A)') ' Structure BLR (dense/bas-rang) exportee : blr_structure_<etape>_rank<N>.csv par processus'
    END IF
  END IF

  PROF_T_START = MPI_WTIME()
  PROF_T_PREV = 0.0D0
  CALL CPU_TIME(PROF_CPU_PREV)
  CALL PROFILE_SAMPLE('debut')

  ! --- 1) assemblage BLR distribue, out-of-core -> disque -------------------
  CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                    T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., &
                                    USE_RSVD=USE_RSVD_FLAG, PROGRESS_CALLBACK=PROFILE_CALLBACK_ASSEMBLY)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (bas-rang / dense equivalent) : ', &
         REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
  END IF
  CALL PROFILE_SAMPLE('assemblage')
  CALL DUMP_BLR_STRUCTURE('assembly')

  ! Cache de A pour CHECK_RESIDUAL : AVANT la factorisation, qui ecrase A
  ! par L/U en place -- GET_BLOCK n'est plus jamais rappelee ensuite.
  IF (CHECK_RESIDUAL_FLAG) CALL CACHE_A_FOR_CHECK(GRID, N, NB)

  ! --- 2) second membre, SEPAREMENT, aussi sur disque ------------------------
  CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
  CALL PROFILE_SAMPLE('second_membre')

  ! Cache de B pour CHECK_RESIDUAL -- FILL_RHS n'est plus jamais rappelee.
  IF (CHECK_RESIDUAL_FLAG) CALL CACHE_B_FOR_CHECK(GRID, N, NB, NRHS, NRHS_BATCH)

  ! --- 3) factorisation LU, avec pivot, out-of-core --------------------------
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                           N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                           SHOW_PROGRESS=.TRUE., PROGRESS_CALLBACK=PROFILE_CALLBACK_FACTOR)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (fin factorisation) : ', &
         REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
  END IF
  CALL PROFILE_SAMPLE('factorisation')
  CALL DUMP_BLR_STRUCTURE('factorization')

  ! --- 4) resolution, out-of-core, par paquets de second membre --------------
  ! Fichier de CE processus, jamais partage -- meme convention de nommage
  ! que a_rank_*/b_rank_* (blr_ooc_mod.f90). STREAM sans POS= explicite :
  ! chaque WRITE reprend a la fin du precedent, comme un flux ordinaire.
  WRITE(RANKSTR,'(I6.6)') GRID%MYID
  OPEN(NEWUNIT=UNIT_X, FILE=TRIM(SCRATCH_DIR)//'/x_rank_'//TRIM(RANKSTR)//'.dat', &
       ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')
  NBATCH_WRITTEN = 0

  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH, SHOW_PROGRESS=.TRUE., &
                                             PROGRESS_CALLBACK=PROFILE_CALLBACK_SOLVE)

  CLOSE(UNIT_X)
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A)') ' Resolution terminee. ', NBATCH_WRITTEN, &
       ' paquet(s) ecrit(s) sur disque par processus (voir SCRATCH_DIR/x_rank_*.dat).'
  CALL PROFILE_SAMPLE('resolution')

  ! --- 5) verification ||A*X-B||/||B|| -- APPEL COLLECTIF (tous les rangs) --
  ! Consomme et supprime x_rank_*.dat et acheck_rank_*.dat au passage ;
  ! bcheck_rank_*.dat survit seul, contenant desormais X (voir la fonction).
  IF (CHECK_RESIDUAL_FLAG) THEN
    RESID = COMPUTE_RESIDUAL(GRID, N, NB, NRHS)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
      WRITE(6,'(A)') ' x_rank_*.dat et acheck_rank_*.dat supprimes ; '// &
           'bcheck_rank_*.dat contient desormais X (economie d''espace disque).'
    END IF
    CALL PROFILE_SAMPLE('verif_residu')
  END IF

  CALL PROFILE_SAMPLE('fin')
  CALL PROFILE_WRITE_CSV()

  CALL OOC_RHS_FREE(RHS)
  CALL MPI_BLR_MATRIX_OOC_FREE(A_BLR)
  CALL MPI_FINALIZE(IERR)

CONTAINS

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

  ! Callbacks PROGRESS_CALLBACK -- voir example_blr_mpi_pivot_ooc.f90 pour
  ! le detail complet (NDONE/NTOTAL de la resolution redemarrent a chaque
  ! PAQUET de second membre).
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
  ! example_blr_mpi_pivot_ooc.f90 pour le detail complet de la technique
  ! (OOC_BLOCK_DESC_T, ownership du bloc diagonal via OFFSET>=0).
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

  ! Voir example_user_matrix.f90 pour le detail du correctif.
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

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

  ! Ecrit CE paquet resolu dans le fichier de CE processus -- AUCUNE
  ! communication MPI (chaque processus ecrit SA PROPRE part locale
  ! X_LOCAL, telle quelle, jamais rassemblee) -- voir le commentaire de
  ! tete. Affiche au passage les premieres valeurs LOCALES du rang 0 pour
  ! le tout premier paquet (BATCH_RSTART=1) -- pas besoin de rassembler
  ! quoi que ce soit pour un simple apercu, le rang (MYROW=0,MYCOL=0)
  ! possede deja localement la ligne globale 1 et la colonne globale 1.
  SUBROUTINE CONSUME_X_BATCH(BATCH_RSTART, BATCH_REND, GSTART, GEND, X_LOCAL, LDX, NRHS_LOC_B)
    INTEGER, INTENT(IN)      :: BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    COMPLEX(WP), INTENT(IN)  :: X_LOCAL(LDX,MAX(NRHS_LOC_B,1))
    INTEGER :: II

    WRITE(UNIT_X) BATCH_RSTART, BATCH_REND, GSTART, GEND, LDX, NRHS_LOC_B
    IF (LDX.GT.0 .AND. NRHS_LOC_B.GT.0) WRITE(UNIT_X) X_LOCAL(1:LDX,1:NRHS_LOC_B)
    NBATCH_WRITTEN = NBATCH_WRITTEN + 1

    IF (MYID0.EQ.0 .AND. BATCH_RSTART.EQ.1 .AND. GSTART.EQ.1) THEN
      WRITE(6,'(A)') ' Premieres valeurs LOCALES du rang 0, X(:,1) :'
      DO II = 1, MIN(5, LDX)
        WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', II, ',1) = ', REAL(X_LOCAL(II,1)), ' + i*', AIMAG(X_LOCAL(II,1))
      END DO
    END IF
  END SUBROUTINE CONSUME_X_BATCH

  ! Ecrit UNE FOIS le cache dense de A (acheck_rank_NNNNNN.dat), par blocs
  ! de lignes bornes par MAX_MB_PER_PROCESS -- meme decoupage que
  ! COMPUTE_RESIDUAL, pour que la relecture s'y cale exactement. Appelee
  ! AVANT la factorisation (qui ecrase A par L/U en place).
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

  ! Ecrit UNE FOIS le cache de B (bcheck_rank_NNNNNN.dat), par (paquet de
  ! RHS x bloc de lignes) -- MEME double decoupage que COMPUTE_RESIDUAL
  ! (colonnes par paquets de NRB, lignes bornees par MAX_MB_PER_PROCESS) :
  ! ce fichier sera relu bloc par bloc puis chaque bloc REECRIT EN PLACE
  ! par X une fois son residu calcule (voir COMPUTE_RESIDUAL).
  SUBROUTINE CACHE_B_FOR_CHECK(GD, NN, NNB, NR, NRB)
    TYPE(GRID_T), INTENT(IN) :: GD
    INTEGER, INTENT(IN)      :: NN, NNB, NR, NRB
    COMPLEX(WP), ALLOCATABLE :: BBLOCK(:,:)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK, CR
    INTEGER :: BSTART, BEND, WIDTH, UNIT_B
    CHARACTER(LEN=16) :: RSTR

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)
    CR = MAX(1, MIN(NNB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(NN,WP)))))

    WRITE(RSTR,'(I6.6)') GD%MYID
    OPEN(NEWUNIT=UNIT_B, FILE=TRIM(SCRATCH_DIR)//'/bcheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')

    DO BSTART = 1, NR, NRB
      BEND = MIN(BSTART+NRB-1, NR)
      WIDTH = BEND - BSTART + 1
      DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CR
        REND_CHUNK = MIN(RSTART_CHUNK+CR-1, REND_MINE)
        NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
        ALLOCATE(BBLOCK(NCHUNK,WIDTH))
        CALL FILL_RHS(RSTART_CHUNK, REND_CHUNK, BSTART, BEND, BBLOCK)
        WRITE(UNIT_B) BSTART, BEND, RSTART_CHUNK, REND_CHUNK, NCHUNK, WIDTH
        WRITE(UNIT_B) BBLOCK(1:NCHUNK,1:WIDTH)
        DEALLOCATE(BBLOCK)
      END DO
    END DO
    CLOSE(UNIT_B)
  END SUBROUTINE CACHE_B_FOR_CHECK

  !===========================================================================
  ! ||A*X-B||/||B|| -- APPEL COLLECTIF, PAR PAQUET (meme decoupage que
  ! l'ecriture de X, voir le commentaire de tete) : chaque rang RELIT son
  ! PROPRE fichier x_rank_NNNNNN.dat de bout en bout (memes enregistrements
  ! qu'ecrits par CONSUME_X_BATCH, dans le meme ordre), et pour CHAQUE
  ! paquet :
  !   1) rassemble ce paquet de X (GATHER_2D_RHS -- EXACTEMENT la meme
  !      technique que B_ORIG_LOCAL -> BG dans example_user_matrix.f90/
  !      COMPUTE_RESIDUAL, X_LOCAL relu ici jouant le role de B_ORIG_LOCAL
  !      la-bas) puis diffuse (MPI_BCAST) a tous les rangs ;
  !   2) relit A (acheck_rank_*.dat, une passe complete par paquet -- ne
  !      depend pas des colonnes mais reste borne en RAM par
  !      MAX_MB_PER_PROCESS via le meme decoupage CHUNK_ROWS) et B
  !      (bcheck_rank_*.dat, ecrit avec exactement ce decoupage par
  !      CACHE_A_FOR_CHECK/CACHE_B_FOR_CHECK) au lieu de rappeler
  !      GET_BLOCK/FILL_RHS, puis calcule A*X par blocs de lignes CGEMM
  !      repartis entre tous les rangs ;
  !   3) une fois le residu de ce bloc calcule, B ne sert plus a rien :
  !      X (deja rassemble dans XG_BATCH) est REECRIT A SA PLACE dans
  !      bcheck_rank_*.dat (POS=, meme forme, comme OOC_WRITE_BLOCK_FIXED)
  !      -- economise le fichier x_rank_*.dat, supprime en fin de fonction
  !      avec acheck_rank_*.dat ; seul bcheck_rank_*.dat survit, avec X.
  ! Memoire bornee : jamais plus d'UN SEUL bloc de lignes de A, ni plus
  ! d'UN SEUL paquet de X pleinement rassemble, a la fois -- coherent
  ! avec le reste de l'algorithme OOC.
  !===========================================================================
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(GD, NN, NNB, NR)
    TYPE(GRID_T), INTENT(IN) :: GD
    INTEGER, INTENT(IN)      :: NN, NNB, NR
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), AXBLOCK(:,:), BBLOCK(:,:)
    COMPLEX(WP), ALLOCATABLE :: X_LOCAL_R(:,:), XG_BATCH(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK, CHUNK_ROWS
    INTEGER :: GI, R, IERR2, UNIT_XR, UNIT_A, UNIT_B, IOS
    INTEGER :: HDR_BRSTART, HDR_BREND, HDR_GSTART, HDR_GEND, HDR_LDX, HDR_NRHSLOC, WIDTH
    INTEGER :: HA_RSTART, HA_REND, HA_NCHUNK
    INTEGER :: HB_BSTART, HB_BEND, HB_RSTART, HB_REND, HB_NCHUNK, HB_WIDTH
    INTEGER(8) :: PAYLOAD_POS
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2
    CHARACTER(LEN=16) :: RSTR

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)
    CHUNK_ROWS = MAX(1, MIN(NNB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(NN,WP)))))

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP

    ! Relit SON PROPRE fichier X (deja ferme en ecriture par le programme
    ! principal avant cet appel) -- memes enregistrements, meme ordre.
    ! bcheck ouvert en LECTURE-ECRITURE : chaque bloc lu sera reecrit par
    ! X juste apres avoir servi au calcul du residu.
    WRITE(RSTR,'(I6.6)') GD%MYID
    OPEN(NEWUNIT=UNIT_XR, FILE=TRIM(SCRATCH_DIR)//'/x_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')
    OPEN(NEWUNIT=UNIT_B, FILE=TRIM(SCRATCH_DIR)//'/bcheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READWRITE')

    DO
      READ(UNIT_XR, IOSTAT=IOS) HDR_BRSTART, HDR_BREND, HDR_GSTART, HDR_GEND, HDR_LDX, HDR_NRHSLOC
      IF (IOS .NE. 0) EXIT   ! fin de fichier -- tous les paquets traites
      ALLOCATE(X_LOCAL_R(MAX(HDR_LDX,1),MAX(HDR_NRHSLOC,1)))
      IF (HDR_LDX.GT.0 .AND. HDR_NRHSLOC.GT.0) READ(UNIT_XR) X_LOCAL_R(1:HDR_LDX,1:HDR_NRHSLOC)

      ! 1) rassemble CE paquet de X (comme B_ORIG_LOCAL -> BG cote RAM),
      ! puis diffuse a tous les rangs (chaque rang doit avoir TOUT X pour
      ! calculer A*X sur sa propre tranche de LIGNES).
      WIDTH = HDR_BREND - HDR_BRSTART + 1
      ALLOCATE(XG_BATCH(NN,WIDTH))
      CALL GATHER_2D_RHS(GD, NN, NNB, WIDTH, X_LOCAL_R, MAX(HDR_LDX,1), HDR_NRHSLOC, 0, XG_BATCH)
      CALL MPI_BCAST(XG_BATCH, NN*WIDTH, MPI_COMPLEX, 0, GD%COMM, IERR2)
      DEALLOCATE(X_LOCAL_R)

      ! 2) A*X et B pour CETTE plage de colonnes uniquement, releus depuis
      ! le cache (acheck rembobine a chaque paquet -- A ne depend pas des
      ! colonnes ; bcheck avance sequentiellement, meme ordre qu'ecrit).
      OPEN(NEWUNIT=UNIT_A, FILE=TRIM(SCRATCH_DIR)//'/acheck_rank_'//TRIM(RSTR)//'.dat', &
           ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')
      DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CHUNK_ROWS
        REND_CHUNK = MIN(RSTART_CHUNK+CHUNK_ROWS-1, REND_MINE)
        NCHUNK = REND_CHUNK - RSTART_CHUNK + 1

        READ(UNIT_A) HA_RSTART, HA_REND, HA_NCHUNK
        ALLOCATE(ABLOCK(NCHUNK,NN))
        READ(UNIT_A) ABLOCK(1:NCHUNK,1:NN)

        READ(UNIT_B) HB_BSTART, HB_BEND, HB_RSTART, HB_REND, HB_NCHUNK, HB_WIDTH
        INQUIRE(UNIT=UNIT_B, POS=PAYLOAD_POS)
        ALLOCATE(BBLOCK(NCHUNK,WIDTH))
        READ(UNIT_B) BBLOCK(1:NCHUNK,1:WIDTH)

        ALLOCATE(AXBLOCK(NCHUNK,WIDTH))
        CALL CGEMM('N','N', NCHUNK, WIDTH, NN, CONE, ABLOCK, NCHUNK, XG_BATCH, NN, CZERO, AXBLOCK, NCHUNK)
        DO R = 1, WIDTH
          DO GI = 1, NCHUNK
            NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(AXBLOCK(GI,R)-BBLOCK(GI,R))**2,WP)
            DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BBLOCK(GI,R))**2,WP)
          END DO
        END DO

        ! 3) B vient de servir : plus utile -- X ecrase sa place sur
        ! disque (meme forme NCHUNK x WIDTH), evite de garder x_rank ET
        ! bcheck simultanement.
        WRITE(UNIT_B, POS=PAYLOAD_POS) XG_BATCH(RSTART_CHUNK:REND_CHUNK, 1:WIDTH)

        DEALLOCATE(ABLOCK, AXBLOCK, BBLOCK)
      END DO
      CLOSE(UNIT_A)
      DEALLOCATE(XG_BATCH)
    END DO
    CLOSE(UNIT_XR, STATUS='DELETE')
    CLOSE(UNIT_B)

    ! A ne sert plus a rien (utilisee ci-dessus pour tous les paquets).
    OPEN(NEWUNIT=UNIT_A, FILE=TRIM(SCRATCH_DIR)//'/acheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')
    CLOSE(UNIT_A, STATUS='DELETE')

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GD%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GD%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

END PROGRAM EXAMPLE_USER_MATRIX_OOC
