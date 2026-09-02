!===============================================================================
! VARIANTE OUT-OF-CORE de example_blr_mpi_pivot.f90 -- MEME algorithme
! (pivot partiel confine au bloc diagonal, MEME matrice de demonstration),
! MEME resultat numerique, pour permettre une comparaison A/B DIRECTE avec
! l'executable existant (example_blr_mpi_pivot, JAMAIS modifie) -- voir
! mpi_blr_matrix_mod.f90, type MPI_BLR_MATRIX_OOC_T, et blr_ooc_mod.f90.
!
! PRINCIPE (demande explicite : ne plus etre limite par la RAM disponible,
! ni pour la matrice, ni pour les facteurs LU apres remplissage de Schur,
! ni pour un second membre massif) -- 4 ETAPES SEPAREES, chacune bornee en
! memoire par un parametre utilisateur :
!   1) assemblage de A       -> ecrit sur un fichier scratch PAR PROCESSUS
!      (jamais partage), borne par MAX_MB_PER_PROCESS (cache LRU).
!   2) assemblage du second membre B -> ECRIT SUR DISQUE LUI AUSSI, par
!      paquets d'au plus NRHS_BATCH colonnes -- jamais tout B en RAM a la
!      fois, meme etape separee de la resolution.
!   3) factorisation LU -> les blocs "finis" (panneau/diagonale de l'etape
!      K) repartent sur disque des qu'ils ne servent plus a la mise a jour
!      de Schur ; la sous-matrice traineante reste bornee par
!      MAX_MB_PER_PROCESS via le meme cache LRU qu'a l'assemblage.
!   4) resolution -> relit L/U depuis le disque, par paquets d'au plus
!      NRHS_BATCH seconds membres (CONSUME_X_BATCH_CB rassemble chaque
!      paquet resolu dans X_GLOBAL ICI, a des fins de verification -- un
!      vrai cas massivement second membre ecrirait plutot chaque paquet
!      directement sur disque au lieu de tout accumuler).
!
! Base sur l'algorithme CLASSIQUE (jamais l'arena -- le cycle allocate/
! deallocate individuel de BLR_BLOCK_T se prete naturellement a l'eviction
! disque, voir la memoire du projet).
!
! CHOISIR MAX_MB_PER_PROCESS (voir le commentaire detaille sur MPI_BLR_
! MATRIX_OOC_T dans mpi_blr_matrix_mod.f90, mesures a l'appui) : NE PAS le
! prendre "le plus petit possible" pour economiser de la RAM -- un budget
! sous l'ensemble de travail reel du processus ne degrade pas gracieusement,
! il fait THRASHER le cache (mesure : 68% du temps de factorisation perdu
! en acces cache/disque a un budget trop juste, /13.5 en portant le budget
! au-dessus de l'ensemble de travail reel, pour une RAM quasi inchangee).
! Choisir large par defaut, et ne descendre que si la RAM disponible
! l'impose vraiment.
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot_ooc PR PC [N] [NB] [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] [EPS] [RATIO] [USE_RSVD] [CHECK_ERROR] [CHECK_RESIDUAL]
!         (defauts : N=4000 NB=200 NRHS=4 TOL=1e-4 MAX_MB_PER_PROCESS=64 NRHS_BATCH=2 SCRATCH_DIR=./blr_ooc_scratch EPS=0.0 RATIO=1/3 USE_RSVD=0 CHECK_ERROR=0 CHECK_RESIDUAL=0)
!
! EPS (BRUIT PLEIN-RANG, POUR EXPERIMENTER FACILEMENT AVEC LA COMPRESSIBI-
! LITE -- IDENTIQUE a example_blr_mpi_pivot.f90, voir son commentaire de
! tete pour le detail) : EPS=0 (defaut) = noyau tres compressible ; EPS
! croissant degrade progressivement la compressibilite, de plus en plus de
! blocs tombant en repli dense a TOL fixe. Passer la MEME valeur d'EPS aux
! deux executables pour comparer classique/OOC sur EXACTEMENT la meme
! matrice.
!
! RATIO (COMPRESS_WORTHWHILE_RATIO -- IDENTIQUE a example_blr_mpi_pivot.f90,
! voir son commentaire de tete pour le detail) : plafond de rang au-dela
! duquel un bloc est repli en DENSE meme si sa SVD trouve un rang valide
! sous TOL. Passer la MEME valeur aux deux executables pour comparer.
!
! USE_RSVD (methode de compression a l'assemblage -- IDENTIQUE a example_
! blr_mpi_pivot.f90, voir son commentaire de tete pour le detail) : 0
! (defaut) = SVD complete (CGESVD) ; 1 = SVD randomisee adaptative.
! Passer la MEME valeur aux deux executables pour comparer.
!
! CHECK_ERROR (verification numerique -- IDENTIQUE a example_blr_mpi_
! pivot.f90, voir son commentaire de tete pour le detail complet) : 0
! (defaut) = comportement INCHANGE (B arbitraire). 1 = B=A*Xexact par
! produit BRUT (JAMAIS via une deuxieme matrice assemblee ou chargee en
! RAM -- l'avantage memoire de l'OOC reste ENTIER, voir plus bas), erreur
! affichee en fin de resolution. Barre de progression + RSS max dedies
! pendant le calcul de B (etape potentiellement longue a grand N, cout
! CALCUL O(N^2*NRHS) mais ZERO cout memoire supplementaire).
!
! CHECK_RESIDUAL (verification SANS solution exacte connue -- INDEPENDANT
! de CHECK_ERROR) : calcule ||A*X-B||/||B|| APRES resolution avec le B
! REELLEMENT utilise -- le cas d'usage GENERAL (A et B fournis par
! l'utilisateur, solution inconnue).
!
! Pour eviter de rappeler GET_BLOCK/FILL_RHS une deuxieme fois (couteux si
! ces callbacks le sont vraiment, cas d'usage reel), A (dense, telle que
! renvoyee par GET_BLOCK) et B (telle qu'ecrite par FILL_RHS, donc deja
! coherente avec CHECK_ERROR) sont mises en cache sur disque UNE SEULE
! FOIS, DES QUE POSSIBLE (juste apres assemblage / apres assemblage du
! RHS, donc AVANT que la factorisation n'ecrase A par L/U en place) --
! acheck_rank_NNNNNN.dat et bcheck_rank_NNNNNN.dat, meme convention
! 1-fichier-par-processus, meme decoupage par blocs de lignes que la
! verification elle-meme (borne par MAX_MB_PER_PROCESS) -- voir
! example_user_matrix_ooc.f90 pour le detail complet de la technique
! (identique ici). Des qu'un bloc de B a servi au calcul du residu, il
! n'est plus utile : X (deja disponible dans X_GLOBAL) est ECRIT A SA
! PLACE, DANS LE MEME FICHIER (reecriture en place par POS=, comme
! OOC_WRITE_BLOCK_FIXED dans blr_ooc_mod.f90) -- a la fin,
! acheck_rank_*.dat est supprime, bcheck_rank_*.dat survit et contient
! desormais X (copie disque de la solution, en plus de X_GLOBAL en RAM).
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_PIVOT_OOC
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, NRHS_BATCH, IERR, MYID0
  REAL(WP) :: TOL, MAX_MB_PER_PROCESS, EPS, RATIO
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=256) :: SCRATCH_DIR
  TYPE(MPI_BLR_MATRIX_OOC_T) :: A_BLR
  TYPE(OOC_RHS_T) :: RHS
  COMPLEX(WP), ALLOCATABLE :: X_GLOBAL(:,:)
  INTEGER :: INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR
  INTEGER :: USE_RSVD_INT, CHECK_ERROR_INT, CHECK_RESIDUAL_INT
  LOGICAL :: USE_RSVD_FLAG, CHECK_ERROR_FLAG, CHECK_RESIDUAL_FLAG
  DOUBLE PRECISION :: RSS_MAX_RHS_MB, RSS_MAX_RESID_MB
  REAL(WP) :: ERR, RESID

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP
  MAX_MB_PER_PROCESS = 64.0_WP; NRHS_BATCH = 2; SCRATCH_DIR = './blr_ooc_scratch'; EPS = 0.0_WP
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 0; CHECK_ERROR_INT = 0; CHECK_RESIDUAL_INT = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./example_blr_mpi_pivot_ooc PR PC [N] [NB] [NRHS] [TOL] '// &
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
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN
    CALL GET_COMMAND_ARGUMENT(10, ARG); READ(ARG,*) EPS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN
    CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) RATIO
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 12) THEN
    CALL GET_COMMAND_ARGUMENT(12, ARG); READ(ARG,*) USE_RSVD_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 13) THEN
    CALL GET_COMMAND_ARGUMENT(13, ARG); READ(ARG,*) CHECK_ERROR_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 14) THEN
    CALL GET_COMMAND_ARGUMENT(14, ARG); READ(ARG,*) CHECK_RESIDUAL_INT
  END IF
  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
  CHECK_ERROR_FLAG = (CHECK_ERROR_INT .NE. 0)
  CHECK_RESIDUAL_FLAG = (CHECK_RESIDUAL_INT .NE. 0)

  ! Voir example_blr_mpi_pivot.f90 pour le detail (SET_COMPRESS_WORTHWHILE_
  ! RATIO, blr_compress_mod.f90).
  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F8.1,A,I0,A,A,A,ES9.2,A,F6.3)') &
         ' [OOC] Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
         ' , MAX_MB_PER_PROCESS=', MAX_MB_PER_PROCESS, ' , NRHS_BATCH=', NRHS_BATCH, &
         ' , SCRATCH_DIR=', TRIM(SCRATCH_DIR), ' , EPS(bruit)=', EPS, ' , RATIO=', RATIO
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
  END IF

  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  ! --- 1) assemblage BLR distribue, out-of-core -> disque ----------------
  CALL MPI_ASSEMBLE_BLR_MATRIX_OOC(GRID, N, NB, TOL, GET_BLOCK, A_BLR, SCRATCH_DIR, MAX_MB_PER_PROCESS, &
                                    T_ASM, N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., &
                                    USE_RSVD=USE_RSVD_FLAG)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    ! Taux de compression REEL, base sur le RANG MOYEN -- voir
    ! REAL_COMPRESSION_PCT ci-dessous (meme correctif que example_blr_mpi_
    ! pivot.f90 : un bloc en repli dense compte NB dans RANK_SUM, PAS 2*NB,
    ! la formule naive donnait a tort 200% quand tout est dense).
    WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (assemblage) : ', REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
  END IF
  CALL REPORT_MEMORY_MB(GRID, 'RAM residente (cache, apres assemblage)', MPI_BLR_MATRIX_OOC_RESIDENT_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)

  ! Cache de A pour CHECK_RESIDUAL : AVANT la factorisation, qui ecrase A
  ! par L/U en place -- GET_BLOCK n'est plus jamais rappelee ensuite.
  IF (CHECK_RESIDUAL_FLAG) CALL CACHE_A_FOR_CHECK(GRID, N, NB)

  ! --- 2) assemblage du second membre, SEPAREMENT, aussi sur disque ------
  RSS_MAX_RHS_MB = 0.0D0
  CALL MPI_ASSEMBLE_RHS_OOC(GRID, A_BLR, NRHS, NRHS_BATCH, FILL_RHS, SCRATCH_DIR, RHS)
  IF (CHECK_ERROR_FLAG .AND. MYID0.EQ.0) &
       WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant B=A*Xexact = ', RSS_MAX_RHS_MB, ' Mo'
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage B)', GET_RSS_KB()*1024_8)

  ! Cache de B pour CHECK_RESIDUAL -- deja coherente avec CHECK_ERROR
  ! puisque FILL_RHS (rappelee ici) applique la meme logique qu'a
  ! l'assemblage reel du second membre.
  IF (CHECK_RESIDUAL_FLAG) CALL CACHE_B_FOR_CHECK(GRID, N, NB, NRHS)

  ! --- 3) factorisation LU, AVEC pivot, out-of-core -----------------------
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC(GRID, A_BLR, TOL, INFO, &
                                           N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                           SHOW_PROGRESS=.TRUE., SHOW_RSS_TRACE=.TRUE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC, INFO=', INFO
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
  CALL REPORT_MEMORY_MB(GRID, 'RAM residente (cache, apres factorisation)', MPI_BLR_MATRIX_OOC_RESIDENT_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)

  ! --- 4) resolution, out-of-core, par paquets de second membre -----------
  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D(GRID, A_BLR, RHS, CONSUME_X_BATCH, SHOW_PROGRESS=.TRUE.)
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres resolution)', GET_RSS_KB()*1024_8)

  ! PAS de verification independante ICI (voir test/test_mpi_blr_matrix_
  ! pivot_ooc.f90 pour la validation rigoureuse, X_exact connu, quel que
  ! soit NRHS_BATCH/MAX_MB_PER_PROCESS -- meme philosophie que les autres
  ! exemples canoniques : demonstration de l'enchainement d'appels, pas un
  ! test de non-regression). NOTE : sur CETTE matrice de demonstration
  ! (chaotique, mal conditionnee), la valeur numerique de X(:,1) peut
  ! varier legerement avec NRHS_BATCH (simple precision + formes d'appels
  ! BLAS differentes amplifiees par le mauvais conditionnement) -- ce n'est
  ! pas un signe d'erreur, voir la memoire du projet.
  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
    BLOCK
      INTEGER :: I
      DO I = 1, MIN(5, N)
        WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
      END DO
    END BLOCK
    IF (CHECK_ERROR_FLAG) THEN
      ERR = COMPUTE_ERR(X_GLOBAL)
      WRITE(6,'(A,ES12.4)') ' Erreur relative ||X-Xexact|| / ||Xexact|| = ', ERR
    END IF
  END IF

  ! APPEL COLLECTIF (MPI_BCAST/MPI_ALLREDUCE a l'interieur) -- TOUS les
  ! rangs doivent l'appeler, pas seulement le rang 0. Consomme et supprime
  ! acheck_rank_*.dat au passage ; bcheck_rank_*.dat survit seul,
  ! contenant desormais X (voir COMPUTE_RESIDUAL).
  IF (CHECK_RESIDUAL_FLAG) THEN
    RSS_MAX_RESID_MB = 0.0D0
    RESID = COMPUTE_RESIDUAL(X_GLOBAL)
    IF (MYID0.EQ.0) THEN
      WRITE(6,'(A,F8.1,A)') ' [RSS trace] MAX observe pendant le residu = ', RSS_MAX_RESID_MB, ' Mo'
      WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
      WRITE(6,'(A)') ' acheck_rank_*.dat supprime ; bcheck_rank_*.dat contient desormais X '// &
           '(economie d''espace disque).'
    END IF
  END IF

  DEALLOCATE(X_GLOBAL)
  CALL OOC_RHS_FREE(RHS)
  CALL MPI_BLR_MATRIX_OOC_FREE(A_BLR)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Voir example_blr_mpi_pivot.f90 pour le detail du correctif.
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  ! Hash deterministe (I,J) -> [-1,1) -- meme recette que example_blr_mpi_
  ! pivot.f90/test/bench_blr_mpi_pivot_compressibility.f90.
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  ! Matrice IDENTIQUE a example_blr_mpi_pivot.f90 (comparaison A/B directe,
  ! y compris le bruit EPS -- passer la MEME valeur d'EPS aux deux
  ! executables). Anti-diagonale CONFINEE au bloc diagonal de GI (voir le
  ! commentaire detaille dans example_blr_mpi_pivot.f90/LOWRANK_A pour la
  ! raison du correctif -- meme correction appliquee ICI a l'identique).
  ! Extrait en fonction independante (au lieu d'etre inline dans GET_BLOCK)
  ! pour etre reutilisable par FILL_RHS en mode CHECK_ERROR.
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

  ! Solution EXACTE deterministe -- IDENTIQUE a example_blr_mpi_pivot.f90,
  ! utilisee UNIQUEMENT si CHECK_ERROR=1.
  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  ! Erreur relative ||X-Xexact||/||Xexact|| -- RANG 0 UNIQUEMENT (X_GLOBAL
  ! y est deja rassemble par CONSUME_X_BATCH, aucune allocation supple-
  ! mentaire de portee globale : l'avantage memoire de l'OOC reste entier).
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

  ! Valeur de B a l'entree (GI,R) -- SOURCE UNIQUE partagee par FILL_RHS
  ! ET COMPUTE_RESIDUAL, voir example_blr_mpi_pivot.f90/B_ENTRY pour le
  ! detail (meme raisonnement, identique ici).
  COMPLEX(WP) FUNCTION B_ENTRY(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    INTEGER :: JJ
    COMPLEX(WP) :: ACC
    COMPLEX(WP), PARAMETER :: CZERO = (0.0_WP,0.0_WP)
    IF (CHECK_ERROR_FLAG) THEN
      ACC = CZERO
      DO JJ = 1, N
        ACC = ACC + LOWRANK_A(GI,JJ) * XEX_FORMULA(JJ,R)
      END DO
      B_ENTRY = ACC
    ELSE
      B_ENTRY = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP), KIND=WP)
    END IF
  END FUNCTION B_ENTRY

  ! Ecrit UNE FOIS le cache dense de A (acheck_rank_NNNNNN.dat), par blocs
  ! de lignes bornes par MAX_MB_PER_PROCESS -- meme decoupage que
  ! COMPUTE_RESIDUAL, pour que la relecture s'y cale exactement. Appelee
  ! AVANT la factorisation (qui ecrase A par L/U en place). Voir
  ! example_user_matrix_ooc.f90/CACHE_A_FOR_CHECK pour le detail complet
  ! (technique identique).
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

  ! Ecrit UNE FOIS le cache de B (bcheck_rank_NNNNNN.dat), par bloc de
  ! lignes -- meme decoupage que COMPUTE_RESIDUAL. B est ici deja LARGE
  ! de NRHS colonnes en une seule passe (pas de sous-paquets : X_GLOBAL
  ! lui-meme n'est jamais borne dans ce fichier de demonstration, voir le
  ! commentaire de COMPUTE_RESIDUAL). Rappelle FILL_RHS -- deja coherente
  ! avec CHECK_ERROR puisque c'est la MEME fonction qui a servi a
  ! l'assemblage reel du second membre.
  SUBROUTINE CACHE_B_FOR_CHECK(GD, NN, NNB, NR)
    TYPE(GRID_T), INTENT(IN) :: GD
    INTEGER, INTENT(IN)      :: NN, NNB, NR
    COMPLEX(WP), ALLOCATABLE :: BBLOCK(:,:)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK, CR
    INTEGER :: UNIT_B
    CHARACTER(LEN=16) :: RSTR

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)
    CR = MAX(1, MIN(NNB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(NN,WP)))))

    WRITE(RSTR,'(I6.6)') GD%MYID
    OPEN(NEWUNIT=UNIT_B, FILE=TRIM(SCRATCH_DIR)//'/bcheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE')

    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CR
      REND_CHUNK = MIN(RSTART_CHUNK+CR-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(BBLOCK(NCHUNK,NR))
      CALL FILL_RHS(RSTART_CHUNK, REND_CHUNK, 1, NR, BBLOCK)
      WRITE(UNIT_B) RSTART_CHUNK, REND_CHUNK, NCHUNK
      WRITE(UNIT_B) BBLOCK(1:NCHUNK,1:NR)
      DEALLOCATE(BBLOCK)
    END DO
    CLOSE(UNIT_B)
  END SUBROUTINE CACHE_B_FOR_CHECK

  ! Residu ||A*X-B||/||B|| -- CAS GENERAL, voir example_blr_mpi_pivot.f90/
  ! COMPUTE_RESIDUAL pour le detail COMPLET du correctif de plafonnage
  ! memoire (identique ici) : BLAS (CGEMM par blocs de lignes) +
  ! DISTRIBUTION MPI sur TOUS les rangs (X diffuse via MPI_BCAST, normes
  ! partielles reduites via MPI_ALLREDUCE) -- APPEL COLLECTIF (voir le
  ! site d'appel).
  !
  ! A et B sont RELUES depuis le cache disque (acheck_rank_*.dat,
  ! bcheck_rank_*.dat -- ecrit par CACHE_A_FOR_CHECK/CACHE_B_FOR_CHECK,
  ! AVANT la factorisation) au lieu de rappeler GET_BLOCK/FILL_RHS/
  ! B_ENTRY -- B etant deja coherente avec CHECK_ERROR (meme FILL_RHS),
  ! XCOMBINED/NCOLS_MULT ne sont plus necessaires : la comparaison se fait
  ! directement contre le B cache, que CHECK_ERROR soit actif ou non.
  ! Des qu'un bloc de B a servi, X (deja disponible dans XG) l'ecrase EN
  ! PLACE dans bcheck_rank_*.dat (POS=, meme forme) -- economie d'espace
  ! disque, voir le commentaire de tete. acheck_rank_*.dat est supprime en
  ! fin de fonction.
  !
  ! LIMITE RESTANTE, HONNETE (prexistante, DOCUMENTEE) : X_GLOBAL(N,NRHS)
  ! reste un tableau COMPLET non borne par MAX_MB_PER_PROCESS -- MEME
  ! limite deja assumee par le commentaire de tete de ce fichier pour
  ! l'affichage des 5 premieres valeurs de X ("un vrai cas massivement
  ! second membre ecrirait plutot chaque paquet directement sur disque").
  ! Cette verification (comme l'affichage) est un outil de DEMONSTRATION,
  ! pas le coeur de l'algorithme OOC (assemblage/factorisation/resolution,
  ! ceux-la restent strictement bornes) -- rester honnete sur cette
  ! distinction plutot que de le cacher.
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG)
    COMPLEX(WP), INTENT(INOUT) :: XG(N,NRHS)
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), AXBLOCK(:,:), BBLOCK(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, NROWS_MINE, NTOTAL_LOCAL
    INTEGER :: RSTART_CHUNK, REND_CHUNK, NCHUNK, GI, R, IERR2, NDONE
    INTEGER :: CHUNK_ROWS, UNIT_A, UNIT_B
    INTEGER :: HA_RSTART, HA_REND, HA_NCHUNK, HB_RSTART, HB_REND, HB_NCHUNK
    INTEGER(8) :: PAYLOAD_POS
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2
    CHARACTER(LEN=16) :: RSTR

    CALL MPI_BCAST(XG, N*NRHS, MPI_COMPLEX, 0, GRID%COMM, IERR2)

    ROWS_PER_PROC = (N + GRID%NPROCS - 1) / GRID%NPROCS
    RSTART_MINE = GRID%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GRID%MYID+1)*ROWS_PER_PROC, N)
    NROWS_MINE = MAX(REND_MINE - RSTART_MINE + 1, 0)
    NTOTAL_LOCAL = MAX(NROWS_MINE, 1)

    ! Taille de bloc de lignes bornee par MAX_MB_PER_PROCESS -- MEME
    ! decoupage que CACHE_A_FOR_CHECK/CACHE_B_FOR_CHECK (relecture calee
    ! exactement dessus).
    CHUNK_ROWS = MAX(1, MIN(NB, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(N,WP)))))

    WRITE(RSTR,'(I6.6)') GRID%MYID
    OPEN(NEWUNIT=UNIT_A, FILE=TRIM(SCRATCH_DIR)//'/acheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')
    OPEN(NEWUNIT=UNIT_B, FILE=TRIM(SCRATCH_DIR)//'/bcheck_rank_'//TRIM(RSTR)//'.dat', &
         ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READWRITE')

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP
    RSS_MAX_RESID_MB = REAL(GET_RSS_KB(),8)/1024.0D0
    NDONE = 0
    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, CHUNK_ROWS
      REND_CHUNK = MIN(RSTART_CHUNK+CHUNK_ROWS-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1

      READ(UNIT_A) HA_RSTART, HA_REND, HA_NCHUNK
      ALLOCATE(ABLOCK(NCHUNK,N))
      READ(UNIT_A) ABLOCK(1:NCHUNK,1:N)

      READ(UNIT_B) HB_RSTART, HB_REND, HB_NCHUNK
      INQUIRE(UNIT=UNIT_B, POS=PAYLOAD_POS)
      ALLOCATE(BBLOCK(NCHUNK,NRHS))
      READ(UNIT_B) BBLOCK(1:NCHUNK,1:NRHS)

      ALLOCATE(AXBLOCK(NCHUNK,NRHS))
      CALL CGEMM('N','N', NCHUNK, NRHS, N, CONE, ABLOCK, NCHUNK, XG, N, CZERO, AXBLOCK, NCHUNK)
      DO R = 1, NRHS
        DO GI = 1, NCHUNK
          NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(AXBLOCK(GI,R)-BBLOCK(GI,R))**2,WP)
          DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BBLOCK(GI,R))**2,WP)
        END DO
      END DO

      ! B vient de servir : plus utile -- X ecrase sa place sur disque.
      WRITE(UNIT_B, POS=PAYLOAD_POS) XG(RSTART_CHUNK:REND_CHUNK, 1:NRHS)

      DEALLOCATE(ABLOCK, AXBLOCK, BBLOCK)
      NDONE = NDONE + NCHUNK
      RSS_MAX_RESID_MB = MAX(RSS_MAX_RESID_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
      IF (GRID%MYID.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul du residu ||A*X-B|| (verification)', NDONE, NTOTAL_LOCAL)
    END DO
    CLOSE(UNIT_B)
    CLOSE(UNIT_A, STATUS='DELETE')

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, RSS_MAX_RESID_MB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

  ! B(GI,GR) : formule arbitraire (defaut, CHECK_ERROR=0, INCHANGEE, O(1)
  ! par entree -- boucle scalaire directe conservee) ou produit BRUT
  ! B=A*Xexact (CHECK_ERROR=1, voir commentaire de tete).
  !
  ! DEMANDE EXPLICITE (2026-09-01) : meme optimisation BLAS que COMPUTE_
  ! RESIDUAL pour le cas CHECK_ERROR=1, ET MEME PLAFONNAGE par
  ! MAX_MB_PER_PROCESS -- [ISTART,IEND] (fourni par MPI_ASSEMBLE_RHS_OOC)
  ! n'est PAS garanti respecter le budget a lui seul (sa largeur depend de
  ! la structure BLR, pas du budget memoire utilisateur), donc le tampon
  ! ABLOCK(CHUNK,N) est encore sous-decoupe si besoin -- EXACTEMENT le
  ! meme raisonnement que le correctif de COMPUTE_RESIDUAL (voir son
  ! commentaire), applique ici pour la meme raison : ne jamais laisser un
  ! tampon BLAS depasser le budget annonce par l'utilisateur.
  SUBROUTINE FILL_RHS(ISTART, IEND, RSTART, REND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, RSTART, REND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, REND-RSTART+1)
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), BBLOCK(:,:), XEXACT_VEC(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: II, RR, NDONE, NCOLS, CHUNK_ROWS, RSTART_CHUNK, REND_CHUNK, NCHUNK
    NDONE = 0
    IF (CHECK_ERROR_FLAG) THEN
      NCOLS = REND - RSTART + 1
      ALLOCATE(XEXACT_VEC(N,NCOLS))
      DO RR = RSTART, REND
        DO II = 1, N
          XEXACT_VEC(II,RR-RSTART+1) = XEX_FORMULA(II,RR)
        END DO
      END DO
      CHUNK_ROWS = MAX(1, MIN(IEND-ISTART+1, INT(MAX_MB_PER_PROCESS*1024.0_WP*1024.0_WP / (8.0_WP*REAL(N,WP)))))
      DO RSTART_CHUNK = ISTART, IEND, CHUNK_ROWS
        REND_CHUNK = MIN(RSTART_CHUNK+CHUNK_ROWS-1, IEND)
        NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
        ALLOCATE(ABLOCK(NCHUNK,N), BBLOCK(NCHUNK,NCOLS))
        CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, N, ABLOCK)
        CALL CGEMM('N','N', NCHUNK, NCOLS, N, CONE, ABLOCK, NCHUNK, XEXACT_VEC, N, CZERO, BBLOCK, NCHUNK)
        BLK(RSTART_CHUNK-ISTART+1:REND_CHUNK-ISTART+1, 1:NCOLS) = BBLOCK
        DEALLOCATE(ABLOCK, BBLOCK)
        NDONE = NDONE + NCHUNK
        RSS_MAX_RHS_MB = MAX(RSS_MAX_RHS_MB, REAL(GET_RSS_KB(),8)/1024.0D0)
        IF (MYID0.EQ.0) CALL PRINT_PROGRESS_BAR('Calcul de B=A*Xexact (verification)', NDONE, IEND-ISTART+1)
      END DO
      DEALLOCATE(XEXACT_VEC)
    ELSE
      DO RR = RSTART, REND
        DO II = ISTART, IEND
          BLK(II-ISTART+1, RR-RSTART+1) = B_ENTRY(II, RR)
        END DO
      END DO
    END IF
  END SUBROUTINE FILL_RHS

  ! Rassemble CE paquet resolu dans X_GLOBAL (verification/demonstration
  ! uniquement -- un vrai cas massivement second membre ecrirait plutot
  ! CE paquet directement sur disque ici, voir commentaire de tete).
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

END PROGRAM EXAMPLE_BLR_MPI_PIVOT_OOC
