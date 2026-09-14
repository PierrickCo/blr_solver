!===============================================================================
! PROTOTYPE SEQUENTIEL -- valide (ou infirme) l'idee d'une permutation des
! blocs BLR guidee par le rang mesure a l'assemblage, sur la matrice
! deliberement heterogene de bench_blr_rank_heterogeneity.f90 (voir son
! commentaire de tete pour la construction : chaque bloc appartient a un
! groupe CACHE assigne par hash SANS RAPPORT avec l'ordre naturel -- blocs du
! meme groupe = "durs" (quasi rang plein), blocs de groupes differents =
! "faciles" (rang EXACTEMENT 1)).
!
! VOLONTAIREMENT SEQUENTIEL (grille MPI 1x1 imposee) : la distribution 2D
! bloc-cyclique attache la propriete d'un bloc a son INDICE (I mod PR, J mod
! PC) -- permuter les indices en MPI distribue demanderait de redistribuer
! des blocs entre processus (un vrai chantier d'ingenierie, voir la
! discussion). En 1x1, tous les blocs sont deja locaux au processus unique :
! permuter revient a un simple reindexage, sans communication -- feu vert
! pour valider l'IDEE avant d'investir dans la version distribuee.
!
! DEROULEMENT :
!   1) Assemble la matrice a l'ordre NATUREL (GET_BLOCK_NATURAL), factorise,
!      mesure (temps, memoire des facteurs, taux de compression).
!   2) Calcule une permutation par un algorithme GLOUTON simple ("chainage
!      plus proche voisin") : place le bloc 1, puis a chaque etape ajoute le
!      bloc NON PLACE dont le rang mesure avec le DERNIER bloc place est le
!      PLUS ELEVE (rang eleve = couplage fort = a garder proche). L'algorithme
!      n'utilise QUE RANK_MAT (mesure a l'etape 1) -- ne regarde JAMAIS le
!      groupe cache reel (CLUSTER_OF, prive a la generation de la matrice) :
!      c'est un algorithme "boite noire", comme il le faudrait sur un
!      probleme reel ou le groupe cache n'existe pas.
!   3) Reassemble la MEME matrice logique mais via GET_BLOCK_PERMUTED (chaque
!      nouveau bloc K lit les donnees de l'ancien bloc ACTIVE_PERM(K) -- meme
!      contenu numerique, ordre different), factorise, mesure.
!   4) MEME chose avec une DEUXIEME strategie : clustering HIERARCHIQUE par
!      bipartition recursive (voir HIER_CLUSTER) -- contrairement au
!      chainage glouton (LOCAL : ne regarde que le dernier bloc place), la
!      bipartition est GLOBALE a chaque niveau (separe tout le groupe
!      courant en 2 sous-groupes de couplage interne maximal), un
!      algorithme plus proche de la nested dissection (mais pondere par le
!      rang mesure, pas par la structure creuse).
!   5) Compare les trois mesures (NATUREL / CHAINAGE / HIERARCHIQUE).
!
! Usage : mpirun -np 1 ./bench_blr_rank_permutation [N] [NB] [TOL] [NCLUSTERS]
!         (defauts : N=4000 NB=200 TOL=1e-4 NCLUSTERS=4 -- IMPOSE PR=PC=1)
!===============================================================================
PROGRAM BENCH_BLR_RANK_PERMUTATION
  USE MPI_BLR_MATRIX_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: N, NB, IERR, MYID0, NCLUSTERS, NBLK, NPROCS_WORLD
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_NAT, A_PERM, A_HIER
  INTEGER :: N_COMPRESSED, N_OFFDIAG, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED2, N_OFFDIAG2, INFO2, N_SWAPS2
  INTEGER :: N_COMPRESSED3, N_OFFDIAG3, INFO3, N_SWAPS3
  INTEGER(8) :: RANK_SUM, RANK_SUM2, RANK_SUM3
  ! MPI_ASSEMBLE_BLR_MATRIX puis MPI_BLR_MATRIX_LU_FACTOR_PIVOT ecrivent
  ! TOUS LES DEUX dans N_COMPRESSED/N_OFFDIAG/RANK_SUM (INTENT(OUT)) --
  ! sauvegarde une copie des valeurs d'ASSEMBLAGE juste apres chaque appel
  ! d'assemblage, AVANT que la factorisation ne les ecrase, pour pouvoir
  ! comparer avant/apres a la fin.
  INTEGER :: NC_ASM, NO_ASM, NC_ASM2, NO_ASM2, NC_ASM3, NO_ASM3
  INTEGER(8) :: RS_ASM, RS_ASM2, RS_ASM3
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_ASM2, T_FACTOR2, T_ASM3, T_FACTOR3
  INTEGER, ALLOCATABLE :: PERM(:), PERM2(:), RANK_MAT(:,:), ACTIVE_PERM(:)
  LOGICAL, ALLOCATABLE :: PLACED(:)
  INTEGER :: I, J, K, BEST, BESTVAL, CUR, NMATCH, OPOS

  N = 4000; NB = 200; TOL = 1.0E-4_WP; NCLUSTERS = 4

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS_WORLD, IERR)
  IF (NPROCS_WORLD .NE. 1) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : ce prototype est volontairement SEQUENTIEL -- lancer avec -np 1.'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NCLUSTERS
  END IF
  IF (MOD(N,NB) .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : ce prototype suppose N multiple de NB (blocs de taille uniforme).'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  NBLK = N / NB

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, 1, 1, IERR)
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,ES9.2,A,I0,A,I0)') &
       ' N=', N, ' , NB=', NB, ' , TOL=', TOL, ' , NCLUSTERS=', NCLUSTERS, ' , NBLK=', NBLK

  ! === 1) Assemblage+factorisation a l'ordre NATUREL (baseline) ===========
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK_NATURAL, A_NAT, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.FALSE.)
  NC_ASM = N_COMPRESSED; NO_ASM = N_OFFDIAG; RS_ASM = RANK_SUM
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' [NATUREL]  Assemblage : ', T_ASM, &
       ' s , compresses ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', &
       100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' % blocs, ', &
       REAL_COMPRESSION_PCT(N_COMPRESSED,N_OFFDIAG,RANK_SUM), ' % stockage reel)'

  ! Extrait la matrice de rangs MESUREE (1 seul processus -- tous les blocs
  ! sont locaux, A_NAT%NBLKR_LOC=A_NAT%NBLKC_LOC=NBLK, indices GLOBAUX=LOCAUX).
  ALLOCATE(RANK_MAT(NBLK,NBLK))
  DO I = 1, NBLK
    DO J = 1, NBLK
      IF (I .EQ. J) THEN
        RANK_MAT(I,J) = 0
      ELSE IF (A_NAT%OFFBLK(I,J)%IS_LOWRANK) THEN
        RANK_MAT(I,J) = A_NAT%OFFBLK(I,J)%RANK
      ELSE
        RANK_MAT(I,J) = NB
      END IF
    END DO
  END DO

  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_NAT, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                       N_SWAPS, T_FACTOR, SHOW_PROGRESS=.FALSE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR factorisation (naturel), INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' [NATUREL]  Factorisation : ', T_FACTOR, ' s , compresses ', &
         N_COMPRESSED, ' / ', N_OFFDIAG, ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), &
         ' % blocs, ', REAL_COMPRESSION_PCT(N_COMPRESSED,N_OFFDIAG,RANK_SUM), ' % stockage reel)'
    WRITE(6,'(A,F10.3,A)') ' [NATUREL]  Memoire facteurs LU : ', &
         REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_NAT),WP)/1048576.0_WP, ' Mo'
  END IF

  ! === 2) Permutation GLOUTONNE a partir de RANK_MAT (boite noire -- ne
  !     regarde JAMAIS CLUSTER_OF) : "chainage plus proche voisin" -- place
  !     le bloc 1, puis a chaque etape le bloc non place dont le rang MESURE
  !     avec le DERNIER bloc place est le plus eleve (couplage fort = a
  !     garder proche dans l'ordre d'elimination). =================
  ALLOCATE(PERM(NBLK), PLACED(NBLK))
  PLACED = .FALSE.
  PERM(1) = 1; PLACED(1) = .TRUE.
  DO K = 2, NBLK
    CUR = PERM(K-1)
    BEST = 0; BESTVAL = -1
    DO J = 1, NBLK
      IF (.NOT. PLACED(J)) THEN
        IF (RANK_MAT(CUR,J) .GT. BESTVAL) THEN
          BESTVAL = RANK_MAT(CUR,J); BEST = J
        END IF
      END IF
    END DO
    PERM(K) = BEST; PLACED(BEST) = .TRUE.
  END DO

  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A)') ' Permutation trouvee (nouvel indice -> ancien indice) :'
    WRITE(6,'(A,20(I3,1X))') '   ', (PERM(K), K=1,NBLK)
    ! Validation A POSTERIORI seulement (n'entre PAS dans l'algorithme
    ! ci-dessus) : combien de paires ADJACENTES dans la permutation trouvee
    ! partagent le MEME groupe cache -- si l'algo a bien fonctionne, ce
    ! nombre doit etre nettement plus eleve qu'a l'ordre naturel.
    NMATCH = 0
    DO K = 1, NBLK-1
      IF (CLUSTER_OF(PERM(K)) .EQ. CLUSTER_OF(PERM(K+1))) NMATCH = NMATCH + 1
    END DO
    WRITE(6,'(A,I0,A,I0)') ' Paires adjacentes de MEME groupe cache apres permutation : ', NMATCH, ' / ', NBLK-1
    NMATCH = 0
    DO K = 1, NBLK-1
      IF (CLUSTER_OF(K) .EQ. CLUSTER_OF(K+1)) NMATCH = NMATCH + 1
    END DO
    WRITE(6,'(A,I0,A,I0)') ' Paires adjacentes de MEME groupe cache a l''ordre NATUREL   : ', NMATCH, ' / ', NBLK-1
  END IF

  ! === 3) Reassemble la MEME matrice logique via GET_BLOCK_PERMUTED (chainage
  !     glouton), refactorise, mesure =====================================
  ALLOCATE(ACTIVE_PERM(NBLK))
  ACTIVE_PERM = PERM
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK_PERMUTED, A_PERM, T_ASM2, &
                                N_COMPRESSED2, N_OFFDIAG2, RANK_SUM2, SHOW_PROGRESS=.FALSE.)
  NC_ASM2 = N_COMPRESSED2; NO_ASM2 = N_OFFDIAG2; RS_ASM2 = RANK_SUM2
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' [CHAINAGE] Assemblage : ', T_ASM2, &
       ' s , compresses ', N_COMPRESSED2, ' / ', N_OFFDIAG2, ' (', &
       100.0_WP*REAL(N_COMPRESSED2,WP)/REAL(MAX(N_OFFDIAG2,1),WP), ' % blocs, ', &
       REAL_COMPRESSION_PCT(N_COMPRESSED2,N_OFFDIAG2,RANK_SUM2), ' % stockage reel)'

  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_PERM, TOL, INFO2, N_COMPRESSED2, N_OFFDIAG2, RANK_SUM2, &
                                       N_SWAPS2, T_FACTOR2, SHOW_PROGRESS=.FALSE.)
  IF (INFO2 .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR factorisation (chainage), INFO=', INFO2
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' [CHAINAGE] Factorisation : ', T_FACTOR2, ' s , compresses ', &
         N_COMPRESSED2, ' / ', N_OFFDIAG2, ' (', 100.0_WP*REAL(N_COMPRESSED2,WP)/REAL(MAX(N_OFFDIAG2,1),WP), &
         ' % blocs, ', REAL_COMPRESSION_PCT(N_COMPRESSED2,N_OFFDIAG2,RANK_SUM2), ' % stockage reel)'
    WRITE(6,'(A,F10.3,A)') ' [CHAINAGE] Memoire facteurs LU : ', &
         REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_PERM),WP)/1048576.0_WP, ' Mo'
  END IF

  ! === 4) Clustering HIERARCHIQUE (bipartition recursive) a partir de la
  !     MEME RANK_MAT -- toujours en boite noire, reassemble+refactorise,
  !     mesure ==============================================================
  ALLOCATE(PERM2(NBLK))
  OPOS = 1
  CALL HIER_CLUSTER((/(K, K=1,NBLK)/), NBLK, PERM2)

  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A)') ' Ordre hierarchique trouve (nouvel indice -> ancien indice) :'
    WRITE(6,'(A,20(I3,1X))') '   ', (PERM2(K), K=1,NBLK)
    NMATCH = 0
    DO K = 1, NBLK-1
      IF (CLUSTER_OF(PERM2(K)) .EQ. CLUSTER_OF(PERM2(K+1))) NMATCH = NMATCH + 1
    END DO
    WRITE(6,'(A,I0,A,I0)') ' Paires adjacentes de MEME groupe cache (hierarchique) : ', NMATCH, ' / ', NBLK-1
  END IF

  ACTIVE_PERM = PERM2
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK_PERMUTED, A_HIER, T_ASM3, &
                                N_COMPRESSED3, N_OFFDIAG3, RANK_SUM3, SHOW_PROGRESS=.FALSE.)
  NC_ASM3 = N_COMPRESSED3; NO_ASM3 = N_OFFDIAG3; RS_ASM3 = RANK_SUM3
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' [HIERARCH] Assemblage : ', T_ASM3, &
       ' s , compresses ', N_COMPRESSED3, ' / ', N_OFFDIAG3, ' (', &
       100.0_WP*REAL(N_COMPRESSED3,WP)/REAL(MAX(N_OFFDIAG3,1),WP), ' % blocs, ', &
       REAL_COMPRESSION_PCT(N_COMPRESSED3,N_OFFDIAG3,RANK_SUM3), ' % stockage reel)'

  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_HIER, TOL, INFO3, N_COMPRESSED3, N_OFFDIAG3, RANK_SUM3, &
                                       N_SWAPS3, T_FACTOR3, SHOW_PROGRESS=.FALSE.)
  IF (INFO3 .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR factorisation (hierarchique), INFO=', INFO3
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A,F6.2,A)') ' [HIERARCH] Factorisation : ', T_FACTOR3, ' s , compresses ', &
         N_COMPRESSED3, ' / ', N_OFFDIAG3, ' (', 100.0_WP*REAL(N_COMPRESSED3,WP)/REAL(MAX(N_OFFDIAG3,1),WP), &
         ' % blocs, ', REAL_COMPRESSION_PCT(N_COMPRESSED3,N_OFFDIAG3,RANK_SUM3), ' % stockage reel)'
    WRITE(6,'(A,F10.3,A)') ' [HIERARCH] Memoire facteurs LU : ', &
         REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_HIER),WP)/1048576.0_WP, ' Mo'
  END IF

  ! === 5) Comparaison ======================================================
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '=== Comparaison NATUREL / CHAINAGE / HIERARCHIQUE ==='
    WRITE(6,'(A,F8.3,A)')          ' Temps factorisation NATUREL    : ', T_FACTOR,  ' s'
    WRITE(6,'(A,F8.3,A,F7.2,A)')   ' Temps factorisation CHAINAGE   : ', T_FACTOR2, ' s  (x', T_FACTOR/MAX(T_FACTOR2,1.0D-9), ')'
    WRITE(6,'(A,F8.3,A,F7.2,A)')   ' Temps factorisation HIERARCHIQ : ', T_FACTOR3, ' s  (x', T_FACTOR/MAX(T_FACTOR3,1.0D-9), ')'
    WRITE(6,'(A,F10.3,A,F10.3,A,F10.3,A)') ' Memoire facteurs LU (Mo)       : ', &
         REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_NAT),WP)/1048576.0_WP, '  /  ', &
         REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_PERM),WP)/1048576.0_WP, '  /  ', &
         REAL(MPI_BLR_MATRIX_MEMORY_BYTES(A_HIER),WP)/1048576.0_WP
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '=== Compression AVANT (assemblage, ordre naturel/renumerote) / APRES (fin de factorisation) ==='
    WRITE(6,'(A)') '                    blocs compresses (bas-rang / total hors-diag)     taux reel (2*NB*rang_moyen/NB^2)'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A,A,I0,A,I0,A,F6.1,A,A,F6.2,A,A,F6.2,A)') &
         ' NATUREL     avant : ', NC_ASM, ' / ', NO_ASM, ' (', 100.0_WP*REAL(NC_ASM,WP)/REAL(MAX(NO_ASM,1),WP), '%)', &
         '   apres : ', N_COMPRESSED, ' / ', N_OFFDIAG, ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), '%)', &
         '   ||  reel avant=', REAL_COMPRESSION_PCT(NC_ASM,NO_ASM,RS_ASM), &
         '%  apres=', REAL_COMPRESSION_PCT(N_COMPRESSED,N_OFFDIAG,RANK_SUM), '%'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A,A,I0,A,I0,A,F6.1,A,A,F6.2,A,A,F6.2,A)') &
         ' CHAINAGE    avant : ', NC_ASM2, ' / ', NO_ASM2, ' (', 100.0_WP*REAL(NC_ASM2,WP)/REAL(MAX(NO_ASM2,1),WP), '%)', &
         '   apres : ', N_COMPRESSED2, ' / ', N_OFFDIAG2, ' (', 100.0_WP*REAL(N_COMPRESSED2,WP)/REAL(MAX(N_OFFDIAG2,1),WP), '%)', &
         '   ||  reel avant=', REAL_COMPRESSION_PCT(NC_ASM2,NO_ASM2,RS_ASM2), &
         '%  apres=', REAL_COMPRESSION_PCT(N_COMPRESSED2,N_OFFDIAG2,RANK_SUM2), '%'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A,A,I0,A,I0,A,F6.1,A,A,F6.2,A,A,F6.2,A)') &
         ' HIERARCHIQ  avant : ', NC_ASM3, ' / ', NO_ASM3, ' (', 100.0_WP*REAL(NC_ASM3,WP)/REAL(MAX(NO_ASM3,1),WP), '%)', &
         '   apres : ', N_COMPRESSED3, ' / ', N_OFFDIAG3, ' (', 100.0_WP*REAL(N_COMPRESSED3,WP)/REAL(MAX(N_OFFDIAG3,1),WP), '%)', &
         '   ||  reel avant=', REAL_COMPRESSION_PCT(NC_ASM3,NO_ASM3,RS_ASM3), &
         '%  apres=', REAL_COMPRESSION_PCT(N_COMPRESSED3,N_OFFDIAG3,RANK_SUM3), '%'
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Taux de compression REEL (stockage bas-rang effectif / stockage dense
  ! equivalent) -- MEME formule qu'example_blr_mpi_pivot.f90/README.md (un
  ! bloc en repli dense compte NB dans RSUM, PAS 2*NB, pour ne pas afficher
  ! a tort 200% quand tout est dense).
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  INTEGER FUNCTION CLUSTER_OF(IBLK)
    INTEGER, INTENT(IN) :: IBLK
    REAL(WP) :: S, FRAC
    S = SIN(REAL(IBLK,WP)*37.719_WP)*29143.671_WP
    FRAC = S - REAL(FLOOR(S),WP)
    CLUSTER_OF = MIN(INT(FRAC*REAL(NCLUSTERS,WP)), NCLUSTERS-1)
  END FUNCTION CLUSTER_OF

  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  COMPLEX(WP) FUNCTION U_FUNC(GI)
    INTEGER, INTENT(IN) :: GI
    U_FUNC = CMPLX(SIN(REAL(GI,WP)*0.013_WP)+2.0_WP, COS(REAL(GI,WP)*0.007_WP), KIND=WP)
  END FUNCTION U_FUNC

  COMPLEX(WP) FUNCTION V_FUNC(GJ)
    INTEGER, INTENT(IN) :: GJ
    V_FUNC = CMPLX(COS(REAL(GJ,WP)*0.011_WP)+2.0_WP, SIN(REAL(GJ,WP)*0.005_WP), KIND=WP)
  END FUNCTION V_FUNC

  ! Formule de la matrice heterogene -- identique a bench_blr_rank_
  ! heterogeneity.f90 (voir son commentaire de tete), indexee directement
  ! par les indices GLOBAUX (naturels) II,JJ.
  SUBROUTINE GET_BLOCK_NATURAL(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ, IBLK, JBLK, CI, CJ
    DO JJ = JSTART, JEND
      JBLK = (JJ-1)/NB + 1
      CJ = CLUSTER_OF(JBLK)
      DO II = ISTART, IEND
        IBLK = (II-1)/NB + 1
        CI = CLUSTER_OF(IBLK)
        IF (CI .EQ. CJ) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(1.0_WP + 0.7_WP*HASH_NOISE(II,JJ), &
                                                 0.4_WP*HASH_NOISE(JJ,II), KIND=WP)
        ELSE
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(0.01_WP,0.0_WP,KIND=WP) * U_FUNC(II) * V_FUNC(JJ)
        END IF
        IF (II .EQ. JJ) BLK(II-ISTART+1, JJ-JSTART+1) = BLK(II-ISTART+1, JJ-JSTART+1) + CMPLX(50.0_WP, 0.0_WP, KIND=WP)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK_NATURAL

  ! GET_BLOCK pour la matrice REORDONNEE : le NOUVEAU bloc K (lignes/colonnes
  ! [ISTART,IEND]/[JSTART,JEND] dans la numerotation NOUVELLE) contient EXACTE-
  ! MENT les donnees de l'ANCIEN bloc PERM(K) -- traduit les indices puis
  ! delegue a GET_BLOCK_NATURAL. Suppose des blocs de taille NB UNIFORME
  ! (N multiple de NB, verifie au demarrage) pour que la traduction
  ! bloc->plage de lignes soit une simple arithmetique.
  ! Utilise ACTIVE_PERM (host, positionne par l'appelant AVANT chaque appel
  ! a MPI_ASSEMBLE_BLR_MATRIX -- PERM pour le chainage glouton, PERM2 pour
  ! le clustering hierarchique) -- permet de reutiliser CE MEME GET_BLOCK
  ! pour les deux strategies sans dupliquer la traduction d'indices.
  SUBROUTINE GET_BLOCK_PERMUTED(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: NEWIBLK, NEWJBLK, OLDIBLK, OLDJBLK
    INTEGER :: OISTART, OIEND, OJSTART, OJEND
    NEWIBLK = (ISTART-1)/NB + 1
    NEWJBLK = (JSTART-1)/NB + 1
    OLDIBLK = ACTIVE_PERM(NEWIBLK)
    OLDJBLK = ACTIVE_PERM(NEWJBLK)
    OISTART = (OLDIBLK-1)*NB + 1; OIEND = OLDIBLK*NB
    OJSTART = (OLDJBLK-1)*NB + 1; OJEND = OLDJBLK*NB
    CALL GET_BLOCK_NATURAL(OISTART, OIEND, OJSTART, OJEND, BLK)
  END SUBROUTINE GET_BLOCK_PERMUTED

  ! Clustering HIERARCHIQUE par bipartition recursive, pondere par RANK_MAT
  ! (boite noire -- ne regarde JAMAIS CLUSTER_OF) : separe le groupe courant
  ! GROUP(1:NG) en 2 sous-groupes A/B en placant chaque bloc du cote dont il
  ! a la coupling MOYENNE la plus forte avec les membres DEJA assignes (sous
  ! contrainte d'equilibre approximatif +-1), puis recurse sur A puis B --
  ! les blocs fortement couples finissent regroupes a TOUS les niveaux de la
  ! recursion, pas seulement localement comme le chainage. Ecrit le resultat
  ! dans ORDER(1:NBLK) au fil des appels recursifs (OPOS = position
  ! d'ecriture courante, INOUT, partagee entre tous les appels).
  RECURSIVE SUBROUTINE HIER_CLUSTER(GROUP, NG, ORDER)
    INTEGER, INTENT(IN)    :: NG
    INTEGER, INTENT(IN)    :: GROUP(NG)
    INTEGER, INTENT(INOUT) :: ORDER(:)
    INTEGER :: SIDE(NG), GA(NG), GB(NG)
    INTEGER :: NA, NBS, KK, JJ, SEED_A, SEED_B, PASS
    REAL(WP) :: SUM_A, SUM_B, RAV_A, RAV_B, BESTVAL

    IF (NG .LE. 1) THEN
      IF (NG .EQ. 1) THEN
        ORDER(OPOS) = GROUP(1); OPOS = OPOS + 1
      END IF
      RETURN
    END IF

    ! Graines : la paire la MOINS couplee de GROUP (2 membres surs de ne
    ! pas finir du meme cote), au lieu des 2 premiers elements arbitraires
    ! -- un mauvais choix de graines biaise toute la bipartition qui suit.
    SEED_A = 1; SEED_B = 2; BESTVAL = REAL(RANK_MAT(GROUP(1),GROUP(2)),WP)
    DO KK = 1, NG-1
      DO JJ = KK+1, NG
        IF (REAL(RANK_MAT(GROUP(KK),GROUP(JJ)),WP) .LT. BESTVAL) THEN
          BESTVAL = REAL(RANK_MAT(GROUP(KK),GROUP(JJ)),WP)
          SEED_A = KK; SEED_B = JJ
        END IF
      END DO
    END DO

    SIDE = -1
    SIDE(SEED_A) = 0; SIDE(SEED_B) = 1
    NA = 1; NBS = 1
    DO KK = 1, NG
      IF (KK.EQ.SEED_A .OR. KK.EQ.SEED_B) CYCLE
      SUM_A = 0.0_WP; SUM_B = 0.0_WP
      DO JJ = 1, NG
        IF (SIDE(JJ) .EQ. 0) THEN
          SUM_A = SUM_A + REAL(RANK_MAT(GROUP(KK),GROUP(JJ)),WP)
        ELSE IF (SIDE(JJ) .EQ. 1) THEN
          SUM_B = SUM_B + REAL(RANK_MAT(GROUP(KK),GROUP(JJ)),WP)
        END IF
      END DO
      RAV_A = SUM_A / MAX(NA,1)
      RAV_B = SUM_B / MAX(NBS,1)
      IF (NA .GE. (NG+1)/2) THEN
        SIDE(KK) = 1; NBS = NBS + 1
      ELSE IF (NBS .GE. NG/2) THEN
        SIDE(KK) = 0; NA = NA + 1
      ELSE IF (RAV_A .GE. RAV_B) THEN
        SIDE(KK) = 0; NA = NA + 1
      ELSE
        SIDE(KK) = 1; NBS = NBS + 1
      END IF
    END DO

    ! Raffinement (Kernighan-Lin simplifie, 2 passes) : l'affectation
    ! initiale ci-dessus est EN FLUX (chaque element ne voit que les
    ! elements DEJA traites) -- ici chaque element est reevalue avec la
    ! composition COMPLETE et A JOUR des 2 cotes, et bascule si ca ameliore
    ! son propre couplage moyen SANS casser l'equilibre.
    DO PASS = 1, 2
      DO KK = 1, NG
        SUM_A = 0.0_WP; SUM_B = 0.0_WP
        DO JJ = 1, NG
          IF (JJ .EQ. KK) CYCLE
          IF (SIDE(JJ) .EQ. 0) THEN
            SUM_A = SUM_A + REAL(RANK_MAT(GROUP(KK),GROUP(JJ)),WP)
          ELSE
            SUM_B = SUM_B + REAL(RANK_MAT(GROUP(KK),GROUP(JJ)),WP)
          END IF
        END DO
        IF (SIDE(KK) .EQ. 0) THEN
          RAV_A = SUM_A / MAX(NA-1,1); RAV_B = SUM_B / MAX(NBS,1)
          IF (RAV_B .GT. RAV_A .AND. NA .GT. (NG+1)/2) THEN
            SIDE(KK) = 1; NA = NA - 1; NBS = NBS + 1
          END IF
        ELSE
          RAV_A = SUM_A / MAX(NA,1); RAV_B = SUM_B / MAX(NBS-1,1)
          IF (RAV_A .GT. RAV_B .AND. NBS .GT. NG/2) THEN
            SIDE(KK) = 0; NBS = NBS - 1; NA = NA + 1
          END IF
        END IF
      END DO
    END DO

    NA = 0; NBS = 0
    DO KK = 1, NG
      IF (SIDE(KK) .EQ. 0) THEN
        NA = NA + 1; GA(NA) = GROUP(KK)
      ELSE
        NBS = NBS + 1; GB(NBS) = GROUP(KK)
      END IF
    END DO

    CALL HIER_CLUSTER(GA(1:NA), NA, ORDER)
    CALL HIER_CLUSTER(GB(1:NBS), NBS, ORDER)
  END SUBROUTINE HIER_CLUSTER

END PROGRAM BENCH_BLR_RANK_PERMUTATION
