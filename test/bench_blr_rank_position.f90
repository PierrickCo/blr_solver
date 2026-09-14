!===============================================================================
! PROTOTYPE SEQUENTIEL -- teste si la POSITION (debut / milieu / fin de l'ordre
! d'elimination 1..NBLK) d'un groupe de blocs fortement couples DEJA COMPACT
! change le temps de factorisation, INDEPENDAMMENT du fait qu'il soit compact.
! Reprend la MEME matrice heterogene deliberee que bench_blr_rank_
! heterogeneity.f90/bench_blr_rank_permutation.f90 (voir leur commentaire de
! tete).
!
! RESTRUCTURE EN LANCEMENTS INDEPENDANTS (un seul VARIANT par appel du
! programme, choisi en argument) -- la premiere version (tout-en-un, un seul
! processus enchainant les 5 variantes) souffrait d'un biais de reechauffement
! systeme : chaque variante suivante beneficiait d'un systeme "plus chaud" que
! la precedente (frequence CPU, cache), produisant un gain monotone factice
! (jusqu'a x3 sur la derniere variante) qui contredisait la mesure, plus
! rigoureuse, faite sur des lancements independants pour NATUREL/CHAINAGE
! (bench_blr_rank_permutation.f90, 6 repetitions). Chaque appel ICI ne fait
! QU'UNE seule paire assemblage+factorisation (celle du VARIANT demande) --
! aucun ordre a l'interieur d'un meme processus ne peut biaiser la mesure. Un
! script externe repete chaque VARIANT plusieurs fois (lancements separes) et
! moyenne, exactement comme pour NATUREL/CHAINAGE precedemment.
!
! VARIANT (obligatoire, argument 6) :
!   natural : ordre 1..NBLK, aucune permutation.
!   chain   : ANCIEN regroupement complet (chainage glouton par plus proche
!             voisin sur le rang MESURE) -- necessite un premier assemblage
!             a l'ordre naturel pour mesurer RANK_MAT avant de permuter (ce
!             premier assemblage n'est PAS chronometre dans T_ASM/T_FACTOR
!             rapportes -- seul le second assemblage, a l'ordre permute, et
!             sa factorisation comptent).
!   start   : UN SEUL groupe cache (TARGET_CLUSTER, connu directement via
!             CLUSTER_OF -- volontaire, voir bench_blr_rank_position.f90
!             version precedente) rendu compact et place TOUT AU DEBUT.
!   mid     : meme groupe compact, insere AU MILIEU des autres blocs (ordre
!             relatif naturel conserve pour les autres).
!   end     : meme groupe compact, place TOUT A LA FIN.
!
! SEPARATION EXPLICITE DES TEMPS (T_PERM/T_ASM/T_FACTOR) ET DE LA COMPRESSION
! (avant assemblage / apres factorisation, comptage de blocs + taux reel base
! sur le rang) -- imprimes par CE processus, a parser par le script externe.
!
! Usage : mpirun -np 1 ./bench_blr_rank_position VARIANT [N] [NB] [TOL] [NCLUSTERS] [TARGET_CLUSTER]
!         (defauts : N=20000 NB=200 TOL=1e-4 NCLUSTERS=8 TARGET_CLUSTER=0)
!===============================================================================
PROGRAM BENCH_BLR_RANK_POSITION
  USE MPI_BLR_MATRIX_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: N, NB, IERR, MYID0, NCLUSTERS, NBLK, NPROCS_WORLD, TARGET_CLUSTER
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=16) :: VARIANT
  TYPE(MPI_BLR_MATRIX_T) :: A_PREP, A_TARGET
  INTEGER, ALLOCATABLE :: PERM(:), ACTIVE_PERM(:), RANK_MAT(:,:)
  LOGICAL, ALLOCATABLE :: PLACED(:)
  INTEGER, ALLOCATABLE :: CLUSTER_MEMBERS(:), OTHER_MEMBERS(:)
  INTEGER :: I, J, K, BEST, BESTVAL, CUR, NCLUST_MEM, NOTHER, HALF
  DOUBLE PRECISION :: TP0, TP1, T_PERM, T_ASM_PREP, T_ASM, T_FACTOR
  INTEGER :: NC, NO, INFO, N_SWAPS, NC_ASM, NO_ASM, NC_FACT, NO_FACT
  INTEGER(8) :: RSUM, RS_ASM, RS_FACT

  N = 20000; NB = 200; TOL = 1.0E-4_WP; NCLUSTERS = 8; TARGET_CLUSTER = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS_WORLD, IERR)
  IF (NPROCS_WORLD .NE. 1) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : ce prototype est volontairement SEQUENTIEL -- lancer avec -np 1.'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF

  IF (COMMAND_ARGUMENT_COUNT() .LT. 1) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') &
         'Usage: mpirun -np 1 ./bench_blr_rank_position VARIANT [N] [NB] [TOL] [NCLUSTERS] [TARGET_CLUSTER]'// &
         '  (VARIANT = natural|chain|start|mid|end)'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, VARIANT)
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NCLUSTERS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TARGET_CLUSTER
  END IF
  IF (MOD(N,NB) .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : ce prototype suppose N multiple de NB.'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  NBLK = N / NB

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, 1, 1, IERR)
  IF (MYID0.EQ.0) WRITE(6,'(A,A,A,I0,A,I0,A,ES9.2,A,I0,A,I0)') &
       ' VARIANT=', TRIM(VARIANT), ' , N=', N, ' , NB=', NB, ' , TOL=', TOL, &
       ' , NCLUSTERS=', NCLUSTERS, ' , NBLK=', NBLK

  ALLOCATE(ACTIVE_PERM(NBLK), PERM(NBLK))
  T_PERM = 0.0D0; T_ASM_PREP = 0.0D0

  SELECT CASE (TRIM(VARIANT))
  CASE ('natural')
    TP0 = MPI_WTIME()
    DO K = 1, NBLK
      PERM(K) = K
    END DO
    TP1 = MPI_WTIME(); T_PERM = TP1 - TP0

  CASE ('chain')
    ! Assemblage PREPARATOIRE a l'ordre naturel, UNIQUEMENT pour mesurer
    ! RANK_MAT -- PAS chronometre dans T_ASM/T_FACTOR rapportes (voir le
    ! commentaire de tete).
    DO K = 1, NBLK
      ACTIVE_PERM(K) = K
    END DO
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK_PERMUTED, A_PREP, T_ASM_PREP, &
                                  NC, NO, RSUM, SHOW_PROGRESS=.FALSE.)
    ALLOCATE(RANK_MAT(NBLK,NBLK))
    DO I = 1, NBLK
      DO J = 1, NBLK
        IF (I .EQ. J) THEN
          RANK_MAT(I,J) = 0
        ELSE IF (A_PREP%OFFBLK(I,J)%IS_LOWRANK) THEN
          RANK_MAT(I,J) = A_PREP%OFFBLK(I,J)%RANK
        ELSE
          RANK_MAT(I,J) = NB
        END IF
      END DO
    END DO

    TP0 = MPI_WTIME()
    ALLOCATE(PLACED(NBLK))
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
    TP1 = MPI_WTIME(); T_PERM = TP1 - TP0

  CASE ('start', 'mid', 'end')
    ALLOCATE(CLUSTER_MEMBERS(NBLK), OTHER_MEMBERS(NBLK))
    NCLUST_MEM = 0; NOTHER = 0
    DO K = 1, NBLK
      IF (CLUSTER_OF(K) .EQ. TARGET_CLUSTER) THEN
        NCLUST_MEM = NCLUST_MEM + 1; CLUSTER_MEMBERS(NCLUST_MEM) = K
      ELSE
        NOTHER = NOTHER + 1; OTHER_MEMBERS(NOTHER) = K
      END IF
    END DO
    IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A)') ' Groupe cible (cluster ', TARGET_CLUSTER, ') : ', NCLUST_MEM, ' blocs'

    TP0 = MPI_WTIME()
    IF (TRIM(VARIANT) .EQ. 'start') THEN
      PERM(1:NCLUST_MEM) = CLUSTER_MEMBERS(1:NCLUST_MEM)
      PERM(NCLUST_MEM+1:NBLK) = OTHER_MEMBERS(1:NOTHER)
    ELSE IF (TRIM(VARIANT) .EQ. 'end') THEN
      PERM(1:NOTHER) = OTHER_MEMBERS(1:NOTHER)
      PERM(NOTHER+1:NBLK) = CLUSTER_MEMBERS(1:NCLUST_MEM)
    ELSE  ! 'mid'
      HALF = NOTHER / 2
      PERM(1:HALF) = OTHER_MEMBERS(1:HALF)
      PERM(HALF+1:HALF+NCLUST_MEM) = CLUSTER_MEMBERS(1:NCLUST_MEM)
      PERM(HALF+NCLUST_MEM+1:NBLK) = OTHER_MEMBERS(HALF+1:NOTHER)
    END IF
    TP1 = MPI_WTIME(); T_PERM = TP1 - TP0

  CASE DEFAULT
    IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : VARIANT inconnu "', TRIM(VARIANT), '" (attendu natural|chain|start|mid|end).'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END SELECT

  ! === Assemblage + factorisation DE LA VARIANTE demandee -- SEULES mesures
  !     rapportees comme T_ASM/T_FACTOR (jamais l'assemblage preparatoire
  !     eventuel de 'chain' ci-dessus). ======================================
  ACTIVE_PERM = PERM
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK_PERMUTED, A_TARGET, T_ASM, &
                                NC, NO, RSUM, SHOW_PROGRESS=.FALSE.)
  NC_ASM = NC; NO_ASM = NO; RS_ASM = RSUM
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_TARGET, TOL, INFO, NC, NO, RSUM, N_SWAPS, T_FACTOR, SHOW_PROGRESS=.FALSE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR factorisation, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  NC_FACT = NC; NO_FACT = NO; RS_FACT = RSUM

  IF (MYID0.EQ.0) THEN
    IF (TRIM(VARIANT).EQ.'chain') WRITE(6,'(A,F8.3,A)') ' [prep, non compte] Assemblage naturel pour mesurer RANK_MAT : ', &
         T_ASM_PREP, ' s'
    WRITE(6,'(A,F8.4,A)') 'RESULT T_PERM   = ', T_PERM, ' s'
    WRITE(6,'(A,F8.3,A)') 'RESULT T_ASM    = ', T_ASM, ' s'
    WRITE(6,'(A,F8.3,A)') 'RESULT T_FACTOR = ', T_FACTOR, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') 'RESULT COMPRESS_ASM  = ', NC_ASM, ' / ', NO_ASM, ' (', &
         100.0_WP*REAL(NC_ASM,WP)/REAL(MAX(NO_ASM,1),WP), ' %)'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') 'RESULT COMPRESS_FACT = ', NC_FACT, ' / ', NO_FACT, ' (', &
         100.0_WP*REAL(NC_FACT,WP)/REAL(MAX(NO_FACT,1),WP), ' %)'
    WRITE(6,'(A,F6.2,A)') 'RESULT REAL_COMPRESS_ASM  = ', REAL_COMPRESSION_PCT(NC_ASM,NO_ASM,RS_ASM), ' %'
    WRITE(6,'(A,F6.2,A)') 'RESULT REAL_COMPRESS_FACT = ', REAL_COMPRESSION_PCT(NC_FACT,NO_FACT,RS_FACT), ' %'
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM_ARG)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM_ARG
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM_ARG,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
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
  ! heterogeneity.f90/bench_blr_rank_permutation.f90, indexee par les
  ! indices GLOBAUX (naturels) II,JJ.
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

  ! ACTIVE_PERM (host, positionne AVANT chaque appel a MPI_ASSEMBLE_BLR_
  ! MATRIX) : le NOUVEAU bloc K contient les donnees de l'ANCIEN bloc
  ! ACTIVE_PERM(K).
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

END PROGRAM BENCH_BLR_RANK_POSITION
