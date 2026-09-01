!===============================================================================
! BENCHMARK -- ordre de grandeur du gain d'une SVD RANDOMISEE (blr_compress_
! mod.f90) face a la SVD complete actuelle (COMPRESS_BLOCK, CGESVD) -- DEUX
! variantes comparees : COMPRESS_BLOCK_RSVD (rang cible connu a l'avance,
! premiere passe) ET COMPRESS_BLOCK_RSVD_ADAPTIVE (rang PAS connu a
! l'avance, decouvert par croissance du sous-espace -- le cas d'usage reel
! en production, ou TOL pilote le rang sans connaissance prealable).
!
! Contexte mesure (bench_blr_mpi_pivot_ooc_compare.f90, N=20000/NB=400) :
! COMPRESS_BLOCK (CGESVD) domine a ~90-95% du temps d'assemblage sur une
! matrice bien compressible -- c'est le cout que ce banc cherche a reduire.
!
! METHODOLOGIE : un SEUL bloc M x N representatif (noyau de Cauchy amorti,
! meme famille que les matrices de demonstration du projet, decroissance
! reglable via DECAY) est COMPRESSE une fois par COMPRESS_BLOCK (rang de
! reference R_TRUE, temps T_FULL), puis :
!   1) compresse NTRIALS fois par COMPRESS_BLOCK_RSVD AU MEME RANG CIBLE
!      R_TRUE (comparaison a QUALITE FINALE comparable, pas a budget de
!      calcul comparable) ;
!   2) compresse NTRIALS fois par COMPRESS_BLOCK_RSVD_ADAPTIVE (TOL seul,
!      AUCUNE connaissance de R_TRUE) -- rapporte en plus le rang trouve
!      (doit converger pres de R_TRUE) et le nombre de doublements de
!      sous-espace (N_GROWTHS, diagnostic du cout du tatonnement).
! Rapporte le temps moyen/pire cas ET l'erreur de reconstruction relative
! de chaque methode (la RSVD etant aleatoire, plusieurs tirages sont
! necessaires pour voir sa variabilite).
!
! Usage : ./bench_rsvd_compress [M] [N] [TOL] [DECAY] [OVERSAMPLE] [NPOWER] [NTRIALS] [RANK_WANT]
!         DECAY = distance (globale) entre les deux blocs -- DOIT etre >> M,N
!         pour un bloc en champ lointain qui compresse (voir plus bas).
!         RANK_WANT > 0 : IGNORE le noyau de Cauchy et construit a la place
!         un bloc a SPECTRE DE VALEURS SINGULIERES PRESCRIT (rang EXACT
!         RANK_WANT par construction, voir plus bas) -- utilise pour tester
!         le comportement de l'adaptatif pile a la limite de LCAP, endroit
!         ou le noyau de Cauchy transitionne de facon trop abrupte (dense
!         puis rang<=8) pour viser un rang precis autour de 100-140.
!         (defauts : M=N=400 TOL=1e-4 DECAY=5000.0 OVERSAMPLE=10 NPOWER=1 NTRIALS=20 RANK_WANT=0)
!===============================================================================
PROGRAM BENCH_RSVD_COMPRESS
  USE BLR_COMPRESS_MOD
  IMPLICIT NONE

  INTEGER :: M, N, OVERSAMPLE, NPOWER, NTRIALS, ITRIAL, RANK_WANT, MN, KK
  REAL(WP) :: TOL, DECAY, GAP
  CHARACTER(LEN=32) :: ARG
  COMPLEX(WP), ALLOCATABLE :: B(:,:), RECON(:,:), UMAT(:,:), VMAT(:,:)
  REAL(WP), ALLOCATABLE :: SIGMA(:)
  REAL :: RTMP1, RTMP2
  TYPE(BLR_BLOCK_T) :: BLK_FULL, BLK_RSVD, BLK_ADAPT
  DOUBLE PRECISION :: T0, T_FULL, T_RSVD, T_RSVD_SUM, T_RSVD_MIN, T_RSVD_MAX
  DOUBLE PRECISION :: T_ADAPT, T_ADAPT_SUM, T_ADAPT_MIN, T_ADAPT_MAX
  REAL(WP) :: NORMB, ERR_FULL, ERR_RSVD, ERR_RSVD_SUM, ERR_RSVD_MIN, ERR_RSVD_MAX
  REAL(WP) :: ERR_ADAPT, ERR_ADAPT_SUM, ERR_ADAPT_MIN, ERR_ADAPT_MAX
  INTEGER :: II, JJ, RANK_ADAPT_SUM, RANK_ADAPT_MIN, RANK_ADAPT_MAX, NGROW, NGROW_SUM, NGROW_MAX
  INTEGER :: N_LOWRANK_ADAPT

  M = 400; N = 400; TOL = 1.0E-4_WP; DECAY = 5000.0_WP
  OVERSAMPLE = 10; NPOWER = 1; NTRIALS = 20; RANK_WANT = 0

  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) M
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) DECAY
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) OVERSAMPLE
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) NPOWER
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) NTRIALS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN
    CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) RANK_WANT
  END IF

  WRITE(6,'(A,I0,A,I0,A,ES9.2,A,F10.1,A,I0,A,I0,A,I0)') ' M=', M, ' , N=', N, ' , TOL=', TOL, &
       ' , DISTANCE=', DECAY, ' , OVERSAMPLE=', OVERSAMPLE, ' , NPOWER=', NPOWER, ' , NTRIALS=', NTRIALS

  ALLOCATE(B(M,N))
  IF (RANK_WANT .GT. 0) THEN
    ! === Bloc a SPECTRE DE VALEURS SINGULIERES PRESCRIT ===
    ! B = UMAT * diag(SIGMA) * VMAT^H avec UMAT (MxMN), VMAT (NxMN) a
    ! colonnes ORTHONORMALES (QR de matrices aleatoires) : c'est LITTERALE-
    ! MENT une SVD reduite valide, donc B a pour valeurs singulieres EXACTES
    ! SIGMA(1:MN) -- rang numerique = RANK_WANT PAR CONSTRUCTION, sans
    ! dependre d'un noyau physique dont la transition rang/DECAY est trop
    ! abrupte pour viser une valeur precise pres de LCAP.
    MN = MIN(M,N)
    IF (RANK_WANT .GT. MN) RANK_WANT = MN
    ALLOCATE(UMAT(M,MN), VMAT(N,MN), SIGMA(MN))
    DO KK = 1, MN
      DO II = 1, M
        CALL RANDOM_NUMBER(RTMP1); CALL RANDOM_NUMBER(RTMP2)
        UMAT(II,KK) = CMPLX(RTMP1-0.5, RTMP2-0.5, KIND=WP)
      END DO
      DO JJ = 1, N
        CALL RANDOM_NUMBER(RTMP1); CALL RANDOM_NUMBER(RTMP2)
        VMAT(JJ,KK) = CMPLX(RTMP1-0.5, RTMP2-0.5, KIND=WP)
      END DO
    END DO
    CALL ORTHONORMALIZE_COLS(UMAT, M, MN)
    CALL ORTHONORMALIZE_COLS(VMAT, N, MN)
    ! Plateau haut = 1 (les RANK_WANT premieres directions), chute NETTE
    ! sous TOL (le reste) -- pas de zone grise intermediaire, pour que le
    ! rang numerique soit EXACTEMENT RANK_WANT a la TOL demandee.
    GAP = TOL * 1.0E-3_WP
    DO KK = 1, MN
      IF (KK .LE. RANK_WANT) THEN
        SIGMA(KK) = 1.0_WP
      ELSE
        SIGMA(KK) = GAP
      END IF
      VMAT(:,KK) = SIGMA(KK) * VMAT(:,KK)
    END DO
    CALL CGEMM('N', 'C', M, N, MN, (1.0_WP,0.0_WP), UMAT, M, VMAT, N, (0.0_WP,0.0_WP), B, M)
    DEALLOCATE(UMAT, VMAT, SIGMA)
    WRITE(6,'(A,I0,A,ES9.2)') ' Bloc a spectre PRESCRIT : RANK_WANT=', RANK_WANT, ' , GAP (valeurs sous le rang)=', GAP
  ELSE
    ! Bloc de demonstration : noyau de Cauchy amorti 1/(1+|DISTANCE+II-JJ|),
    ! DISTANCE >> M,N -- represente un bloc hors-diagonale EN CHAMP LOINTAIN
    ! (deux blocs de la matrice globale suffisamment eloignes l'un de
    ! l'autre), meme situation que les blocs qui compressent bien dans le
    ! reste du projet : la fonction varie LENTEMENT sur toute la fenetre
    ! locale (DISTANCE domine II-JJ), d'ou un rang numerique faible. Un
    ! noyau centre sur II=JJ (decroissance locale) resterait au contraire
    ! quasi plein rang -- ce n'est PAS le regime des blocs hors-diagonale
    ! BLR reels, qui sont toujours suffisamment loin de la diagonale globale
    ! pour tomber en champ lointain (DECAY, ici, est cette DISTANCE).
    DO JJ = 1, N
      DO II = 1, M
        B(II,JJ) = CMPLX(1.0_WP/(1.0_WP+ABS(DECAY+REAL(II-JJ,WP))), &
                          0.3_WP/(1.0_WP+ABS(DECAY+REAL(II-JJ,WP))), KIND=WP)
      END DO
    END DO
  END IF
  NORMB = SQRT(REAL(SUM(ABS(B)**2),WP))

  ! === Reference : SVD complete (CGESVD) ===
  T0 = OMP_GET_WTIME_LOCAL()
  CALL COMPRESS_BLOCK(M, N, B, M, TOL, BLK_FULL)
  T_FULL = OMP_GET_WTIME_LOCAL() - T0
  ALLOCATE(RECON(M,N))

  IF (.NOT. BLK_FULL%IS_LOWRANK) THEN
    ! Bloc NON compressible a ce DECAY/TOL (reste dense cote reference) --
    ! on saute la comparaison a rang cible fixe (pas de R_TRUE a copier),
    ! mais on teste QUAND MEME l'ADAPTATIF plus bas : c'est justement le
    ! cas interessant a verifier (l'adaptatif doit lui aussi conclure
    ! "dense", pas inventer un faux rang).
    WRITE(6,'(A)') ' Ce bloc ne compresse pas a ce DECAY/TOL (reference = DENSE).'
    ERR_FULL = 0.0_WP
  ELSE
    CALL EXPAND_BLOCK(BLK_FULL, RECON, M)
    ERR_FULL = SQRT(REAL(SUM(ABS(B-RECON)**2),WP)) / NORMB

    WRITE(6,'(A,I0)') ' Rang trouve par COMPRESS_BLOCK (reference, R_TRUE) = ', BLK_FULL%RANK
    WRITE(6,'(A,F10.6,A)') ' Temps COMPRESS_BLOCK (CGESVD complet)     = ', T_FULL, ' s'
    WRITE(6,'(A,ES12.4)')  ' Erreur relative ||B-UV||/||B|| (reference) = ', ERR_FULL

    ! === COMPRESS_BLOCK_RSVD au MEME rang cible, NTRIALS tirages ===
    T_RSVD_SUM = 0.0D0; T_RSVD_MIN = HUGE(T_RSVD_MIN); T_RSVD_MAX = 0.0D0
    ERR_RSVD_SUM = 0.0_WP; ERR_RSVD_MIN = HUGE(ERR_RSVD_MIN); ERR_RSVD_MAX = 0.0_WP
    DO ITRIAL = 1, NTRIALS
      T0 = OMP_GET_WTIME_LOCAL()
      CALL COMPRESS_BLOCK_RSVD(M, N, B, M, BLK_FULL%RANK, OVERSAMPLE, NPOWER, BLK_RSVD)
      T_RSVD = OMP_GET_WTIME_LOCAL() - T0
      T_RSVD_SUM = T_RSVD_SUM + T_RSVD
      T_RSVD_MIN = MIN(T_RSVD_MIN, T_RSVD)
      T_RSVD_MAX = MAX(T_RSVD_MAX, T_RSVD)

      CALL EXPAND_BLOCK(BLK_RSVD, RECON, M)
      ERR_RSVD = SQRT(REAL(SUM(ABS(B-RECON)**2),WP)) / NORMB
      ERR_RSVD_SUM = ERR_RSVD_SUM + ERR_RSVD
      ERR_RSVD_MIN = MIN(ERR_RSVD_MIN, ERR_RSVD)
      ERR_RSVD_MAX = MAX(ERR_RSVD_MAX, ERR_RSVD)

      DEALLOCATE(BLK_RSVD%UFAC, BLK_RSVD%VFAC)
    END DO

    WRITE(6,'(A)') ' '
    WRITE(6,'(A,F10.6,A,F10.6,A,F10.6,A)') ' Temps COMPRESS_BLOCK_RSVD (moy/min/max) = ', &
         T_RSVD_SUM/REAL(NTRIALS,8), ' / ', T_RSVD_MIN, ' / ', T_RSVD_MAX, ' s'
    WRITE(6,'(A,F8.2)') ' Ratio de temps (CGESVD complet / RSVD, moyenne) = ', T_FULL/(T_RSVD_SUM/REAL(NTRIALS,8))
    WRITE(6,'(A,ES12.4,A,ES12.4,A,ES12.4)') ' Erreur relative RSVD (moy/min/max) = ', &
         ERR_RSVD_SUM/REAL(NTRIALS,WP), ' / ', ERR_RSVD_MIN, ' / ', ERR_RSVD_MAX
    WRITE(6,'(A,F8.2)') ' Erreur RSVD / erreur reference (moyenne, x fois plus imprecis) = ', &
         (ERR_RSVD_SUM/REAL(NTRIALS,WP)) / MAX(ERR_FULL,TINY(1.0_WP))
  END IF

  ! === COMPRESS_BLOCK_RSVD_ADAPTIVE -- TOL seul, AUCUN rang connu a
  ! l'avance (le cas d'usage reel), NTRIALS tirages ===
  T_ADAPT_SUM = 0.0D0; T_ADAPT_MIN = HUGE(T_ADAPT_MIN); T_ADAPT_MAX = 0.0D0
  ERR_ADAPT_SUM = 0.0_WP; ERR_ADAPT_MIN = HUGE(ERR_ADAPT_MIN); ERR_ADAPT_MAX = 0.0_WP
  RANK_ADAPT_SUM = 0; RANK_ADAPT_MIN = HUGE(RANK_ADAPT_MIN); RANK_ADAPT_MAX = 0
  NGROW_SUM = 0; NGROW_MAX = 0; N_LOWRANK_ADAPT = 0
  DO ITRIAL = 1, NTRIALS
    T0 = OMP_GET_WTIME_LOCAL()
    CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(M, N, B, M, TOL, BLK_ADAPT, OVERSAMPLE, NPOWER, N_GROWTHS=NGROW)
    T_ADAPT = OMP_GET_WTIME_LOCAL() - T0
    T_ADAPT_SUM = T_ADAPT_SUM + T_ADAPT
    T_ADAPT_MIN = MIN(T_ADAPT_MIN, T_ADAPT)
    T_ADAPT_MAX = MAX(T_ADAPT_MAX, T_ADAPT)
    NGROW_SUM = NGROW_SUM + NGROW
    NGROW_MAX = MAX(NGROW_MAX, NGROW)

    IF (BLK_ADAPT%IS_LOWRANK) THEN
      N_LOWRANK_ADAPT = N_LOWRANK_ADAPT + 1
      RANK_ADAPT_SUM = RANK_ADAPT_SUM + BLK_ADAPT%RANK
      RANK_ADAPT_MIN = MIN(RANK_ADAPT_MIN, BLK_ADAPT%RANK)
      RANK_ADAPT_MAX = MAX(RANK_ADAPT_MAX, BLK_ADAPT%RANK)
      CALL EXPAND_BLOCK(BLK_ADAPT, RECON, M)
      ERR_ADAPT = SQRT(REAL(SUM(ABS(B-RECON)**2),WP)) / NORMB
      DEALLOCATE(BLK_ADAPT%UFAC, BLK_ADAPT%VFAC)
    ELSE
      ! Replie dense (n'a pas trouve de rang satisfaisant sous LCAP) --
      ! erreur nulle par definition (reconstruction exacte), mais signale
      ! un DESACCORD avec la reference si celle-ci, elle, a compresse.
      ERR_ADAPT = 0.0_WP
      DEALLOCATE(BLK_ADAPT%DENSE)
    END IF
    ERR_ADAPT_SUM = ERR_ADAPT_SUM + ERR_ADAPT
    ERR_ADAPT_MIN = MIN(ERR_ADAPT_MIN, ERR_ADAPT)
    ERR_ADAPT_MAX = MAX(ERR_ADAPT_MAX, ERR_ADAPT)
  END DO

  WRITE(6,'(A)') ' '
  WRITE(6,'(A,I0,A,I0,A)') ' [ADAPTATIF] Decisions : ', N_LOWRANK_ADAPT, '/', NTRIALS, ' tirages BASSE-RANG (le reste DENSE)'
  IF (N_LOWRANK_ADAPT .GT. 0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0)') ' [ADAPTATIF] Rang trouve quand basse-rang (moy/min/max) = ', &
         RANK_ADAPT_SUM/N_LOWRANK_ADAPT, ' / ', RANK_ADAPT_MIN, ' / ', RANK_ADAPT_MAX
  END IF
  IF (BLK_FULL%IS_LOWRANK) THEN
    WRITE(6,'(A,I0,A)') ' [ADAPTATIF] Reference (COMPRESS_BLOCK) a trouve R_TRUE = ', BLK_FULL%RANK, &
         ' -- comparer aux valeurs ci-dessus'
  ELSE
    WRITE(6,'(A)') ' [ADAPTATIF] Reference (COMPRESS_BLOCK) = DENSE -- coherent SEULEMENT si 0/'// &
         'NTRIALS ci-dessus (sinon l''adaptatif invente a tort un rang la ou la reference exacte n''en trouve pas)'
  END IF
  WRITE(6,'(A,F5.2,A,I0)') ' [ADAPTATIF] Doublements de sous-espace (moy/max) = ', &
       REAL(NGROW_SUM,WP)/REAL(NTRIALS,WP), ' / ', NGROW_MAX
  WRITE(6,'(A,F10.6,A,F10.6,A,F10.6,A)') ' [ADAPTATIF] Temps (moy/min/max) = ', &
       T_ADAPT_SUM/REAL(NTRIALS,8), ' / ', T_ADAPT_MIN, ' / ', T_ADAPT_MAX, ' s'
  WRITE(6,'(A,F8.2)') ' [ADAPTATIF] Ratio de temps (CGESVD complet / adaptatif, moyenne) = ', &
       T_FULL/(T_ADAPT_SUM/REAL(NTRIALS,8))
  WRITE(6,'(A,ES12.4,A,ES12.4,A,ES12.4)') ' [ADAPTATIF] Erreur relative (moy/min/max) = ', &
       ERR_ADAPT_SUM/REAL(NTRIALS,WP), ' / ', ERR_ADAPT_MIN, ' / ', ERR_ADAPT_MAX

  DEALLOCATE(B, RECON)
  IF (ALLOCATED(BLK_FULL%UFAC)) DEALLOCATE(BLK_FULL%UFAC)
  IF (ALLOCATED(BLK_FULL%VFAC)) DEALLOCATE(BLK_FULL%VFAC)

CONTAINS

  ! Horloge murale portable (pas de MPI ici, donc pas de MPI_WTIME) --
  ! SYSTEM_CLOCK intrinseque, independant d'OpenMP/MPI.
  DOUBLE PRECISION FUNCTION OMP_GET_WTIME_LOCAL()
    INTEGER(8) :: COUNT, COUNT_RATE
    CALL SYSTEM_CLOCK(COUNT, COUNT_RATE)
    OMP_GET_WTIME_LOCAL = REAL(COUNT,8) / REAL(COUNT_RATE,8)
  END FUNCTION OMP_GET_WTIME_LOCAL

END PROGRAM BENCH_RSVD_COMPRESS
