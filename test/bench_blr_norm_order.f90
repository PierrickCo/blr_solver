!===============================================================================
! PROTOTYPE SEQUENTIEL -- teste si reordonner les blocs PAR NORME (au lieu de
! par rang, voir bench_blr_rank_permutation.f90/bench_blr_rank_position.f90)
! affecte la PRECISION de la solution calculee, plutot que la vitesse.
!
! RAPPEL (voir la discussion) : le conditionnement de A lui-meme (rapport des
! valeurs singulieres extremes) est INVARIANT sous toute permutation lignes/
! colonnes -- une permutation est une matrice orthogonale, elle ne change
! jamais les valeurs singulieres. Ce qui PEUT reellement varier avec l'ordre :
! la facon dont l'erreur de troncature BLR (compression bas-rang, <= TOL par
! bloc) se propage au fil des mises a jour de Schur -- un bloc de GRANDE
! norme compresse/traite tot influence potentiellement plus de blocs qu'un
! bloc de meme norme traite tard. Teste ICI en mesurant l'erreur relative
! ||X-Xexact||/||Xexact|| REELLEMENT obtenue, pas une estimation theorique.
!
! MATRICE : noyau decroissant simple (1/(1+|I-J|)), compressible PARTOUT
! (contrairement a bench_blr_rank_heterogeneity.f90 -- ICI on ne veut PAS
! confondre la variable "rang" avec la variable "norme", donc la
! compressibilite reste uniforme). Seule variation : chaque BLOC recoit un
! facteur d'echelle NORM_OF(bloc) tire d'un hash deterministe, etale sur
! NORM_DECADES decades (defaut 4 : facteur 10^4 entre le bloc le plus fort et
! le plus faible) -- applique en sqrt(NORM_OF(I)*NORM_OF(J)) a chaque entree
! (I,J), pour que la norme d'un bloc soit bien une propriete DU bloc, pas
! d'une seule entree.
!
! COMPRESSION : RSVD adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE, voir
! blr_compress_mod.f90) au lieu de la SVD complete par defaut -- USE_RSVD=
! .TRUE. transmis a MPI_ASSEMBLE_BLR_MATRIX, x34-54 plus rapide sur un bloc
! bien compressible (mesure ailleurs dans le projet), permet de tester a une
! taille de probleme raisonnable sans attendre des minutes par variante.
!
! VARIANT (argument 1, obligatoire) :
!   natural   : ordre 1..NBLK, aucune permutation.
!   norm_desc : blocs elimines du plus FORT au plus FAIBLE (norme decroissante).
!   norm_asc  : blocs elimines du plus FAIBLE au plus FORT (norme croissante).
!
! Contrairement a bench_blr_rank_position.f90, PAS besoin de repeter chaque
! variante plusieurs fois : l'erreur relative est une quantite deterministe
! (memes operations flottantes, meme resultat a chaque lancement) -- aucun
! bruit de mesure systeme a moyenner, contrairement au temps d'horloge.
!
! Usage : mpirun -np 1 ./bench_blr_norm_order VARIANT [N] [NB] [TOL] [NRHS] [NORM_DECADES]
!         (defauts : N=12000 NB=200 TOL=1e-4 NRHS=4 NORM_DECADES=4.0)
!===============================================================================
PROGRAM BENCH_BLR_NORM_ORDER
  USE MPI_BLR_MATRIX_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: N, NB, NRHS, IERR, MYID0, NPROCS_WORLD, NBLK
  REAL(WP) :: TOL, NORM_DECADES
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=16) :: VARIANT
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  INTEGER, ALLOCATABLE :: PERM(:), ACTIVE_PERM(:)
  REAL(WP), ALLOCATABLE :: NORMS(:)
  LOGICAL, ALLOCATABLE :: PLACED(:)
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), XEXACT_NEW(:,:), ABLOCK(:,:), BTMP(:,:)
  COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
  INTEGER :: I, J, K, KK, L, BEST, INFO, N_SWAPS, MLOC, RS, OLDBLK
  INTEGER :: N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR
  REAL(WP) :: ERR, BESTNORM, NUM2, DEN2

  N = 12000; NB = 200; TOL = 1.0E-4_WP; NRHS = 4; NORM_DECADES = 4.0_WP

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
         'Usage: mpirun -np 1 ./bench_blr_norm_order VARIANT [N] [NB] [TOL] [NRHS] [NORM_DECADES]'// &
         '  (VARIANT = natural|norm_desc|norm_asc)'
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
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) NORM_DECADES
  END IF
  IF (MOD(N,NB) .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A)') ' ERREUR : ce prototype suppose N multiple de NB.'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  NBLK = N / NB

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, 1, 1, IERR)
  IF (MYID0.EQ.0) WRITE(6,'(A,A,A,I0,A,I0,A,ES9.2,A,I0,A,F4.1)') &
       ' VARIANT=', TRIM(VARIANT), ' , N=', N, ' , NB=', NB, ' , TOL=', TOL, &
       ' , NRHS=', NRHS, ' , NORM_DECADES=', NORM_DECADES

  ALLOCATE(NORMS(NBLK), PERM(NBLK), ACTIVE_PERM(NBLK), PLACED(NBLK))
  DO K = 1, NBLK
    NORMS(K) = NORM_OF(K)
  END DO
  IF (MYID0.EQ.0) WRITE(6,'(A,ES10.3,A,ES10.3)') ' Norme de bloc : min=', MINVAL(NORMS), '  max=', MAXVAL(NORMS)

  SELECT CASE (TRIM(VARIANT))
  CASE ('natural')
    DO K = 1, NBLK
      PERM(K) = K
    END DO
  CASE ('norm_desc', 'norm_asc')
    ! Tri simple par selection (NBLK petit, O(NBLK^2) largement suffisant).
    PLACED = .FALSE.
    DO K = 1, NBLK
      BEST = 0
      IF (TRIM(VARIANT) .EQ. 'norm_desc') THEN
        BESTNORM = -HUGE(BESTNORM)
        DO J = 1, NBLK
          IF (.NOT.PLACED(J) .AND. NORMS(J) .GT. BESTNORM) THEN
            BESTNORM = NORMS(J); BEST = J
          END IF
        END DO
      ELSE
        BESTNORM = HUGE(BESTNORM)
        DO J = 1, NBLK
          IF (.NOT.PLACED(J) .AND. NORMS(J) .LT. BESTNORM) THEN
            BESTNORM = NORMS(J); BEST = J
          END IF
        END DO
      END IF
      PERM(K) = BEST; PLACED(BEST) = .TRUE.
    END DO
  CASE DEFAULT
    IF (MYID0.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : VARIANT inconnu "', TRIM(VARIANT), '" (attendu natural|norm_desc|norm_asc).'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END SELECT

  ! === Assemblage (RSVD) + factorisation de la matrice PERMUTEE ===========
  ACTIVE_PERM = PERM
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK_PERMUTED, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.FALSE., USE_RSVD=.TRUE.)
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A)') ' Assemblage (RSVD) : ', T_ASM, ' s , compresses ', &
       N_COMPRESSED, ' / ', N_OFFDIAG, ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'

  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                       N_SWAPS, T_FACTOR, SHOW_PROGRESS=.FALSE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,'(A,I0)') ' ERREUR factorisation, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A,I0,A,I0,A,F6.1,A)') ' Factorisation : ', T_FACTOR, ' s , compresses ', &
       N_COMPRESSED, ' / ', N_OFFDIAG, ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'

  ! === Second membre B = A_permutee * Xexact_permute, resolution, erreur ===
  ! Xexact_permute(nouvelle_ligne) = Xexact(ancienne_ligne) -- meme
  ! traduction bloc->bloc que GET_BLOCK_PERMUTED, pour comparer le X calcule
  ! (dans la numerotation PERMUTEE, jamais reconverti) a la BONNE reference.
  MLOC = N
  ALLOCATE(B_LOCAL(MLOC,NRHS), XEXACT_NEW(N,NRHS))
  DO K = 1, NBLK
    OLDBLK = PERM(K)
    DO KK = 1, NB
      DO J = 1, NRHS
        XEXACT_NEW((K-1)*NB+KK, J) = XEX_FORMULA((OLDBLK-1)*NB+KK, J)
      END DO
    END DO
  END DO

  RS = 1
  DO K = 1, NBLK
    ALLOCATE(ABLOCK(NB,N), BTMP(NB,NRHS))
    ! GET_BLOCK_PERMUTED traduit UN SEUL bloc a la fois (NEWJBLK derive de
    ! JSTART seul) -- convention valide pour MPI_ASSEMBLE_BLR_MATRIX (appelle
    ! toujours bloc par bloc), mais PAS pour un appel direct sur toute la
    ! largeur N : il faut boucler explicitement sur les blocs-colonnes L.
    DO L = 1, NBLK
      CALL GET_BLOCK_PERMUTED((K-1)*NB+1, K*NB, (L-1)*NB+1, L*NB, ABLOCK(1:NB,(L-1)*NB+1:L*NB))
    END DO
    ! CGEMM (BLAS externe, sans interface explicite) ecrit son resultat en
    ! supposant un pas memoire EGAL au LDC fourni -- passer directement une
    ! tranche de lignes B_LOCAL(RS:RS+NB-1,:) donnerait un LDC=NB alors que
    ! le vrai pas memoire de B_LOCAL est MLOC : CGEMM ecrirait alors au
    ! mauvais endroit. Tampon contigu BTMP, copie normale (Fortran, pas
    ! BLAS) vers B_LOCAL -- l'affectation respecte le vrai pas memoire de
    ! B_LOCAL, elle.
    CALL CGEMM('N','N', NB, NRHS, N, CONE, ABLOCK, NB, XEXACT_NEW, N, CZERO, BTMP, NB)
    B_LOCAL(RS:RS+NB-1,1:NRHS) = BTMP
    DEALLOCATE(ABLOCK, BTMP)
    RS = RS + NB
  END DO

  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS, B_LOCAL, MLOC, SHOW_PROGRESS=.FALSE.)

  NUM2 = SUM(ABS(B_LOCAL - XEXACT_NEW)**2)
  DEN2 = SUM(ABS(XEXACT_NEW)**2)
  ERR = SQRT(NUM2/DEN2)

  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A)') ''
    WRITE(6,'(A,ES14.6)') 'RESULT ERR_RELATIVE = ', ERR
    WRITE(6,'(A,F8.3,A)') 'RESULT T_ASM        = ', T_ASM, ' s'
    WRITE(6,'(A,F8.3,A)') 'RESULT T_FACTOR      = ', T_FACTOR, ' s'
    WRITE(6,'(A,I0)')     'RESULT N_SWAPS       = ', N_SWAPS
  END IF

  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Facteur d'echelle par bloc -- hash deterministe SANS RAPPORT avec l'ordre
  ! naturel, etale sur NORM_DECADES decades (ex : 4 -> facteur 1e4 entre le
  ! bloc le plus fort et le plus faible). Constantes DIFFERENTES de tout
  ! autre hash du projet (CLUSTER_OF, HASH_NOISE) pour eviter toute
  ! correlation involontaire.
  REAL(WP) FUNCTION NORM_OF(IBLK)
    INTEGER, INTENT(IN) :: IBLK
    REAL(WP) :: S, FRAC
    S = SIN(REAL(IBLK,WP)*53.471_WP)*19349.331_WP
    FRAC = S - REAL(FLOOR(S),WP)
    NORM_OF = 10.0_WP**(NORM_DECADES*(FRAC - 0.5_WP))
  END FUNCTION NORM_OF

  ! Noyau decroissant simple, compressible PARTOUT (contrairement a bench_
  ! blr_rank_heterogeneity.f90) -- volontaire, pour ne pas confondre la
  ! variable "rang" avec la variable "norme" testee ici.
  COMPLEX(WP) FUNCTION BASE_VALUE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: DF
    DF = 1.0_WP / (1.0_WP + ABS(REAL(GI-GJ,WP)))
    BASE_VALUE = CMPLX(DF, 0.3_WP*DF, KIND=WP)
  END FUNCTION BASE_VALUE

  COMPLEX(WP) FUNCTION XEX_FORMULA(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    XEX_FORMULA = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION XEX_FORMULA

  ! A(GI,GJ) = sqrt(NORM_OF(bloc GI)*NORM_OF(bloc GJ)) * noyau_decroissant --
  ! la norme d'un bloc est bien une propriete DU BLOC (facteur commun a
  ! toutes ses entrees), pas d'une entree isolee. Diagonale : boost fixe
  ! +50, NON mis a l'echelle (stabilite/inversibilite, comme les autres
  ! exemples du projet).
  SUBROUTINE GET_BLOCK_NATURAL(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ, IBLK, JBLK
    REAL(WP) :: SCALE
    DO JJ = JSTART, JEND
      JBLK = (JJ-1)/NB + 1
      DO II = ISTART, IEND
        IBLK = (II-1)/NB + 1
        SCALE = SQRT(NORM_OF(IBLK)*NORM_OF(JBLK))
        BLK(II-ISTART+1, JJ-JSTART+1) = SCALE * BASE_VALUE(II, JJ)
        IF (II .EQ. JJ) BLK(II-ISTART+1, JJ-JSTART+1) = BLK(II-ISTART+1, JJ-JSTART+1) + CMPLX(50.0_WP, 0.0_WP, KIND=WP)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK_NATURAL

  ! ACTIVE_PERM (host, positionne AVANT l'assemblage) : le NOUVEAU bloc K
  ! contient les donnees de l'ANCIEN bloc ACTIVE_PERM(K).
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

END PROGRAM BENCH_BLR_NORM_ORDER
