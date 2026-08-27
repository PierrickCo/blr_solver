!===============================================================================
! ETAPE 2 (partie 2/2) : factorisation LU par blocs AVEC compression BLR
! des blocs hors-diagonale, complexe simple precision.
!
! CHOIX DE PORTEE DELIBERE : la compression n'intervient QUE dans le calcul
! de la mise a jour de Schur (A22 <- A22 - L21*U12), pas dans le STOCKAGE
! final des facteurs -- L et U restent stockes DENSES dans A, exactement
! comme a l'etape 1, et BLOCK_LU_SOLVE (module BLOCK_LU_MOD, etape 1,
! INCHANGE) reste directement utilisable pour la resolution. Ceci
! correspond exactement a l'option ICNTL(35)=3 de MUMPS / au mode
! "BLR factorization but full-rank solve" de STRUMPACK (voir la memoire du
! projet STRUMPACK) : le gain vise ici est la REDUCTION DU COUT DE CALCUL
! de la factorisation (le GEMM de mise a jour de Schur devient basse-rang
! quand les blocs le permettent), pas encore la reduction memoire des
! facteurs stockes -- ce sera une extension naturelle pour plus tard
! (stocker L21/U12 compresses eux-memes, et adapter le solve en
! consequence).
!
! CONSEQUENCE NUMERIQUE ATTENDUE : la factorisation devient APPROCHEE (plus
! seulement a la precision machine comme a l'etape 1) -- l'erreur est
! bornee par la tolerance de compression TOL, exactement comme dans un vrai
! solveur BLR. Le test associe (test/test_blr_lu.f90) verifie que le residu
! final reste de l'ordre de TOL, pas qu'il est nul.
!
! ETAPE 5 : PARALLELISATION OPENMP.
!
! OU PARALLELISER : PAS le facteur diagonal (CGETRF/CTRSM) -- ce sont des
! appels UNIQUES par panneau, deja parallelises EN INTERNE par MKL
! (threading BLAS/LAPACK standard, "gratuit" avec la MKL threadee deja
! utilisee dans ce projet). La ou NOUS avons du travail a repartir, c'est
! la ou l'etape 2 a transforme UN GROS CGEMM (etape 1) en BEAUCOUP DE
! PETITS appels independants : les boucles de compression (une par bloc de
! colonnes/lignes trainantes) et la double boucle de mise a jour de Schur
! (une par PAIRE de blocs trainants (I,J)) -- exactement le genre de
! "beaucoup de petites taches independantes" qu'OpenMP fait bien. Chaque
! iteration ecrit dans une zone de A OU un element de tableau DISJOINT des
! autres (voir les commentaires sur chaque boucle) -- pas de section
! critique necessaire.
!
! SUR-SOUSCRIPTION MKL/OPENMP (piege classique, gere explicitement) : si on
! laisse MKL utiliser plusieurs threads EN INTERNE pour chaque petit CGEMM
! DANS notre propre boucle OpenMP, on obtient un threading imbrique
! (N_OMP x N_OMP threads en concurrence sur le CPU) qui ralentit au lieu
! d'accelerer. On force donc MKL a 1 thread par appel PENDANT les boucles
! paralleles OpenMP (MKL_SET_NUM_THREADS(1)), et on le restaure a N_OMP
! juste apres pour que les GROS appels sequentiels (CGETRF/CTRSM du
! panneau suivant) profitent eux du threading interne de MKL. Cet appel
! utilise l'interface F77 classique de MKL (comme CGETRF/CGEMM, pas de
! module MKL_SERVICE necessaire).
!===============================================================================
MODULE BLR_LU_MOD
  USE BLR_COMPRESS_MOD
  USE OMP_LIB
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! Factorisation LU par blocs avec compression BLR des mises a jour de
  ! Schur. Memes conventions de sortie (A, IPIV) que BLOCK_LU_FACTOR de
  ! l'etape 1 -- utiliser BLOCK_LU_SOLVE (module BLOCK_LU_MOD) pour
  ! resoudre. TOL : tolerance relative de compression (0 = desactive la
  ! compression, redonne exactement l'etape 1). Statistiques de
  ! compression renvoyees pour diagnostic/tests.
  !-----------------------------------------------------------------------
  ! N_SWAPS/MIN_PIVOT (optionnels) : memes diagnostics de pivot que
  ! BLOCK_LU_FACTOR (module BLOCK_LU_MOD) -- voir la note associee la-bas.
  ! T_COMPRESS (optionnel) : temps CUMULE (horloge murale, OMP_GET_WTIME)
  ! passe DANS les deux boucles de compression (etapes 4+5 ci-dessous),
  ! mesure separement du reste (facteur diagonal dense + mise a jour de
  ! Schur) -- demande explicite de l'utilisateur pour voir OU va le temps,
  ! pas seulement le total. L'appelant, qui chronometre deja l'appel
  ! entier, peut en deduire le temps de factorisation "pure" par
  ! soustraction (T_TOTAL - T_COMPRESS) sans qu'on ait besoin de le
  ! renvoyer separement.
  SUBROUTINE BLR_LU_FACTOR(N, NB, A, LDA, IPIV, TOL, INFO, &
                            N_COMPRESSED, N_OFFDIAG_BLOCKS, RANK_SUM, &
                            N_SWAPS, MIN_PIVOT, T_COMPRESS)
    INTEGER, INTENT(IN)          :: N, NB, LDA
    COMPLEX(WP), INTENT(INOUT)   :: A(LDA,N)
    INTEGER, INTENT(OUT)         :: IPIV(N)
    REAL(WP), INTENT(IN)         :: TOL
    INTEGER, INTENT(OUT)         :: INFO
    INTEGER, INTENT(OUT)         :: N_COMPRESSED, N_OFFDIAG_BLOCKS
    INTEGER(8), INTENT(OUT)      :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL  :: N_SWAPS
    REAL(WP), INTENT(OUT), OPTIONAL :: MIN_PIVOT
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_COMPRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP, 0.0_WP)
    INTEGER :: NBLK, K, I, J, KSTART, KEND, KW, IINFO, MI, NJ, II, NS
    INTEGER :: NTHREADS, TID
    REAL(WP) :: MP
    DOUBLE PRECISION :: TC, TC0
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(BLR_BLOCK_T), ALLOCATABLE :: LBLK(:), UBLK(:)
    TYPE(BLR_WORKSPACE_T), ALLOCATABLE :: WSARR(:)

    CALL BLOCK_PARTITION(N, NB, NBLK, BSTART, BEND)
    ALLOCATE(LBLK(NBLK), UBLK(NBLK))

    ! UN workspace SVD PAR THREAD OpenMP (jamais partage -- voir la note de
    ! thread-safety de BLR_WORKSPACE_T) pour les deux boucles paralleles de
    ! compression ci-dessous (etapes 4 et 5) : chaque thread lit/ecrit
    ! exclusivement WSARR(OMP_GET_THREAD_NUM()+1), sans aucune contention.
    NTHREADS = OMP_GET_MAX_THREADS()
    ALLOCATE(WSARR(NTHREADS))
    DO TID = 1, NTHREADS
      CALL BLR_WORKSPACE_ALLOC(WSARR(TID), NB)
    END DO

    INFO = 0
    N_COMPRESSED = 0
    N_OFFDIAG_BLOCKS = 0
    RANK_SUM = 0_8
    NS = 0
    MP = HUGE(1.0_WP)
    TC = 0.0D0

    DO K = 1, NBLK
      KSTART = BSTART(K); KEND = BEND(K); KW = KEND - KSTART + 1

      ! 1) Factoriser le panneau diagonal (dense, pivot partiel -- comme
      !    l'etape 1)
      CALL CGETRF(N-KSTART+1, KW, A(KSTART,KSTART), LDA, IPIV(KSTART), IINFO)
      IF (IINFO .NE. 0 .AND. INFO .EQ. 0) INFO = KSTART - 1 + IINFO

      DO II = 1, KW
        IF (IPIV(KSTART-1+II) .NE. II) NS = NS + 1
        MP = MIN(MP, ABS(A(KSTART-1+II,KSTART-1+II)))
      END DO

      ! 2) Propager les pivots (comme l'etape 1)
      IF (KSTART .GT. 1) THEN
        CALL CLASWP(KSTART-1, A(KSTART,1), LDA, 1, KW, IPIV(KSTART), 1)
      END IF

      IF (KEND .LT. N) THEN
        CALL CLASWP(N-KEND, A(KSTART,KEND+1), LDA, 1, KW, IPIV(KSTART), 1)

        ! 3) U12 = L11^-1 * A12 (dense, comme l'etape 1)
        CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', &
                    KW, N-KEND, CONE, A(KSTART,KSTART), LDA, &
                    A(KSTART,KEND+1), LDA)

        ! MKL a 1 thread/appel pour la duree des boucles paralleles
        ! OpenMP ci-dessous (voir note de tete de fichier).
        CALL MKL_SET_NUM_THREADS(1)

        ! 4) Compresser U12, bloc de colonnes trainantes par bloc de
        !    colonnes trainantes -- chaque iteration J ecrit dans un
        !    element DIFFERENT de UBLK (tableau), aucun recouvrement.
        TC0 = OMP_GET_WTIME()
        !$OMP PARALLEL DO PRIVATE(J, NJ) &
        !$OMP   REDUCTION(+:N_OFFDIAG_BLOCKS, N_COMPRESSED, RANK_SUM) &
        !$OMP   SCHEDULE(DYNAMIC)
        DO J = K+1, NBLK
          NJ = BEND(J) - BSTART(J) + 1
          CALL COMPRESS_BLOCK(KW, NJ, A(KSTART,BSTART(J)), LDA, TOL, UBLK(J), &
                               WSARR(OMP_GET_THREAD_NUM()+1))
          N_OFFDIAG_BLOCKS = N_OFFDIAG_BLOCKS + 1
          IF (UBLK(J)%IS_LOWRANK) THEN
            N_COMPRESSED = N_COMPRESSED + 1
            RANK_SUM = RANK_SUM + UBLK(J)%RANK
          ELSE
            RANK_SUM = RANK_SUM + MIN(KW,NJ)
          END IF
        END DO
        !$OMP END PARALLEL DO

        ! 5) Compresser L21 (deja present dans A, issu de CGETRF), bloc de
        !    lignes trainantes par bloc de lignes trainantes -- chaque
        !    iteration I ecrit dans un element DIFFERENT de LBLK.
        !$OMP PARALLEL DO PRIVATE(I, MI) &
        !$OMP   REDUCTION(+:N_OFFDIAG_BLOCKS, N_COMPRESSED, RANK_SUM) &
        !$OMP   SCHEDULE(DYNAMIC)
        DO I = K+1, NBLK
          MI = BEND(I) - BSTART(I) + 1
          CALL COMPRESS_BLOCK(MI, KW, A(BSTART(I),KSTART), LDA, TOL, LBLK(I), &
                               WSARR(OMP_GET_THREAD_NUM()+1))
          N_OFFDIAG_BLOCKS = N_OFFDIAG_BLOCKS + 1
          IF (LBLK(I)%IS_LOWRANK) THEN
            N_COMPRESSED = N_COMPRESSED + 1
            RANK_SUM = RANK_SUM + LBLK(I)%RANK
          ELSE
            RANK_SUM = RANK_SUM + MIN(MI,KW)
          END IF
        END DO
        !$OMP END PARALLEL DO
        TC = TC + (OMP_GET_WTIME() - TC0)

        ! 6) Mise a jour de Schur, bloc (I,J) par bloc (I,J), en tirant
        !    parti des representations basse-rang qui viennent d'etre
        !    calculees -- chaque paire (I,J) ecrit dans une zone de A
        !    DISJOINTE de toutes les autres paires (blocs differents de
        !    la sous-matrice trainante), donc parallelisable sans section
        !    critique. COLLAPSE(2) fusionne les deux boucles pour mieux
        !    repartir la charge (NBLK^2 taches, tailles heterogenes selon
        !    quels blocs sont compresses ou non -- SCHEDULE(DYNAMIC)
        !    plutot que STATIC pour cette raison).
        !$OMP PARALLEL DO PRIVATE(I, J, MI, NJ) COLLAPSE(2) SCHEDULE(DYNAMIC)
        DO J = K+1, NBLK
          DO I = K+1, NBLK
            MI = BEND(I) - BSTART(I) + 1
            NJ = BEND(J) - BSTART(J) + 1
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LBLK(I), UBLK(J), &
                                     A(BSTART(I),BSTART(J)), LDA)
          END DO
        END DO
        !$OMP END PARALLEL DO

        ! Restaure MKL multi-thread pour le CGETRF/CTRSM du panneau suivant.
        CALL MKL_SET_NUM_THREADS(OMP_GET_MAX_THREADS())
      END IF
    END DO

    DO TID = 1, NTHREADS
      CALL BLR_WORKSPACE_FREE(WSARR(TID))
    END DO
    DEALLOCATE(WSARR)
    DEALLOCATE(BSTART, BEND, LBLK, UBLK)

    IF (PRESENT(N_SWAPS))   N_SWAPS   = NS
    IF (PRESENT(MIN_PIVOT)) MIN_PIVOT = MP
    IF (PRESENT(T_COMPRESS)) T_COMPRESS = TC
  END SUBROUTINE BLR_LU_FACTOR

END MODULE BLR_LU_MOD
