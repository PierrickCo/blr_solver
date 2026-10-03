!===============================================================================
! Pipeline "matrice compressee de bout en bout", en 5 etapes clairement
! separees :
!   1) Assemblage bloc par bloc, chaque bloc compresse immediatement apres
!      assemblage (pas de dense intermediaire garde en memoire) -- mesure
!      combinee assemblage+compression, et taux de compression de A.
!   2) Assemblage du second membre (mesure separee).
!   3) Factorisation LU directement sur la matrice compressee (mesure du
!      temps, taux de compression des facteurs L/U finaux).
!   4) Resolution.
!   5) Verification en appliquant la matrice compressee (pas une reference
!      dense recalculee a part).
!
! La matrice entiere ne passe jamais par une forme dense intermediaire
! hors des blocs individuels en cours de traitement -- compressee des
! l'assemblage, elle reste compressee pendant toute la factorisation.
!
! CHOIX DE PORTEE :
!   - Sans pivot (matrices a diagonale dominante). Ajouter le pivot sur
!     une representation compressee est une extension notable (voir
!     mpi_blr_lu_pivot_mod.f90 pour le meme saut de complexite deja fait
!     cote MPI).
!   - Blocs diagonaux restent toujours denses (pratique standard BLR,
!     necessaire pour la factorisation dense du pivot).
!   - Mise a jour de Schur d'un bloc hors-diagonale compresse : EXPAND
!     (dense transitoire) -> SCHUR_UPDATE_BLOCK -> COMPRESS_BLOCK a
!     nouveau (re-troncature). Une combinaison directe "basse-rang moins
!     basse-rang, re-tronque" sans repasser par le dense serait plus
!     rapide mais beaucoup plus delicate a implementer correctement (c'est
!     ce qui rend les vrais solveurs BLR/H-matrices difficiles) -- choix
!     delibere de simplicite/fiabilite.
!   - Transformation des panneaux L21/U12 eux-memes compresses : pas
!     d'expansion necessaire ici -- resoudre CTRSM sur le seul facteur
!     UFAC (ou VFAC) d'un bloc basse-rang donne directement le bon
!     resultat (L11^-1*(UFAC*VFAC) = (L11^-1*UFAC)*VFAC), un vrai gain
!     (CTRSM sur une matrice KWxRANG au lieu de KWxN) -- pas une
!     simplification, un benefice reel de garder la donnee compressee.
!
! GENERALISATION (librairie d'algebre lineaire bloc full/low-rank) : le
! facteur diagonal sans pivot, les TRSM appliques a un facteur low-rank
! seul, le decoupage en blocs et le calcul de taille memoire d'un bloc
! sont les primitives partagees de blr_compress_mod.f90 (BLR_LU_NOPIV_
! UNBLOCKED, BLR_TRSM_LEFT_LOWER_UNIT/BLR_TRSM_RIGHT_UPPER_NONUNIT,
! BLOCK_PARTITION, BLR_BLOCK_BYTES) -- ce module ne les redefinit pas
! localement.
!===============================================================================
MODULE BLR_MATRIX_MOD
  USE BLR_COMPRESS_MOD
  USE PROGRESS_BAR_MOD
  USE THREAD_CONTROL_MOD, ONLY: SET_NUM_THREADS
  USE OMP_LIB
  USE MKL_INTERFACES_MOD, ONLY: CGEMM, CTRSM
  IMPLICIT NONE

  TYPE :: DENSE_PANEL_T
    COMPLEX(WP), ALLOCATABLE :: D(:,:)
  END TYPE DENSE_PANEL_T

  TYPE :: BLR_MATRIX_T
    INTEGER :: N = 0, NB = 0, NBLK = 0
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    TYPE(DENSE_PANEL_T), ALLOCATABLE :: DIAGBLK(:)   ! toujours dense
    TYPE(BLR_BLOCK_T), ALLOCATABLE   :: OFFBLK(:,:)  ! hors-diagonale, compresse des l'assemblage
  END TYPE BLR_MATRIX_T

  ABSTRACT INTERFACE
    SUBROUTINE GET_BLOCK_IFACE(ISTART, IEND, JSTART, JEND, BLK)
      IMPORT :: WP
      INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
      COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    END SUBROUTINE GET_BLOCK_IFACE
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------
  ! Assemblage bloc par bloc, avec compression immediate de chaque bloc
  ! hors-diagonale (GET_BLOCK : callback fourni par l'appelant, ne calcule
  ! qu'un bloc a la fois -- jamais de tableau N x N dense complet en
  ! memoire). T_ASSEMBLE_COMPRESS : temps combine (assemblage+compression).
  !-----------------------------------------------------------------------
  SUBROUTINE ASSEMBLE_BLR_MATRIX(N, NB, TOL, GET_BLOCK, A_BLR, &
                                  T_ASSEMBLE_COMPRESS, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                  SHOW_PROGRESS, USE_RSVD)
    INTEGER, INTENT(IN)          :: N, NB
    REAL(WP), INTENT(IN)         :: TOL
    PROCEDURE(GET_BLOCK_IFACE)   :: GET_BLOCK
    TYPE(BLR_MATRIX_T), INTENT(OUT) :: A_BLR
    DOUBLE PRECISION, INTENT(OUT) :: T_ASSEMBLE_COMPRESS
    INTEGER, INTENT(OUT)         :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)      :: RANK_SUM
    LOGICAL, INTENT(IN), OPTIONAL :: SHOW_PROGRESS
    ! USE_RSVD (defaut .TRUE.) : SVD randomisee adaptative ; .FALSE. =
    ! CGESVD avec workspace par thread. RANDOM_NUMBER est thread-safe
    ! (ifx/gfortran), suites differentes d'un thread a l'autre.
    LOGICAL, INTENT(IN), OPTIONAL :: USE_RSVD

    DOUBLE PRECISION :: T0
    LOGICAL :: RSVD_FLAG
    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: I, J, MI, NJ, NDONE, NTOTAL, NTHREADS, TID, PBAR_STATE
    LOGICAL :: PROGRESS
    TYPE(BLR_WORKSPACE_T), ALLOCATABLE :: WSARR(:)

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS
    RSVD_FLAG = .TRUE.
    IF (PRESENT(USE_RSVD)) RSVD_FLAG = USE_RSVD

    ! La double boucle I,J est independante (chaque (I,J) ecrit dans une
    ! case disjointe de DIAGBLK/OFFBLK) -- un workspace SVD par thread
    ! (jamais partage, voir la note de thread-safety de BLR_WORKSPACE_T).
    ! Suppose GET_BLOCK thread-safe (formule pure, sans etat partage).
    NTHREADS = OMP_GET_MAX_THREADS()
    ALLOCATE(WSARR(NTHREADS))
    DO TID = 1, NTHREADS
      CALL BLR_WORKSPACE_ALLOC(WSARR(TID), NB)
    END DO

    T0 = OMP_GET_WTIME()
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%N = N
    A_BLR%NB = NB
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(A_BLR%NBLK, A_BLR%NBLK))

    N_COMPRESSED = 0
    N_OFFDIAG = 0
    RANK_SUM = 0_8
    NDONE = 0
    PBAR_STATE = -1
    NTOTAL = A_BLR%NBLK * A_BLR%NBLK

    ! NDONE reste PARTAGE (pas en REDUCTION) : un compteur de progression en
    ! REDUCTION serait prive par thread pendant la boucle et ne refleterait
    ! pas le total reel avant la toute fin -- ici il est incremente et lu
    ! DANS la meme section critique que l'affichage, pour un pourcentage
    ! exact malgre le parallelisme.
    !$OMP PARALLEL DO PRIVATE(I, J, MI, NJ, TMP) SHARED(RSVD_FLAG) COLLAPSE(2) &
    !$OMP   REDUCTION(+:N_OFFDIAG, N_COMPRESSED, RANK_SUM) SCHEDULE(DYNAMIC)
    DO I = 1, A_BLR%NBLK
      DO J = 1, A_BLR%NBLK
        MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
        NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
        IF (I .EQ. J) THEN
          ALLOCATE(A_BLR%DIAGBLK(I)%D(MI,MI))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), A_BLR%DIAGBLK(I)%D)
        ELSE
          ALLOCATE(TMP(MI,NJ))
          CALL GET_BLOCK(A_BLR%BSTART(I), A_BLR%BEND(I), A_BLR%BSTART(J), A_BLR%BEND(J), TMP)
          IF (RSVD_FLAG) THEN
            CALL COMPRESS_BLOCK_RSVD_ADAPTIVE(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(I,J))
          ELSE
            CALL COMPRESS_BLOCK(MI, NJ, TMP, MI, TOL, A_BLR%OFFBLK(I,J), WSARR(OMP_GET_THREAD_NUM()+1))
          END IF
          N_OFFDIAG = N_OFFDIAG + 1
          IF (A_BLR%OFFBLK(I,J)%IS_LOWRANK) THEN
            N_COMPRESSED = N_COMPRESSED + 1
            RANK_SUM = RANK_SUM + A_BLR%OFFBLK(I,J)%RANK
          ELSE
            RANK_SUM = RANK_SUM + MIN(MI,NJ)
          END IF
          DEALLOCATE(TMP)
        END IF
        !$OMP CRITICAL
        NDONE = NDONE + 1
        IF (PROGRESS) CALL PRINT_PROGRESS_BAR('1) Assemblage+compression', NDONE, NTOTAL, PBAR_STATE)
        !$OMP END CRITICAL
      END DO
    END DO
    !$OMP END PARALLEL DO
    T_ASSEMBLE_COMPRESS = OMP_GET_WTIME() - T0

    DO TID = 1, NTHREADS
      CALL BLR_WORKSPACE_FREE(WSARR(TID))
    END DO
    DEALLOCATE(WSARR)
  END SUBROUTINE ASSEMBLE_BLR_MATRIX

  !-----------------------------------------------------------------------
  ! Empreinte memoire reelle de A_BLR (octets), pour caracteriser le taux
  ! de compression de la matrice assemblee.
  !-----------------------------------------------------------------------
  FUNCTION BLR_MATRIX_MEMORY_BYTES(A_BLR) RESULT(NBYTES)
    TYPE(BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER(8) :: NBYTES
    INTEGER :: I, J
    NBYTES = 0_8
    DO I = 1, A_BLR%NBLK
      NBYTES = NBYTES + 8_8*INT(SIZE(A_BLR%DIAGBLK(I)%D),8)
      DO J = 1, A_BLR%NBLK
        IF (I .NE. J) NBYTES = NBYTES + BLR_BLOCK_BYTES(A_BLR%OFFBLK(I,J))
      END DO
    END DO
  END FUNCTION BLR_MATRIX_MEMORY_BYTES

  !-----------------------------------------------------------------------
  ! Factorisation LU en place sur A_BLR, qui reste compresse tout du long
  ! (sauf expansion transitoire d'un bloc hors-diagonale, le temps de sa
  ! propre mise a jour de Schur -- voir note de tete de fichier). Sans
  ! pivot. N_COMPRESSED/N_OFFDIAG/RANK_SUM en sortie caracterisent les
  ! facteurs L/U finaux (releves au moment ou chaque bloc atteint sa forme
  ! definitive, pas les etats intermediaires).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_MATRIX_LU_FACTOR(A_BLR, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, T_FACTOR, &
                                   SHOW_PROGRESS)
    TYPE(BLR_MATRIX_T), INTENT(INOUT) :: A_BLR
    REAL(WP), INTENT(IN)              :: TOL
    INTEGER, INTENT(OUT)              :: INFO
    INTEGER, INTENT(OUT)              :: N_COMPRESSED, N_OFFDIAG
    INTEGER(8), INTENT(OUT)           :: RANK_SUM
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR
    LOGICAL, INTENT(IN), OPTIONAL     :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    INTEGER :: NBLK, K, I, J, KW, MI, NJ, IINFO, NTHREADS, TID, PBAR_STATE
    DOUBLE PRECISION :: T0
    LOGICAL :: PROGRESS
    TYPE(BLR_WORKSPACE_T), ALLOCATABLE :: WSARR(:)

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS
    PBAR_STATE = -1

    ! Un workspace QR+SVD par thread (mise a jour de Schur parallelisee
    ! ci-dessous) -- NB suffit comme taille : MI,NJ<=NB et le rang
    ! recompresse RTOT<=RMAX<NB.
    NTHREADS = OMP_GET_MAX_THREADS()
    ALLOCATE(WSARR(NTHREADS))
    DO TID = 1, NTHREADS
      CALL BLR_WORKSPACE_ALLOC(WSARR(TID), A_BLR%NB)
    END DO

    T0 = OMP_GET_WTIME()
    NBLK = A_BLR%NBLK
    INFO = 0
    N_COMPRESSED = 0
    N_OFFDIAG = 0
    RANK_SUM = 0_8

    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1

      ! 1) facteur diagonal, dense, sans pivot -- appel unique et
      !    sequentiel (une elimination sur un seul bloc KWxKW ne se
      !    decoupe pas en taches independantes), mais MKL tourne a son
      !    nombre de threads max pendant cet appel (le seul moyen de le
      !    paralleliser, une elimination scalaire n'offrant pas de
      !    parallelisme de tache) -- restaure a 1 juste avant les boucles
      !    OpenMP ci-dessous pour eviter le threading imbrique (threads
      !    OpenMP x threads MKL internes).
      CALL SET_NUM_THREADS(NTHREADS)
      CALL BLR_LU_NOPIV_UNBLOCKED(KW, A_BLR%DIAGBLK(K)%D, KW, IINFO)
      IF (IINFO .NE. 0 .AND. INFO .EQ. 0) INFO = A_BLR%BSTART(K) - 1 + IINFO
      CALL SET_NUM_THREADS(1)

      ! 2) U12(K,J) = L11^-1 * A(K,J) -- transforme UFAC (ou DENSE) SEUL,
      !    jamais besoin d'expanser un bloc basse-rang pour cette etape --
      !    PARALLELISE (chaque J ecrit un bloc disjoint).
      !$OMP PARALLEL DO PRIVATE(J) REDUCTION(+:N_OFFDIAG, N_COMPRESSED, RANK_SUM) SCHEDULE(DYNAMIC)
      DO J = K+1, NBLK
        CALL BLR_TRSM_LEFT_LOWER_UNIT(A_BLR%DIAGBLK(K)%D, KW, A_BLR%OFFBLK(K,J))
        N_OFFDIAG = N_OFFDIAG + 1
        IF (A_BLR%OFFBLK(K,J)%IS_LOWRANK) THEN
          N_COMPRESSED = N_COMPRESSED + 1
          RANK_SUM = RANK_SUM + A_BLR%OFFBLK(K,J)%RANK
        ELSE
          RANK_SUM = RANK_SUM + MIN(A_BLR%OFFBLK(K,J)%M, A_BLR%OFFBLK(K,J)%N)
        END IF
      END DO
      !$OMP END PARALLEL DO

      ! 3) L21(I,K) = A(I,K) * U11^-1 -- transforme VFAC (ou DENSE) SEUL --
      !    PARALLELISE (chaque I ecrit un bloc disjoint).
      !$OMP PARALLEL DO PRIVATE(I) REDUCTION(+:N_OFFDIAG, N_COMPRESSED, RANK_SUM) SCHEDULE(DYNAMIC)
      DO I = K+1, NBLK
        CALL BLR_TRSM_RIGHT_UPPER_NONUNIT(A_BLR%OFFBLK(I,K), A_BLR%DIAGBLK(K)%D, KW)
        N_OFFDIAG = N_OFFDIAG + 1
        IF (A_BLR%OFFBLK(I,K)%IS_LOWRANK) THEN
          N_COMPRESSED = N_COMPRESSED + 1
          RANK_SUM = RANK_SUM + A_BLR%OFFBLK(I,K)%RANK
        ELSE
          RANK_SUM = RANK_SUM + MIN(A_BLR%OFFBLK(I,K)%M, A_BLR%OFFBLK(I,K)%N)
        END IF
      END DO
      !$OMP END PARALLEL DO

      ! 4) mise a jour de Schur de la sous-matrice trainante -- PARALLELISE
      !    (chaque paire (I,J) est disjointe), un workspace SVD/QR PAR
      !    THREAD pour le chemin low-rank (FAST_RECOMPRESS, via
      !    SCHUR_UPDATE_BLR_BLOCK) qui garde le bloc compresse aussi
      !    souvent que possible -- voir le commentaire de tete de fichier.
      !$OMP PARALLEL DO PRIVATE(I, J, MI, NJ) COLLAPSE(2) SCHEDULE(DYNAMIC)
      DO J = K+1, NBLK
        DO I = K+1, NBLK
          MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
          NJ = A_BLR%BEND(J) - A_BLR%BSTART(J) + 1
          IF (I .EQ. J) THEN
            ! futur bloc diagonal : toujours dense, mise a jour directe
            CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, A_BLR%OFFBLK(I,K), A_BLR%OFFBLK(K,J), &
                                     A_BLR%DIAGBLK(I)%D, MI)
          ELSE
            ! bloc hors-diagonale : reste compresse si possible (voir
            ! SCHUR_UPDATE_BLR_BLOCK, blr_compress_mod.f90) -- evite les
            ! ~NBLK^3/3 SVD MIxNJ completes qui dominaient le cout de
            ! factorisation (835s a N=10000/NB=100 avant cette optimisation)
            CALL SCHUR_UPDATE_BLR_BLOCK(A_BLR%OFFBLK(I,J), MI, NJ, KW, &
                                         A_BLR%OFFBLK(I,K), A_BLR%OFFBLK(K,J), TOL, &
                                         WSARR(OMP_GET_THREAD_NUM()+1))
          END IF
        END DO
      END DO
      !$OMP END PARALLEL DO

      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('3) Factorisation LU', K, NBLK, PBAR_STATE)
    END DO

    CALL SET_NUM_THREADS(NTHREADS)
    IF (PRESENT(T_FACTOR)) T_FACTOR = OMP_GET_WTIME() - T0

    DO TID = 1, NTHREADS
      CALL BLR_WORKSPACE_FREE(WSARR(TID))
    END DO
    DEALLOCATE(WSARR)
  END SUBROUTINE BLR_MATRIX_LU_FACTOR

  !-----------------------------------------------------------------------
  ! Resolution A*X=B par substitution avant/arriere par blocs, directement
  ! sur les facteurs L21/U12 compresses (BLR_BLOCK_APPLY -- jamais
  ! d'expansion dense necessaire ici). B ecrase par X.
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_MATRIX_LU_SOLVE(A_BLR, NRHS, B, LDB, SHOW_PROGRESS)
    TYPE(BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)            :: NRHS, LDB
    COMPLEX(WP), INTENT(INOUT)     :: B(LDB,NRHS)
    LOGICAL, INTENT(IN), OPTIONAL  :: SHOW_PROGRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: NBLK, K, I, J, KW, NDONE, NTOTAL, NTHREADS, PBAR_STATE
    LOGICAL :: PROGRESS
    COMPLEX(WP), ALLOCATABLE :: YACC(:,:)

    PROGRESS = .FALSE.
    IF (PRESENT(SHOW_PROGRESS)) PROGRESS = SHOW_PROGRESS

    NBLK = A_BLR%NBLK
    NDONE = 0
    PBAR_STATE = -1
    NTOTAL = 2*NBLK
    NTHREADS = OMP_GET_MAX_THREADS()

    ! L*Y = B (substitution avant, L unitaire, diagonale implicite). La
    ! boucle externe sur K est sequentielle (K depend du resultat de K-1) ;
    ! la boucle interne sur I, elle, est parallelisee : les contributions
    ! des blocs deja resolus s'accumulent dans la meme case B(K,:), donc
    ! REDUCTION(+:YACC) plutot qu'ecrire directement dans B depuis
    ! plusieurs threads (chaque thread accumule dans sa copie privee de
    ! YACC, combinees par + a la fin du parallel do). BLR_BLOCK_APPLY
    ! exploite deja le low-rank en interne (2 petits GEMM via UFAC/VFAC) --
    ! la parallelisation en multiplie le benefice sur plusieurs coeurs. Ne
    ! vaut le cout de la region parallele que si K-1 est assez grand.
    DO K = 1, NBLK
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      IF (K .GT. 1) THEN
        ALLOCATE(YACC(KW,NRHS))
        YACC = CZERO
        ! MKL a 1 thread/appel pendant la boucle parallele (le CGEMM dense
        ! de BLR_BLOCK_APPLY, cas bloc-plein, est assez gros pour que MKL
        ! le threade lui-meme en interne -- meme garde-fou que
        ! BLR_MATRIX_LU_FACTOR) ; restaure juste avant le CTRSM sequentiel.
        CALL SET_NUM_THREADS(1)
        !$OMP PARALLEL DO PRIVATE(I) REDUCTION(+:YACC) SCHEDULE(DYNAMIC)
        DO I = 1, K-1
          CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(K,I), NRHS, -CONE, B(A_BLR%BSTART(I),1), LDB, &
                                CONE, YACC, KW)
        END DO
        !$OMP END PARALLEL DO
        CALL SET_NUM_THREADS(NTHREADS)
        B(A_BLR%BSTART(K):A_BLR%BEND(K),1:NRHS) = B(A_BLR%BSTART(K):A_BLR%BEND(K),1:NRHS) + YACC
        DEALLOCATE(YACC)
      END IF
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS, &
                  CONE, A_BLR%DIAGBLK(K)%D, KW, B(A_BLR%BSTART(K),1), LDB)
      NDONE = NDONE + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('4) Resolution', NDONE, NTOTAL, PBAR_STATE)
    END DO

    ! U*X = Y (substitution arriere) -- meme principe, boucle interne sur J.
    DO K = NBLK, 1, -1
      KW = A_BLR%BEND(K) - A_BLR%BSTART(K) + 1
      IF (K .LT. NBLK) THEN
        ALLOCATE(YACC(KW,NRHS))
        YACC = CZERO
        CALL SET_NUM_THREADS(1)
        !$OMP PARALLEL DO PRIVATE(J) REDUCTION(+:YACC) SCHEDULE(DYNAMIC)
        DO J = K+1, NBLK
          CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(K,J), NRHS, -CONE, B(A_BLR%BSTART(J),1), LDB, &
                                CONE, YACC, KW)
        END DO
        !$OMP END PARALLEL DO
        CALL SET_NUM_THREADS(NTHREADS)
        B(A_BLR%BSTART(K):A_BLR%BEND(K),1:NRHS) = B(A_BLR%BSTART(K):A_BLR%BEND(K),1:NRHS) + YACC
        DEALLOCATE(YACC)
      END IF
      CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS, &
                  CONE, A_BLR%DIAGBLK(K)%D, KW, B(A_BLR%BSTART(K),1), LDB)
      NDONE = NDONE + 1
      IF (PROGRESS) CALL PRINT_PROGRESS_BAR('4) Resolution', NDONE, NTOTAL, PBAR_STATE)
    END DO
  END SUBROUTINE BLR_MATRIX_LU_SOLVE

  !-----------------------------------------------------------------------
  ! Verification : Y = A_BLR * X, en utilisant directement la matrice
  ! compressee (avant factorisation -- l'appelant doit garder une copie de
  ! A_BLR non factorisee pour ceci, peu couteux puisque compresse).
  ! Reutilise BLR_BLOCK_APPLY pour les blocs hors-diagonale, CGEMM dense
  ! pour les blocs diagonaux.
  !-----------------------------------------------------------------------
  SUBROUTINE APPLY_BLR_MATRIX(A_BLR, NRHS, X, LDX, Y, LDY)
    TYPE(BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN)            :: NRHS, LDX, LDY
    COMPLEX(WP), INTENT(IN)        :: X(LDX,NRHS)
    COMPLEX(WP), INTENT(OUT)       :: Y(LDY,NRHS)

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: NBLK, I, J, MI

    NBLK = A_BLR%NBLK
    DO I = 1, NBLK
      MI = A_BLR%BEND(I) - A_BLR%BSTART(I) + 1
      CALL CGEMM('No transpose', 'No transpose', MI, NRHS, MI, CONE, A_BLR%DIAGBLK(I)%D, MI, &
                  X(A_BLR%BSTART(I),1), LDX, CZERO, Y(A_BLR%BSTART(I),1), LDY)
      DO J = 1, NBLK
        IF (J .NE. I) THEN
          CALL BLR_BLOCK_APPLY(A_BLR%OFFBLK(I,J), NRHS, CONE, X(A_BLR%BSTART(J),1), LDX, &
                                CONE, Y(A_BLR%BSTART(I),1), LDY)
        END IF
      END DO
    END DO
  END SUBROUTINE APPLY_BLR_MATRIX

END MODULE BLR_MATRIX_MOD
