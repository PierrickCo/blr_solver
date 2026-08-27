!===============================================================================
! PROPOSITION 2 (post-etape 6) : stocker L21/U12 COMPRESSES (pas seulement
! accelerer la mise a jour de Schur pendant la factorisation, comme aux
! etapes 2/6e) -- objectif : reduire la memoire du RESULTAT stocke, pas
! seulement le cout de calcul de la factorisation.
!
! POURQUOI C'EST PLUS SIMPLE/SUR QUE LA PROPOSITION 1 (LDLt+BLR, abandonnee
! -- voir blr_ldlt_mod.f90) : le pivot de la LU est DEJA entierement
! resolu PENDANT la factorisation, via la propagation CLASWP deja validee
! (etapes 1/2/6e) -- pas de schema de mise a jour DIFFEREE a deviner comme
! CLASYF. On peut donc se permettre une strategie en DEUX PHASES bien plus
! sure :
!   PHASE 1 : appeler BLR_LU_FACTOR (etape 2, INCHANGEE) telle quelle --
!     produit A dense + IPIV, EXACTEMENT comme avant, aucune modification
!     de l'algorithme de factorisation lui-meme.
!   PHASE 2 : une fois A ENTIEREMENT factorisee (donc CONSISTENTE -- plus
!     aucune propagation de pivot future ne viendra la modifier), on la
!     PARCOURT bloc par bloc et on COMPRESSE chaque bloc L21/U12 dans une
!     structure PERSISTANTE (BLR_FACTORED_LU_T). Ceci evite un piege
!     rencontre en essayant la meme idee pour LDLt : compresser un bloc
!     "au fil de l'eau" PENDANT la factorisation risquerait de le rendre
!     perime des qu'un pivot ULTERIEUR est propage en arriere dans ses
!     colonnes -- en compressant seulement APRES la fin complete de la
!     factorisation, ce risque n'existe pas.
!
! CONSEQUENCE : la memoire DENSE (A, taille N x N) reste necessaire
! PENDANT la factorisation (comme etapes 2/6e) -- le gain memoire de cette
! proposition s'applique au STOCKAGE FINAL (BLR_FACTORED_LU_T, que
! l'appelant garde a la place de A une fois la Phase 2 terminee), pas au
! pic memoire de calcul. C'est la meme logique que MUMPS/STRUMPACK en
! mode "BLR partout" (facto avec espace de travail dense, seul le
! stockage final est compresse).
!
! SOLVE PAR BLOCS : substitution avant/arriere bloc par bloc (pas un seul
! CTRSM sur toute la matrice comme BLOCK_LU_SOLVE, puisqu'il n'y a plus de
! tableau N x N dense a passer a CTRSM) -- chaque bloc diagonal reste
! dense (CTRSM classique), chaque contribution hors-diagonale utilise
! BLR_BLOCK_APPLY (blr_compress_mod.f90, PROPOSITION 2 aussi) qui evite de
! reconstruire un bloc basse-rang en dense pour un simple produit
! matrice-vecteur -- c'est ICI que la compression profite aussi au TEMPS
! DE RESOLUTION, pas seulement a la memoire.
!===============================================================================
MODULE BLR_LU_STORAGE_MOD
  USE BLOCK_LU_MOD
  USE BLR_LU_MOD
  IMPLICIT NONE

  TYPE :: DENSE_PANEL_T
    COMPLEX(WP), ALLOCATABLE :: D(:,:)
  END TYPE DENSE_PANEL_T

  TYPE :: BLR_FACTORED_LU_T
    INTEGER :: N = 0, NB = 0, NBLK = 0
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    TYPE(DENSE_PANEL_T), ALLOCATABLE :: DIAGBLK(:)      ! DIAGBLK(K)%D : KWxKW, format CGETRF
    TYPE(BLR_BLOCK_T), ALLOCATABLE   :: LBLK(:,:)       ! LBLK(I,K), I>K : bloc L(I,K)
    TYPE(BLR_BLOCK_T), ALLOCATABLE   :: UBLK(:,:)       ! UBLK(K,J), J>K : bloc U(K,J)
  END TYPE BLR_FACTORED_LU_T

CONTAINS

  !-----------------------------------------------------------------------
  ! Phase 1 (BLR_LU_FACTOR, inchangee) + Phase 2 (compression persistante
  ! bloc par bloc du resultat final). A est ECRASE par la factorisation
  ! dense intermediaire (comme BLR_LU_FACTOR) -- l'appelant peut la
  ! liberer une fois cette routine terminee, seul FACT est necessaire pour
  ! resoudre ensuite (BLR_LU_SOLVE_COMPRESSED).
  !-----------------------------------------------------------------------
  ! T_FACTOR/T_COMPRESS (optionnels) : decoupage du temps total en "calcul
  ! dense pur" (CGETRF/CTRSM/CGEMM, Phase 1 hors compression) et
  ! "compression basse-rang" (SVD tronquee, Phase 1 -- mise a jour de
  ! Schur -- ET Phase 2 -- passe persistante -- CUMULEES) -- demande
  ! explicite de l'utilisateur pour voir ou va le temps.
  SUBROUTINE BLR_LU_FACTOR_COMPRESSED(N, NB, A, LDA, TOL, FACT, INFO, &
                                       N_COMPRESSED, N_OFFDIAG_BLOCKS, RANK_SUM, &
                                       N_SWAPS, MIN_PIVOT, T_FACTOR, T_COMPRESS)
    INTEGER, INTENT(IN)             :: N, NB, LDA
    COMPLEX(WP), INTENT(INOUT)      :: A(LDA,N)
    REAL(WP), INTENT(IN)            :: TOL
    TYPE(BLR_FACTORED_LU_T), INTENT(OUT) :: FACT
    INTEGER, INTENT(OUT)            :: INFO
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_OFFDIAG_BLOCKS
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL  :: N_SWAPS
    REAL(WP), INTENT(OUT), OPTIONAL :: MIN_PIVOT
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_FACTOR, T_COMPRESS

    INTEGER :: K, I, J, KW, MI, NJ
    DOUBLE PRECISION :: T0, T1, TC_PHASE1, TC_PHASE2
    TYPE(BLR_WORKSPACE_T) :: WS

    FACT%N = N
    FACT%NB = NB
    ALLOCATE(FACT%IPIV(N))

    ! --- Phase 1 : factorisation dense, deja validee, INCHANGEE ---
    T0 = OMP_GET_WTIME()
    CALL BLR_LU_FACTOR(N, NB, A, LDA, FACT%IPIV, TOL, INFO, &
                        N_COMPRESSED, N_OFFDIAG_BLOCKS, RANK_SUM, N_SWAPS, MIN_PIVOT, &
                        T_COMPRESS=TC_PHASE1)
    T1 = OMP_GET_WTIME()
    IF (PRESENT(T_FACTOR)) T_FACTOR = (T1-T0) - TC_PHASE1
    IF (INFO .NE. 0) RETURN

    ! --- Phase 2 : compression persistante bloc par bloc du resultat final ---
    CALL BLOCK_PARTITION(N, NB, FACT%NBLK, FACT%BSTART, FACT%BEND)
    ALLOCATE(FACT%DIAGBLK(FACT%NBLK))
    ALLOCATE(FACT%LBLK(FACT%NBLK,FACT%NBLK), FACT%UBLK(FACT%NBLK,FACT%NBLK))

    ! Workspace SVD reutilise pour tous les COMPRESS_BLOCK de cette boucle
    ! (sequentielle) au lieu d'un ALLOCATE/DEALLOCATE a chaque bloc L/U.
    CALL BLR_WORKSPACE_ALLOC(WS, NB)

    T0 = OMP_GET_WTIME()
    DO K = 1, FACT%NBLK
      KW = FACT%BEND(K) - FACT%BSTART(K) + 1
      ALLOCATE(FACT%DIAGBLK(K)%D(KW,KW))
      FACT%DIAGBLK(K)%D = A(FACT%BSTART(K):FACT%BEND(K), FACT%BSTART(K):FACT%BEND(K))

      DO I = K+1, FACT%NBLK
        MI = FACT%BEND(I) - FACT%BSTART(I) + 1
        CALL COMPRESS_BLOCK(MI, KW, A(FACT%BSTART(I),FACT%BSTART(K)), LDA, TOL, FACT%LBLK(I,K), WS)
      END DO
      DO J = K+1, FACT%NBLK
        NJ = FACT%BEND(J) - FACT%BSTART(J) + 1
        CALL COMPRESS_BLOCK(KW, NJ, A(FACT%BSTART(K),FACT%BSTART(J)), LDA, TOL, FACT%UBLK(K,J), WS)
      END DO
    END DO
    TC_PHASE2 = OMP_GET_WTIME() - T0
    CALL BLR_WORKSPACE_FREE(WS)

    IF (PRESENT(T_COMPRESS)) T_COMPRESS = TC_PHASE1 + TC_PHASE2
  END SUBROUTINE BLR_LU_FACTOR_COMPRESSED

  !-----------------------------------------------------------------------
  ! Resolution A*X=B a partir du stockage COMPRESSE (BLR_FACTORED_LU_T) --
  ! substitution par blocs, chaque contribution hors-diagonale via
  ! BLR_BLOCK_APPLY (evite l'expansion dense d'un bloc basse-rang).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_LU_SOLVE_COMPRESSED(FACT, NRHS, B, LDB)
    TYPE(BLR_FACTORED_LU_T), INTENT(IN) :: FACT
    INTEGER, INTENT(IN)                 :: NRHS, LDB
    COMPLEX(WP), INTENT(INOUT)          :: B(LDB,NRHS)

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    INTEGER :: K, I, J, KW

    ! --- P*B : memes echanges de lignes que la factorisation, panneau par
    !     panneau (meme convention/format IPIV que BLOCK_LU_SOLVE) ---
    DO K = 1, FACT%NBLK
      KW = FACT%BEND(K) - FACT%BSTART(K) + 1
      CALL CLASWP(NRHS, B(FACT%BSTART(K),1), LDB, 1, KW, FACT%IPIV(FACT%BSTART(K)), 1)
    END DO

    ! --- L*Y = P*B (substitution avant, par blocs) ---
    DO K = 1, FACT%NBLK
      KW = FACT%BEND(K) - FACT%BSTART(K) + 1
      DO I = 1, K-1
        CALL BLR_BLOCK_APPLY(FACT%LBLK(K,I), NRHS, -CONE, B(FACT%BSTART(I),1), LDB, &
                              CONE, B(FACT%BSTART(K),1), LDB)
      END DO
      CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NRHS, &
                  CONE, FACT%DIAGBLK(K)%D, KW, B(FACT%BSTART(K),1), LDB)
    END DO

    ! --- U*X = Y (substitution arriere, par blocs) ---
    DO K = FACT%NBLK, 1, -1
      KW = FACT%BEND(K) - FACT%BSTART(K) + 1
      DO J = K+1, FACT%NBLK
        CALL BLR_BLOCK_APPLY(FACT%UBLK(K,J), NRHS, -CONE, B(FACT%BSTART(J),1), LDB, &
                              CONE, B(FACT%BSTART(K),1), LDB)
      END DO
      CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', KW, NRHS, &
                  CONE, FACT%DIAGBLK(K)%D, KW, B(FACT%BSTART(K),1), LDB)
    END DO
  END SUBROUTINE BLR_LU_SOLVE_COMPRESSED

  !-----------------------------------------------------------------------
  ! Diagnostic : empreinte memoire REELLE (octets) du stockage compresse
  ! -- pour PROUVER le gain annonce par cette proposition, pas seulement
  ! l'esperer (meme discipline que les mesures de temps a l'etape 5).
  ! Complexe simple precision = 8 octets/element (2 x REAL 4 octets).
  !-----------------------------------------------------------------------
  FUNCTION BLR_FACTORED_LU_MEMORY_BYTES(FACT) RESULT(NBYTES)
    TYPE(BLR_FACTORED_LU_T), INTENT(IN) :: FACT
    INTEGER(8) :: NBYTES
    INTEGER :: K, I, J

    NBYTES = 0_8
    DO K = 1, FACT%NBLK
      NBYTES = NBYTES + 8_8 * INT(SIZE(FACT%DIAGBLK(K)%D), 8)
      DO I = K+1, FACT%NBLK
        NBYTES = NBYTES + BLR_BLOCK_BYTES(FACT%LBLK(I,K))
      END DO
      DO J = K+1, FACT%NBLK
        NBYTES = NBYTES + BLR_BLOCK_BYTES(FACT%UBLK(K,J))
      END DO
    END DO
  END FUNCTION BLR_FACTORED_LU_MEMORY_BYTES

END MODULE BLR_LU_STORAGE_MOD
