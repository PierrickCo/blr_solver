!===============================================================================
! ETAPE 1 du solveur BLR "from scratch" : factorisation LU PAR BLOCS, SANS
! compression (elle arrive a l'etape 2), complexe SIMPLE PRECISION.
!
! PRINCIPE : c'est exactement l'algorithme "right-looking blocked LU" que
! LAPACK utilise en interne dans sa version bloquee de CGETRF -- on le
! reecrit nous-memes ici (au lieu d'appeler CGETRF directement sur la
! matrice entiere) car c'est le squelette sur lequel on va greffer, a
! l'etape suivante, la compression BLR des blocs hors-diagonale dans la
! mise a jour de Schur (etape 5 ci-dessous) : on ne peut pas intercepter ce
! qui se passe DANS un CGETRF, donc on reimplemente l'orchestration par
! blocs nous-memes, en s'appuyant sur les noyaux denses de LAPACK/MKL
! (CGETRF, CLASWP, CTRSM, CGEMM) pour ce qui se passe A L'INTERIEUR de
! chaque bloc/panneau -- comme le font MUMPS et STRUMPACK (ils n'ecrivent
! pas leur propre BLAS non plus).
!
! DECOUPAGE EN BLOCS : taille FIXE NB (derniere colonne de blocs
! eventuellement plus petite si N n'est pas multiple de NB).
!
! PIVOT : partiel, CONFINE AU PANNEAU COURANT (colonnes du bloc K, TOUTES
! les lignes restantes K..N) -- c'est le pivot standard d'une LU bloquee,
! PAS un pivot global sur toute la matrice a chaque etape scalaire (ce qui
! serait couteux et inutile), mais PAS NON PLUS confine a un tout petit
! bloc NBxNB isole (ce qui degraderait la stabilite numerique sans raison,
! puisqu'on a acces a toutes les lignes du panneau). C'est le meme choix
! que fait MUMPS/STRUMPACK au niveau d'un front. Les echanges de lignes
! sont ensuite propages a TOUTE la matrice (blocs deja factorises a
! gauche ET blocs pas encore touches a droite), exactement comme le fait
! CGETRF en interne.
!
! ALGORITHME (bloc K = colonnes [kstart:kend], kw = kend-kstart+1) :
!   1) Factoriser le panneau [kstart:N, kstart:kend] avec CGETRF (pivot
!      partiel sur ce panneau).
!   2) Appliquer les echanges de lignes obtenus aux colonnes A GAUCHE
!      (1:kstart-1, deja factorisees) ET A DROITE (kend+1:N, pas encore
!      touchees) du panneau, via CLASWP.
!   3) Resoudre L11*U12 = A12 (CTRSM) pour obtenir le bloc-ligne U de la
!      partie a droite du panneau.
!   4) Mettre a jour la sous-matrice restante : A22 <- A22 - L21*U12
!      (CGEMM) -- C'EST CETTE ETAPE que l'etape 2 remplacera par une
!      version basse-rang quand les blocs L21/U12 seront compresses.
!
! A la sortie, A contient L (strictement inferieure, diagonale unite
! implicite) et U (triangulaire superieure, diagonale incluse) -- comme
! CGETRF. IPIV(1:N) est au format LAPACK standard (IPIV(i) : a l'etape i,
! la ligne i a ete echangee avec la ligne IPIV(i)) -- directement
! compatible avec CLASWP/CGETRS.
!
! PAS DE DIRECTIVE OPENMP ICI (etape 5) : contrairement a blr_lu_mod.f90
! (etape 2), la mise a jour de Schur de CE module est UN SEUL GROS CGEMM
! par panneau (pas de compression -> pas de multiples petits blocs
! independants a repartir nous-memes) -- ce gros appel est DEJA parallelise
! EN INTERNE par MKL (threading BLAS standard). Rien a paralleliser a notre
! niveau ; le parallelisme vient entierement de la MKL threadee.
!===============================================================================
MODULE BLOCK_LU_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! Factorisation LU par blocs, sans compression.
  !
  ! N_SWAPS/MIN_PIVOT (optionnels) : diagnostics de pivot -- N_SWAPS = nombre
  ! d'echanges de lignes REELLEMENT effectues (IPIV(i)/=i, au sens local du
  ! panneau -- voir BLOCK_LU_SOLVE plus bas pour cette convention), MIN_PIVOT
  ! = plus petit |pivot| (element diagonal de U) rencontre sur toute la
  ! factorisation. Sert a VERIFIER qu'un test sollicite reellement le pivot
  ! partiel (matrice a diagonale dominante = souvent N_SWAPS~0, peu
  ! probant) plutot que de le supposer -- cf. test/test_pivoting.f90.
  !-----------------------------------------------------------------------
  SUBROUTINE BLOCK_LU_FACTOR(N, NB, A, LDA, IPIV, INFO, N_SWAPS, MIN_PIVOT)
    INTEGER, INTENT(IN)             :: N, NB, LDA
    COMPLEX(WP), INTENT(INOUT)      :: A(LDA,N)
    INTEGER, INTENT(OUT)            :: IPIV(N)
    INTEGER, INTENT(OUT)            :: INFO
    INTEGER, INTENT(OUT), OPTIONAL  :: N_SWAPS
    REAL(WP), INTENT(OUT), OPTIONAL :: MIN_PIVOT

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP, 0.0_WP)
    INTEGER :: KSTART, KEND, KW, IINFO, II, NS
    REAL(WP) :: MP

    INFO = 0
    NS = 0
    MP = HUGE(1.0_WP)
    KSTART = 1
    DO WHILE (KSTART .LE. N)
      KEND = MIN(KSTART + NB - 1, N)
      KW   = KEND - KSTART + 1

      ! 1) Factoriser le panneau [kstart:N, kstart:kend]
      CALL CGETRF(N-KSTART+1, KW, A(KSTART,KSTART), LDA, IPIV(KSTART), IINFO)
      IF (IINFO .NE. 0 .AND. INFO .EQ. 0) INFO = KSTART - 1 + IINFO

      ! Statistiques de pivot pour ce panneau (indices locaux, cf. IPIV(i')
      ! vaut i' si aucun echange n'a ete necessaire pour la ligne locale i')
      DO II = 1, KW
        IF (IPIV(KSTART-1+II) .NE. II) NS = NS + 1
        MP = MIN(MP, ABS(A(KSTART-1+II,KSTART-1+II)))
      END DO

      ! 2) Propager les echanges de lignes aux colonnes A GAUCHE (deja
      !    factorisees) du panneau
      IF (KSTART .GT. 1) THEN
        CALL CLASWP(KSTART-1, A(KSTART,1), LDA, 1, KW, IPIV(KSTART), 1)
      END IF

      IF (KEND .LT. N) THEN
        ! 2bis) Propager les echanges de lignes aux colonnes A DROITE
        !       (pas encore touchees)
        CALL CLASWP(N-KEND, A(KSTART,KEND+1), LDA, 1, KW, IPIV(KSTART), 1)

        ! 3) Resoudre L11 * U12 = A12
        CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', &
                    KW, N-KEND, CONE, A(KSTART,KSTART), LDA, &
                    A(KSTART,KEND+1), LDA)

        ! 4) Mise a jour de Schur : A22 <- A22 - L21*U12
        CALL CGEMM('No transpose', 'No transpose', N-KEND, N-KEND, KW, &
                    -CONE, A(KEND+1,KSTART), LDA, A(KSTART,KEND+1), LDA, &
                    CONE, A(KEND+1,KEND+1), LDA)
      END IF

      KSTART = KEND + 1
    END DO

    IF (PRESENT(N_SWAPS))   N_SWAPS   = NS
    IF (PRESENT(MIN_PIVOT)) MIN_PIVOT = MP
  END SUBROUTINE BLOCK_LU_FACTOR

  !-----------------------------------------------------------------------
  ! Resolution A*X = B a partir des facteurs L/U et IPIV calcules par
  ! BLOCK_LU_FACTOR. B est ecrase par X.
  !
  ! IMPORTANT (bug corrige ici) : les indices de pivot stockes dans
  ! IPIV(KSTART:KEND) par BLOCK_LU_FACTOR sont LOCAUX A CHAQUE PANNEAU
  ! (CGETRF, appele sur le sous-tableau A(KSTART:N,...), n'a aucune notion
  ! de KSTART et renvoie des indices relatifs a SA PROPRE ligne 1 = notre
  ! ligne globale KSTART). Ca fonctionne directement pour le premier
  ! panneau (KSTART=1, donc local=global) mais PAS pour les suivants -- un
  ! unique CLASWP(NRHS,B,LDB,1,N,IPIV,1) sur tout B interpreterait a tort
  ! les valeurs locales des panneaux K>1 comme des indices globaux. On
  ! applique donc les pivots PANNEAU PAR PANNEAU, dans le meme ordre et
  ! avec la meme convention (sous-tableau demarrant a KSTART) que
  ! BLOCK_LU_FACTOR -- exactement les memes appels CLASWP que ceux deja
  ! valides lors de la factorisation.
  !-----------------------------------------------------------------------
  SUBROUTINE BLOCK_LU_SOLVE(N, NB, NRHS, A, LDA, IPIV, B, LDB)
    INTEGER, INTENT(IN)        :: N, NB, NRHS, LDA, LDB
    COMPLEX(WP), INTENT(IN)    :: A(LDA,N)
    INTEGER, INTENT(IN)        :: IPIV(N)
    COMPLEX(WP), INTENT(INOUT) :: B(LDB,NRHS)

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP, 0.0_WP)
    INTEGER :: KSTART, KEND, KW

    ! P*B : memes echanges de lignes que lors de la factorisation,
    ! panneau par panneau (voir note ci-dessus).
    KSTART = 1
    DO WHILE (KSTART .LE. N)
      KEND = MIN(KSTART + NB - 1, N)
      KW   = KEND - KSTART + 1
      CALL CLASWP(NRHS, B(KSTART,1), LDB, 1, KW, IPIV(KSTART), 1)
      KSTART = KEND + 1
    END DO

    ! L*Y = P*B  (L unitaire inferieure)
    CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', &
                N, NRHS, CONE, A, LDA, B, LDB)
    ! U*X = Y
    CALL CTRSM('Left', 'Upper', 'No transpose', 'Non-unit', &
                N, NRHS, CONE, A, LDA, B, LDB)
  END SUBROUTINE BLOCK_LU_SOLVE

END MODULE BLOCK_LU_MOD
