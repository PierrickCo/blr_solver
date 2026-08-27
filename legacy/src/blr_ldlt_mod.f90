!===============================================================================
! *** NON FONCTIONNEL -- CONSERVE UNIQUEMENT POUR MEMOIRE, NE PAS UTILISER ***
!
! Echec confirme par experimentation (voir memoire du projet, section
! "Proposition 1") : CLASYF utilise un schema de mise a jour DIFFEREE
! (accumule le rang-KB du panneau et l'applique EN UNE SEULE FOIS a la fin,
! comme les autres routines bloquees de LAPACK) -- les valeurs
! intermediaires qu'on peut lire dans une copie apres l'appel NE
! permettent PAS de reconstruire, par une simple permutation composee +
! une seule mise a jour L21*D*L21^T, ce que CLASYF calcule reellement des
! qu'un panneau contient PLUSIEURS echanges de pivot (verifie a la main :
! ni l'ordre "avant" ni l'ordre "apres" pour composer les echanges ne
! reproduit le resultat de CLASYF au-dela d'un seul echange isole).
! Correctement exposer L21 avant la mise a jour demanderait de reimplanter
! a la main le schema bloque de Bunch-Kaufman lui-meme (accumulation
! differee incluse), pas seulement d'appeler CLASYF intelligemment --
! exactement le risque que le commentaire de tete de block_ldlt_mod.f90
! signalait des l'etape 4. Utilisateur informe, a choisi de passer aux
! propositions 2 et 3 plutot que de poursuivre cette voie.
!===============================================================================
!===============================================================================
! PROPOSITION 1 (post-etape 6) : compression BLR pour LDLt -- extension
! explicitement deferee a l'etape 4 (voir block_ldlt_mod.f90) car CLASYF
! fusionne factorisation de panneau ET mise a jour de Schur en un seul
! appel "boite noire", empechant d'intercepter L21 avant la mise a jour
! (contrairement a CGETRF, utilise SEUL pour le panneau dans block_lu_mod
! /blr_lu_mod, laissant l'orchestration -- et donc l'injection de
! compression -- a notre charge).
!
! TECHNIQUE POUR CONTOURNER LA BOITE NOIRE (sans reimplementer Bunch-
! Kaufman a la main -- trop risque, cf. la note de block_ldlt_mod.f90) :
! on appelle CLASYF sur une COPIE JETABLE de la sous-matrice active
! [KSTART:N,KSTART:N] (ACOPY). CLASYF y calcule -- exactement comme
! d'habitude -- le panneau factorise (L21, D par blocs 1x1/2x2 de Bunch-
! Kaufman, choix de pivot FIABLE car c'est LAPACK qui le fait) ET met a
! jour la sous-matrice restante DANS LA COPIE -- mise a jour qu'on IGNORE
! et RECALCULE nous-memes (avec compression BLR) sur le VRAI tableau A.
! On recupere ainsi le meilleur des deux mondes : le choix de pivot
! Bunch-Kaufman reste 100% LAPACK (aucun risque de bug de selection 1x1/
! 2x2), et on regagne le controle sur la mise a jour de Schur necessaire
! pour injecter la compression -- au prix d'une copie O(NREMAIN^2)
! supplementaire par panneau (cout reel, documente : total O(N^3/NB) sur
! toute la factorisation, meme ordre de grandeur que l'algorithme lui-meme
! mais avec un facteur constant en plus -- accepte pour la simplicite et
! la fiabilite de cette premiere version).
!
! RECONSTRUCTION DE LA MISE A JOUR : A22 <- A22 - L21*D*L21^T. On calcule
! DL21T = D*L21^T (KB x MTR) a la main (D est PAR BLOCS 1x1/2x2, lus
! directement dans ACOPY -- convention de stockage CSYTRF standard :
! ACOPY(j,j) = D(j,j) ; pour un bloc 2x2 en (j,j+1), ACOPY(j+1,j) = D(j+1,j)
! l'element hors-diagonale, PAS conjugue puisque LDLt ici est SYMETRIQUE,
! pas hermitienne, meme convention que partout ailleurs dans ce projet),
! puis on COMPRESSE L21 et DL21T separement (COMPRESS_BLOCK, etape 2) et on
! appelle SCHUR_UPDATE_BLOCK(MTR,MTR,KB,LBLK=L21,UBLK=DL21T,...) -- EXACTEMENT
! le meme dispatch dense/basse-rang que la LU, sans aucune modification :
! A22 -= L21 * DL21T = L21 * D * L21^T, la formule voulue.
!
! GRANULARITE -- MEME SIMPLIFICATION DELIBEREE QU'A L'ETAPE 6e (MPI) : on
! compresse L21 et DL21T COMME UN SEUL BLOC chacun (pas sous-decoupes en
! blocs NB x NB comme le fait blr_lu_mod.f90 pour la LU) -- plus simple,
! le rang numerique trouve reste borne par le rang global reel de la
! matrice quel que soit le decoupage. Consequence : la mise a jour est
! appliquee au CARRE COMPLET de la sous-matrice trainante (SCHUR_UPDATE_
! BLOCK n'a pas de variante "symetrique, triangle inferieur seulement"),
! alors que CSYTRF/CLASYF ne mettent a jour QUE le triangle inferieur --
! du travail perdu sur le triangle superieur (~2x flops sur cette etape),
! mais SANS CONSEQUENCE de justesse puisque BLOCK_LDLT_SOLVE (CSYTRS) ne
! lit jamais que le triangle inferieur en aval. Documente comme piste
! d'optimisation future (comme la granularite de compression en 6e).
!===============================================================================
MODULE BLR_LDLT_MOD
  USE BLR_COMPRESS_MOD
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! Factorisation A = L*D*L^T par blocs, AVEC compression BLR de la mise a
  ! jour de Schur. Meme convention de sortie (A, IPIV) que BLOCK_LDLT_FACTOR
  ! (etape 4) -- BLOCK_LDLT_SOLVE (module BLOCK_LDLT_MOD, INCHANGE) reste
  ! directement utilisable pour resoudre. TOL=0 doit redonner (a l'ordre
  ! de l'arrondi flottant pres) exactement l'etape 4.
  !
  ! N_2X2_PIVOTS (optionnel) : nombre de pivots 2x2 de Bunch-Kaufman
  ! rencontres -- diagnostic pour PROUVER que la conversion D bloc-2x2 ci-
  ! dessus est reellement exercisee par un test, pas seulement le cas 1x1
  ! (meme discipline que N_SWAPS/MIN_PIVOT ailleurs dans ce projet).
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_LDLT_FACTOR(N, NB, A, LDA, IPIV, TOL, INFO, &
                              N_COMPRESSED, N_PANELS, RANK_SUM, N_2X2_PIVOTS)
    INTEGER, INTENT(IN)             :: N, NB, LDA
    COMPLEX(WP), INTENT(INOUT)      :: A(LDA,N)
    INTEGER, INTENT(OUT)            :: IPIV(N)
    REAL(WP), INTENT(IN)            :: TOL
    INTEGER, INTENT(OUT)            :: INFO
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_PANELS
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL  :: N_2X2_PIVOTS

    COMPLEX(WP), ALLOCATABLE :: WORK(:,:), ACOPY(:,:), DL21T(:,:)
    TYPE(BLR_BLOCK_T) :: LBLK, UBLK
    COMPLEX(WP) :: D11, D21, D22
    INTEGER :: KSTART, NREMAIN, KB, IINFO, II, J, MTR, N2X2

    ALLOCATE(WORK(N,NB))
    INFO = 0
    N_COMPRESSED = 0
    N_PANELS = 0
    RANK_SUM = 0_8
    N2X2 = 0

    KSTART = 1
    DO WHILE (KSTART .LE. N)
      NREMAIN = N - KSTART + 1

      ! --- copie jetable : CLASYF y choisit le pivot Bunch-Kaufman et y
      !     calcule L21/D, mais sa propre mise a jour de la sous-matrice
      !     trainante (aussi ecrite dans ACOPY) est IGNOREE -- voir note
      !     de tete de fichier ---
      ALLOCATE(ACOPY(NREMAIN,NREMAIN))
      ACOPY = A(KSTART:N,KSTART:N)
      CALL CLASYF('Lower', NREMAIN, MIN(NB,NREMAIN), KB, &
                   ACOPY, NREMAIN, IPIV(KSTART), WORK, N, IINFO)
      IF (IINFO .NE. 0 .AND. INFO .EQ. 0) INFO = KSTART - 1 + IINFO

      ! Conversion locale -> globale des indices de pivot (meme convention
      ! signee que block_ldlt_mod.f90, voir sa note detaillee).
      DO II = KSTART, KSTART+KB-1
        IPIV(II) = IPIV(II) + SIGN(KSTART-1, IPIV(II))
      END DO

      ! --- ecrire le panneau reellement factorise (L21+D) dans le VRAI A ---
      A(KSTART:N, KSTART:KSTART+KB-1) = ACOPY(1:NREMAIN,1:KB)

      MTR = NREMAIN - KB
      IF (MTR .GT. 0) THEN
        ! --- construire DL21T = D*L21^T a la main, bloc 1x1 ou 2x2 par
        !     bloc 1x1 ou 2x2 (lus dans ACOPY, convention CSYTRF standard) ---
        ALLOCATE(DL21T(KB,MTR))
        J = 1
        DO WHILE (J .LE. KB)
          IF (IPIV(KSTART-1+J) .GT. 0) THEN
            ! pivot 1x1
            D11 = ACOPY(J,J)
            DL21T(J,1:MTR) = D11 * ACOPY(KB+1:NREMAIN,J)
            J = J + 1
          ELSE
            ! pivot 2x2, colonnes (J,J+1) -- IPIV(KSTART-1+J) et
            ! IPIV(KSTART-1+J+1) sont egaux et negatifs (convention Bunch-
            ! Kaufman standard, cf. doc CLASYF)
            N2X2 = N2X2 + 1
            D11 = ACOPY(J,J)
            D21 = ACOPY(J+1,J)
            D22 = ACOPY(J+1,J+1)
            DL21T(J,  1:MTR) = D11*ACOPY(KB+1:NREMAIN,J) + D21*ACOPY(KB+1:NREMAIN,J+1)
            DL21T(J+1,1:MTR) = D21*ACOPY(KB+1:NREMAIN,J) + D22*ACOPY(KB+1:NREMAIN,J+1)
            J = J + 2
          END IF
        END DO

        CALL COMPRESS_BLOCK(MTR, KB, ACOPY(KB+1,1), NREMAIN, TOL, LBLK)
        N_PANELS = N_PANELS + 1
        IF (LBLK%IS_LOWRANK) THEN
          N_COMPRESSED = N_COMPRESSED + 1
          RANK_SUM = RANK_SUM + LBLK%RANK
        ELSE
          RANK_SUM = RANK_SUM + MIN(MTR,KB)
        END IF

        CALL COMPRESS_BLOCK(KB, MTR, DL21T, KB, TOL, UBLK)
        N_PANELS = N_PANELS + 1
        IF (UBLK%IS_LOWRANK) THEN
          N_COMPRESSED = N_COMPRESSED + 1
          RANK_SUM = RANK_SUM + UBLK%RANK
        ELSE
          RANK_SUM = RANK_SUM + MIN(KB,MTR)
        END IF

        ! A22 <- A22 - L21*DL21T = A22 - L21*D*L21^T (carre complet, voir
        ! note de tete de fichier sur le triangle superieur "gaspille")
        CALL SCHUR_UPDATE_BLOCK(MTR, MTR, KB, LBLK, UBLK, A(KSTART+KB,KSTART+KB), LDA)

        DEALLOCATE(DL21T)
      END IF

      DEALLOCATE(ACOPY)
      KSTART = KSTART + KB
    END DO

    DEALLOCATE(WORK)
    IF (PRESENT(N_2X2_PIVOTS)) N_2X2_PIVOTS = N2X2
  END SUBROUTINE BLR_LDLT_FACTOR

END MODULE BLR_LDLT_MOD
