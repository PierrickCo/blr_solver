!===============================================================================
! ETAPE 4 (partie 1) : variante LDLt (SYMETRIQUE complexe, PAS hermitienne
! -- A(i,j)=A(j,i), meme convention que CMUMPS SYM=2 utilise partout
! ailleurs dans ce projet), sans compression pour l'instant (la compression
! BLR pour LDLt est une extension future -- voir la note de portee en fin
! de fichier).
!
! DIFFERENCE STRUCTURELLE AVEC LA LU (block_lu_mod.f90) : pour une matrice
! symetrique indefinie, le pivot doit rester symetrique (echanger la ligne
! i ET la colonne i ensemble) -- un simple pivot 1x1 par ligne/colonne
! N'EST PAS TOUJOURS STABLE (ex. matrice avec un bloc [[0,1],[1,0]] : aucun
! pivot diagonal 1x1 n'est utilisable). La methode de Bunch-Kaufman, que
! nous N'AVONS PAS reimplementee nous-memes, choisit donc entre un pivot
! 1x1 (cas normal) ou un pivot 2x2 (quand necessaire pour la stabilite),
! donnant une factorisation A = L*D*L^T avec D diagonale PAR BLOCS (blocs
! 1x1 ou 2x2).
!
! POURQUOI ON UTILISE CLASYF (LAPACK) DIRECTEMENT, PAS UNE REIMPLEMENTATION
! MAISON COMME POUR LA LU : contrairement a CGETRF (qui ne fait QUE la
! factorisation d'un panneau, laissant l'orchestration bloc/mise a jour de
! Schur a notre charge), CLASYF est DEJA le "primitif de blocage" que
! LAPACK utilise en interne pour sa PROPRE version bloquee de CSYTRF -- un
! seul appel factorise jusqu'a NB colonnes (pivot de Bunch-Kaufman confine
! au sous-tableau ACTIF qu'on lui passe -- meme philosophie "confine au
! panneau, pas global" que pour la LU) ET applique la mise a jour de Schur
! resultante a la sous-matrice restante, en une seule fois. Reimplementer
! ce fused panel-factor + trailing-update nous-memes (comme on l'a fait a
! la main pour CGETRF+CLASWP+CTRSM+CGEMM) demanderait de re-deriver
! nous-memes la logique de selection 1x1/2x2 de Bunch-Kaufman -- l'un des
! algorithmes les plus delicats de l'algebre lineaire dense, avec un vrai
! risque de bug subtil. On accepte ici de perdre le controle fin necessaire
! pour injecter la compression BLR (contrairement a la LU) en echange
! d'une factorisation LDLt correcte et fiable des le depart.
!
! SOLVE : notre factorisation (obtenue bloc par bloc via CLASYF) produit
! EXACTEMENT le meme format (A, IPIV) qu'un appel direct a CSYTRF sur la
! matrice entiere -- CLASYF EST le mecanisme de blocage interne de CSYTRF,
! donc l'enchainement de plusieurs appels CLASYF (un par panneau) est
! numeriquement equivalent a UN appel CSYTRF bloque. On reutilise donc
! CSYTRS de LAPACK telle quelle pour resoudre, sans le reimplementer
! (contrairement a BLOCK_LU_SOLVE, qui devait gerer manuellement la
! convention "pivots locaux au panneau" -- ici CSYTRS sait deja faire).
!
! PAS DE DIRECTIVE OPENMP NON PLUS (etape 5) : CLASYF est UN SEUL appel par
! panneau qui fait factorisation+mise a jour en une fois (boite noire, cf.
! ci-dessus) -- comme pour block_lu_mod.f90, pas de boucle de blocs
! independants a notre niveau a paralleliser ; le parallelisme vient du
! threading interne de MKL sur cet appel.
!
! PORTEE : PAS DE COMPRESSION BLR ICI (contrairement a blr_lu_mod.f90).
! CLASYF etant une boite noire qui fait factorisation+mise a jour en un
! seul appel, on n'a pas acces aux blocs L21 individuels pour les
! compresser avant la mise a jour de Schur, contrairement a la LU ou on
! ecrivait nous-memes le CGEMM et pouvait donc le remplacer. Une extension
! future exploiterait le tableau de travail W de CLASYF (qui contient des
! valeurs intermediaires exploitables) ou reimplementerait Bunch-Kaufman a
! la main -- hors perimetre de cette etape.
!===============================================================================
MODULE BLOCK_LDLT_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! Factorisation A = L*D*L^T par blocs (UPLO='Lower' -- seule la
  ! partie triangulaire inferieure de A, diagonale incluse, est lue/ecrite ;
  ! la partie superieure n'est ni lue ni modifiee, comme CSYTRF).
  !-----------------------------------------------------------------------
  SUBROUTINE BLOCK_LDLT_FACTOR(N, NB, A, LDA, IPIV, INFO)
    INTEGER, INTENT(IN)        :: N, NB, LDA
    COMPLEX(WP), INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(OUT)       :: IPIV(N)
    INTEGER, INTENT(OUT)       :: INFO

    COMPLEX(WP), ALLOCATABLE :: WORK(:,:)
    INTEGER :: KSTART, NREMAIN, KB, IINFO, II

    ALLOCATE(WORK(N,NB))
    INFO = 0
    KSTART = 1
    DO WHILE (KSTART .LE. N)
      NREMAIN = N - KSTART + 1
      ! CLASYF factorise jusqu'a NB colonnes du sous-tableau ACTIF
      ! [KSTART:N,KSTART:N] (KB<=NB colonnes reellement traitees, KB peut
      ! valoir NB-1 si un pivot 2x2 tombe juste sur la frontiere du bloc),
      ! et met a jour la sous-matrice restante EN PLACE.
      CALL CLASYF('Lower', NREMAIN, MIN(NB,NREMAIN), KB, &
                   A(KSTART,KSTART), LDA, IPIV(KSTART), WORK, N, IINFO)
      IF (IINFO .NE. 0 .AND. INFO .EQ. 0) INFO = KSTART - 1 + IINFO

      ! IMPORTANT (meme categorie de bug que BLOCK_LU_FACTOR a l'etape 1,
      ! mais avec une subtilite en plus ici) : CLASYF, comme CGETRF, ne
      ! connait pas KSTART -- les indices qu'elle renvoie dans IPIV sont
      ! LOCAUX au sous-tableau [KSTART:N] qu'on lui a passe (valeur 1 =
      ! notre ligne/colonne globale KSTART). Contrairement a CGETRF ou ces
      ! indices locaux restaient directement utilisables (on reappliquait
      ! nous-memes les echanges via CLASWP sur des sous-tableaux demarrant
      ! a KSTART, cf. block_lu_mod.f90), ici on doit les convertir en
      ! indices GLOBAUX avant de les stocker, car on delegue ensuite TOUTE
      ! la resolution a CSYTRS (qui, lui, attend des indices globaux,
      ! exactement comme s'ils venaient d'un unique appel CSYTRF sur la
      ! matrice entiere). Convention Bunch-Kaufman (voir doc CLASYF) :
      ! IPIV(k)>0 = pivot 1x1, valeur = ligne/colonne LOCALE echangee ;
      ! IPIV(k)=IPIV(k+1)<0 = pivot 2x2, |valeur| = ligne/colonne LOCALE
      ! echangee. Le signe porte l'information "1x1 ou 2x2" et doit donc
      ! etre PRESERVE pendant la conversion local -> global : on decale
      ! uniquement la MAGNITUDE de (KSTART-1), jamais le signe.
      DO II = KSTART, KSTART+KB-1
        IPIV(II) = IPIV(II) + SIGN(KSTART-1, IPIV(II))
      END DO

      KSTART = KSTART + KB
    END DO
    DEALLOCATE(WORK)
  END SUBROUTINE BLOCK_LDLT_FACTOR

  !-----------------------------------------------------------------------
  ! Resolution A*X = B a partir des facteurs calcules par
  ! BLOCK_LDLT_FACTOR. Simple relais vers CSYTRS -- voir la note en tete
  ! de fichier sur pourquoi c'est correct malgre une factorisation faite
  ! panneau par panneau ici plutot qu'en un seul appel CSYTRF.
  !-----------------------------------------------------------------------
  SUBROUTINE BLOCK_LDLT_SOLVE(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    INTEGER, INTENT(IN)        :: N, NRHS, LDA, LDB
    COMPLEX(WP), INTENT(IN)    :: A(LDA,N)
    INTEGER, INTENT(IN)        :: IPIV(N)
    COMPLEX(WP), INTENT(INOUT) :: B(LDB,NRHS)
    INTEGER, INTENT(OUT)       :: INFO

    CALL CSYTRS('Lower', N, NRHS, A, LDA, IPIV, B, LDB, INFO)
  END SUBROUTINE BLOCK_LDLT_SOLVE

END MODULE BLOCK_LDLT_MOD
