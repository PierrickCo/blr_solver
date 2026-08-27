!===============================================================================
! Parametre de precision partage par TOUS les modules du solveur -- doit
! etre le SEUL endroit ou WP est defini (sinon "USE"-associer deux modules
! qui redefinissent chacun leur propre WP, meme a la meme valeur, est une
! reference ambigue pour le compilateur -- rencontre en pratique en
! ecrivant les tests de l'etape 2).
!===============================================================================
MODULE BLR_KINDS_MOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: WP = KIND(1.0)   ! simple precision (reel et complexe)
END MODULE BLR_KINDS_MOD
