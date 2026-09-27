!===============================================================================
! Parametre de precision partage par l'ensemble du solveur. Doit rester le
! SEUL endroit ou WP est defini : deux modules redefinissant chacun leur
! propre WP, meme a la meme valeur, produisent une reference ambigue a la
! compilation des qu'ils sont USE-associes ensemble.
!===============================================================================
MODULE BLR_KINDS_MOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: WP = KIND(1.0)   ! simple precision (reel et complexe)
END MODULE BLR_KINDS_MOD
