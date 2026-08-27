!===============================================================================
! Implementation THREAD_CONTROL_MOD de repli : ne fait RIEN. Utilisee quand
! aucune bibliotheque BLAS/LAPACK avec controle de threads n'est disponible
! (ex: BLAS de reference Netlib, intrinsequement mono-thread -- il n'y a
! alors rien a limiter/restaurer, les boucles !$OMP de blr_matrix_mod.f90
! restent correctes, simplement sans le garde-fou anti-sursouscription qui
! n'a pas lieu d'etre dans ce cas).
!===============================================================================
MODULE THREAD_CONTROL_MOD
  IMPLICIT NONE
CONTAINS
  SUBROUTINE SET_NUM_THREADS(N)
    INTEGER, INTENT(IN) :: N
    IF (N .LT. 0) RETURN  ! evite un avertissement "argument inutilise"
  END SUBROUTINE SET_NUM_THREADS
END MODULE THREAD_CONTROL_MOD
