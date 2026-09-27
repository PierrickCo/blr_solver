!===============================================================================
! Implementation THREAD_CONTROL_MOD pour la chaine Intel (ifx+MKL) -- appelle
! directement MKL_SET_NUM_THREADS (extension Intel, liee implicitement comme
! toute routine BLAS/LAPACK dans ce projet, sans USE MKL_SERVICE).
!
! Selectionnee au build par le Makefile (TOOLCHAIN=intel, defaut) -- SEUL ce
! fichier ou thread_control_openblas.f90/thread_control_noop.f90 doit etre
! compile, jamais plusieurs a la fois (meme nom de module).
!===============================================================================
MODULE THREAD_CONTROL_MOD
  USE MKL_INTERFACES_MOD, ONLY: MKL_SET_NUM_THREADS
  IMPLICIT NONE
CONTAINS
  SUBROUTINE SET_NUM_THREADS(N)
    INTEGER, INTENT(IN) :: N
    CALL MKL_SET_NUM_THREADS(N)
  END SUBROUTINE SET_NUM_THREADS
END MODULE THREAD_CONTROL_MOD
