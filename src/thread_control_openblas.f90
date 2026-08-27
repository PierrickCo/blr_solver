!===============================================================================
! Implementation THREAD_CONTROL_MOD pour la chaine gcc+OpenBLAS -- appelle
! openblas_set_num_threads (fonction C d'OpenBLAS) via BIND(C), l'equivalent
! d'OpenBLAS a MKL_SET_NUM_THREADS.
!
! Selectionnee au build par le Makefile (TOOLCHAIN=gcc, si OpenBLAS est
! installe -- voir INSTALL.md). Necessite -lopenblas a l'edition de liens.
!===============================================================================
MODULE THREAD_CONTROL_MOD
  USE ISO_C_BINDING, ONLY: C_INT
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE OPENBLAS_SET_NUM_THREADS_C(N) BIND(C, NAME="openblas_set_num_threads")
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: N
    END SUBROUTINE OPENBLAS_SET_NUM_THREADS_C
  END INTERFACE

CONTAINS
  SUBROUTINE SET_NUM_THREADS(N)
    INTEGER, INTENT(IN) :: N
    CALL OPENBLAS_SET_NUM_THREADS_C(INT(N, C_INT))
  END SUBROUTINE SET_NUM_THREADS
END MODULE THREAD_CONTROL_MOD
