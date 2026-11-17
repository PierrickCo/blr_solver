!===============================================================================
! Definitions partagees par les solveurs MPI LU et LDLt (extrait de
! src/blr_matrix_mod.f90, sans le solveur sequentiel) :
!   DENSE_PANEL_T   : bloc dense (blocs diagonaux)
!   GET_BLOCK_IFACE : interface de la routine qui fournit A(ISTART:IEND,
!                     JSTART:JEND) a l'assemblage (ici : PCMATGEN)
!===============================================================================
MODULE BLR_MATRIX_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE

  TYPE :: DENSE_PANEL_T
    COMPLEX(WP), ALLOCATABLE :: D(:,:)
  END TYPE DENSE_PANEL_T

  ABSTRACT INTERFACE
    SUBROUTINE GET_BLOCK_IFACE(ISTART, IEND, JSTART, JEND, BLK)
      IMPORT :: WP
      INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
      COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    END SUBROUTINE GET_BLOCK_IFACE
  END INTERFACE

END MODULE BLR_MATRIX_MOD
