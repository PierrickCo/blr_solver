!===============================================================================
! Verification ||A*X-B||/||B|| SANS solution exacte connue -- fonctionne
! avec n'importe quel B (formule, fichier, etc.), tant qu'on lui fournit
! GET_BLOCK (interface GET_BLOCK_IFACE, meme convention que la
! bibliotheque : renvoie un bloc de A, jamais A entiere). Recoit GET_BLOCK
! en argument plutot que de dependre de user_matrix_mod.f90 directement --
! reste ainsi reutilisable avec n'importe quelle matrice/GET_BLOCK
! conforme a cette interface, pas seulement USER_GET_BLOCK.
!===============================================================================
MODULE RESIDUAL_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: GRID_T, WP
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BLR_MATRIX_MOD, ONLY: GET_BLOCK_IFACE
  USE MKL_INTERFACES_MOD, ONLY: CGEMM
  USE MPI
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ACCUMULATE_RESIDUAL

CONTAINS

  ! Accumule ||A*X-B||^2 (NUM2) et ||B||^2 (DEN2) pour CE lot dans
  ! NUM2_ACC/DEN2_ACC (deja des sommes globales, reduites ci-dessous) -- la
  ! norme de Frobenius totale se recompose exactement comme la somme des
  ! normes par lot (colonnes independantes), donc le residu final
  ! SQRT(somme des NUM2)/SQRT(somme des DEN2), calcule une seule fois a la
  ! toute fin (voir le site d'appel), reste EXACT.
  !
  ! Appel COLLECTIF (tous les rangs). A*X est calcule par blocs de NB
  ! lignes via GET_BLOCK + CGEMM (BLAS), lignes reparties entre tous les
  ! rangs -- jamais de matrice N x N complete en RAM. BORIG_LOCAL (le b
  ! reel, sauvegarde avant que la resolution ne l'ecrase en place) est
  ! rassemble ici -- ca marche donc que B vienne d'une formule, d'un
  ! fichier, ou de n'importe quelle autre source.
  SUBROUTINE ACCUMULATE_RESIDUAL(GRID, GET_BLOCK, N, NB, MLOC, XG, BORIG_LOCAL, NR, NRL, NUM2_ACC, DEN2_ACC)
    TYPE(GRID_T), INTENT(IN)   :: GRID
    PROCEDURE(GET_BLOCK_IFACE) :: GET_BLOCK
    INTEGER, INTENT(IN)        :: N, NB, MLOC, NR, NRL
    COMPLEX(WP), INTENT(INOUT) :: XG(N,NR)
    COMPLEX(WP), INTENT(IN)    :: BORIG_LOCAL(MAX(MLOC,1),MAX(NRL,1))
    REAL(WP), INTENT(INOUT)    :: NUM2_ACC, DEN2_ACC
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), AXBLOCK(:,:), BG(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK
    INTEGER :: GI, R, IERR2
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2

    ALLOCATE(BG(N,NR))
    CALL GATHER_2D_RHS(GRID, N, NB, NR, BORIG_LOCAL, MAX(MLOC,1), NRL, 0, BG)
    CALL MPI_BCAST(XG, N*NR, MPI_COMPLEX, 0, GRID%COMM, IERR2)
    CALL MPI_BCAST(BG, N*NR, MPI_COMPLEX, 0, GRID%COMM, IERR2)

    ROWS_PER_PROC = (N + GRID%NPROCS - 1) / GRID%NPROCS
    RSTART_MINE = GRID%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GRID%MYID+1)*ROWS_PER_PROC, N)

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP
    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, NB
      REND_CHUNK = MIN(RSTART_CHUNK+NB-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(ABLOCK(NCHUNK,N), AXBLOCK(NCHUNK,NR))
      CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, N, ABLOCK)
      CALL CGEMM('N','N', NCHUNK, NR, N, CONE, ABLOCK, NCHUNK, XG, N, CZERO, AXBLOCK, NCHUNK)
      DO R = 1, NR
        DO GI = 1, NCHUNK
          NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(AXBLOCK(GI,R)-BG(RSTART_CHUNK+GI-1,R))**2,WP)
          DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BG(RSTART_CHUNK+GI-1,R))**2,WP)
        END DO
      END DO
      DEALLOCATE(ABLOCK, AXBLOCK)
    END DO
    DEALLOCATE(BG)

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GRID%COMM, IERR2)
    NUM2_ACC = NUM2_ACC + NUM2
    DEN2_ACC = DEN2_ACC + DEN2
  END SUBROUTINE ACCUMULATE_RESIDUAL

END MODULE RESIDUAL_MOD
