!===============================================================================
! COMPLETE LA GRILLE dense/BLR x OpenMP/MPI-sans-pivot/MPI-avec-pivot :
! seule case manquante apres les etapes 6a-6e et les 3 propositions --
! compression BLR sur le squelette DISTRIBUE SANS PIVOT (mpi_lu_nopiv_mod.f90,
! etape 6c), qui n'avait jamais ete recombine avec la compression (6e ne
! recombine que la version AVEC pivot, 6d).
!
! CHOIX D'IMPLEMENTATION : contrairement a 6e (qui diffuse LBUF/UBUF DENSE
! puis compresse cote recepteur -- pattern initial, corrige plus tard par
! la proposition 3), ce module compresse DIRECTEMENT a la source AVANT
! diffusion (BCAST_BLR_BLOCK, proposition 3) des le depart : desormais
! qu'on sait que c'est strictement meilleur (moins de calcul redondant,
! moins d'octets sur les configurations realistes NB<<N -- voir
! mpi_blr_lu_pivot_bcast_mod.f90), inutile de reproduire ici le pattern
! moins bon de 6e puis de le corriger plus tard.
!
! Le reste (rassemblement/diffusion du bloc diagonal, panneaux U/L par
! CTRSM) est repris tel quel de 6c (mpi_lu_nopiv_mod.f90, deja valide) --
! SANS PIVOT, donc matrices a diagonale dominante uniquement, exactement
! comme 6c.
!===============================================================================
MODULE MPI_BLR_LU_NOPIV_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE MPI_LU_NOPIV_MOD, ONLY: LU_NOPIV_UNBLOCKED
  USE MPI_BLR_LU_PIVOT_BCAST_MOD, ONLY: BCAST_BLR_BLOCK
  IMPLICIT NONE

CONTAINS

  SUBROUTINE MPI_BLR_LU_NOPIV_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, TOL, INFO, &
                                      N_COMPRESSED, N_TOTAL_PANELS, RANK_SUM, T_COMPRESS)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)     :: GRID
    INTEGER, INTENT(IN)          :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(INOUT)   :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    REAL(WP), INTENT(IN)         :: TOL
    INTEGER, INTENT(OUT)         :: INFO
    INTEGER, INTENT(OUT)         :: N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8), INTENT(OUT)      :: RANK_SUM
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_COMPRESS

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: DKK(:,:)
    TYPE(BLR_BLOCK_T) :: LBLK, UBLK
    INTEGER :: NBLK, K, KSTART, KEND, KW
    INTEGER :: PROW_K, PCOL_K, OWNER_RANK, RS, CS
    INTEGER :: RTRAIL0, CTRAIL0, NROWTRAIL, NCOLTRAIL
    INTEGER :: INFO_STEP, IERR
    INTEGER :: NC_LOCAL, NP_LOCAL
    INTEGER(8) :: RSUM_LOCAL
    DOUBLE PRECISION :: TC_LOCAL, TC0

    INFO = 0
    NC_LOCAL = 0
    NP_LOCAL = 0
    RSUM_LOCAL = 0_8
    TC_LOCAL = 0.0D0
    NBLK = (N + NB - 1) / NB

    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1

      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)
      OWNER_RANK = PROW_K*GRID%PC + PCOL_K

      RTRAIL0 = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      CTRAIL0 = NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      NROWTRAIL = MLOC - RTRAIL0
      NCOLTRAIL = NLOC - CTRAIL0

      ! --- 1) Factorisation locale du bloc diagonal (proprietaire seul,
      !     identique a 6c) ---
      INFO_STEP = 0
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        CS = INDXG2L(KSTART, NB, GRID%PC)
        CALL LU_NOPIV_UNBLOCKED(KW, A_LOCAL(RS,CS), MLOC, INFO_STEP)
        IF (INFO_STEP .NE. 0) INFO_STEP = INFO_STEP + KSTART - 1
      END IF
      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        RETURN
      END IF

      ! --- 2) Diffusion du bloc diagonal factorise (identique a 6c) ---
      ALLOCATE(DKK(KW,KW))
      IF (GRID%MYROW.EQ.PROW_K .AND. GRID%MYCOL.EQ.PCOL_K) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        CS = INDXG2L(KSTART, NB, GRID%PC)
        DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
      END IF
      IF (GRID%MYROW .EQ. PROW_K) THEN
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
      END IF
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
      END IF

      ! --- 3) Panneau U (CTRSM, identique a 6c) ---
      IF (GRID%MYROW.EQ.PROW_K .AND. NCOLTRAIL.GT.0) THEN
        RS = INDXG2L(KSTART, NB, GRID%PR)
        CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NCOLTRAIL, &
                    CONE, DKK, KW, A_LOCAL(RS,CTRAIL0+1), MLOC)
      END IF

      ! --- 4) Panneau L (CTRSM, identique a 6c) ---
      IF (GRID%MYCOL.EQ.PCOL_K .AND. NROWTRAIL.GT.0) THEN
        CS = INDXG2L(KSTART, NB, GRID%PC)
        CALL CTRSM('Right', 'Upper', 'No transpose', 'Non-unit', NROWTRAIL, KW, &
                    CONE, DKK, KW, A_LOCAL(RTRAIL0+1,CS), MLOC)
      END IF
      DEALLOCATE(DKK)

      ! --- 5) Compresser cote SOURCE puis diffuser compresse (nouveaute,
      !     technique de la proposition 3, PAS le pattern initial de 6e) ---
      IF (NCOLTRAIL .GT. 0) THEN
        IF (GRID%MYROW .EQ. PROW_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          TC0 = MPI_WTIME()
          CALL COMPRESS_BLOCK(KW, NCOLTRAIL, A_LOCAL(RS,CTRAIL0+1), MLOC, TOL, UBLK)
          TC_LOCAL = TC_LOCAL + (MPI_WTIME() - TC0)
          NP_LOCAL = NP_LOCAL + 1
          IF (UBLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RSUM_LOCAL = RSUM_LOCAL + UBLK%RANK
          ELSE
            RSUM_LOCAL = RSUM_LOCAL + MIN(KW,NCOLTRAIL)
          END IF
        END IF
        CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, NCOLTRAIL, UBLK)
      END IF

      IF (NROWTRAIL .GT. 0) THEN
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          TC0 = MPI_WTIME()
          CALL COMPRESS_BLOCK(NROWTRAIL, KW, A_LOCAL(RTRAIL0+1,CS), MLOC, TOL, LBLK)
          TC_LOCAL = TC_LOCAL + (MPI_WTIME() - TC0)
          NP_LOCAL = NP_LOCAL + 1
          IF (LBLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RSUM_LOCAL = RSUM_LOCAL + LBLK%RANK
          ELSE
            RSUM_LOCAL = RSUM_LOCAL + MIN(NROWTRAIL,KW)
          END IF
        END IF
        CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, NROWTRAIL, KW, LBLK)
      END IF

      ! --- 6) Mise a jour de Schur, basse-rang quand possible ---
      IF (NROWTRAIL.GT.0 .AND. NCOLTRAIL.GT.0) THEN
        CALL SCHUR_UPDATE_BLOCK(NROWTRAIL, NCOLTRAIL, KW, LBLK, UBLK, &
                                 A_LOCAL(RTRAIL0+1,CTRAIL0+1), MLOC)
      END IF
    END DO

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NP_LOCAL, N_TOTAL_PANELS, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RSUM_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(T_COMPRESS)) THEN
      CALL MPI_ALLREDUCE(TC_LOCAL, T_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
  END SUBROUTINE MPI_BLR_LU_NOPIV_FACTOR

END MODULE MPI_BLR_LU_NOPIV_MOD
