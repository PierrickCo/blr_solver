!===============================================================================
! ETAPE 6e : recombine 6d (factorisation LU distribuee, pivot partiel) avec
! la compression BLR de l'etape 2 -- dernier maillon du chantier MPI.
!
! MEME CHOIX DE PORTEE QU'A L'ETAPE 2 (blr_lu_mod.f90), applique ici sans
! changement de philosophie : la compression n'intervient QUE dans la mise
! a jour de Schur (le CGEMM final de 6d devient une dispatch basse-rang) --
! les facteurs stockes dans A_LOCAL restent DENSES, donc MPI_LU_PIVOT_FACTOR
! (6d) et sa methode de validation (rassembler + reutiliser BLOCK_LU_SOLVE)
! restent directement valables ici aussi.
!
! GRANULARITE DE COMPRESSION -- SIMPLIFICATION DELIBEREE (a documenter
! clairement) : en sequentiel (blr_lu_mod.f90), chaque PAIRE de blocs
! NB x NB de la sous-matrice trainante est compressee INDIVIDUELLEMENT
! (boucle sur tous les blocs-lignes/colonnes trainants). Ici, on compresse
! le panneau L LOCAL D'UN PROCESSUS (LBUF, toutes ses lignes trainantes
! locales, potentiellement plusieurs blocs-lignes globaux cycliques
! regroupes) comme UN SEUL bloc, et de meme pour le panneau U local (UBUF)
! -- une granularite plus grossiere que le sequentiel (moins d'occasions de
! compression fine si les blocs globaux ont des rangs tres differents),
! mais bien plus simple a mettre en oeuvre correctement (pas de
! sous-decoupage NB x NB a gerer en plus du decoupage bloc-cyclique
! process-par-process deja complexe depuis 6c/6d). Le rang NUMERIQUE d'un
! panneau reste borne par le rang global de la partie basse-rang de la
! matrice (independant du decoupage choisi), donc la compression reste
! reelle et mesurable avec cette granularite -- voir test_mpi_blr_lu_pivot
! qui verifie explicitement qu'une fraction significative des panneaux est
! compressee. Un affinage futur (compression NB x NB a l'interieur du
! panneau local) est possible sans changer l'architecture globale, si
! necessaire plus tard pour des matrices ou le rang varie beaucoup d'un
! bloc global a l'autre.
!
! COMMUNICATION : LBUF/UBUF restent diffuses DENSES (comme en 6d) avant
! compression -- on ne cherche PAS ici a reduire le volume de
! communication en diffusant des facteurs UFAC/VFAC compresses a la place
! (ce serait un gain reel et propre a la version distribuee, sans
! equivalent sequentiel, mais ajoute une complexite de diffusion a taille
! variable -- laisse en extension future explicite, pas fait ici pour
! garder ce premier passage BLR+MPI aussi simple que possible).
!===============================================================================
MODULE MPI_BLR_LU_PIVOT_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE MPI_LU_PIVOT_MOD, ONLY: SWAP_ROWS_RANGE
  IMPLICIT NONE

CONTAINS

  ! T_BCAST (optionnel, NOUVEAU) : temps CUMULE (max/processus) passe DANS
  ! les 2 MPI_BCAST denses de LBUF/UBUF -- ajoute pour comparer directement
  ! au T_BCAST equivalent de la proposition 3 (mpi_blr_lu_pivot_bcast_mod.f90)
  ! et verifier l'hypothese de serialisation discutee avec l'utilisateur :
  ! ici, la diffusion demarre IMMEDIATEMENT (rien a calculer avant), et
  ! CHAQUE processus recompresse ENSUITE en parallele (T_COMPRESS, deja
  ! existant) -- pas de dependance source-doit-finir-avant-que-les-autres-
  ! demarrent comme en P3.
  SUBROUTINE MPI_BLR_LU_PIVOT_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, TOL, &
                                      IPIV_GLOBAL, INFO, N_COMPRESSED, N_TOTAL_PANELS, &
                                      RANK_SUM, N_SWAPS, MIN_PIVOT, T_COMPRESS, T_BCAST)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(INOUT)      :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    REAL(WP), INTENT(IN)            :: TOL
    INTEGER, INTENT(OUT)            :: IPIV_GLOBAL(N)
    INTEGER, INTENT(OUT)            :: INFO
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    INTEGER, INTENT(OUT), OPTIONAL  :: N_SWAPS
    REAL(WP), INTENT(OUT), OPTIONAL :: MIN_PIVOT
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_COMPRESS, T_BCAST

    COMPLEX(WP), ALLOCATABLE :: PANEL(:,:), RBUF(:,:), SBUF(:,:), DKK(:,:), LBUF(:,:), UBUF(:,:)
    TYPE(BLR_BLOCK_T) :: LBLK, UBLK
    INTEGER :: IPIV_LOCAL(NB), IPIV_PANEL_GLOBAL(NB)
    INTEGER :: NBLK, K, KSTART, KEND, KW, MPANEL
    INTEGER :: PROW_K, PCOL_K, OWNER_RANK, RS, CS, RSTART0, NROWS_MINE
    INTEGER :: R, RMLOC_R, RSTART0_R, RNROWS_R, I, J, GI
    INTEGER :: LEFT_COUNT, RIGHT_COUNT, RTRAIL0, CTRAIL0, NROWTRAIL, NCOLTRAIL
    INTEGER :: R1, R2, PROW2, IINFO, INFO_STEP, IERR, NS
    INTEGER :: NC_LOCAL, NP_LOCAL
    INTEGER(8) :: RSUM_LOCAL
    INTEGER :: STATUS(MPI_STATUS_SIZE)
    REAL(WP) :: MP
    DOUBLE PRECISION :: TC_LOCAL, TB_LOCAL, TC0

    INFO = 0
    NS = 0
    MP = HUGE(1.0_WP)
    NC_LOCAL = 0
    NP_LOCAL = 0
    RSUM_LOCAL = 0_8
    TC_LOCAL = 0.0D0
    TB_LOCAL = 0.0D0
    NBLK = (N + NB - 1) / NB

    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      MPANEL = N - KSTART + 1

      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)
      OWNER_RANK = PROW_K*GRID%PC + PCOL_K

      INFO_STEP = 0

      ! === A+B+C : rassembler le panneau, factoriser (CGETRF), renvoyer ===
      ! (identique a 6d -- voir mpi_lu_pivot_mod.f90 pour les commentaires
      ! detailles de cette partie, inchangee)
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        RSTART0    = NUMROC(KSTART-1, NB, GRID%MYROW, GRID%PR)
        NROWS_MINE = MLOC - RSTART0
        CS         = INDXG2L(KSTART, NB, GRID%PC)

        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(PANEL(MPANEL,KW))
          IF (NROWS_MINE .GT. 0) THEN
            DO J = 1, KW
              DO I = 1, NROWS_MINE
                GI = INDXL2G(RSTART0+I, NB, GRID%MYROW, GRID%PR)
                PANEL(GI-KSTART+1,J) = A_LOCAL(RSTART0+I, CS+J-1)
              END DO
            END DO
          END IF
          DO R = 0, GRID%PR-1
            IF (R .EQ. PROW_K) CYCLE
            RMLOC_R   = NUMROC(N, NB, R, GRID%PR)
            RSTART0_R = NUMROC(KSTART-1, NB, R, GRID%PR)
            RNROWS_R  = RMLOC_R - RSTART0_R
            IF (RNROWS_R .LE. 0) CYCLE
            ALLOCATE(RBUF(RNROWS_R,KW))
            CALL MPI_RECV(RBUF, RNROWS_R*KW, MPI_COMPLEX, R, 10, GRID%COL_COMM, STATUS, IERR)
            DO J = 1, KW
              DO I = 1, RNROWS_R
                GI = INDXL2G(RSTART0_R+I, NB, R, GRID%PR)
                PANEL(GI-KSTART+1,J) = RBUF(I,J)
              END DO
            END DO
            DEALLOCATE(RBUF)
          END DO

          CALL CGETRF(MPANEL, KW, PANEL, MPANEL, IPIV_LOCAL, IINFO)
          IF (IINFO .NE. 0) INFO_STEP = KSTART - 1 + IINFO
          DO J = 1, KW
            IPIV_PANEL_GLOBAL(J) = IPIV_LOCAL(J)
            MP = MIN(MP, ABS(PANEL(J,J)))
          END DO

          IF (NROWS_MINE .GT. 0) THEN
            DO J = 1, KW
              DO I = 1, NROWS_MINE
                GI = INDXL2G(RSTART0+I, NB, GRID%MYROW, GRID%PR)
                A_LOCAL(RSTART0+I, CS+J-1) = PANEL(GI-KSTART+1,J)
              END DO
            END DO
          END IF
          DO R = 0, GRID%PR-1
            IF (R .EQ. PROW_K) CYCLE
            RMLOC_R   = NUMROC(N, NB, R, GRID%PR)
            RSTART0_R = NUMROC(KSTART-1, NB, R, GRID%PR)
            RNROWS_R  = RMLOC_R - RSTART0_R
            IF (RNROWS_R .LE. 0) CYCLE
            ALLOCATE(SBUF(RNROWS_R,KW))
            DO J = 1, KW
              DO I = 1, RNROWS_R
                GI = INDXL2G(RSTART0_R+I, NB, R, GRID%PR)
                SBUF(I,J) = PANEL(GI-KSTART+1,J)
              END DO
            END DO
            CALL MPI_SEND(SBUF, RNROWS_R*KW, MPI_COMPLEX, R, 11, GRID%COL_COMM, IERR)
            DEALLOCATE(SBUF)
          END DO
          DEALLOCATE(PANEL)
        ELSE
          IF (NROWS_MINE .GT. 0) THEN
            ALLOCATE(SBUF(NROWS_MINE,KW))
            SBUF = A_LOCAL(RSTART0+1:MLOC, CS:CS+KW-1)
            CALL MPI_SEND(SBUF, NROWS_MINE*KW, MPI_COMPLEX, PROW_K, 10, GRID%COL_COMM, IERR)
            DEALLOCATE(SBUF)

            ALLOCATE(RBUF(NROWS_MINE,KW))
            CALL MPI_RECV(RBUF, NROWS_MINE*KW, MPI_COMPLEX, PROW_K, 11, GRID%COL_COMM, STATUS, IERR)
            A_LOCAL(RSTART0+1:MLOC, CS:CS+KW-1) = RBUF
            DEALLOCATE(RBUF)
          END IF
        END IF
      END IF

      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        RETURN
      END IF
      CALL MPI_BCAST(IPIV_PANEL_GLOBAL, KW, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IPIV_GLOBAL(KSTART:KEND) = IPIV_PANEL_GLOBAL(1:KW)

      ! === E : propager les echanges de lignes (identique a 6d) ===
      LEFT_COUNT  = NUMROC(KSTART-1, NB, GRID%MYCOL, GRID%PC)
      RIGHT_COUNT = NLOC - NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      DO I = 1, KW
        R1 = KSTART - 1 + I
        R2 = KSTART - 1 + IPIV_PANEL_GLOBAL(I)
        IF (R1 .NE. R2) THEN
          NS = NS + 1
          PROW2 = INDXG2P(R2, NB, GRID%PR)
          IF (LEFT_COUNT .GT. 0) THEN
            CALL SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW_K, PROW2, &
                                  1, LEFT_COUNT)
          END IF
          IF (RIGHT_COUNT .GT. 0) THEN
            CALL SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW_K, PROW2, &
                                  NLOC-RIGHT_COUNT+1, NLOC)
          END IF
        END IF
      END DO

      ! === F : U12 par CTRSM (identique a 6d), puis Schur update BLR ===
      RTRAIL0   = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      CTRAIL0   = NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      NROWTRAIL = MLOC - RTRAIL0
      NCOLTRAIL = NLOC - CTRAIL0

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NCOLTRAIL .GT. 0) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NCOLTRAIL, &
                      (1.0_WP,0.0_WP), DKK, KW, A_LOCAL(RS,CTRAIL0+1), MLOC)
        END IF
        DEALLOCATE(DKK)
      END IF

      IF (NROWTRAIL .GT. 0) THEN
        ALLOCATE(LBUF(NROWTRAIL,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          CS = INDXG2L(KSTART, NB, GRID%PC)
          LBUF = A_LOCAL(RTRAIL0+1:MLOC, CS:CS+KW-1)
        END IF
        TC0 = MPI_WTIME()
        CALL MPI_BCAST(LBUF, NROWTRAIL*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        TB_LOCAL = TB_LOCAL + (MPI_WTIME() - TC0)
      END IF

      IF (NCOLTRAIL .GT. 0) THEN
        ALLOCATE(UBUF(KW,NCOLTRAIL))
        IF (GRID%MYROW .EQ. PROW_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          UBUF = A_LOCAL(RS:RS+KW-1, CTRAIL0+1:NLOC)
        END IF
        TC0 = MPI_WTIME()
        CALL MPI_BCAST(UBUF, KW*NCOLTRAIL, MPI_COMPLEX, PROW_K, GRID%COL_COMM, IERR)
        TB_LOCAL = TB_LOCAL + (MPI_WTIME() - TC0)
      END IF

      ! --- mise a jour de Schur : compression LOCALE (par processus) des
      !     panneaux L/U juste recus, puis dispatch basse-rang (au lieu du
      !     CGEMM dense de 6c/6d) ---
      IF (NROWTRAIL.GT.0 .AND. NCOLTRAIL.GT.0) THEN
        TC0 = MPI_WTIME()
        CALL COMPRESS_BLOCK(NROWTRAIL, KW, LBUF, NROWTRAIL, TOL, LBLK)
        NP_LOCAL = NP_LOCAL + 1
        IF (LBLK%IS_LOWRANK) THEN
          NC_LOCAL = NC_LOCAL + 1
          RSUM_LOCAL = RSUM_LOCAL + LBLK%RANK
        ELSE
          RSUM_LOCAL = RSUM_LOCAL + MIN(NROWTRAIL,KW)
        END IF

        CALL COMPRESS_BLOCK(KW, NCOLTRAIL, UBUF, KW, TOL, UBLK)
        NP_LOCAL = NP_LOCAL + 1
        IF (UBLK%IS_LOWRANK) THEN
          NC_LOCAL = NC_LOCAL + 1
          RSUM_LOCAL = RSUM_LOCAL + UBLK%RANK
        ELSE
          RSUM_LOCAL = RSUM_LOCAL + MIN(KW,NCOLTRAIL)
        END IF
        TC_LOCAL = TC_LOCAL + (MPI_WTIME() - TC0)

        CALL SCHUR_UPDATE_BLOCK(NROWTRAIL, NCOLTRAIL, KW, LBLK, UBLK, &
                                 A_LOCAL(RTRAIL0+1,CTRAIL0+1), MLOC)
      END IF

      IF (NROWTRAIL .GT. 0) DEALLOCATE(LBUF)
      IF (NCOLTRAIL .GT. 0) DEALLOCATE(UBUF)
    END DO

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NP_LOCAL, N_TOTAL_PANELS, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RSUM_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS
    IF (PRESENT(MIN_PIVOT)) THEN
      CALL MPI_ALLREDUCE(MP, MIN_PIVOT, 1, MPI_REAL, MPI_MIN, GRID%COMM, IERR)
    END IF
    IF (PRESENT(T_COMPRESS)) THEN
      ! max sur tous les processus : cout de compression sur le chemin
      ! critique, pas une simple moyenne.
      CALL MPI_ALLREDUCE(TC_LOCAL, T_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
    IF (PRESENT(T_BCAST)) THEN
      CALL MPI_ALLREDUCE(TB_LOCAL, T_BCAST, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
  END SUBROUTINE MPI_BLR_LU_PIVOT_FACTOR

END MODULE MPI_BLR_LU_PIVOT_MOD
