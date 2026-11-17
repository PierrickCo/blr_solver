!===============================================================================
! Sauvegarde/reprise de A_BLR et du RHS/solution -- un fichier PAR
! PROCESSUS dans CHECKPOINT_DIR, STREAM non formate. Routines sans etat
! cache (GRID/A_BLR/N/NB/PR/PC/MYID0 passes explicitement) -- reutilisables
! par n'importe quel appelant. Deux familles :
!
!   - SAVE_BLR_MATRIX/LOAD_BLR_MATRIX, OPEN_RHS_CHECKPOINT_WRITE/READ,
!     SAVE_RHS_BATCH/LOAD_RHS_BATCH : appelees PAR CHAQUE RANG d'un run MPI
!     en cours (le sien uniquement), jamais de communication a l'ecriture/
!     la lecture sauf une reduction MPI_ALLREDUCE si une incompatibilite
!     d'en-tete doit arreter tous les rangs proprement.
!
!   - LOAD_SOLUTION/LOAD_RHS : post-traitement, SEQUENTIELLES (aucun GRID
!     ni communication MPI), appelables independamment et apres coup par
!     un seul appelant -- reconstruisent elles-memes la sous-matrice
!     GLOBALE demandee en parcourant tous les PR*PC fichiers.
!===============================================================================
MODULE CHECKPOINT_IO_MOD
  USE MPI_BLR_MATRIX_MOD
  USE MPI
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SAVE_BLR_MATRIX, LOAD_BLR_MATRIX
  PUBLIC :: OPEN_RHS_CHECKPOINT_WRITE, SAVE_RHS_BATCH
  PUBLIC :: OPEN_RHS_CHECKPOINT_READ, LOAD_RHS_BATCH
  PUBLIC :: LOAD_SOLUTION
  PUBLIC :: RECOMPUTE_COMPRESSION_STATS, REAL_COMPRESSION_PCT

CONTAINS

  ! Ecrit A_BLR (SAVE_MATRIX/SAVE_FACTORED) -- un fichier par processus,
  ! STREAM non formate, sans communication MPI. En-tete d'abord (verifie
  ! par LOAD_BLR_MATRIX), puis les blocs diagonaux localement possedes,
  ! puis tous les blocs hors-diagonale locaux (toujours "possedes" par
  ! construction du decoupage 2D block-cyclique). WITH_IPIV n'ecrit IPIV
  ! que si elle est allouee (apres factorisation).
  SUBROUTINE SAVE_BLR_MATRIX(GRID, A_BLR, N, NB, PR, PC, MYID0, PREFIX, WITH_IPIV, OK)
    TYPE(GRID_T), INTENT(IN) :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(IN) :: N, NB, PR, PC, MYID0
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: U, li2, lj2, I2, J2, ND, MI, IOS, IERR2
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    OK_LOCAL = (IOS .EQ. 0)
    IF (.NOT. OK_LOCAL) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible de creer ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- le repertoire CHECKPOINT_DIR doit deja exister.'
    END IF
    ! Un seul rang qui ne peut pas ecrire doit arreter tous les rangs --
    ! sinon les autres continueraient seuls vers un appel MPI collectif
    ! que ce rang n'atteindra jamais.
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U, STATUS='DELETE')
      RETURN
    END IF

    ! Marqueur de format (N > 0 dans l'ancien format, sans marqueur) +
    ! distribution des blocs-colonnes (CMAP, voir BLR_COL_OWNER).
    WRITE(U) -7411, A_BLR%CMAP
    WRITE(U) N, NB, A_BLR%NBLK, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, &
         A_BLR%NBLKR_LOC, A_BLR%NBLKC_LOC

    ND = 0
    DO I2 = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(I2)%D)) ND = ND + 1
    END DO
    WRITE(U) ND
    DO I2 = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(I2)%D)) THEN
        MI = SIZE(A_BLR%DIAGBLK(I2)%D, 1)
        WRITE(U) I2, MI
        WRITE(U) A_BLR%DIAGBLK(I2)%D
      END IF
    END DO

    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        ! Position diagonale d'OFFBLK : jamais compressee (I.EQ.J), rien
        ! a ecrire.
        IF (I2 .EQ. J2) CYCLE
        WRITE(U) I2, J2, A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK, A_BLR%OFFBLK(li2,lj2)%M, &
             A_BLR%OFFBLK(li2,lj2)%N, A_BLR%OFFBLK(li2,lj2)%RANK
        IF (A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK) THEN
          WRITE(U) A_BLR%OFFBLK(li2,lj2)%UFAC
          WRITE(U) A_BLR%OFFBLK(li2,lj2)%VFAC
        ELSE
          WRITE(U) A_BLR%OFFBLK(li2,lj2)%DENSE
        END IF
      END DO
    END DO

    IF (WITH_IPIV) WRITE(U) A_BLR%IPIV
    CLOSE(U)
  END SUBROUTINE SAVE_BLR_MATRIX

  ! Relit A_BLR ecrite par SAVE_BLR_MATRIX -- reconstruit BSTART/BEND/
  ! NBLKR_LOC/NBLKC_LOC/DIAGBLK/OFFBLK comme le ferait MPI_ASSEMBLE_BLR_
  ! MATRIX, donc la factorisation qui suit ne voit aucune difference.
  ! OK=.FALSE. (tous les rangs s'arretent) si le fichier est absent ou si
  ! son en-tete ne correspond pas au run courant.
  SUBROUTINE LOAD_BLR_MATRIX(GRID, A_BLR, N, NB, PR, PC, MYID0, PREFIX, WITH_IPIV, OK)
    TYPE(GRID_T), INTENT(IN) :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(OUT) :: A_BLR
    INTEGER, INTENT(IN) :: N, NB, PR, PC, MYID0
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    LOGICAL, INTENT(IN) :: WITH_IPIV
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: U, li2, lj2, I2, J2, ND, K2, MI, IOS, IERR2
    INTEGER :: F_CMAP, F_N, F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
    INTEGER :: R_I, R_J, R_M, R_N, R_RANK
    LOGICAL :: R_LOWRANK, OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    OK_LOCAL = .TRUE.
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- avez-vous sauvegarde avec SAVE_MATRIX/SAVE_FACTORED=1 auparavant ?'
      OK_LOCAL = .FALSE.
    ELSE
      ! Nouveau format : marqueur -7411 puis CMAP ; ancien format (premier
      ! entier = N > 0) : distribution cyclique pure (CMAP=0).
      READ(U) F_N
      F_CMAP = 0
      IF (F_N .EQ. -7411) READ(U) F_CMAP, F_N
      READ(U) F_NB, F_NBLK, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_NBLKR, F_NBLKC
      IF (F_N.NE.N .OR. F_NB.NE.NB .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NB/PR/PC ou grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    ! Un seul rang incompatible doit arreter tous les rangs -- sinon les
    ! rangs corrects continueraient seuls vers un appel MPI collectif :
    ! blocage garanti plutot qu'un arret net.
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(U)
      RETURN
    END IF

    A_BLR%N = N; A_BLR%NB = NB
    CALL BLOCK_PARTITION(N, NB, A_BLR%NBLK, A_BLR%BSTART, A_BLR%BEND)
    A_BLR%NBLKR_LOC = BLOCK_NUMROC(A_BLR%NBLK, GRID%PR, GRID%MYROW)
    A_BLR%CMAP = F_CMAP
    A_BLR%NBLKC_LOC = BLR_NBLKC_LOC(A_BLR%CMAP, A_BLR%NBLK, GRID)
    ALLOCATE(A_BLR%DIAGBLK(A_BLR%NBLK))
    ALLOCATE(A_BLR%OFFBLK(MAX(A_BLR%NBLKR_LOC,1), MAX(A_BLR%NBLKC_LOC,1)))

    READ(U) ND
    DO K2 = 1, ND
      READ(U) I2, MI
      ALLOCATE(A_BLR%DIAGBLK(I2)%D(MI,MI))
      READ(U) A_BLR%DIAGBLK(I2)%D
    END DO

    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        ! Meme position diagonale sautee a l'ecriture -- rien a relire,
        ! l'element reste dans son etat par defaut.
        IF (I2 .EQ. J2) CYCLE
        READ(U) R_I, R_J, R_LOWRANK, R_M, R_N, R_RANK
        A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK = R_LOWRANK
        A_BLR%OFFBLK(li2,lj2)%M = R_M
        A_BLR%OFFBLK(li2,lj2)%N = R_N
        A_BLR%OFFBLK(li2,lj2)%RANK = R_RANK
        IF (R_LOWRANK) THEN
          ALLOCATE(A_BLR%OFFBLK(li2,lj2)%UFAC(R_M,R_RANK), A_BLR%OFFBLK(li2,lj2)%VFAC(R_RANK,R_N))
          READ(U) A_BLR%OFFBLK(li2,lj2)%UFAC
          READ(U) A_BLR%OFFBLK(li2,lj2)%VFAC
        ELSE
          ALLOCATE(A_BLR%OFFBLK(li2,lj2)%DENSE(R_M,R_N))
          READ(U) A_BLR%OFFBLK(li2,lj2)%DENSE
        END IF
      END DO
    END DO

    IF (WITH_IPIV) THEN
      ALLOCATE(A_BLR%IPIV(NB, A_BLR%NBLK))
      READ(U) A_BLR%IPIV
    END IF
    CLOSE(U)
  END SUBROUTINE LOAD_BLR_MATRIX

  ! Ouvre un checkpoint RHS/solution en ECRITURE et ecrit son en-tete une
  ! seule fois (N, NRHS, grille, MLOC, NRHS_BATCH_HDR, NBATCHES_HDR) --
  ! reste ouvert pendant toute la boucle de lots, ferme par l'appelant.
  ! NRHS_BATCH_HDR/NBATCHES_HDR permettent a LOAD_RHS de reproduire
  ! exactement le meme decoupage, sans dependre du NRHS_BATCH du run qui
  ! relit.
  SUBROUTINE OPEN_RHS_CHECKPOINT_WRITE(GRID, N, NRHS, PR, PC, MYID0, PREFIX, UNIT, &
                                        NBATCHES_HDR, NRHS_BATCH_HDR, MLOC_HDR, OK)
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER, INTENT(IN) :: N, NRHS, PR, PC, MYID0
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER, INTENT(OUT) :: UNIT
    INTEGER, INTENT(IN) :: NBATCHES_HDR, NRHS_BATCH_HDR, MLOC_HDR
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: IOS, IERR2
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=UNIT, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    OK_LOCAL = (IOS .EQ. 0)
    IF (.NOT. OK_LOCAL) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible de creer ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- le repertoire CHECKPOINT_DIR doit deja exister.'
    END IF
    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(UNIT, STATUS='DELETE')
      RETURN
    END IF
    WRITE(UNIT) N, NRHS, PR, PC, GRID%MYROW, GRID%MYCOL, MYID0, MLOC_HDR, NRHS_BATCH_HDR, NBATCHES_HDR
  END SUBROUTINE OPEN_RHS_CHECKPOINT_WRITE

  ! Ecrit un lot (NLOC_B colonnes locales de BDATA) -- appelee pour
  ! SAVE_RHS (apres USER_FILL_B) et SAVE_SOLUTION (apres la resolution) --
  ! meme routine, deux fichiers/unites distincts.
  SUBROUTINE SAVE_RHS_BATCH(UNIT, NLOC_B, BDATA, ML)
    INTEGER, INTENT(IN) :: UNIT, NLOC_B, ML
    COMPLEX(WP), INTENT(IN) :: BDATA(MAX(ML,1),MAX(NLOC_B,1))
    WRITE(UNIT) NLOC_B
    IF (ML.GT.0 .AND. NLOC_B.GT.0) WRITE(UNIT) BDATA(1:ML,1:NLOC_B)
  END SUBROUTINE SAVE_RHS_BATCH

  ! Ouvre un checkpoint RHS en LECTURE, verifie son en-tete contre le run
  ! courant, et retourne NBATCHES_OUT/NRHS_BATCH_OUT lus depuis le fichier
  ! -- c'est ce decoupage, pas le NRHS_BATCH de ce run, qui pilote la
  ! boucle de lots. Meme logique OK/MPI_ALLREDUCE que LOAD_BLR_MATRIX.
  SUBROUTINE OPEN_RHS_CHECKPOINT_READ(GRID, N, NRHS, PR, PC, MYID0, MLOC, PREFIX, UNIT, &
                                       NBATCHES_OUT, NRHS_BATCH_OUT, OK)
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER, INTENT(IN) :: N, NRHS, PR, PC, MYID0, MLOC
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX
    INTEGER, INTENT(OUT) :: UNIT, NBATCHES_OUT, NRHS_BATCH_OUT
    LOGICAL, INTENT(OUT) :: OK
    INTEGER :: IOS, IERR2
    INTEGER :: F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC
    LOGICAL :: OK_LOCAL
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME
    OK_LOCAL = .TRUE.
    NBATCHES_OUT = 0; NRHS_BATCH_OUT = 0
    WRITE(RSTR,'(I6.6)') MYID0
    FNAME = TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
    OPEN(NEWUNIT=UNIT, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A,A,A,I0,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
           ' (rang ', MYID0, ') -- avez-vous sauvegarde avec SAVE_RHS=1 auparavant ?'
      OK_LOCAL = .FALSE.
    ELSE
      READ(UNIT) F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, NRHS_BATCH_OUT, NBATCHES_OUT
      IF (F_N.NE.N .OR. F_NRHS.NE.NRHS .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
          F_MYROW.NE.GRID%MYROW .OR. F_MYCOL.NE.GRID%MYCOL .OR. F_MYID.NE.MYID0 .OR. F_MLOC.NE.MLOC) THEN
        WRITE(6,'(A,A,A,I0,A)') ' ERREUR : checkpoint ', TRIM(FNAME), &
             ' incompatible avec le run courant (rang ', MYID0, &
             ') -- N/NRHS/PR/PC/grille different de la sauvegarde.'
        OK_LOCAL = .FALSE.
      END IF
    END IF

    CALL MPI_ALLREDUCE(OK_LOCAL, OK, 1, MPI_LOGICAL, MPI_LAND, GRID%COMM, IERR2)
    IF (.NOT. OK) THEN
      IF (OK_LOCAL) CLOSE(UNIT)
    END IF
  END SUBROUTINE OPEN_RHS_CHECKPOINT_READ

  ! Relit un lot (NLOC_B colonnes locales attendues, connues de l'appelant
  ! via le meme decoupage RHS_LOCAL_RANGE_MBM que la sauvegarde) dans
  ! BDATA.
  SUBROUTINE LOAD_RHS_BATCH(UNIT, NLOC_B, BDATA, ML)
    INTEGER, INTENT(IN) :: UNIT, NLOC_B, ML
    COMPLEX(WP), INTENT(OUT) :: BDATA(MAX(ML,1),MAX(NLOC_B,1))
    INTEGER :: NLOC_B_FILE
    READ(UNIT) NLOC_B_FILE
    IF (ML.GT.0 .AND. NLOC_B_FILE.GT.0) READ(UNIT) BDATA(1:ML,1:NLOC_B_FILE)
  END SUBROUTINE LOAD_RHS_BATCH

  ! Relit une sauvegarde --save-solution et en reconstruit la sous-matrice
  ! GLOBALE X(IMIN:IMAX,NMIN:NMAX) demandee -- simple appel de
  ! LOAD_RHS_OR_SOLUTION_RANGE avec le prefixe 'saved_solution' (voir son
  ! en-tete pour le fonctionnement et les limites communes aux deux).
  SUBROUTINE LOAD_SOLUTION(CHECKPOINT_DIR, N, NRHS, PR, PC, NB, IMIN, IMAX, NMIN, NMAX, XOUT, OK)
    CHARACTER(LEN=*), INTENT(IN) :: CHECKPOINT_DIR
    INTEGER, INTENT(IN)          :: N, NRHS, PR, PC, NB, IMIN, IMAX, NMIN, NMAX
    COMPLEX(WP), INTENT(OUT)     :: XOUT(IMIN:IMAX, NMIN:NMAX)
    LOGICAL, INTENT(OUT)         :: OK
    CALL LOAD_RHS_OR_SOLUTION_RANGE('saved_solution', CHECKPOINT_DIR, N, NRHS, PR, PC, NB, &
                                     IMIN, IMAX, NMIN, NMAX, XOUT, OK)
  END SUBROUTINE LOAD_SOLUTION

  ! Coeur de LOAD_SOLUTION : reconstruit la sous-matrice
  ! GLOBALE XOUT(IMIN:IMAX,NMIN:NMAX) d'une sauvegarde par lot (--save-rhs
  ! ou --save-solution, memes en-tetes/format, seul PREFIX distingue les
  ! deux jeux de fichiers), en parcourant tous les PR*PC fichiers
  ! <PREFIX>_rank*.dat. SEQUENTIEL (aucun GRID/A_BLR requis, aucune
  ! communication MPI) -- pensee comme outil de post-traitement, appelable
  ! independamment (et apres coup) du run MPI qui a produit la sauvegarde
  ! -- typiquement par le seul rang 0, ou par un petit executable separe.
  ! Reutilise EXACTEMENT les memes fonctions de decoupage
  ! (BLOCK_PARTITION/BLOCK_NUMROC pour les lignes bloc-cycliques,
  ! RHS_LOCAL_RANGE_MBM pour les colonnes par lot) que le run d'origine,
  ! pour garantir la meme correspondance globale/locale qu'a l'ecriture.
  !
  ! N/NRHS/PR/PC doivent correspondre au run d'origine (verifies contre
  ! l'en-tete de chaque fichier, comme OPEN_RHS_CHECKPOINT_READ). NB n'est
  ! PAS stocke dans l'en-tete : s'il ne correspond pas au run d'origine, le
  ! MLOC recalcule ici ne collera plus a MLOC_HDR lu dans le fichier et OK
  ! sera mis a .FALSE. -- mais la valeur de NB elle-meme reste a la charge
  ! de l'appelant (aucun moyen de la deviner depuis les fichiers seuls).
  SUBROUTINE LOAD_RHS_OR_SOLUTION_RANGE(PREFIX, CHECKPOINT_DIR, N, NRHS, PR, PC, NB, &
                                         IMIN, IMAX, NMIN, NMAX, XOUT, OK)
    CHARACTER(LEN=*), INTENT(IN) :: PREFIX, CHECKPOINT_DIR
    INTEGER, INTENT(IN)          :: N, NRHS, PR, PC, NB, IMIN, IMAX, NMIN, NMAX
    COMPLEX(WP), INTENT(OUT)     :: XOUT(IMIN:IMAX, NMIN:NMAX)
    LOGICAL, INTENT(OUT)         :: OK
    INTEGER, ALLOCATABLE :: BSTART(:), BEND(:)
    INTEGER :: NBLK, MYROW, MYCOL, RANKID, IOS, UNIT
    INTEGER :: NBLKR_LOC, MLOC_EXPECT, li, I, MI, RS
    INTEGER :: F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC
    INTEGER :: NRHS_BATCH_EFF, NBATCHES, B, RSTART_B, REND_B, WIDTH_B
    INTEGER :: LOCAL_OFFSET_B, NRHS_LOC_B, GSTART_B, GEND_B, NLOC_B_FILE
    INTEGER :: IROW_LO, IROW_HI, JCOL_LO, JCOL_HI
    LOGICAL :: ROW_OVERLAP
    COMPLEX(WP), ALLOCATABLE :: BUF(:,:)
    CHARACTER(LEN=16) :: RSTR
    CHARACTER(LEN=160) :: FNAME

    OK = .TRUE.
    IF (IMIN.LT.1 .OR. IMAX.GT.N .OR. IMIN.GT.IMAX .OR. &
        NMIN.LT.1 .OR. NMAX.GT.NRHS .OR. NMIN.GT.NMAX) THEN
      WRITE(6,'(A,A,A)') ' ERREUR (', TRIM(PREFIX), &
           ') : plage [IMIN,IMAX]x[NMIN,NMAX] invalide (hors de [1,N]x[1,NRHS]).'
      OK = .FALSE.
      RETURN
    END IF

    CALL BLOCK_PARTITION(N, NB, NBLK, BSTART, BEND)

    DO MYROW = 0, PR-1
      NBLKR_LOC = BLOCK_NUMROC(NBLK, PR, MYROW)
      MLOC_EXPECT = 0
      ROW_OVERLAP = .FALSE.
      DO li = 1, NBLKR_LOC
        I = (li-1)*PR + MYROW + 1
        MLOC_EXPECT = MLOC_EXPECT + (BEND(I)-BSTART(I)+1)
        IF (BEND(I).GE.IMIN .AND. BSTART(I).LE.IMAX) ROW_OVERLAP = .TRUE.
      END DO
      ! Ce processus-ligne ne possede aucune ligne demandee -- inutile
      ! d'ouvrir ses PC fichiers.
      IF (.NOT. ROW_OVERLAP) CYCLE

      DO MYCOL = 0, PC-1
        RANKID = MYROW*PC + MYCOL
        WRITE(RSTR,'(I6.6)') RANKID
        FNAME = TRIM(CHECKPOINT_DIR)//'/'//TRIM(PREFIX)//'_rank'//TRIM(RSTR)//'.dat'
        OPEN(NEWUNIT=UNIT, FILE=TRIM(FNAME), ACCESS='STREAM', FORM='UNFORMATTED', &
             STATUS='OLD', ACTION='READ', IOSTAT=IOS)
        IF (IOS .NE. 0) THEN
          WRITE(6,'(A,A,A)') ' ERREUR : impossible d''ouvrir ', TRIM(FNAME), &
               ' -- cette sauvegarde a-t-elle bien ete faite au prealable (option --save-... correspondante) ?'
          OK = .FALSE.
          RETURN
        END IF

        READ(UNIT) F_N, F_NRHS, F_PR, F_PC, F_MYROW, F_MYCOL, F_MYID, F_MLOC, NRHS_BATCH_EFF, NBATCHES
        IF (F_N.NE.N .OR. F_NRHS.NE.NRHS .OR. F_PR.NE.PR .OR. F_PC.NE.PC .OR. &
            F_MYROW.NE.MYROW .OR. F_MYCOL.NE.MYCOL .OR. F_MYID.NE.RANKID .OR. F_MLOC.NE.MLOC_EXPECT) THEN
          WRITE(6,'(A,A,A)') ' ERREUR : en-tete incoherent dans ', TRIM(FNAME), &
               ' -- N/NRHS/PR/PC/NB ne correspondent pas au run d''origine.'
          CLOSE(UNIT)
          OK = .FALSE.
          RETURN
        END IF

        DO B = 1, NBATCHES
          RSTART_B = (B-1)*NRHS_BATCH_EFF + 1
          REND_B   = MIN(B*NRHS_BATCH_EFF, NRHS)
          WIDTH_B  = REND_B - RSTART_B + 1
          CALL RHS_LOCAL_RANGE_MBM(WIDTH_B, PC, MYCOL, LOCAL_OFFSET_B, NRHS_LOC_B)
          GSTART_B = RSTART_B + LOCAL_OFFSET_B - 1
          GEND_B   = GSTART_B + NRHS_LOC_B - 1

          READ(UNIT) NLOC_B_FILE
          IF (MLOC_EXPECT.GT.0 .AND. NLOC_B_FILE.GT.0) THEN
            ! Toujours lire le lot en entier (flux STREAM sans marqueurs
            ! d'enregistrement a sauter) -- seule l'intersection demandee
            ! est ensuite recopiee dans XOUT.
            ALLOCATE(BUF(MLOC_EXPECT, NLOC_B_FILE))
            READ(UNIT) BUF
            JCOL_LO = MAX(GSTART_B, NMIN)
            JCOL_HI = MIN(GEND_B, NMAX)
            IF (JCOL_HI .GE. JCOL_LO) THEN
              RS = 0
              DO li = 1, NBLKR_LOC
                I = (li-1)*PR + MYROW + 1
                MI = BEND(I) - BSTART(I) + 1
                IROW_LO = MAX(BSTART(I), IMIN)
                IROW_HI = MIN(BEND(I), IMAX)
                IF (IROW_HI .GE. IROW_LO) THEN
                  XOUT(IROW_LO:IROW_HI, JCOL_LO:JCOL_HI) = &
                       BUF(RS+(IROW_LO-BSTART(I))+1 : RS+(IROW_HI-BSTART(I))+1, &
                           JCOL_LO-GSTART_B+1 : JCOL_HI-GSTART_B+1)
                END IF
                RS = RS + MI
              END DO
            END IF
            DEALLOCATE(BUF)
          END IF
        END DO
        CLOSE(UNIT)
      END DO
    END DO
  END SUBROUTINE LOAD_RHS_OR_SOLUTION_RANGE

  ! Recalcule N_COMPRESSED/N_OFFDIAG/RANK_SUM depuis A_BLR%OFFBLK (meme
  ! calcul qu'en interne dans MPI_ASSEMBLE_BLR_MATRIX/LU_FACTOR_PIVOT) --
  ! necessaire apres LOAD_MATRIX/LOAD_FACTORED, qui sautent justement
  ! l'appel qui produirait normalement ces stats.
  SUBROUTINE RECOMPUTE_COMPRESSION_STATS(GRID, A_BLR, NC_OUT, NO_OUT, RS_OUT)
    TYPE(GRID_T), INTENT(IN) :: GRID
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: A_BLR
    INTEGER, INTENT(OUT) :: NC_OUT, NO_OUT
    INTEGER(8), INTENT(OUT) :: RS_OUT
    INTEGER :: li2, lj2, I2, J2, NC_LOCAL, NO_LOCAL, IERR2
    INTEGER(8) :: RS_LOCAL
    NC_LOCAL = 0; NO_LOCAL = 0; RS_LOCAL = 0_8
    DO li2 = 1, A_BLR%NBLKR_LOC
      I2 = (li2-1)*GRID%PR + GRID%MYROW + 1
      DO lj2 = 1, A_BLR%NBLKC_LOC
        J2 = BLR_COL_GLOBAL(A_BLR%CMAP, lj2, GRID)
        IF (I2 .EQ. J2) CYCLE
        NO_LOCAL = NO_LOCAL + 1
        IF (A_BLR%OFFBLK(li2,lj2)%IS_LOWRANK) THEN
          NC_LOCAL = NC_LOCAL + 1
          RS_LOCAL = RS_LOCAL + A_BLR%OFFBLK(li2,lj2)%RANK
        ELSE
          RS_LOCAL = RS_LOCAL + MIN(A_BLR%OFFBLK(li2,lj2)%M, A_BLR%OFFBLK(li2,lj2)%N)
        END IF
      END DO
    END DO
    CALL MPI_ALLREDUCE(NC_LOCAL, NC_OUT, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(NO_LOCAL, NO_OUT, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR2)
    CALL MPI_ALLREDUCE(RS_LOCAL, RS_OUT, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR2)
  END SUBROUTINE RECOMPUTE_COMPRESSION_STATS

  ! Taux de compression REEL (stockage bas-rang effectif / dense
  ! equivalent) -- pas juste 2*NB*RANK_SUM/(NB^2*N_OFFDIAG) : RANK_SUM
  ! melange le rang reel des blocs compresses et un "rang virtuel"=NB pour
  ! chaque bloc dense (convention de comptage interne, pas un vrai
  ! stockage UV), qui occupe NB^2 elements et non 2*NB^2 (une formule
  ! naive donnerait a tort 200% quand tout est dense). La contribution des
  ! blocs denses est donc retiree de RANK_SUM avant de les recompter a
  ! leur cout reel NB^2.
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NB, NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NB, NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

END MODULE CHECKPOINT_IO_MOD
