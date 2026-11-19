!===============================================================================
! Relecture des fichiers de solution ecrits par solveur_blr --save-solution
! (un fichier par processus : <DIR>/saved_solution_rank<NNNNNN>.dat).
!
! FORMAT d'un fichier (Fortran ACCESS='STREAM', sans marqueurs) :
!   en-tete : 10 entiers 32 bits
!             N, NRHS, PR, PC, MYROW, MYCOL, RANG, MLOC, NRHS_BATCH, NBATCHES
!   puis, pour chaque lot b = 1..NBATCHES :
!             NLOC (entier 32 bits) : colonnes du lot possedees par ce processus
!             MLOC x NLOC complexes simple precision (8 octets), par colonnes
!
! REPARTITION (identique au solveur) :
!   lignes   : blocs de NB lignes ; le bloc I appartient a la ligne de grille
!              MOD(I-1,PR) ; les blocs locaux sont concatenes dans l'ordre.
!   colonnes : chaque lot [RSTART_B, REND_B] est decoupe en PC tranches
!              contigues (les MOD(largeur,PC) premieres ont une colonne de plus).
!   fichier  : rang = MYROW*PC + MYCOL.
!
! Module autonome : ni MPI, ni bibliotheque du solveur. Seules les colonnes
! demandees sont lues (positionnement direct dans le fichier, POS=).
!===============================================================================
MODULE SOLUTION_READER_MOD
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SOLUTION_HEADER_T, READ_SOLUTION_HEADER, READ_SOLUTION_RANGE

  INTEGER, PARAMETER :: SP = KIND(1.0)

  TYPE :: SOLUTION_HEADER_T
    INTEGER :: N = 0, NRHS = 0, PR = 0, PC = 0, NRHS_BATCH = 0, NBATCHES = 0
  END TYPE SOLUTION_HEADER_T

CONTAINS

  FUNCTION RANK_FILE(DIR, RANKID) RESULT(FNAME)
    CHARACTER(LEN=*), INTENT(IN) :: DIR
    INTEGER, INTENT(IN) :: RANKID
    CHARACTER(LEN=1024) :: FNAME
    CHARACTER(LEN=6) :: RSTR
    WRITE(RSTR,'(I6.6)') RANKID
    FNAME = TRIM(DIR)//'/saved_solution_rank'//RSTR//'.dat'
  END FUNCTION RANK_FILE

  ! Lit l'en-tete du fichier du rang 0 (N, NRHS, grille, decoupage en lots).
  SUBROUTINE READ_SOLUTION_HEADER(DIR, H, OK)
    CHARACTER(LEN=*), INTENT(IN)          :: DIR
    TYPE(SOLUTION_HEADER_T), INTENT(OUT)  :: H
    LOGICAL, INTENT(OUT)                  :: OK
    INTEGER :: U, IOS, HDR(10)
    OPEN(NEWUNIT=U, FILE=TRIM(RANK_FILE(DIR,0)), ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    OK = (IOS .EQ. 0)
    IF (.NOT. OK) THEN
      WRITE(0,'(A,A)') ' ERREUR : impossible d''ouvrir ', TRIM(RANK_FILE(DIR,0))
      RETURN
    END IF
    READ(U, IOSTAT=IOS) HDR
    CLOSE(U)
    OK = (IOS .EQ. 0)
    IF (.NOT. OK) THEN
      WRITE(0,'(A)') ' ERREUR : en-tete illisible dans le fichier du rang 0.'
      RETURN
    END IF
    H%N = HDR(1); H%NRHS = HDR(2); H%PR = HDR(3); H%PC = HDR(4)
    H%NRHS_BATCH = HDR(9); H%NBATCHES = HDR(10)
  END SUBROUTINE READ_SOLUTION_HEADER

  ! Colonnes du lot (largeur W) possedees par la colonne de grille C :
  ! premiere colonne relative RSTART (1-based) et nombre RLOC.
  SUBROUTINE RHS_RANGE(W, PC, C, RSTART, RLOC)
    INTEGER, INTENT(IN)  :: W, PC, C
    INTEGER, INTENT(OUT) :: RSTART, RLOC
    INTEGER :: BASE, REM
    BASE = W / PC
    REM  = MOD(W, PC)
    IF (C .LT. REM) THEN
      RLOC = BASE + 1; RSTART = C*(BASE+1) + 1
    ELSE
      RLOC = BASE; RSTART = REM*(BASE+1) + (C-REM)*BASE + 1
    END IF
  END SUBROUTINE RHS_RANGE

  ! X(IMIN:IMAX, JMIN:JMAX) <- solution globale. NB doit etre celui du
  ! calcul : il n'est pas stocke dans les fichiers, et un NB faux n'est
  ! detecte que s'il change MLOC (nombre de lignes par processus) -- par
  ! exemple N=1000, PR=2 donne MLOC=500 avec NB=50 comme avec NB=100.
  ! OK = .FALSE. en cas d'erreur (message sur la sortie d'erreur).
  SUBROUTINE READ_SOLUTION_RANGE(DIR, H, NB, IMIN, IMAX, JMIN, JMAX, X, OK)
    CHARACTER(LEN=*), INTENT(IN)        :: DIR
    TYPE(SOLUTION_HEADER_T), INTENT(IN) :: H
    INTEGER, INTENT(IN)                 :: NB, IMIN, IMAX, JMIN, JMAX
    COMPLEX(SP), INTENT(OUT)            :: X(IMIN:IMAX, JMIN:JMAX)
    LOGICAL, INTENT(OUT)                :: OK
    INTEGER :: NBLK, MYROW, MYCOL, RANKID, U, IOS, HDR(10)
    INTEGER :: MLOC, I, MI, B, RSTART_B, REND_B, W, OFFC, NLOC, NLOC_F, GS, J, RS
    INTEGER :: I1, I2
    INTEGER(8) :: POS_BATCH, POS
    COMPLEX(SP), ALLOCATABLE :: COL(:)
    LOGICAL :: ROWS_HIT

    OK = .FALSE.
    IF (NB .LE. 0) THEN
      WRITE(0,'(A)') ' ERREUR : NB doit etre strictement positif.'
      RETURN
    END IF
    NBLK = (H%N + NB - 1) / NB

    DO MYROW = 0, H%PR-1
      ! MLOC et presence d'au moins une ligne demandee pour cette ligne de grille
      MLOC = 0
      ROWS_HIT = .FALSE.
      DO I = MYROW+1, NBLK, H%PR
        MI = MIN(I*NB, H%N) - (I-1)*NB
        MLOC = MLOC + MI
        IF (MIN(I*NB,H%N) .GE. IMIN .AND. (I-1)*NB+1 .LE. IMAX) ROWS_HIT = .TRUE.
      END DO
      IF (.NOT. ROWS_HIT) CYCLE
      ALLOCATE(COL(MAX(MLOC,1)))

      DO MYCOL = 0, H%PC-1
        RANKID = MYROW*H%PC + MYCOL
        OPEN(NEWUNIT=U, FILE=TRIM(RANK_FILE(DIR,RANKID)), ACCESS='STREAM', FORM='UNFORMATTED', &
             STATUS='OLD', ACTION='READ', IOSTAT=IOS)
        IF (IOS .NE. 0) THEN
          WRITE(0,'(A,A)') ' ERREUR : impossible d''ouvrir ', TRIM(RANK_FILE(DIR,RANKID))
          RETURN
        END IF
        READ(U, IOSTAT=IOS) HDR
        IF (IOS.NE.0 .OR. HDR(1).NE.H%N .OR. HDR(2).NE.H%NRHS .OR. HDR(3).NE.H%PR .OR. HDR(4).NE.H%PC &
            .OR. HDR(5).NE.MYROW .OR. HDR(6).NE.MYCOL .OR. HDR(7).NE.RANKID) THEN
          WRITE(0,'(A,A)') ' ERREUR : en-tete incoherent dans ', TRIM(RANK_FILE(DIR,RANKID))
          CLOSE(U); RETURN
        END IF
        IF (HDR(8) .NE. MLOC) THEN
          WRITE(0,'(A,I0,A)') ' ERREUR : nombre de lignes locales incoherent -- NB=', NB, &
               ' est-il bien celui du calcul ?'
          CLOSE(U); RETURN
        END IF

        ! Parcours des lots : en-tete de 40 octets, puis par lot 4 + 8*MLOC*NLOC octets.
        POS_BATCH = 41_8
        DO B = 1, H%NBATCHES
          RSTART_B = (B-1)*H%NRHS_BATCH + 1
          REND_B = MIN(B*H%NRHS_BATCH, H%NRHS)
          W = REND_B - RSTART_B + 1
          CALL RHS_RANGE(W, H%PC, MYCOL, OFFC, NLOC)
          GS = RSTART_B + OFFC - 1           ! premiere colonne globale de ce processus
          READ(U, POS=POS_BATCH, IOSTAT=IOS) NLOC_F
          IF (IOS.NE.0 .OR. NLOC_F.NE.NLOC) THEN
            WRITE(0,'(A,A)') ' ERREUR : lot illisible ou incoherent dans ', TRIM(RANK_FILE(DIR,RANKID))
            CLOSE(U); RETURN
          END IF
          DO J = MAX(GS,JMIN), MIN(GS+NLOC-1,JMAX)
            POS = POS_BATCH + 4_8 + 8_8*INT(MLOC,8)*INT(J-GS,8)
            IF (MLOC .GT. 0) READ(U, POS=POS, IOSTAT=IOS) COL(1:MLOC)
            IF (IOS .NE. 0) THEN
              WRITE(0,'(A,A)') ' ERREUR : donnees tronquees dans ', TRIM(RANK_FILE(DIR,RANKID))
              CLOSE(U); RETURN
            END IF
            RS = 0
            DO I = MYROW+1, NBLK, H%PR
              MI = MIN(I*NB, H%N) - (I-1)*NB
              I1 = MAX((I-1)*NB+1, IMIN); I2 = MIN(I*NB, H%N, IMAX)
              IF (I2 .GE. I1) X(I1:I2, J) = COL(RS + I1-(I-1)*NB : RS + I2-(I-1)*NB)
              RS = RS + MI
            END DO
          END DO
          POS_BATCH = POS_BATCH + 4_8 + 8_8*INT(MLOC,8)*INT(NLOC,8)
        END DO
        CLOSE(U)
      END DO
      DEALLOCATE(COL)
    END DO
    OK = .TRUE.
  END SUBROUTINE READ_SOLUTION_RANGE

END MODULE SOLUTION_READER_MOD
