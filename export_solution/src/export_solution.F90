!===============================================================================
! export_solution : relit la solution X ecrite par
!   solveur_blr ... --save-solution [--checkpoint-dir=DIR]
! (un fichier par processus) et l'ecrit dans un seul fichier, au(x)
! format(s) choisi(s) : texte, binaire et/ou NetCDF.
!
! Usage : ./export_solution --nb=NB --format=txt,bin,netcdf [options]
!         (--help pour le detail). Programme sequentiel : pas de mpirun.
!
! La solution est lue et ecrite par paquets de colonnes (--chunk) : la
! memoire utilisee est de l'ordre de (lignes exportees) x (chunk) complexes,
! quel que soit NRHS.
!===============================================================================
PROGRAM EXPORT_SOLUTION
  USE SOLUTION_READER_MOD
#ifdef WITH_NETCDF
  USE NETCDF_C_MOD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
#endif
  IMPLICIT NONE

  INTEGER, PARAMETER :: SP = KIND(1.0)
  CHARACTER(LEN=1024) :: DIR, OUTP, ARG, VAL
  CHARACTER(LEN=64)   :: KEY, FMTS
  INTEGER :: NB, IMIN, IMAX, JMIN, JMAX, CHUNK, IARG, IOS, J0, J1, I, J
  LOGICAL :: DO_TXT, DO_BIN, DO_NC, OK, ROWS_SET, COLS_SET
  TYPE(SOLUTION_HEADER_T) :: H
  COMPLEX(SP), ALLOCATABLE :: XC(:,:)
  INTEGER :: UTXT, UBIN, NROWS, NCOLS
#ifdef WITH_NETCDF
  INTEGER(C_INT) :: NCID, DROW, DCOL, VRE, VIM, VROW, VCOL
  INTEGER(C_INT), ALLOCATABLE :: IDX(:)
#endif

  DIR = '.'; OUTP = 'solution'; FMTS = ''
  NB = 0; CHUNK = 16
  ROWS_SET = .FALSE.; COLS_SET = .FALSE.

  ! --- options ---
  IF (COMMAND_ARGUMENT_COUNT() .EQ. 0) THEN
    CALL USAGE(); STOP
  END IF
  DO IARG = 1, COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(IARG, ARG)
    IF (TRIM(ARG).EQ.'--help' .OR. TRIM(ARG).EQ.'-h') THEN
      CALL USAGE(); STOP
    END IF
    IF (ARG(1:2) .NE. '--' .OR. INDEX(ARG,'=') .EQ. 0) CALL FAIL('argument invalide "'//TRIM(ARG)//'" (forme --cle=valeur, voir --help)')
    KEY = ARG(3:INDEX(ARG,'=')-1)
    VAL = ARG(INDEX(ARG,'=')+1:)
    IOS = 0
    SELECT CASE (TRIM(KEY))
    CASE ('dir');    DIR = VAL
    CASE ('out');    OUTP = VAL
    CASE ('format'); FMTS = VAL
    CASE ('nb');     READ(VAL,*,IOSTAT=IOS) NB
    CASE ('chunk');  READ(VAL,*,IOSTAT=IOS) CHUNK
    CASE ('rows');   CALL READ_RANGE(VAL, IMIN, IMAX, IOS); ROWS_SET = .TRUE.
    CASE ('cols');   CALL READ_RANGE(VAL, JMIN, JMAX, IOS); COLS_SET = .TRUE.
    CASE DEFAULT;    CALL FAIL('option inconnue --'//TRIM(KEY)//' (voir --help)')
    END SELECT
    IF (IOS .NE. 0) CALL FAIL('valeur invalide pour --'//TRIM(KEY))
  END DO

  IF (NB .LE. 0) CALL FAIL('--nb=NB est obligatoire (meme NB que le calcul)')
  IF (CHUNK .LE. 0) CALL FAIL('--chunk doit etre strictement positif')
  IF (LEN_TRIM(FMTS) .EQ. 0) CALL FAIL('--format est obligatoire : txt, bin, netcdf (separes par des virgules) ou all')
  DO_TXT = HAS_WORD(FMTS, 'txt')    .OR. HAS_WORD(FMTS, 'all')
  DO_BIN = HAS_WORD(FMTS, 'bin')    .OR. HAS_WORD(FMTS, 'all')
  DO_NC  = HAS_WORD(FMTS, 'netcdf') .OR. HAS_WORD(FMTS, 'nc') .OR. HAS_WORD(FMTS, 'all')
  IF (.NOT. (DO_TXT .OR. DO_BIN .OR. DO_NC)) CALL FAIL('--format : valeurs reconnues txt, bin, netcdf, all')
#ifndef WITH_NETCDF
  IF (DO_NC) CALL FAIL('format netcdf indisponible : programme compile sans NetCDF (make NETCDF=1)')
#endif

  ! --- en-tete et plage exportee ---
  CALL READ_SOLUTION_HEADER(DIR, H, OK)
  IF (.NOT. OK) STOP 1
  IF (.NOT. ROWS_SET) THEN
    IMIN = 1; IMAX = H%N
  END IF
  IF (.NOT. COLS_SET) THEN
    JMIN = 1; JMAX = H%NRHS
  END IF
  IF (IMIN.LT.1 .OR. IMAX.GT.H%N .OR. IMIN.GT.IMAX) CALL FAIL('--rows hors de [1,N]')
  IF (JMIN.LT.1 .OR. JMAX.GT.H%NRHS .OR. JMIN.GT.JMAX) CALL FAIL('--cols hors de [1,NRHS]')
  NROWS = IMAX - IMIN + 1
  NCOLS = JMAX - JMIN + 1

  WRITE(6,'(A,A)') ' Solution lue dans : ', TRIM(DIR)
  WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' N=', H%N, ' NRHS=', H%NRHS, ' grille ', H%PR, 'x', H%PC, ' NB=', NB
  WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A)') ' Export de X(', IMIN, ':', IMAX, ', ', JMIN, ':', JMAX, ')'

  ! --- ouverture des sorties ---
  IF (DO_TXT) THEN
    OPEN(NEWUNIT=UTXT, FILE=TRIM(OUTP)//'.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    IF (IOS .NE. 0) CALL FAIL('impossible de creer '//TRIM(OUTP)//'.txt')
    WRITE(UTXT,'(A)') '# Solution X du solveur BLR (complexe simple precision)'
    WRITE(UTXT,'(A,I0,A,I0)') '# N = ', H%N, ' , NRHS = ', H%NRHS
    WRITE(UTXT,'(A,I0,A,I0,A,I0,A,I0)') '# lignes ', IMIN, ' a ', IMAX, ' , colonnes ', JMIN, ' a ', JMAX
    WRITE(UTXT,'(A)') '# colonnes du fichier : i  j  Re(X(i,j))  Im(X(i,j))   (ordre : colonne par colonne)'
  END IF
  IF (DO_BIN) THEN
    OPEN(NEWUNIT=UBIN, FILE=TRIM(OUTP)//'.bin', ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
    IF (IOS .NE. 0) CALL FAIL('impossible de creer '//TRIM(OUTP)//'.bin')
    WRITE(UBIN) NROWS, NCOLS, IMIN, JMIN
  END IF
#ifdef WITH_NETCDF
  IF (DO_NC) THEN
    OK = NC_CHECK(NC_CREATE_FILE(TRIM(OUTP)//'.nc', NCID), 'creation de '//TRIM(OUTP)//'.nc')
    IF (OK) OK = NC_CHECK(NC_DEF_DIM_F(NCID, 'row', NROWS, DROW), 'dimension row')
    IF (OK) OK = NC_CHECK(NC_DEF_DIM_F(NCID, 'col', NCOLS, DCOL), 'dimension col')
    ! convention C : (col, row) <=> X(row, col) en Fortran
    IF (OK) OK = NC_CHECK(NC_DEF_VAR_F(NCID, 'x_re', NC_FLOAT, [DCOL, DROW], VRE), 'variable x_re')
    IF (OK) OK = NC_CHECK(NC_DEF_VAR_F(NCID, 'x_im', NC_FLOAT, [DCOL, DROW], VIM), 'variable x_im')
    IF (OK) OK = NC_CHECK(NC_DEF_VAR_F(NCID, 'row', NC_INT, [DROW], VROW), 'variable row')
    IF (OK) OK = NC_CHECK(NC_DEF_VAR_F(NCID, 'col', NC_INT, [DCOL], VCOL), 'variable col')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_TEXT_F(NCID, VRE, 'long_name', 'partie reelle de X'), 'attribut')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_TEXT_F(NCID, VIM, 'long_name', 'partie imaginaire de X'), 'attribut')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_TEXT_F(NCID, VROW, 'long_name', 'indice global de ligne (1-based)'), 'attribut')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_TEXT_F(NCID, VCOL, 'long_name', 'indice global de second membre (1-based)'), 'attribut')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_TEXT_F(NCID, NC_GLOBAL, 'title', 'Solution X du solveur BLR, X = x_re + i*x_im'), 'attribut')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_INT_F(NCID, NC_GLOBAL, 'N', H%N), 'attribut')
    IF (OK) OK = NC_CHECK(NC_PUT_ATT_INT_F(NCID, NC_GLOBAL, 'NRHS', H%NRHS), 'attribut')
    IF (OK) OK = NC_CHECK(NC_ENDDEF_F(NCID), 'fin de definition')
    IF (OK) THEN
      ALLOCATE(IDX(MAX(NROWS,NCOLS)))
      IDX(1:NROWS) = [(I, I = IMIN, IMAX)]
      OK = NC_CHECK(NC_PUT_VARA_INT_F(NCID, VROW, [0], [NROWS], IDX), 'ecriture de row')
      IDX(1:NCOLS) = [(J, J = JMIN, JMAX)]
      IF (OK) OK = NC_CHECK(NC_PUT_VARA_INT_F(NCID, VCOL, [0], [NCOLS], IDX), 'ecriture de col')
      DEALLOCATE(IDX)
    END IF
    IF (.NOT. OK) STOP 1
  END IF
#endif

  ! --- lecture et ecriture par paquets de colonnes ---
  DO J0 = JMIN, JMAX, CHUNK
    J1 = MIN(J0+CHUNK-1, JMAX)
    ALLOCATE(XC(IMIN:IMAX, J0:J1))
    CALL READ_SOLUTION_RANGE(DIR, H, NB, IMIN, IMAX, J0, J1, XC, OK)
    IF (.NOT. OK) STOP 1
    IF (DO_TXT) THEN
      DO J = J0, J1
        DO I = IMIN, IMAX
          WRITE(UTXT,'(I0,1X,I0,2(1X,ES16.8))') I, J, REAL(XC(I,J)), AIMAG(XC(I,J))
        END DO
      END DO
    END IF
    IF (DO_BIN) WRITE(UBIN) XC
#ifdef WITH_NETCDF
    IF (DO_NC) THEN
      OK = NC_CHECK(NC_PUT_VARA_FLOAT_F(NCID, VRE, [J0-JMIN, 0], [J1-J0+1, NROWS], REAL(XC, C_FLOAT)), 'ecriture de x_re')
      IF (OK) OK = NC_CHECK(NC_PUT_VARA_FLOAT_F(NCID, VIM, [J0-JMIN, 0], [J1-J0+1, NROWS], &
                                                REAL(AIMAG(XC), C_FLOAT)), 'ecriture de x_im')
      IF (.NOT. OK) STOP 1
    END IF
#endif
    DEALLOCATE(XC)
  END DO

  IF (DO_TXT) THEN
    CLOSE(UTXT); WRITE(6,'(A,A)') ' Ecrit : ', TRIM(OUTP)//'.txt'
  END IF
  IF (DO_BIN) THEN
    CLOSE(UBIN); WRITE(6,'(A,A)') ' Ecrit : ', TRIM(OUTP)//'.bin'
  END IF
#ifdef WITH_NETCDF
  IF (DO_NC) THEN
    IF (.NOT. NC_CHECK(NC_CLOSE_F(NCID), 'fermeture')) STOP 1
    WRITE(6,'(A,A)') ' Ecrit : ', TRIM(OUTP)//'.nc'
  END IF
#endif

CONTAINS

  SUBROUTINE FAIL(MSG)
    CHARACTER(LEN=*), INTENT(IN) :: MSG
    WRITE(0,'(A,A)') ' ERREUR : ', MSG
    ERROR STOP 1
  END SUBROUTINE FAIL

  ! "I1:I2" ou "I" (une seule valeur).
  SUBROUTINE READ_RANGE(S, A, B, IOS)
    CHARACTER(LEN=*), INTENT(IN) :: S
    INTEGER, INTENT(OUT) :: A, B, IOS
    INTEGER :: P
    P = INDEX(S, ':')
    IF (P .EQ. 0) THEN
      READ(S,*,IOSTAT=IOS) A
      B = A
    ELSE
      READ(S(1:P-1),*,IOSTAT=IOS) A
      IF (IOS .EQ. 0) READ(S(P+1:),*,IOSTAT=IOS) B
    END IF
  END SUBROUTINE READ_RANGE

  ! WORD figure dans la liste LIST separee par des virgules.
  LOGICAL FUNCTION HAS_WORD(LIST, WORD)
    CHARACTER(LEN=*), INTENT(IN) :: LIST, WORD
    HAS_WORD = (INDEX(','//TRIM(LIST)//',', ','//WORD//',') .GT. 0)
  END FUNCTION HAS_WORD

  SUBROUTINE USAGE()
    WRITE(6,'(A)') 'Usage : ./export_solution --nb=NB --format=FORMATS [options]'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Relit la solution ecrite par  solveur_blr ... --save-solution  et l''ecrit dans'
    WRITE(6,'(A)') 'un seul fichier par format. Programme sequentiel (pas de mpirun).'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Obligatoires :'
    WRITE(6,'(A)') '  --nb=NB            taille de bloc NB utilisee pour le calcul. ATTENTION : NB n''est'
    WRITE(6,'(A)') '                     pas stocke dans les fichiers ; un NB faux n''est detecte que'
    WRITE(6,'(A)') '                     s''il change le nombre de lignes par processus, sinon X est'
    WRITE(6,'(A)') '                     exporte dans le desordre SANS message d''erreur.'
    WRITE(6,'(A)') '  --format=LISTE     txt, bin, netcdf, separes par des virgules, ou all'
    WRITE(6,'(A)') '                       txt    -> PREFIXE.txt : "i j Re Im" par ligne'
    WRITE(6,'(A)') '                       bin    -> PREFIXE.bin : 4 entiers 32 bits (NROWS NCOLS IMIN JMIN)'
    WRITE(6,'(A)') '                                 puis X en complexes 2x32 bits, colonne par colonne'
    WRITE(6,'(A)') '                       netcdf -> PREFIXE.nc  : variables x_re, x_im (col,row), row, col'
    WRITE(6,'(A)') 'Facultatives :'
    WRITE(6,'(A)') '  --dir=CHEMIN       dossier des fichiers saved_solution_rank*.dat (defaut .)'
    WRITE(6,'(A)') '                     = le --checkpoint-dir du calcul'
    WRITE(6,'(A)') '  --out=PREFIXE      nom des fichiers produits, sans extension (defaut solution)'
    WRITE(6,'(A)') '  --rows=I1:I2       lignes exportees (defaut toutes, 1:N)'
    WRITE(6,'(A)') '  --cols=J1:J2       seconds membres exportes (defaut tous, 1:NRHS)'
    WRITE(6,'(A)') '  --chunk=K          colonnes lues a la fois (defaut 16) : borne la memoire'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'N, NRHS et la grille PRxPC sont lus dans les fichiers.'
    WRITE(6,'(A)') 'Exemple :'
    WRITE(6,'(A)') '  ./export_solution --dir=resultats --nb=200 --format=txt,netcdf --out=X'
  END SUBROUTINE USAGE

END PROGRAM EXPORT_SOLUTION
