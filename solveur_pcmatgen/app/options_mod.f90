!===============================================================================
! Options de la ligne de commande du solveur (voir PRINT_USAGE, ou --help).
!
!   mpirun -np PR*PC ./solveur_blr PR PC [N] [NB] [NRHS] [TOL] [--cle=valeur ...]
!
! PR PC obligatoires ; N NB NRHS TOL positionnels optionnels dans cet ordre ;
! toutes les autres options sont nommees (--cle=valeur, ou --cle seul pour
! --cle=1), dans n'importe quel ordre.
!===============================================================================
MODULE OPTIONS_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: OPTIONS_T, PARSE_OPTIONS, PRINT_USAGE, SET_ENV

  TYPE :: OPTIONS_T
    ! Solveur : .FALSE. = LU avec pivot (A generale), .TRUE. = LDLt (A symetrique)
    LOGICAL  :: LDLT = .FALSE.
    INTEGER  :: PR = 0, PC = 0
    INTEGER  :: N = 4000, NB = 200, NRHS = 4
    REAL(WP) :: TOL = 1.0E-4_WP
    ! Compression
    LOGICAL  :: USE_RSVD = .TRUE.        ! SVD randomisee adaptative (sinon CGESVD complete)
    LOGICAL  :: ORTHO_RECOMP = .TRUE.    ! recompression incrementale sur bases orthonormees
    ! Factorisation / resolution
    LOGICAL  :: DIAG_SPREAD = .TRUE.     ! distribution decalee des blocs-colonnes
    LOGICAL  :: LOOKAHEAD = .TRUE.       ! anticipation d'une etape
    LOGICAL  :: PANEL_PACK = .TRUE.      ! diffusion groupee des panneaux
    INTEGER  :: DESCENT_MODE = 0         ! 0 auto, 1 colonne, 2 distribuee
    INTEGER  :: NRHS_BATCH = 0           ! 0 = tous les seconds membres en une fois
    ! Verification / diagnostic
    LOGICAL  :: CHECK_RESIDUAL = .FALSE.
    LOGICAL  :: PROFILE_USAGE = .FALSE.
    LOGICAL  :: TRACE_FACT = .FALSE.
    LOGICAL  :: DUMP_WORK = .FALSE.
    LOGICAL  :: DUMP_STRUCTURE = .FALSE.
    LOGICAL  :: MALLOC_TRIM = .TRUE.
    ! Sauvegarde / reprise
    LOGICAL  :: SAVE_MATRIX = .FALSE., SAVE_RHS = .FALSE., SAVE_FACTORED = .FALSE., SAVE_SOLUTION = .FALSE.
    LOGICAL  :: LOAD_MATRIX = .FALSE., LOAD_RHS = .FALSE., LOAD_FACTORED = .FALSE.
    LOGICAL  :: RELOAD_SOLUTION = .FALSE.
    CHARACTER(LEN=256) :: CHECKPOINT_DIR = '.'
  END TYPE OPTIONS_T

CONTAINS

  ! STAT : 0 = continuer, 1 = arret normal (--help), 2 = erreur.
  ! Messages affiches par le rang 0 seulement.
  SUBROUTINE PARSE_OPTIONS(MYID, OPT, STAT)
    INTEGER, INTENT(IN)          :: MYID
    TYPE(OPTIONS_T), INTENT(OUT) :: OPT
    INTEGER, INTENT(OUT)         :: STAT
    CHARACTER(LEN=256) :: ARG, VAL
    CHARACTER(LEN=64)  :: KEY
    INTEGER :: IARG, NARG, POS, IOS

    STAT = 0
    NARG = COMMAND_ARGUMENT_COUNT()
    DO IARG = 1, NARG
      CALL GET_COMMAND_ARGUMENT(IARG, ARG)
      IF (TRIM(ARG).EQ.'--help' .OR. TRIM(ARG).EQ.'-h') THEN
        IF (MYID.EQ.0) CALL PRINT_USAGE()
        STAT = 1
        RETURN
      END IF
    END DO
    IF (NARG .LT. 2) THEN
      IF (MYID.EQ.0) CALL PRINT_USAGE()
      STAT = 2
      RETURN
    END IF
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*,IOSTAT=IOS) OPT%PR
    IF (IOS .NE. 0) GOTO 900
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*,IOSTAT=IOS) OPT%PC
    IF (IOS .NE. 0) GOTO 900

    ! Positionnels optionnels N NB NRHS TOL, jusqu'au premier "--".
    POS = 0
    IARG = 3
    DO WHILE (IARG .LE. NARG)
      CALL GET_COMMAND_ARGUMENT(IARG, ARG)
      IF (ARG(1:2) .EQ. '--') EXIT
      POS = POS + 1
      SELECT CASE (POS)
      CASE (1); READ(ARG,*,IOSTAT=IOS) OPT%N
      CASE (2); READ(ARG,*,IOSTAT=IOS) OPT%NB
      CASE (3); READ(ARG,*,IOSTAT=IOS) OPT%NRHS
      CASE (4); READ(ARG,*,IOSTAT=IOS) OPT%TOL
      CASE DEFAULT
        IF (MYID.EQ.0) WRITE(6,'(A)') ' ERREUR : trop d''arguments positionnels (apres PR PC N NB NRHS TOL, '// &
             'utiliser --cle=valeur ; voir --help).'
        STAT = 2
        RETURN
      END SELECT
      IF (IOS .NE. 0) GOTO 900
      IARG = IARG + 1
    END DO

    ! Options nommees.
    DO WHILE (IARG .LE. NARG)
      CALL GET_COMMAND_ARGUMENT(IARG, ARG)
      IF (ARG(1:2) .NE. '--') THEN
        IF (MYID.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : argument inattendu "', TRIM(ARG), &
             '" -- les options doivent commencer par -- (voir --help).'
        STAT = 2
        RETURN
      END IF
      CALL SPLIT_KV(ARG, KEY, VAL)
      IOS = 0
      SELECT CASE (TRIM(KEY))
      CASE ('solver')
        SELECT CASE (TRIM(VAL))
        CASE ('lu', 'LU');     OPT%LDLT = .FALSE.
        CASE ('ldlt', 'LDLT', 'ldlt^t'); OPT%LDLT = .TRUE.
        CASE DEFAULT
          IF (MYID.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : --solver="', TRIM(VAL), '" inconnu (lu ou ldlt).'
          STAT = 2
          RETURN
        END SELECT
      CASE ('n');               READ(VAL,*,IOSTAT=IOS) OPT%N
      CASE ('nb');              READ(VAL,*,IOSTAT=IOS) OPT%NB
      CASE ('nrhs');            READ(VAL,*,IOSTAT=IOS) OPT%NRHS
      CASE ('tol');             READ(VAL,*,IOSTAT=IOS) OPT%TOL
      CASE ('use-rsvd');        CALL READ_FLAG(VAL, OPT%USE_RSVD, IOS)
      CASE ('ortho-recomp');    CALL READ_FLAG(VAL, OPT%ORTHO_RECOMP, IOS)
      CASE ('diag-spread');     CALL READ_FLAG(VAL, OPT%DIAG_SPREAD, IOS)
      CASE ('lookahead');       CALL READ_FLAG(VAL, OPT%LOOKAHEAD, IOS)
      CASE ('panel-pack');      CALL READ_FLAG(VAL, OPT%PANEL_PACK, IOS)
      CASE ('descent-mode')
        READ(VAL,*,IOSTAT=IOS) OPT%DESCENT_MODE
        IF (IOS.EQ.0 .AND. (OPT%DESCENT_MODE.LT.0 .OR. OPT%DESCENT_MODE.GT.2)) IOS = 1
      CASE ('nrhs-batch');      READ(VAL,*,IOSTAT=IOS) OPT%NRHS_BATCH
      CASE ('check-residual');  CALL READ_FLAG(VAL, OPT%CHECK_RESIDUAL, IOS)
      CASE ('profile-usage');   CALL READ_FLAG(VAL, OPT%PROFILE_USAGE, IOS)
      CASE ('trace-fact');      CALL READ_FLAG(VAL, OPT%TRACE_FACT, IOS)
      CASE ('dump-work');       CALL READ_FLAG(VAL, OPT%DUMP_WORK, IOS)
      CASE ('dump-structure');  CALL READ_FLAG(VAL, OPT%DUMP_STRUCTURE, IOS)
      CASE ('malloc-trim');     CALL READ_FLAG(VAL, OPT%MALLOC_TRIM, IOS)
      CASE ('save-matrix');     CALL READ_FLAG(VAL, OPT%SAVE_MATRIX, IOS)
      CASE ('save-rhs');        CALL READ_FLAG(VAL, OPT%SAVE_RHS, IOS)
      CASE ('save-factored');   CALL READ_FLAG(VAL, OPT%SAVE_FACTORED, IOS)
      CASE ('save-solution');   CALL READ_FLAG(VAL, OPT%SAVE_SOLUTION, IOS)
      CASE ('load-matrix');     CALL READ_FLAG(VAL, OPT%LOAD_MATRIX, IOS)
      CASE ('load-rhs');        CALL READ_FLAG(VAL, OPT%LOAD_RHS, IOS)
      CASE ('load-factored');   CALL READ_FLAG(VAL, OPT%LOAD_FACTORED, IOS)
      CASE ('reload-solution'); CALL READ_FLAG(VAL, OPT%RELOAD_SOLUTION, IOS)
      CASE ('checkpoint-dir');  OPT%CHECKPOINT_DIR = VAL
      CASE DEFAULT
        IF (MYID.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : option inconnue "--', TRIM(KEY), '" (voir --help).'
        STAT = 2
        RETURN
      END SELECT
      IF (IOS .NE. 0) THEN
        IF (MYID.EQ.0) WRITE(6,'(A,A,A,A,A)') ' ERREUR : valeur invalide "', TRIM(VAL), &
             '" pour --', TRIM(KEY), ' (voir --help).'
        STAT = 2
        RETURN
      END IF
      IARG = IARG + 1
    END DO

    ! --reload-solution n'a rien a relire sans sauvegarde de la solution.
    IF (OPT%RELOAD_SOLUTION) OPT%SAVE_SOLUTION = .TRUE.
    IF (OPT%N.LE.0 .OR. OPT%NB.LE.0 .OR. OPT%NRHS.LE.0 .OR. OPT%PR.LE.0 .OR. OPT%PC.LE.0) THEN
      IF (MYID.EQ.0) WRITE(6,'(A)') ' ERREUR : PR, PC, N, NB et NRHS doivent etre strictement positifs.'
      STAT = 2
    END IF
    RETURN

900 IF (MYID.EQ.0) WRITE(6,'(A,A,A)') ' ERREUR : argument positionnel invalide "', TRIM(ARG), '" (voir --help).'
    STAT = 2
  END SUBROUTINE PARSE_OPTIONS

  SUBROUTINE SPLIT_KV(ARGSTR, KEYOUT, VALOUT)
    CHARACTER(LEN=*), INTENT(IN)  :: ARGSTR
    CHARACTER(LEN=*), INTENT(OUT) :: KEYOUT, VALOUT
    INTEGER :: EQPOS
    EQPOS = INDEX(ARGSTR, '=')
    IF (EQPOS .EQ. 0) THEN
      KEYOUT = ARGSTR(3:)
      VALOUT = '1'
    ELSE
      KEYOUT = ARGSTR(3:EQPOS-1)
      VALOUT = ARGSTR(EQPOS+1:)
    END IF
  END SUBROUTINE SPLIT_KV

  ! 0|1 (ou tout entier : non nul = vrai).
  SUBROUTINE READ_FLAG(VAL, FLAG, IOS)
    CHARACTER(LEN=*), INTENT(IN) :: VAL
    LOGICAL, INTENT(INOUT)       :: FLAG
    INTEGER, INTENT(OUT)         :: IOS
    INTEGER :: IV
    READ(VAL,*,IOSTAT=IOS) IV
    IF (IOS .EQ. 0) FLAG = (IV .NE. 0)
  END SUBROUTINE READ_FLAG

  ! Positionne une variable d'environnement du processus (setenv de la libc) :
  ! certaines options de la bibliotheque (BLR_DUMP_WORK, BLR_MALLOC_TRIM) ne
  ! sont lues que dans l'environnement.
  SUBROUTINE SET_ENV(NAME, VALUE)
    USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR
    CHARACTER(LEN=*), INTENT(IN) :: NAME, VALUE
    INTERFACE
      INTEGER(C_INT) FUNCTION C_SETENV(N, V, OVERWRITE) BIND(C, NAME='setenv')
        IMPORT :: C_INT, C_CHAR
        CHARACTER(KIND=C_CHAR), INTENT(IN) :: N(*), V(*)
        INTEGER(C_INT), VALUE :: OVERWRITE
      END FUNCTION C_SETENV
    END INTERFACE
    INTEGER(C_INT) :: IR
    IR = C_SETENV(TRIM(NAME)//C_NULL_CHAR, TRIM(VALUE)//C_NULL_CHAR, 1_C_INT)
  END SUBROUTINE SET_ENV

  SUBROUTINE PRINT_USAGE()
    WRITE(6,'(A)') 'Usage : mpirun -np PR*PC ./solveur_blr PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  PR PC          grille de processus MPI (obligatoire, PR*PC = nombre de processus)'
    WRITE(6,'(A)') '  N NB NRHS TOL  optionnels, dans cet ordre (defauts 4000 200 4 1e-4)'
    WRITE(6,'(A)') '                 N : taille de A, NB : taille des blocs, NRHS : seconds membres,'
    WRITE(6,'(A)') '                 TOL : tolerance de compression bas-rang'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'A et B sont generees par PCMATGEN (generateur ScaLAPACK), graines fixes.'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Options (--cle=valeur, ou --cle seul = --cle=1 ; dans n''importe quel ordre) :'
    WRITE(6,'(A)') ' Solveur'
    WRITE(6,'(A)') '  --solver=lu|ldlt        (defaut lu) LU avec pivot (A generale) ou'
    WRITE(6,'(A)') '                          LDLt Bunch-Kaufman (A symetrique, PCMATGEN ''S'')'
    WRITE(6,'(A)') '  --n= --nb= --nrhs= --tol=   equivalents des positionnels'
    WRITE(6,'(A)') ' Compression'
    WRITE(6,'(A)') '  --use-rsvd=0|1          (defaut 1) SVD randomisee adaptative ; 0 = CGESVD complete'
    WRITE(6,'(A)') '  --ortho-recomp=0|1      (defaut 1) recompression incrementale du Schur (bases orthonormees)'
    WRITE(6,'(A)') ' Factorisation et resolution'
    WRITE(6,'(A)') '  --diag-spread=0|1       (defaut 1) distribution decalee des blocs-colonnes'
    WRITE(6,'(A)') '  --lookahead=0|1         (defaut 1) anticipation d''une etape de factorisation'
    WRITE(6,'(A)') '  --panel-pack=0|1        (defaut 1) diffusion groupee des panneaux'
    WRITE(6,'(A)') '  --descent-mode=0|1|2    (defaut 0) resolution : 0 auto, 1 colonne, 2 distribuee'
    WRITE(6,'(A)') '  --nrhs-batch=ENTIER     (defaut 0 = tout) seconds membres traites par lots'
    WRITE(6,'(A)') ' Verification et diagnostic'
    WRITE(6,'(A)') '  --check-residual=0|1    (defaut 0) affiche ||A*X-B||/||B||'
    WRITE(6,'(A)') '  --profile-usage=0|1     (defaut 0) RSS/CPU par etape (usage_profile_rank<N>.csv)'
    WRITE(6,'(A)') '  --trace-fact=0|1        (defaut 0) trace et compteurs de la factorisation'
    WRITE(6,'(A)') '  --dump-work=0|1         (defaut 0) releve du travail de Schur (blr_work_rank<N>.bin)'
    WRITE(6,'(A)') '  --dump-structure=0|1    (defaut 0) structure BLR (blr_structure_<etape>_rank<N>.csv)'
    WRITE(6,'(A)') '  --malloc-trim=0|1       (defaut 1) rend la memoire liberee au systeme a chaque etape'
    WRITE(6,'(A)') ' Sauvegarde et reprise (un fichier par processus dans --checkpoint-dir)'
    WRITE(6,'(A)') '  --save-matrix --save-factored --save-rhs --save-solution     (=0|1, defaut 0)'
    WRITE(6,'(A)') '  --load-matrix --load-factored --load-rhs                     (=0|1, defaut 0)'
    WRITE(6,'(A)') '                          --load-factored saute assemblage et factorisation'
    WRITE(6,'(A)') '  --checkpoint-dir=CHEMIN (defaut .) repertoire existant des sauvegardes'
    WRITE(6,'(A)') '  --reload-solution=0|1   (defaut 0) en fin de run, libere A puis relit X (rang 0) ;'
    WRITE(6,'(A)') '                          implique --save-solution'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Exemple :'
    WRITE(6,'(A)') '  mpirun -np 4 ./solveur_blr 2 2 4000 200 8 1e-4 --solver=ldlt --check-residual'
  END SUBROUTINE PRINT_USAGE

END MODULE OPTIONS_MOD
