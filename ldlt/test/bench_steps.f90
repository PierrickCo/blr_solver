!===============================================================================
! Benchmark par etapes d'UN solveur (LU ou LDL^T BLR+MPI) sur la matrice
! complexe symetrique de sym_testmat_mod.f90 : temps, RSS et taux de
! compression a chaque etape. Un seul solveur par execution, pour que la RSS
! mesuree ne soit pas polluee par l'autre.
!
! Etapes : (1) assemblage + compression, (2) factorisation, puis pour
! chacun des NBATCH paquets de NRHS seconds membres : construction du second
! membre (hors solveur, pour information) et resolution. Un seul paquet de
! seconds membres en memoire a la fois.
! RSS par etape (MAX par processus et SOMME sur la grille) : RSS en fin
! d'etape et PIC pendant l'etape (VmHWM, remis a zero avant chaque etape par
! /proc/self/clear_refs ; repli sur le pic depuis le debut si le noyau ne le
! permet pas -- signale dans la sortie).
!
! Second membre : NCHECK colonnes de controle PAR PAQUET (reparties)
! valent exactement A*Xexact (calcule par tuiles NB x NB, sans tableau
! pleine ligne), les autres sont pseudo-aleatoires -- le temps de resolution
! ne depend pas des valeurs. L'erreur ||X-Xexact||/||Xexact|| est mesuree
! sur les colonnes de controle.
!
! Taux de compression : stockage BLR (blocs diagonaux + hors-diagonale) /
! stockage dense du meme solveur (LU : N^2 ; LDLt : N(N+1)/2, triangle
! inferieur), apres assemblage et apres factorisation ; plus la part de
! blocs hors-diagonale compresses et leur rang moyen.
!
! Usage : mpirun -np PR*PC ./bench_steps SOLVER PR PC N NB NRHS NBATCH TOL K [NCHECK] [STOP_AFTER]
!   SOLVER : LU ou LDLT ; NRHS : seconds membres par paquet ; NBATCH :
!   nombre de paquets ; K : nombre d'onde du noyau de Helmholtz ;
!   NCHECK (defaut 4, par paquet) ; STOP_AFTER=1 : s'arrete apres l'assemblage
!   (calibrage de la compressibilite).
! Options : RSVD, distribution decalee (DIAG_SPREAD) ; LDLt : anticipation
! et resolution automatique (DESCENT_MODE=0, distribuee si NRHS grand).
!===============================================================================
PROGRAM BENCH_STEPS
  USE MPI_BLR_LDLT_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: MPI_BLR_MATRIX_T, MPI_ASSEMBLE_BLR_MATRIX, MPI_BLR_MATRIX_LU_FACTOR_PIVOT, &
                                MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D, MPI_BLR_MATRIX_MEMORY_BYTES
  USE SYM_TESTMAT_MOD
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  TYPE(MPI_BLR_MATRIX_T) :: AL
  TYPE(MPI_BLR_SYM_T) :: AS
  INTEGER :: PR, PC, N, NB, NRHS, NBATCH, IB, ROFF, NCHECK, STOP_AFTER, IERR, MYID0, MLOC, NRHS_LOC, RSTART_MINE
  INTEGER :: NC, NO, INFO, NSW, li, I, NBLKR_LOC
  INTEGER(8) :: RS
  REAL(WP) :: TOL, KWAVE
  LOGICAL :: IS_LU, HWM_RESET_OK
  CHARACTER(LEN=32) :: ARG, SOLVER
  COMPLEX(WP), ALLOCATABLE :: BL(:,:)
  INTEGER, ALLOCATABLE :: BSTART_(:), BEND_(:)
  ! etapes : 1 assemblage, 2 factorisation, puis par paquet IB :
  ! 1+2*IB second membre, 2+2*IB resolution
  DOUBLE PRECISION, ALLOCATABLE :: T(:), RSS_END(:,:), RSS_PEAK(:,:)
  REAL(WP), ALLOCATABLE :: ERRB(:)
  CHARACTER(LEN=26), ALLOCATABLE :: STEP_NAMES(:)
  DOUBLE PRECISION :: T0, RSS_BASE(2)
  INTEGER :: NSTEPS_TOT
  ! compression : (1) apres assemblage, (2) apres factorisation
  INTEGER :: NCMP(2), NOFF(2)
  INTEGER(8) :: RSUM(2), MEMB(2), MEMB_MAX(2)
  REAL(WP) :: ERR

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)
  IF (COMMAND_ARGUMENT_COUNT() .LT. 9) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./bench_steps LU|LDLT PR PC N NB NRHS NBATCH TOL K [NCHECK] [STOP_AFTER]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, SOLVER)
  IS_LU = (TRIM(SOLVER) .EQ. 'LU' .OR. TRIM(SOLVER) .EQ. 'lu')
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) PC
  CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) NRHS
  CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) NBATCH
  CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) TOL
  CALL GET_COMMAND_ARGUMENT(9, ARG); READ(ARG,*) KWAVE
  NCHECK = 4; STOP_AFTER = 0
  IF (COMMAND_ARGUMENT_COUNT() .GE. 10) THEN; CALL GET_COMMAND_ARGUMENT(10, ARG); READ(ARG,*) NCHECK; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 11) THEN; CALL GET_COMMAND_ARGUMENT(11, ARG); READ(ARG,*) STOP_AFTER; END IF
  NCHECK = MAX(1, MIN(NCHECK, NRHS))
  NBATCH = MAX(1, NBATCH)
  NSTEPS_TOT = 2 + 2*NBATCH
  ALLOCATE(T(NSTEPS_TOT), RSS_END(NSTEPS_TOT,2), RSS_PEAK(NSTEPS_TOT,2), STEP_NAMES(NSTEPS_TOT), ERRB(NBATCH))
  STEP_NAMES(1) = 'assemblage + compression'
  STEP_NAMES(2) = 'factorisation'
  DO IB = 1, NBATCH
    WRITE(STEP_NAMES(1+2*IB),'(A,I0,A)') 'paquet ', IB, ' : second membre'
    WRITE(STEP_NAMES(2+2*IB),'(A,I0,A)') 'paquet ', IB, ' : resolution'
  END DO

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL SYM_TESTMAT_SETUP(N, NB, WAVENUMBER=KWAVE)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES8.1,A,F7.1,A,I0)') ' Solveur ', TRIM(SOLVER), ' , grille ', PR, 'x', PC, &
         ' , N=', N, ' , NB=', NB, ' , ', NBATCH, ' paquets de NRHS=', NRHS, ' , TOL=', TOL, ' , K=', KWAVE, &
         ' , NCHECK/paquet=', NCHECK
    IF (IS_LU) THEN
      WRITE(6,'(A)') ' Options : RSVD, DIAG_SPREAD'
    ELSE
      WRITE(6,'(A)') ' Options : RSVD, DIAG_SPREAD ; LDLt : anticipation, resolution automatique'
    END IF
  END IF

  T = 0.0D0; RSS_END = 0.0D0; RSS_PEAK = 0.0D0; ERRB = 0.0_WP
  MEMB = 0_8; MEMB_MAX = 0_8; RSUM = 0_8; NCMP = 0; NOFF = 0
  CALL RSS_STATS(GET_RSS_KB(), RSS_BASE)
  HWM_RESET_OK = RESET_PEAK()

  ! ---- 1) assemblage + compression
  CALL MPI_BARRIER(GRID%COMM, IERR)
  T0 = MPI_WTIME()
  IF (IS_LU) THEN
    CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, SYM_GET_BLOCK, AL, T(1), NCMP(1), NOFF(1), RSUM(1), &
                                 SHOW_PROGRESS=.TRUE., USE_RSVD=.TRUE., DIAG_SPREAD=.TRUE.)
    MEMB(1) = MPI_BLR_MATRIX_MEMORY_BYTES(AL)
  ELSE
    CALL MPI_ASSEMBLE_BLR_SYM(GRID, N, NB, TOL, SYM_GET_BLOCK, AS, T(1), NCMP(1), NOFF(1), RSUM(1), &
                              SHOW_PROGRESS=.TRUE., USE_RSVD=.TRUE., DIAG_SPREAD=.TRUE.)
    MEMB(1) = MPI_BLR_SYM_MEMORY_BYTES(AS)
  END IF
  CALL END_STEP(1, T0)
  CALL REDUCE_BYTES(MEMB(1), MEMB_MAX(1))
  IF (STOP_AFTER .EQ. 1) THEN
    CALL REPORT(1)
    CALL MPI_FINALIZE(IERR)
    STOP
  END IF

  ! ---- 2) factorisation
  HWM_RESET_OK = RESET_PEAK() .AND. HWM_RESET_OK
  CALL MPI_BARRIER(GRID%COMM, IERR)
  T0 = MPI_WTIME()
  IF (IS_LU) THEN
    CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, AL, TOL, INFO, NCMP(2), NOFF(2), RSUM(2), NSW, SHOW_PROGRESS=.TRUE.)
    MEMB(2) = MPI_BLR_MATRIX_MEMORY_BYTES(AL)
  ELSE
    CALL MPI_BLR_SYM_LDLT_FACTOR(GRID, AS, TOL, INFO, NCMP(2), NOFF(2), RSUM(2), SHOW_PROGRESS=.TRUE.)
    MEMB(2) = MPI_BLR_SYM_MEMORY_BYTES(AS)
  END IF
  CALL END_STEP(2, T0)
  CALL REDUCE_BYTES(MEMB(2), MEMB_MAX(2))
  IF (INFO .NE. 0 .AND. MYID0.EQ.0) WRITE(6,*) 'ERREUR factorisation, INFO=', INFO

  ! ---- 3) second membre (hors solveur)
  IF (IS_LU) THEN
    NBLKR_LOC = AL%NBLKR_LOC
    ALLOCATE(BSTART_(AL%NBLK), BEND_(AL%NBLK)); BSTART_ = AL%BSTART; BEND_ = AL%BEND
  ELSE
    NBLKR_LOC = AS%NBLKR_LOC
    ALLOCATE(BSTART_(AS%NBLK), BEND_(AS%NBLK)); BSTART_ = AS%BSTART; BEND_ = AS%BEND
  END IF
  MLOC = 0
  DO li = 1, NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + BEND_(I) - BSTART_(I) + 1
  END DO
  ALLOCATE(BL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  DO IB = 1, NBATCH
    ROFF = (IB-1)*NRHS                    ! decalage global des colonnes du paquet
    HWM_RESET_OK = RESET_PEAK() .AND. HWM_RESET_OK
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = MPI_WTIME()
    CALL FILL_B()
    CALL END_STEP(1+2*IB, T0)
    HWM_RESET_OK = RESET_PEAK() .AND. HWM_RESET_OK
    CALL MPI_BARRIER(GRID%COMM, IERR)
    T0 = MPI_WTIME()
    IF (IS_LU) THEN
      CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, AL, NRHS, NRHS_LOC, BL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE.)
    ELSE
      CALL MPI_BLR_SYM_LDLT_SOLVE_2D(GRID, AS, NRHS, NRHS_LOC, BL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE.)
    END IF
    CALL END_STEP(2+2*IB, T0)
    ERRB(IB) = CHECK_ERROR()
  END DO
  ERR = MAXVAL(ERRB)

  CALL REPORT(NSTEPS_TOT)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! (MAX par processus, SOMME) en Mo d'une valeur locale en Ko.
  SUBROUTINE RSS_STATS(KB, OUT)
    INTEGER(8), INTENT(IN)        :: KB
    DOUBLE PRECISION, INTENT(OUT) :: OUT(2)
    DOUBLE PRECISION :: V
    INTEGER :: IE
    V = DBLE(KB)/1024.0D0
    CALL MPI_ALLREDUCE(V, OUT(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IE)
    CALL MPI_ALLREDUCE(V, OUT(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, GRID%COMM, IE)
  END SUBROUTINE RSS_STATS

  ! Temps (MAX sur la grille), RSS de fin et pic de l'etape K.
  SUBROUTINE END_STEP(K, TSTART)
    INTEGER, INTENT(IN)          :: K
    DOUBLE PRECISION, INTENT(IN) :: TSTART
    DOUBLE PRECISION :: TL
    INTEGER :: IE
    TL = MPI_WTIME() - TSTART
    CALL MPI_ALLREDUCE(TL, T(K), 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IE)
    CALL RSS_STATS(GET_RSS_KB(), RSS_END(K,:))
    CALL RSS_STATS(GET_HWM_KB(), RSS_PEAK(K,:))
  END SUBROUTINE END_STEP

  SUBROUTINE REDUCE_BYTES(B, BMAX)
    INTEGER(8), INTENT(INOUT) :: B
    INTEGER(8), INTENT(OUT)   :: BMAX
    INTEGER(8) :: BS
    INTEGER :: IE
    CALL MPI_ALLREDUCE(B, BMAX, 1, MPI_INTEGER8, MPI_MAX, GRID%COMM, IE)
    CALL MPI_ALLREDUCE(B, BS, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IE)
    B = BS
  END SUBROUTINE REDUCE_BYTES

  ! Remet a zero le pic de RSS (VmHWM) du processus : ecrire "5" dans
  ! /proc/self/clear_refs (Linux >= 4.0). .FALSE. si impossible.
  LOGICAL FUNCTION RESET_PEAK()
    INTEGER :: U, IOS
    RESET_PEAK = .FALSE.
    OPEN(NEWUNIT=U, FILE='/proc/self/clear_refs', STATUS='OLD', ACTION='WRITE', IOSTAT=IOS)
    IF (IOS .NE. 0) RETURN
    WRITE(U, '(A)', IOSTAT=IOS) '5'
    CLOSE(U)
    RESET_PEAK = (IOS .EQ. 0)
  END FUNCTION RESET_PEAK

  ! Pic de RSS (VmHWM) en Ko, -1 si illisible.
  INTEGER(8) FUNCTION GET_HWM_KB()
    INTEGER :: U, IOS
    CHARACTER(LEN=256) :: LINE
    GET_HWM_KB = -1_8
    OPEN(NEWUNIT=U, FILE='/proc/self/status', STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) RETURN
    DO
      READ(U, '(A)', IOSTAT=IOS) LINE
      IF (IOS .NE. 0) EXIT
      IF (LINE(1:6) .EQ. 'VmHWM:') THEN
        READ(LINE(7:), *, IOSTAT=IOS) GET_HWM_KB
        EXIT
      END IF
    END DO
    CLOSE(U)
  END FUNCTION GET_HWM_KB

  ! Colonne globale R est-elle une colonne de controle ?
  LOGICAL FUNCTION IS_CHECK(R)
    INTEGER, INTENT(IN) :: R
    INTEGER :: STRIDE
    STRIDE = MAX(1, NRHS/NCHECK)
    IS_CHECK = (MOD(R-1, STRIDE) .EQ. 0 .AND. (R-1)/STRIDE .LT. NCHECK)
  END FUNCTION IS_CHECK

  ! Valeur pseudo-aleatoire deterministe des colonnes hors controle.
  COMPLEX(WP) FUNCTION FILLER(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    REAL(WP) :: S1, S2
    S1 = SIN(REAL(GI,WP)*12.9898_WP + REAL(R,WP)*78.233_WP) * 43758.5453_WP
    S2 = SIN(REAL(GI,WP)*39.3468_WP + REAL(R,WP)*11.135_WP) * 24634.6345_WP
    FILLER = CMPLX(S1 - REAL(FLOOR(S1),WP) - 0.5_WP, S2 - REAL(FLOOR(S2),WP) - 0.5_WP, KIND=WP)
  END FUNCTION FILLER

  ! BL : colonnes de controle = A*Xexact (tuiles NB x NB de MES lignes),
  ! autres = FILLER.
  SUBROUTINE FILL_B()
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: TILE(:,:), XC(:,:), BC(:,:)
    INTEGER, ALLOCATABLE :: CLOC(:)
    INTEGER :: NCL, C, R, GI, RSB, J, BSZ, NJ, NBLK_
    NCL = 0
    ALLOCATE(CLOC(MAX(NRHS_LOC,1)))
    DO C = 1, NRHS_LOC
      R = RSTART_MINE + C - 1                 ! indice dans le paquet
      IF (IS_CHECK(R)) THEN
        NCL = NCL + 1; CLOC(NCL) = C
      ELSE
        RSB = 1
        DO li = 1, NBLKR_LOC
          I = (li-1)*GRID%PR + GRID%MYROW + 1
          DO GI = BSTART_(I), BEND_(I)
            BL(RSB + GI - BSTART_(I), C) = FILLER(GI, ROFF + R)
          END DO
          RSB = RSB + BEND_(I) - BSTART_(I) + 1
        END DO
      END IF
    END DO
    IF (NCL .EQ. 0) THEN
      DEALLOCATE(CLOC)
      RETURN
    END IF
    NBLK_ = SIZE(BSTART_)
    ALLOCATE(XC(N,NCL))
    DO C = 1, NCL
      DO GI = 1, N
        XC(GI,C) = SYM_XEX(GI, ROFF + RSTART_MINE + CLOC(C) - 1)
      END DO
    END DO
    RSB = 1
    DO li = 1, NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      BSZ = BEND_(I) - BSTART_(I) + 1
      ALLOCATE(BC(BSZ,NCL))
      BC = (0.0_WP,0.0_WP)
      DO J = 1, NBLK_
        NJ = BEND_(J) - BSTART_(J) + 1
        ALLOCATE(TILE(BSZ,NJ))
        CALL SYM_GET_BLOCK(BSTART_(I), BEND_(I), BSTART_(J), BEND_(J), TILE)
        CALL CGEMM('N', 'N', BSZ, NCL, NJ, CONE, TILE, BSZ, XC(BSTART_(J),1), N, CONE, BC, BSZ)
        DEALLOCATE(TILE)
      END DO
      DO C = 1, NCL
        BL(RSB:RSB+BSZ-1, CLOC(C)) = BC(:,C)
      END DO
      DEALLOCATE(BC)
      RSB = RSB + BSZ
    END DO
    DEALLOCATE(XC, CLOC)
  END SUBROUTINE FILL_B

  ! ||X-Xexact|| / ||Xexact|| sur les colonnes de controle (sommes locales
  ! reduites, aucun rassemblement de X).
  REAL(WP) FUNCTION CHECK_ERROR()
    DOUBLE PRECISION :: LOC(2), GLB(2)
    INTEGER :: C, R, GI, RSB, IE
    LOC = 0.0D0
    DO C = 1, NRHS_LOC
      R = RSTART_MINE + C - 1
      IF (.NOT. IS_CHECK(R)) CYCLE
      RSB = 1
      DO li = 1, NBLKR_LOC
        I = (li-1)*GRID%PR + GRID%MYROW + 1
        DO GI = BSTART_(I), BEND_(I)
          LOC(1) = LOC(1) + ABS(BL(RSB + GI - BSTART_(I), C) - SYM_XEX(GI, ROFF + R))**2
          LOC(2) = LOC(2) + ABS(SYM_XEX(GI, ROFF + R))**2
        END DO
        RSB = RSB + BEND_(I) - BSTART_(I) + 1
      END DO
    END DO
    CALL MPI_ALLREDUCE(LOC, GLB, 2, MPI_DOUBLE_PRECISION, MPI_SUM, GRID%COMM, IE)
    CHECK_ERROR = REAL(SQRT(GLB(1)/MAX(GLB(2),1.0D-300)), WP)
  END FUNCTION CHECK_ERROR

  ! Tableau final ; NSTEP = derniere etape faite. Ligne STEPS_RESULT
  ! (cle=valeur) pour l'agregation par script.
  SUBROUTINE REPORT(NSTEP)
    INTEGER, INTENT(IN) :: NSTEP
    DOUBLE PRECISION :: DENSE_B, TAUX(2), RMOY(2)
    INTEGER :: K, S
    IF (MYID0 .NE. 0) RETURN
    IF (IS_LU) THEN
      DENSE_B = 8.0D0*DBLE(N)*DBLE(N)
    ELSE
      DENSE_B = 8.0D0*DBLE(N)*DBLE(N+1)/2.0D0
    END IF
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') ' ============================================================================================='
    WRITE(6,'(A,A)') '  Solveur ', TRIM(SOLVER)
    WRITE(6,'(A,F9.1,A,F9.1,A)') '  RSS de depart : max/processus ', RSS_BASE(1), ' Mo , total ', RSS_BASE(2), ' Mo'
    IF (.NOT. HWM_RESET_OK) WRITE(6,'(A)') '  (pic de RSS : remise a zero impossible -> pic DEPUIS LE DEBUT)'
    WRITE(6,'(A)') '  Etape                         temps (s)   RSS fin max/tot (Mo)     pic max/tot (Mo)     erreur'
    DO K = 1, NSTEP
      IF (K .GE. 4 .AND. MOD(K,2) .EQ. 0) THEN
        WRITE(6,'(2X,A26,F14.3,F12.1,F10.1,F12.1,F10.1,ES11.3)') STEP_NAMES(K), T(K), RSS_END(K,1), RSS_END(K,2), &
             RSS_PEAK(K,1), RSS_PEAK(K,2), ERRB((K-2)/2)
      ELSE
        WRITE(6,'(2X,A26,F14.3,F12.1,F10.1,F12.1,F10.1)') STEP_NAMES(K), T(K), RSS_END(K,1), RSS_END(K,2), &
             RSS_PEAK(K,1), RSS_PEAK(K,2)
      END IF
    END DO
    IF (NSTEP .GE. 4) WRITE(6,'(2X,A26,F14.3)') 'total solveur hors 2nd m. ', T(1)+T(2)+SUM(T(4:NSTEP:2))
    WRITE(6,'(A)') '  Compression                 stockage tot (Mo)  max/proc (Mo)   taux (% dense)  blocs compr.  rang moyen'
    DO S = 1, MIN(2, MAX(1, NSTEP-0))
      IF (NSTEP .EQ. 1 .AND. S .EQ. 2) EXIT
      TAUX(S) = 100.0D0*DBLE(MEMB(S))/DENSE_B
      RMOY(S) = DBLE(RSUM(S))/DBLE(MAX(NOFF(S),1))
      WRITE(6,'(2X,A24,F16.1,F15.1,F16.2,F11.1,A,F12.1)') MERGE('apres assemblage        ', 'apres factorisation     ', S.EQ.1), &
           DBLE(MEMB(S))/1048576.0D0, DBLE(MEMB_MAX(S))/1048576.0D0, TAUX(S), &
           100.0D0*DBLE(NCMP(S))/DBLE(MAX(NOFF(S),1)), ' %', RMOY(S)
    END DO
    IF (NSTEP .GE. 4) WRITE(6,'(A,ES11.4,A,I0,A)') '  ||X-Xexact||/||Xexact|| max sur les paquets = ', ERR, '  (', NCHECK, &
         ' colonnes de controle par paquet)'
    WRITE(6,'(A)') ' ============================================================================================='
    WRITE(6,'(A,A,A,F7.1,99(A,G0.4))') 'STEPS_RESULT solver=', TRIM(SOLVER), ' k=', KWAVE, &
         ' t_asm=', T(1), ' t_fac=', T(2), ' t_sol_tot=', SUM(T(4:MAX(NSTEP,4):2)), &
         ' rss_asm=', RSS_END(1,1), ' rss_fac=', RSS_END(2,1), ' rss_sol=', RSS_END(NSTEP,1), &
         ' pic_asm=', RSS_PEAK(1,1), ' pic_fac=', RSS_PEAK(2,1), ' pic_sol=', MAXVAL(RSS_PEAK(4:MAX(NSTEP,4):2,1)), &
         ' taux_asm=', 100.0D0*DBLE(MEMB(1))/DENSE_B, ' taux_fac=', 100.0D0*DBLE(MEMB(2))/DENSE_B, &
         ' rang_asm=', DBLE(RSUM(1))/DBLE(MAX(NOFF(1),1)), ' err=', ERR
  END SUBROUTINE REPORT

END PROGRAM BENCH_STEPS
