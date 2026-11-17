!===============================================================================
! Memoire reellement allouee par l'OS au processus courant (RSS, "Resident
! Set Size"), lue dans /proc/self/status (Linux uniquement). Complement aux
! compteurs analytiques de taille de structure de donnees : la RSS englobe
! aussi les tableaux de travail transitoires et l'allocateur runtime, et
! les deux mesures sont rapprochees dans les benchmarks.
!===============================================================================
MODULE BENCH_UTILS_MOD
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! RSS courant en kilooctets (-1 si /proc/self/status illisible, ex. hors
  ! Linux).
  !-----------------------------------------------------------------------
  FUNCTION GET_RSS_KB() RESULT(RSS)
    INTEGER(8) :: RSS
    INTEGER :: U, IOS
    CHARACTER(LEN=256) :: LINE

    RSS = -1_8
    OPEN(NEWUNIT=U, FILE='/proc/self/status', STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) RETURN

    DO
      READ(U, '(A)', IOSTAT=IOS) LINE
      IF (IOS .NE. 0) EXIT
      IF (LINE(1:6) .EQ. 'VmRSS:') THEN
        READ(LINE(7:), *, IOSTAT=IOS) RSS
        IF (IOS .NE. 0) RSS = -1_8
        EXIT
      END IF
    END DO
    CLOSE(U)
  END FUNCTION GET_RSS_KB

  !-----------------------------------------------------------------------
  ! Rend au systeme la memoire liberee mais retenue par malloc (glibc,
  ! malloc_trim(0)). Appelee a la fin de chaque etape de factorisation :
  ! les recompressions reallouent U/V a des rangs croissants et fragmentent
  ! le tas (mesure, LDLt N=40000 NB=500 K=5000 grille 2x2, compilation
  ! OPENMP=0 : RSS 519 -> 450 Mo par processus, temps +1 %). Sans effet
  ! quand les ALLOCATE ne passent pas par malloc -- cas d'ifx des que
  ! libiomp5 est lie (MKL threadee), d'ou la compilation OPENMP=0 (voir
  ! CLUSTER.md). Desactivable par la variable d'environnement
  ! BLR_MALLOC_TRIM=0. Linux/glibc.
  !-----------------------------------------------------------------------
  SUBROUTINE BLR_TRIM_HEAP()
    USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T
    INTERFACE
      INTEGER(C_INT) FUNCTION MALLOC_TRIM(PAD) BIND(C, NAME='malloc_trim')
        IMPORT :: C_INT, C_SIZE_T
        INTEGER(C_SIZE_T), VALUE :: PAD
      END FUNCTION MALLOC_TRIM
    END INTERFACE
    LOGICAL, SAVE :: INIT = .FALSE., ENABLED = .TRUE.
    CHARACTER(LEN=8) :: V
    INTEGER :: L, ST
    INTEGER(C_INT) :: IR
    IF (.NOT. INIT) THEN
      CALL GET_ENVIRONMENT_VARIABLE('BLR_MALLOC_TRIM', V, L, ST)
      IF (ST .EQ. 0 .AND. L .GT. 0) ENABLED = (INDEX('0nNfF', V(1:1)) .EQ. 0)
      INIT = .TRUE.
    END IF
    IF (ENABLED) IR = MALLOC_TRIM(0_C_SIZE_T)
  END SUBROUTINE BLR_TRIM_HEAP

END MODULE BENCH_UTILS_MOD
