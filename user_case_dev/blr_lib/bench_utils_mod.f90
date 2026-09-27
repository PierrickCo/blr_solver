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

END MODULE BENCH_UTILS_MOD
