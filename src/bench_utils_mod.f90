!===============================================================================
! Utilitaire de benchmark : memoire REELLE (RSS, "Resident Set Size") du
! processus courant, lue directement dans /proc/self/status (specifique
! Linux/WSL2, l'environnement de developpement de ce projet). Complement
! aux compteurs ANALYTIQUES deja existants (BLR_FACTORED_LU_MEMORY_BYTES,
! proposition 2) : RSS mesure la memoire REELLEMENT allouee par l'OS pour
! CE processus a l'instant de l'appel (tableaux de travail, allocateur
! runtime inclus), pas seulement la taille theorique d'une structure de
! donnees precise -- les deux sont rapportes cote a cote dans les
! benchmarks pour se recouper.
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
