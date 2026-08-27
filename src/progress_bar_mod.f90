!===============================================================================
! Petite barre de progression texte (carriage-return, une seule ligne mise a
! jour en place). Utilisee par blr_matrix_mod.f90 pour suivre l'avancement
! de chacune des 3 etapes couteuses (assemblage+compression, factorisation,
! resolution) -- UNE barre distincte par etape, chacune appelee avec son
! propre LABEL. Desactivee par defaut (SHOW_PROGRESS=.FALSE. cote appelant) :
! pensee pour l'observation manuelle sur de petits cas test, pas pour tourner
! en continu sur les gros benchmarks (des dizaines de milliers de rafraichis-
! sements de ligne n'apportent rien une fois le calcul lance en tache de
! fond).
!===============================================================================
MODULE PROGRESS_BAR_MOD
  IMPLICIT NONE
CONTAINS

  SUBROUTINE PRINT_PROGRESS_BAR(LABEL, CURRENT, TOTAL)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(IN)          :: CURRENT, TOTAL

    INTEGER, PARAMETER :: BARWIDTH = 30
    INTEGER :: NFILL, PCT, I
    CHARACTER(LEN=BARWIDTH) :: BAR
    CHARACTER(LEN=256) :: LINE

    IF (TOTAL .LE. 0) RETURN
    PCT   = MIN(100, MAX(0, NINT(100.0*REAL(CURRENT)/REAL(TOTAL))))
    NFILL = MIN(BARWIDTH, MAX(0, (BARWIDTH*CURRENT)/TOTAL))
    DO I = 1, BARWIDTH
      IF (I .LE. NFILL) THEN
        BAR(I:I) = '='
      ELSE
        BAR(I:I) = ' '
      END IF
    END DO

    WRITE(LINE,'(A,A,A,A,A,I4,A)') '  ', TRIM(LABEL), ' [', BAR, '] ', PCT, ' %'
    WRITE(6,'(A,A)',ADVANCE='NO') CHAR(13), TRIM(LINE)
    FLUSH(6)
    IF (CURRENT .GE. TOTAL) WRITE(6,*)
  END SUBROUTINE PRINT_PROGRESS_BAR

END MODULE PROGRESS_BAR_MOD
