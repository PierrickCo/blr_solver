!===============================================================================
! Deux variantes d'affichage de progression, desactivees par defaut
! (SHOW_PROGRESS=.FALSE. cote appelant) -- pensees pour l'observation
! manuelle sur de petits cas test, pas pour tourner en continu sur les gros
! benchmarks :
!
! PRINT_PROGRESS_BAR : barre texte a rafraichissement en place (carriage
! return, CHAR(13), sans saut de ligne) -- FONCTIONNE seulement attachee a
! un vrai terminal interactif (utilisee par blr_matrix_mod.f90, processus
! sequentiel unique). NE PAS utiliser sous MPI/mpirun : le lanceur relaie
! le stdout de chaque processus via un pipe, pas un pseudo-terminal -- le
! retour chariot n'est alors PAS interprete comme "ecraser la ligne
! precedente" et chaque appel finit sur SA PROPRE ligne (bug constate et
! corrige : `mpirun ./example_blr_mpi_pivot` polluait la sortie d'une
! ligne par mise a jour au lieu d'une barre unique).
!
! PRINT_PROGRESS_MILESTONE : affichage par PALIERS (tous les 5% par
! defaut), une ligne complete normale (avec saut de ligne) par palier --
! robuste quel que soit le lanceur (mpirun/srun) ou le terminal, puisqu'il
! ne repose sur aucun caractere de controle. Utilisee par
! mpi_blr_matrix_mod.f90. L'appelant fournit et conserve LAST_PCT (un
! entier local, initialise a une valeur < 0 avant la boucle) : c'est ce
! qui permet le throttling SANS etat cache dans le module (chaque appelant
! garde son propre compteur, pas de risque d'interference entre plusieurs
! barres actives en meme temps dans le meme programme).
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

  SUBROUTINE PRINT_PROGRESS_MILESTONE(LABEL, CURRENT, TOTAL, LAST_PCT, STEP_PCT)
    CHARACTER(LEN=*), INTENT(IN)    :: LABEL
    INTEGER, INTENT(IN)             :: CURRENT, TOTAL
    INTEGER, INTENT(INOUT)          :: LAST_PCT
    INTEGER, INTENT(IN), OPTIONAL   :: STEP_PCT

    INTEGER :: PCT, STEP

    IF (TOTAL .LE. 0) RETURN
    STEP = 5
    IF (PRESENT(STEP_PCT)) STEP = MAX(1, STEP_PCT)
    PCT = MIN(100, MAX(0, (100*CURRENT)/TOTAL))
    IF (PCT .GE. LAST_PCT + STEP .OR. (CURRENT .GE. TOTAL .AND. PCT .GT. LAST_PCT)) THEN
      WRITE(6,'(A,A,A,I4,A)') '  ', TRIM(LABEL), ' : ', PCT, ' %'
      FLUSH(6)
      LAST_PCT = PCT
    END IF
  END SUBROUTINE PRINT_PROGRESS_MILESTONE

END MODULE PROGRESS_BAR_MOD
