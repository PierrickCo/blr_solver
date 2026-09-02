!===============================================================================
! Deux variantes d'affichage de progression, desactivees par defaut
! (SHOW_PROGRESS=.FALSE. cote appelant) -- pensees pour l'observation
! manuelle sur de petits cas test, pas pour tourner en continu sur les gros
! benchmarks :
!
! PRINT_PROGRESS_BAR : barre texte sur UNE SEULE LIGNE qui se REMPLIT au
! fur et a mesure (demande explicite : plus de reecriture complete de la
! ligne a chaque appel). AUCUN caractere de controle (ni CHAR(13), ni
! saut de ligne avant la fin) : chaque appel se contente d'AJOUTER les
! nouveaux caracteres '=' devenus dus depuis le dernier appel, avec
! ADVANCE='NO' -- la ligne grandit normalement, comme n'importe quelle
! sortie texte. Robuste PARTOUT (terminal interactif, sortie redirigee
! vers un fichier, mpirun/srun relayant le stdout par pipe) puisqu'aucun
! rendu special n'est requis : c'est ce qui a motive l'abandon du
! rafraichissement par CHAR(13) (qui ne s'affichait correctement que sur
! un vrai terminal, cf. l'historique dans la memoire du projet).
!
! ETAT PORTE PAR L'APPELANT (NFILL_STATE, INTENT(INOUT)) : un entier
! local, initialise a une valeur < 0 AVANT la boucle -- meme convention
! que LAST_PCT dans PRINT_PROGRESS_MILESTONE ci-dessous (aucun etat cache
! dans le module, donc aucun risque d'interference entre plusieurs barres
! actives en meme temps dans le meme programme). Le premier appel (etat
! < 0) ouvre la ligne (etiquette + crochet) ; le dernier appel
! (CURRENT>=TOTAL) referme le crochet et saute une ligne.
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

  SUBROUTINE PRINT_PROGRESS_BAR(LABEL, CURRENT, TOTAL, NFILL_STATE)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(IN)          :: CURRENT, TOTAL
    INTEGER, INTENT(INOUT)       :: NFILL_STATE

    INTEGER, PARAMETER :: BARWIDTH = 30
    INTEGER :: NFILL, PCT, I

    IF (TOTAL .LE. 0) RETURN
    NFILL = MIN(BARWIDTH, MAX(0, (BARWIDTH*CURRENT)/TOTAL))

    IF (NFILL_STATE .LT. 0) THEN
      WRITE(6,'(A,A,A)',ADVANCE='NO') '  ', TRIM(LABEL), ' ['
      NFILL_STATE = 0
    END IF

    DO I = NFILL_STATE+1, NFILL
      WRITE(6,'(A)',ADVANCE='NO') '='
    END DO
    NFILL_STATE = NFILL
    FLUSH(6)

    IF (CURRENT .GE. TOTAL) THEN
      PCT = MIN(100, MAX(0, NINT(100.0*REAL(CURRENT)/REAL(TOTAL))))
      WRITE(6,'(A,I4,A)') '] ', PCT, ' %'
      FLUSH(6)
    END IF
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
