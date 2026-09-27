!===============================================================================
! Deux mecanismes d'affichage de progression, desactives par defaut cote
! appelant (SHOW_PROGRESS=.FALSE.) -- destines a l'observation manuelle sur
! des cas de taille reduite, pas a tourner en continu sur des benchmarks
! longs.
!
! PRINT_PROGRESS_BAR : barre texte sur une seule ligne qui se remplit
! progressivement. N'emet aucun caractere de controle (ni retour chariot,
! ni saut de ligne avant la fin) : chaque appel ajoute uniquement les
! caracteres '=' devenus dus depuis le dernier appel (ADVANCE='NO'), si
! bien que la ligne grandit comme une sortie texte ordinaire. Rendu
! identique sur un terminal interactif, une sortie redirigee vers un
! fichier, ou un lanceur MPI relayant le flux par un pipe.
!
! Etat porte par l'appelant (NFILL_STATE, INTENT(INOUT)) : un entier local
! initialise a une valeur negative avant la boucle -- aucun etat n'est
! conserve dans le module, donc plusieurs barres actives simultanement
! dans le meme programme n'interferent jamais entre elles. Le premier
! appel (etat negatif) ouvre la ligne (etiquette puis crochet) ; le
! dernier appel (CURRENT>=TOTAL) referme le crochet et passe a la ligne
! suivante.
!
! PRINT_PROGRESS_MILESTONE : affichage par paliers (tous les 5% par
! defaut), une ligne complete par palier -- ne repose sur aucun caractere
! de controle, donc robuste quel que soit le lanceur ou le terminal.
! Utilisee par mpi_blr_matrix_mod.f90. Meme convention d'etat porte par
! l'appelant que ci-dessus (LAST_PCT, INTENT(INOUT), initialise a une
! valeur negative avant la boucle).
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
