!===============================================================================
! MODULE A PERSONNALISER -- deux subroutines a remplacer par VOTRE matrice A
! et VOTRE second membre B (voir example_blr_mpi_pivot_user.f90, dans ce
! meme dossier, pour le programme qui les appelle). Rien d'autre dans le
! solveur n'a besoin d'etre modifie.
!
! USER_MATRIX_INIT(N_GLOBAL) : a appeler UNE FOIS avant l'assemblage, pour
! transmettre la taille globale N a USER_GET_BLOCK ci-dessous -- son
! interface est imposee par la bibliotheque (GET_BLOCK_IFACE, sans
! parametre supplementaire possible), N est donc stocke ici dans le module.
! Si votre matrice a besoin d'autres parametres (fichier a lire, constante
! physique...), ajoutez-les de la meme facon (variable de module + routine
! d'initialisation).
!
! USER_GET_BLOCK(ISTART,IEND,JSTART,JEND,BLK) : renvoie BLK=A(ISTART:IEND,
! JSTART:JEND). Appelee PLUSIEURS FOIS par la bibliotheque pour des
! sous-blocs, JAMAIS pour A entiere d'un coup -- c'est ce qui permet a
! l'assemblage de compresser au fur et a mesure, sans jamais materialiser A
! dense en RAM.
!
! USER_FILL_B(AB,GD,NRL,RSTART,BL,ML) : renvoie VOTRE second membre b pour
! les lignes GLOBALES possedees LOCALEMENT par ce processus (deja calculees
! via AB%BSTART/BEND, ne pas y toucher) et les colonnes RSTART..RSTART+NRL-1.
!
! Le placeholder ci-dessous (noyau 1/(1+|i-j|), diagonale=N, b=vecteur de 1)
! est un exemple RUNNABLE pour vous permettre de verifier que tout compile
! et tourne avant de brancher vos vraies donnees -- aucune signification
! physique.
!
! COMPRESSIBILITE : la compression BLR ne fait gagner du temps/memoire QUE
! si les blocs hors-diagonale de VOTRE matrice ont un rang numerique
! faible (noyau a decroissance geometrique, methode des elements de
! frontiere, etc.). Sinon les blocs tombent simplement en repli DENSE
! (comportement SUR, jamais de plantage), sans gain -- ajustez --ratio/TOL
! si besoin (voir example_blr_mpi_pivot_user.f90 / example_blr_mpi_pivot.f90
! pour le detail de ces parametres).
!===============================================================================
MODULE USER_MATRIX_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: MPI_BLR_MATRIX_T, GRID_T, WP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: USER_MATRIX_INIT, USER_GET_BLOCK, USER_FILL_B

  INTEGER, SAVE :: N_GLOBAL_SAVED = 0

CONTAINS

  SUBROUTINE USER_MATRIX_INIT(N_GLOBAL)
    INTEGER, INTENT(IN) :: N_GLOBAL
    N_GLOBAL_SAVED = N_GLOBAL
  END SUBROUTINE USER_MATRIX_INIT

  !=============================================================================
  ! >>> REMPLACEZ LE CORPS DE CETTE SUBROUTINE par VOTRE matrice A <<<
  !=============================================================================
  SUBROUTINE USER_GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        IF (II .EQ. JJ) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(REAL(N_GLOBAL_SAVED,WP), 0.0_WP, KIND=WP)
        ELSE
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(1.0_WP/(1.0_WP+ABS(REAL(II-JJ,WP))), 0.0_WP, KIND=WP)
        END IF
      END DO
    END DO
  END SUBROUTINE USER_GET_BLOCK

  !=============================================================================
  ! >>> REMPLACEZ LE CORPS DE CETTE SUBROUTINE par VOTRE second membre B <<<
  !=============================================================================
  SUBROUTINE USER_FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: BL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          BL(RS + (GII-AB%BSTART(II)), RR) = CMPLX(1.0_WP, 0.0_WP, KIND=WP)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE USER_FILL_B

END MODULE USER_MATRIX_MOD
