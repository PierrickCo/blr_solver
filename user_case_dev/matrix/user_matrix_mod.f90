!===============================================================================
! MODULE A PERSONNALISER -- USER_GET_BLOCK (matrice A) et USER_FILL_B
! (second membre B) sont les deux subroutines a remplacer par vos
! propres donnees ; rien d'autre dans le solveur n'a besoin de changer.
!
! Version actuelle : A generee par un vrai appel a PCMATGEN (scalapack/
! pcmatgen.f, copie non modifiee de scalapack_ooc_ref/), B=A*X ou X est
! une seconde matrice PCMATGEN (graine differente) -- meme recette que
! scalapack_ooc_ref/testcdriver.f ("Generate rhs B <- A*X"). PCMATGEN
! s'appuie sur PSRAND/XJUMPM/JUMPIT/SETRAN, reimplementes dans
! scalapack/scalapack_randgen.f90 (absents de MKL et de
! scalapack_ooc_ref/) -- voir son en-tete pour la limite de fidelite.
!
! COMPRESSIBILITE : matrice PLEINE (diagonale dominante juste pour rester
! resoluble) -- la compression BLR n'y trouve rien a compresser, tous les
! blocs tombent en repli dense (sans plantage, sans gain). Utile pour
! tester la robustesse du solveur sur une matrice non structuree, pas
! pour demontrer l'avantage BLR (qui demande un noyau a rang faible).
!===============================================================================
MODULE USER_MATRIX_MOD
  USE MPI_BLR_MATRIX_MOD, ONLY: MPI_BLR_MATRIX_T, GRID_T, WP
  USE MKL_INTERFACES_MOD, ONLY: CGEMM
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: USER_MATRIX_INIT, USER_GET_BLOCK, USER_FILL_B

  INTEGER, SAVE :: N_GLOBAL_SAVED = 0
  INTEGER, SAVE :: NB_GLOBAL_SAVED = 0
  INTEGER, SAVE :: NRHS_GLOBAL_SAVED = 0
  ! Graines PCMATGEN -- fixes, arbitraires (n'importe quelle valeur donne
  ! une matrice/solution differente mais tout aussi valide/reproductible).
  ! DEUX graines DISTINCTES (meme convention que testcdriver.f, ISEED=13
  ! pour A / 999999 pour X) : sans ca, X serait correle a A (memes
  ! valeurs), ce qui n'a aucune raison physique d'etre le cas pour une
  ! vraie solution.
  INTEGER, PARAMETER :: PCMATGEN_SEED_A = 1
  INTEGER, PARAMETER :: PCMATGEN_SEED_X = 999999

CONTAINS

  SUBROUTINE USER_MATRIX_INIT(N_GLOBAL, NB_GLOBAL, NRHS_GLOBAL)
    INTEGER, INTENT(IN) :: N_GLOBAL, NB_GLOBAL, NRHS_GLOBAL
    N_GLOBAL_SAVED = N_GLOBAL
    NB_GLOBAL_SAVED = NB_GLOBAL
    NRHS_GLOBAL_SAVED = NRHS_GLOBAL
  END SUBROUTINE USER_MATRIX_INIT

  !=============================================================================
  ! >>> REMPLACEZ LE CORPS DE CETTE SUBROUTINE par VOTRE matrice A <<<
  ! Grille virtuelle 1x1 : ce processus se presente a PCMATGEN comme s'il
  ! possedait la matrice entiere, et lui demande de n'en generer que la
  ! sous-plage voulue via IROFF/IRNUM/ICOFF/ICNUM. PCMATGEN exige ces
  ! offsets multiples de NB -- toujours vrai pour l'assemblage BLR, mais
  ! pas pour un appelant generique (ACCUMULATE_RESIDUAL) qui decoupe
  ! independamment de NB : on boucle donc sur les blocs NB alignes qui
  ! intersectent la plage demandee, et on ne copie que l'intersection.
  !=============================================================================
  SUBROUTINE USER_GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    EXTERNAL :: PCMATGEN
    COMPLEX(WP), ALLOCATABLE :: TMPBLK(:,:)
    INTEGER :: RB0, CB0, RALIGNED, CALIGNED, RCHUNKEND, CCHUNKEND
    INTEGER :: RLO, RHI, CLO, CHI

    RB0 = (ISTART-1)/NB_GLOBAL_SAVED
    DO WHILE (RB0*NB_GLOBAL_SAVED .LT. IEND)
      RALIGNED = RB0*NB_GLOBAL_SAVED + 1
      RCHUNKEND = MIN(RALIGNED+NB_GLOBAL_SAVED-1, N_GLOBAL_SAVED)
      RLO = MAX(RALIGNED, ISTART)
      RHI = MIN(RCHUNKEND, IEND)

      CB0 = (JSTART-1)/NB_GLOBAL_SAVED
      DO WHILE (CB0*NB_GLOBAL_SAVED .LT. JEND)
        CALIGNED = CB0*NB_GLOBAL_SAVED + 1
        CCHUNKEND = MIN(CALIGNED+NB_GLOBAL_SAVED-1, N_GLOBAL_SAVED)
        CLO = MAX(CALIGNED, JSTART)
        CHI = MIN(CCHUNKEND, JEND)

        ALLOCATE(TMPBLK(RCHUNKEND-RALIGNED+1, CCHUNKEND-CALIGNED+1))
        CALL PCMATGEN(0, 'N', 'D', N_GLOBAL_SAVED, N_GLOBAL_SAVED, NB_GLOBAL_SAVED, NB_GLOBAL_SAVED, &
                      TMPBLK, RCHUNKEND-RALIGNED+1, 0, 0, PCMATGEN_SEED_A, &
                      RALIGNED-1, RCHUNKEND-RALIGNED+1, CALIGNED-1, CCHUNKEND-CALIGNED+1, &
                      0, 0, 1, 1)
        BLK(RLO-ISTART+1:RHI-ISTART+1, CLO-JSTART+1:CHI-JSTART+1) = &
             TMPBLK(RLO-RALIGNED+1:RHI-RALIGNED+1, CLO-CALIGNED+1:CHI-CALIGNED+1)
        DEALLOCATE(TMPBLK)

        CB0 = CB0 + 1
      END DO
      RB0 = RB0 + 1
    END DO
  END SUBROUTINE USER_GET_BLOCK

  ! Genere XBLK=X(ISTART:IEND,JSTART:JEND), la matrice "solution" utilisee
  ! par USER_FILL_B pour calculer B=A*X -- meme technique que
  ! USER_GET_BLOCK, avec une graine differente (decorreler X de A) et un
  ! blocking colonne de 1 (NRHS est petit, jamais aligne sur NB : un
  ! blocking de 1 rend n'importe quel ICOFF valide sans boucle colonne).
  SUBROUTINE GEN_X_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    EXTERNAL :: PCMATGEN
    COMPLEX(WP), ALLOCATABLE :: TMPBLK(:,:)
    INTEGER :: RB0, RALIGNED, RCHUNKEND, RLO, RHI

    RB0 = (ISTART-1)/NB_GLOBAL_SAVED
    DO WHILE (RB0*NB_GLOBAL_SAVED .LT. IEND)
      RALIGNED = RB0*NB_GLOBAL_SAVED + 1
      RCHUNKEND = MIN(RALIGNED+NB_GLOBAL_SAVED-1, N_GLOBAL_SAVED)
      RLO = MAX(RALIGNED, ISTART)
      RHI = MIN(RCHUNKEND, IEND)

      ALLOCATE(TMPBLK(RCHUNKEND-RALIGNED+1, JEND-JSTART+1))
      CALL PCMATGEN(0, 'N', 'N', N_GLOBAL_SAVED, NRHS_GLOBAL_SAVED, NB_GLOBAL_SAVED, 1, &
                    TMPBLK, RCHUNKEND-RALIGNED+1, 0, 0, PCMATGEN_SEED_X, &
                    RALIGNED-1, RCHUNKEND-RALIGNED+1, JSTART-1, JEND-JSTART+1, &
                    0, 0, 1, 1)
      BLK(RLO-ISTART+1:RHI-ISTART+1, :) = TMPBLK(RLO-RALIGNED+1:RHI-RALIGNED+1, :)
      DEALLOCATE(TMPBLK)

      RB0 = RB0 + 1
    END DO
  END SUBROUTINE GEN_X_BLOCK

  !=============================================================================
  ! >>> REMPLACEZ LE CORPS DE CETTE SUBROUTINE par VOTRE second membre B <<<
  ! B = A*X, calcule par bloc de lignes BLR (jamais A ni X entieres en
  ! RAM) : materialise A(bloc,1:N) et X(1:N,RSTART:RSTART+NRL-1) puis
  ! CGEMM. BTMP (tampon contigu) est necessaire car BL(RS:RS+BSZ-1,:) est
  ! une tranche non contigue (BL a pour dimension ML, differente de BSZ)
  ! -- la passer directement a CGEMM corromprait silencieusement la
  ! memoire.
  !=============================================================================
  SUBROUTINE USER_FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: BL(MAX(ML,1),MAX(NRL,1))
    COMPLEX(WP), ALLOCATABLE :: AROWBLK(:,:), XBLK(:,:), BTMP(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: lli, II, BSZ, RS

    IF (NRL .LE. 0) RETURN

    ALLOCATE(XBLK(N_GLOBAL_SAVED, NRL))
    CALL GEN_X_BLOCK(1, N_GLOBAL_SAVED, RSTART, RSTART+NRL-1, XBLK)

    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      BSZ = AB%BEND(II) - AB%BSTART(II) + 1
      ALLOCATE(AROWBLK(BSZ, N_GLOBAL_SAVED), BTMP(BSZ, NRL))
      CALL USER_GET_BLOCK(AB%BSTART(II), AB%BEND(II), 1, N_GLOBAL_SAVED, AROWBLK)
      CALL CGEMM('N','N', BSZ, NRL, N_GLOBAL_SAVED, CONE, AROWBLK, BSZ, XBLK, N_GLOBAL_SAVED, CZERO, BTMP, BSZ)
      BL(RS:RS+BSZ-1, 1:NRL) = BTMP
      DEALLOCATE(AROWBLK, BTMP)
      RS = RS + BSZ
    END DO
    DEALLOCATE(XBLK)
  END SUBROUTINE USER_FILL_B

END MODULE USER_MATRIX_MOD
