!===============================================================================
! Generation de la matrice A et du second membre B par PCMATGEN
! (scalapack/pcmatgen.f, generateur de test ScaLAPACK), a la demande et par
! morceaux : jamais de matrice N x N ni N x NRHS complete en memoire.
!
!   A : PCMATGEN, diagonale dominante ('D'), graine SEED_A.
!       - solveur LU   : matrice generale        (AFORM = 'N')
!       - solveur LDLt : matrice SYMETRIQUE A^T=A (AFORM = 'S')
!   B : PCMATGEN, matrice generale N x NRHS, graine SEED_B -- generee
!       DIRECTEMENT (pas de B = A*X) : seul le residu ||A*X-B||/||B||
!       permet de verifier la solution (--check-residual).
!
! Chaque processus se presente a PCMATGEN comme une grille 1x1 (contexte 0,
! MYROW=MYCOL=0, NPROW=NPCOL=1) et ne demande que la sous-plage voulue via
! IROFF/IRNUM/ICOFF/ICNUM. Le flux aleatoire de PCMATGEN etant indexe par
! la position GLOBALE de l'element, le resultat ne depend ni de la grille
! MPI ni du decoupage. PCMATGEN exige IROFF multiple de MB et ICOFF multiple
! de NB : on boucle donc sur les morceaux alignes et on ne garde que
! l'intersection demandee.
!===============================================================================
MODULE MATRIX_GEN_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GEN_INIT, GEN_GET_A_BLOCK, GEN_FILL_B

  INTEGER, SAVE :: N_G = 0, NB_G = 0, NRHS_G = 0
  CHARACTER(LEN=1), SAVE :: AFORM_G = 'N'
  ! Graines fixes (meme convention que scalapack_ooc_ref/testcdriver.f :
  ! deux graines distinctes pour A et pour le second membre).
  INTEGER, PARAMETER :: SEED_A = 1
  INTEGER, PARAMETER :: SEED_B = 999999
  EXTERNAL :: PCMATGEN

CONTAINS

  ! SYMMETRIC = .TRUE. pour le LDLt (A symetrique), .FALSE. pour le LU.
  SUBROUTINE GEN_INIT(N, NB, NRHS, SYMMETRIC)
    INTEGER, INTENT(IN) :: N, NB, NRHS
    LOGICAL, INTENT(IN) :: SYMMETRIC
    N_G = N; NB_G = NB; NRHS_G = NRHS
    AFORM_G = MERGE('S', 'N', SYMMETRIC)
  END SUBROUTINE GEN_INIT

  ! BLK = A(ISTART:IEND, JSTART:JEND) -- interface GET_BLOCK_IFACE, appelee
  ! par l'assemblage BLR et par le calcul du residu.
  SUBROUTINE GEN_GET_A_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    COMPLEX(WP), ALLOCATABLE :: TMP(:,:)
    INTEGER :: RB, CB, RA, CA, RE, CE, RLO, RHI, CLO, CHI

    RB = (ISTART-1)/NB_G
    DO WHILE (RB*NB_G .LT. IEND)
      RA = RB*NB_G + 1
      RE = MIN(RA+NB_G-1, N_G)
      RLO = MAX(RA, ISTART); RHI = MIN(RE, IEND)
      CB = (JSTART-1)/NB_G
      DO WHILE (CB*NB_G .LT. JEND)
        CA = CB*NB_G + 1
        CE = MIN(CA+NB_G-1, N_G)
        CLO = MAX(CA, JSTART); CHI = MIN(CE, JEND)
        ALLOCATE(TMP(RE-RA+1, CE-CA+1))
        CALL PCMATGEN(0, AFORM_G, 'D', N_G, N_G, NB_G, NB_G, TMP, RE-RA+1, 0, 0, SEED_A, &
                      RA-1, RE-RA+1, CA-1, CE-CA+1, 0, 0, 1, 1)
        BLK(RLO-ISTART+1:RHI-ISTART+1, CLO-JSTART+1:CHI-JSTART+1) = &
             TMP(RLO-RA+1:RHI-RA+1, CLO-CA+1:CHI-CA+1)
        DEALLOCATE(TMP)
        CB = CB + 1
      END DO
      RB = RB + 1
    END DO
  END SUBROUTINE GEN_GET_A_BLOCK

  ! Remplit BL(1:ML, 1:NRL) = B(lignes locales, GSTART:GSTART+NRL-1).
  ! Les lignes locales sont les blocs-lignes I = (li-1)*PR + MYROW + 1,
  ! li = 1..NBLKR_LOC, concatenes dans l'ordre (meme convention que le
  ! solveur). Les debuts de blocs BSTART(I) = (I-1)*NB+1 sont alignes sur
  ! NB ; en colonne, un blocking de 1 rend tout GSTART valide.
  SUBROUTINE GEN_FILL_B(BSTART, BEND, NBLKR_LOC, PR, MYROW, NRL, GSTART, BL, ML)
    INTEGER, INTENT(IN)      :: BSTART(:), BEND(:), NBLKR_LOC, PR, MYROW, NRL, GSTART, ML
    COMPLEX(WP), INTENT(OUT) :: BL(MAX(ML,1), MAX(NRL,1))
    INTEGER :: li, I, MI, RS

    IF (NRL .LE. 0 .OR. ML .LE. 0) RETURN
    RS = 1
    DO li = 1, NBLKR_LOC
      I = (li-1)*PR + MYROW + 1
      MI = BEND(I) - BSTART(I) + 1
      CALL PCMATGEN(0, 'N', 'N', N_G, NRHS_G, NB_G, 1, BL(RS,1), MAX(ML,1), 0, 0, SEED_B, &
                    BSTART(I)-1, MI, GSTART-1, NRL, 0, 0, 1, 1)
      RS = RS + MI
    END DO
  END SUBROUTINE GEN_FILL_B

END MODULE MATRIX_GEN_MOD
