!===============================================================================
! Tests unitaires des primitives de compression (module BLR_COMPRESS_MOD),
! INDEPENDANTS de la factorisation LU -- a executer avant de faire
! confiance a BLR_LU_FACTOR (etape 2, partie 2).
!
! Cas testes :
!   1) Bloc EXACTEMENT de rang R (construit comme U*V) : la compression
!      doit trouver le bon rang et reconstruire quasi-exactement.
!   2) Bloc PLEIN RANG (aleatoire) : doit basculer en mode dense (pas de
!      compression), reconstruction EXACTE (recopie, pas d'approximation).
!   3) SCHUR_UPDATE_BLOCK, les 4 combinaisons dense/basse-rang, comparees a
!      un GEMM dense direct sur les blocs ORIGINAUX (blocs construits
!      EXACTEMENT basse-rang pour ce test, donc l'update basse-rang doit
!      egaler l'update dense a la precision machine pres, pas juste "a TOL
!      pres").
!
! Usage : ./test_blr_compress   (pas d'argument)
!===============================================================================
PROGRAM TEST_BLR_COMPRESS
  USE BLR_COMPRESS_MOD
  IMPLICIT NONE

  LOGICAL :: ALL_OK
  ALL_OK = .TRUE.

  CALL TEST_EXACT_LOW_RANK(ALL_OK)
  CALL TEST_FULL_RANK_FALLBACK(ALL_OK)
  CALL TEST_SCHUR_UPDATE_ALL_CASES(ALL_OK)

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TOUS LES TESTS BLR_COMPRESS : OK ==='
  ELSE
    WRITE(6,*) '=== AU MOINS UN TEST BLR_COMPRESS A ECHOUE ==='
    STOP 1
  END IF

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE RANDOM_COMPLEX(M, N, X)
    INTEGER, INTENT(IN) :: M, N
    COMPLEX(WP), INTENT(OUT) :: X(M,N)
    REAL(WP), ALLOCATABLE :: RP(:,:), IP(:,:)
    ALLOCATE(RP(M,N), IP(M,N))
    CALL RANDOM_NUMBER(RP)
    CALL RANDOM_NUMBER(IP)
    X = CMPLX(RP-0.5_WP, IP-0.5_WP, KIND=WP)
    DEALLOCATE(RP, IP)
  END SUBROUTINE RANDOM_COMPLEX

  !-----------------------------------------------------------------------
  SUBROUTINE TEST_EXACT_LOW_RANK(OK)
    LOGICAL, INTENT(INOUT) :: OK
    INTEGER, PARAMETER :: M = 60, N = 45, R = 5
    COMPLEX(WP) :: U(M,R), V(R,N), B(M,N), REC(M,N)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    TYPE(BLR_BLOCK_T) :: BLK
    REAL(WP) :: ERR, NORMB

    CALL RANDOM_COMPLEX(M, R, U)
    CALL RANDOM_COMPLEX(R, N, V)
    CALL CGEMM('N','N', M, N, R, CONE, U, M, V, R, CZERO, B, M)

    CALL COMPRESS_BLOCK(M, N, B, M, 1.0E-5_WP, BLK)
    CALL EXPAND_BLOCK(BLK, REC, M)

    NORMB = SQRT(SUM(ABS(B)**2))
    ERR = SQRT(SUM(ABS(REC-B)**2)) / NORMB

    WRITE(6,'(A,I0,A,L1,A,I0,A,ES10.3)') &
      ' [rang exact] taille attendue=', R, ', IS_LOWRANK=', BLK%IS_LOWRANK, &
      ', rang trouve=', BLK%RANK, ', erreur reconstruction=', ERR

    IF (.NOT.BLK%IS_LOWRANK .OR. BLK%RANK .GT. R+1 .OR. ERR .GT. 1.0E-3_WP) THEN
      WRITE(6,*) '  -> ECHEC (rang ou reconstruction incorrects)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
  END SUBROUTINE TEST_EXACT_LOW_RANK

  !-----------------------------------------------------------------------
  SUBROUTINE TEST_FULL_RANK_FALLBACK(OK)
    LOGICAL, INTENT(INOUT) :: OK
    INTEGER, PARAMETER :: M = 40, N = 40
    COMPLEX(WP) :: B(M,N), REC(M,N)
    TYPE(BLR_BLOCK_T) :: BLK
    REAL(WP) :: ERR
    INTEGER :: II

    CALL RANDOM_COMPLEX(M, N, B)
    ! diagonale dominante -> plein rang numerique, aucune raison d'etre
    ! compressible a une tolerance serree
    DO II = 1, MIN(M,N)
      B(II,II) = B(II,II) + CMPLX(5.0_WP, 0.0_WP, KIND=WP)
    END DO

    CALL COMPRESS_BLOCK(M, N, B, M, 1.0E-6_WP, BLK)
    CALL EXPAND_BLOCK(BLK, REC, M)
    ERR = SQRT(SUM(ABS(REC-B)**2)) / SQRT(SUM(ABS(B)**2))

    WRITE(6,'(A,L1,A,ES10.3)') &
      ' [plein rang]  IS_LOWRANK=', BLK%IS_LOWRANK, ', erreur reconstruction=', ERR

    IF (BLK%IS_LOWRANK .OR. ERR .GT. 1.0E-5_WP) THEN
      WRITE(6,*) '  -> ECHEC (aurait du rester dense, exact)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '  -> OK'
    END IF
  END SUBROUTINE TEST_FULL_RANK_FALLBACK

  !-----------------------------------------------------------------------
  SUBROUTINE TEST_SCHUR_UPDATE_ALL_CASES(OK)
    LOGICAL, INTENT(INOUT) :: OK
    INTEGER, PARAMETER :: MI = 30, NJ = 25, KW = 20, RL = 3, RU = 4
    COMPLEX(WP) :: UL(MI,RL), VL(RL,KW), UU(KW,RU), VU(RU,NJ)
    COMPLEX(WP) :: LDENSE(MI,KW), UDENSE(KW,NJ)
    COMPLEX(WP) :: A_REF(MI,NJ), A_TEST(MI,NJ)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    TYPE(BLR_BLOCK_T) :: LBLK_LR, UBLK_LR, LBLK_D, UBLK_D
    REAL(WP) :: ERR
    INTEGER :: CASE_ID

    ! Construire L et U EXACTEMENT basse-rang (rangs RL, RU) pour que la
    ! version compressee egale la version dense a la precision machine.
    CALL RANDOM_COMPLEX(MI, RL, UL); CALL RANDOM_COMPLEX(RL, KW, VL)
    CALL CGEMM('N','N', MI, KW, RL, CONE, UL, MI, VL, RL, CZERO, LDENSE, MI)
    CALL RANDOM_COMPLEX(KW, RU, UU); CALL RANDOM_COMPLEX(RU, NJ, VU)
    CALL CGEMM('N','N', KW, NJ, RU, CONE, UU, KW, VU, RU, CZERO, UDENSE, KW)

    ! Representations : LBLK_LR/UBLK_LR basse-rang (tolerance large pour
    ! forcer/accepter la compression), LBLK_D/UBLK_D dense forcee (TOL=0)
    CALL COMPRESS_BLOCK(MI, KW, LDENSE, MI, 1.0E-4_WP, LBLK_LR)
    CALL COMPRESS_BLOCK(KW, NJ, UDENSE, KW, 1.0E-4_WP, UBLK_LR)
    CALL COMPRESS_BLOCK(MI, KW, LDENSE, MI, 0.0_WP, LBLK_D)
    CALL COMPRESS_BLOCK(KW, NJ, UDENSE, KW, 0.0_WP, UBLK_D)

    DO CASE_ID = 1, 4
      CALL RANDOM_COMPLEX(MI, NJ, A_REF)
      A_TEST = A_REF

      ! Reference : GEMM dense direct sur les blocs ORIGINAUX
      CALL CGEMM('N','N', MI, NJ, KW, -CONE, LDENSE, MI, UDENSE, KW, &
                  CONE, A_REF, MI)

      SELECT CASE (CASE_ID)
      CASE (1)
        CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LBLK_D,  UBLK_D,  A_TEST, MI)
      CASE (2)
        CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LBLK_LR, UBLK_D,  A_TEST, MI)
      CASE (3)
        CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LBLK_D,  UBLK_LR, A_TEST, MI)
      CASE (4)
        CALL SCHUR_UPDATE_BLOCK(MI, NJ, KW, LBLK_LR, UBLK_LR, A_TEST, MI)
      END SELECT

      ERR = SQRT(SUM(ABS(A_TEST-A_REF)**2)) / SQRT(SUM(ABS(A_REF)**2))
      WRITE(6,'(A,I0,A,L1,A,L1,A,ES10.3)') ' [Schur cas ', CASE_ID, &
        '] L_lowrank=', MERGE(LBLK_LR%IS_LOWRANK,.FALSE.,CASE_ID==2.OR.CASE_ID==4), &
        ', U_lowrank=', MERGE(UBLK_LR%IS_LOWRANK,.FALSE.,CASE_ID==3.OR.CASE_ID==4), &
        ', erreur = ', ERR

      IF (ERR .GT. 1.0E-3_WP) THEN
        WRITE(6,*) '  -> ECHEC'
        OK = .FALSE.
      ELSE
        WRITE(6,*) '  -> OK'
      END IF
    END DO
  END SUBROUTINE TEST_SCHUR_UPDATE_ALL_CASES

END PROGRAM TEST_BLR_COMPRESS
