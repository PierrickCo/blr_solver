!===============================================================================
! Matrice de test COMPLEXE SYMETRIQUE (A^T = A, non hermitienne), partagee
! par le test, l'exemple et le benchmark du LDL^T BLR+MPI.
!
!   A(i,j) = STRONG                   si i et j sont dans le MEME bloc NB et
!                                     j = miroir de i dans ce bloc
!                                     (anti-diagonale LOCALE, symetrique)
!          = DIAG_WEAK                si i = j (et pas sur l'anti-diagonale)
!          = H(r_ij) + EPS*bruit(i,j) sinon
!
! H(r) = exp(i*K*r)/(1+r) : noyau de Helmholtz regularise entre N points
! equirepartis sur le cercle unite (r = corde = 2|sin(pi(i-j)/N)|) --
! indices proches = champ proche, indices lointains = champ lointain, donc
! blocs hors-diagonale genuinement compressibles (meme noyau que
! test/bench_blr_mpi_pivot_ooc_geometric.f90). Bruit symetrique (hash de
! (min,max)) pour degrader la compressibilite a volonte.
!
! Diagonale FAIBLE + anti-diagonale FORTE confinee au bloc : Bunch-Kaufman
! est force de pivoter, avec des pivots 2x2 -- mais l'element fort reste
! toujours visible du pivot CONFINE au bloc diagonal (meme principe que la
! matrice du LU pivote, voir examples/example_blr_mpi_pivot.f90).
!===============================================================================
MODULE SYM_TESTMAT_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE
  INTEGER, SAVE  :: SM_N = 0, SM_NB = 1
  REAL(WP), SAVE :: SM_K = 8.0_WP, SM_EPS = 0.0_WP, SM_STRONG = 50.0_WP
  COMPLEX(WP), SAVE :: SM_DIAG = (0.01_WP, 0.005_WP)
  REAL(WP), PARAMETER, PRIVATE :: PI = 3.14159265358979_WP

CONTAINS

  ! STRONG=0 supprime l'anti-diagonale ; DIAG (reel) fixe la diagonale :
  ! grand (ex. 3) -> pivots 1x1 sans echange, intermediaire (ex. 0.6) ->
  ! melange 1x1/2x2, faible (defaut 0.01, avec STRONG) -> pivots 2x2.
  SUBROUTINE SYM_TESTMAT_SETUP(N, NB, WAVENUMBER, EPS, STRONG, DIAG)
    INTEGER, INTENT(IN)            :: N, NB
    REAL(WP), INTENT(IN), OPTIONAL :: WAVENUMBER, EPS, STRONG, DIAG
    SM_N = N; SM_NB = NB
    IF (PRESENT(WAVENUMBER)) SM_K = WAVENUMBER
    IF (PRESENT(EPS)) SM_EPS = EPS
    IF (PRESENT(STRONG)) SM_STRONG = STRONG
    IF (PRESENT(DIAG)) SM_DIAG = CMPLX(DIAG, 0.5_WP*DIAG, KIND=WP)
  END SUBROUTINE SYM_TESTMAT_SETUP

  REAL(WP) FUNCTION SYM_HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(MIN(GI,GJ),WP)*12.9898_WP + REAL(MAX(GI,GJ),WP)*78.233_WP) * 43758.5453_WP
    SYM_HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION SYM_HASH_NOISE

  COMPLEX(WP) FUNCTION SYM_A(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    INTEGER :: BSTART0, BEND0
    REAL(WP) :: R
    BSTART0 = ((GI-1)/SM_NB)*SM_NB + 1
    BEND0 = MIN(BSTART0+SM_NB-1, SM_N)
    IF (GJ.GE.BSTART0 .AND. GJ.LE.BEND0 .AND. GJ.EQ.BEND0-(GI-BSTART0) .AND. SM_STRONG.NE.0.0_WP) THEN
      SYM_A = CMPLX(SM_STRONG, 0.0_WP, KIND=WP)
    ELSE IF (GI .EQ. GJ) THEN
      SYM_A = SM_DIAG
    ELSE
      R = 2.0_WP*ABS(SIN(PI*REAL(GI-GJ,WP)/REAL(SM_N,WP)))
      SYM_A = EXP(CMPLX(0.0_WP, SM_K*R, KIND=WP)) / CMPLX(1.0_WP+R, 0.0_WP, KIND=WP)
      IF (SM_EPS .NE. 0.0_WP) SYM_A = SYM_A + SM_EPS*SYM_HASH_NOISE(GI,GJ)
    END IF
  END FUNCTION SYM_A

  ! Interface GET_BLOCK_IFACE (blr_matrix_mod.f90).
  SUBROUTINE SYM_GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        BLK(II-ISTART+1, JJ-JSTART+1) = SYM_A(II, JJ)
      END DO
    END DO
  END SUBROUTINE SYM_GET_BLOCK

  ! Solution exacte deterministe (verifications).
  COMPLEX(WP) FUNCTION SYM_XEX(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    SYM_XEX = CMPLX(REAL(GI,WP)*0.01_WP, -0.3_WP*REAL(R,WP)+0.01_WP*REAL(GI,WP), KIND=WP)
  END FUNCTION SYM_XEX

END MODULE SYM_TESTMAT_MOD
