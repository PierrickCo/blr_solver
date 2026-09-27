!===============================================================================
! Reimplementation de PSRAND/XJUMPM/JUMPIT/SETRAN, le generateur "a saut"
! interne dont PCMATGEN a besoin -- absent de MKL (routines de test
! ScaLAPACK, jamais expediees dans une bibliotheque de calcul) et de
! scalapack_ooc_ref/ (qui suppose une vraie bibliotheque ScaLAPACK/BLACS
! compilee, introuvable ici).
!
! FIDELITE : meme interface et memes constantes que pcmatgen.f leur passe
! deja (MULT0=20077, MULT1=16838, IADD0=12345, IADD1=0) -- pcmatgen.f
! tourne donc SANS MODIFICATION. L'arithmetique interne est une
! reimplementation propre en INTEGER(8) (congruence lineaire mod 2^32,
! saut-en-avant par doublement), PAS une reconstruction bit-a-bit du code
! Fortran77 original : le resultat n'est donc pas bit-identique a un vrai
! ScaLAPACK avec la meme graine, mais reste deterministe et reproductible
! independamment du decoupage (meme garantie que pcmatgen.f exploite).
!===============================================================================

! Combine un handle 2-elements (HI,LO in [0,65536)) en une valeur 32 bits
! portee par un INTEGER(8), pour l'arithmetique modulo 2^32 ci-dessous.
FUNCTION RANGEN_COMBINE(X) RESULT(V)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: X(2)
  INTEGER(8) :: V
  V = INT(X(1),8)*65536_8 + INT(X(2),8)
END FUNCTION RANGEN_COMBINE

! Operation inverse.
SUBROUTINE RANGEN_SPLIT(V, X)
  IMPLICIT NONE
  INTEGER(8), INTENT(IN) :: V
  INTEGER, INTENT(OUT) :: X(2)
  INTEGER(8) :: VM
  VM = IAND(V, 4294967295_8)
  X(1) = INT(VM/65536_8)
  X(2) = INT(MOD(VM,65536_8))
END SUBROUTINE RANGEN_SPLIT

! (A*B) mod 2^32 -- IAND est bit a bit, insensible a un depassement
! signe de INTEGER(8) pendant le produit intermediaire (comportement de
! wrap standard x86-64/ifx) : les 32 bits de poids faible restent
! corrects. Masque = 4294967295 en litteral DECIMAL explicitement
! INTEGER(8) -- PAS Z'FFFFFFFF' converti via INT(...,8), qui serait
! d'abord lu comme l'entier 32 bits signe -1 puis etendu en signe vers
! -1_8 (masque inoperant).
FUNCTION RANGEN_MULMOD(A, B) RESULT(C)
  IMPLICIT NONE
  INTEGER(8), INTENT(IN) :: A, B
  INTEGER(8) :: C
  C = IAND(A*B, 4294967295_8)
END FUNCTION RANGEN_MULMOD

! (A+B) mod 2^32 -- meme raisonnement que RANGEN_MULMOD.
FUNCTION RANGEN_ADDMOD(A, B) RESULT(C)
  IMPLICIT NONE
  INTEGER(8), INTENT(IN) :: A, B
  INTEGER(8) :: C
  C = IAND(A+B, 4294967295_8)
END FUNCTION RANGEN_ADDMOD

! Saut-en-avant par doublement : pour x_{n+1}=A0*x_n+C0 (mod 2^32),
! calcule (AJ,CJ) tel qu'appliquer (AJ,CJ) UNE fois equivaut a appliquer
! (A0,C0) NJUMP fois. Deux transformations affines T1(x)=a1*x+c1,
! T2(x)=a2*x+c2 se composent en T2(T1(x))=(a1*a2)*x+(a2*c1+c2) ; comme
! toutes les transformations combinees ici sont des puissances de la
! MEME base, l'ordre de composition n'a pas d'importance.
SUBROUTINE RANGEN_JUMPAHEAD(NJUMP, A0, C0, AJ, CJ)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: NJUMP
  INTEGER(8), INTENT(IN)  :: A0, C0
  INTEGER(8), INTENT(OUT) :: AJ, CJ
  INTEGER(8) :: A, C, RA, RC
  INTEGER(8) :: RANGEN_MULMOD, RANGEN_ADDMOD
  EXTERNAL   :: RANGEN_MULMOD, RANGEN_ADDMOD
  INTEGER :: K
  A = A0; C = C0
  RA = 1_8; RC = 0_8
  K = NJUMP
  DO WHILE (K .GT. 0)
    IF (MOD(K,2) .EQ. 1) THEN
      RC = RANGEN_ADDMOD(RANGEN_MULMOD(A,RC), C)
      RA = RANGEN_MULMOD(A,RA)
    END IF
    C = RANGEN_ADDMOD(RANGEN_MULMOD(A,C), C)
    A = RANGEN_MULMOD(A,A)
    K = K/2
  END DO
  AJ = RA; CJ = RC
END SUBROUTINE RANGEN_JUMPAHEAD

! XJUMPM(NJUMP, MULT, IADD, IRAN, IRANN, IAJUMP, ICJUMP) : etant donne le
! "pas unitaire" (MULT,IADD), calcule l'etat apres NJUMP applications a
! partir de IRAN (-> IRANN), et le couple (IAJUMP,ICJUMP) qui realise ce
! meme saut en un seul appel a JUMPIT.
SUBROUTINE XJUMPM(NJUMP, MULT, IADD, IRAN, IRANN, IAJUMP, ICJUMP)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NJUMP, MULT(2), IADD(2), IRAN(2)
  INTEGER, INTENT(OUT) :: IRANN(2), IAJUMP(2), ICJUMP(2)
  INTEGER(8) :: RANGEN_COMBINE, RANGEN_MULMOD, RANGEN_ADDMOD
  EXTERNAL   :: RANGEN_COMBINE, RANGEN_MULMOD, RANGEN_ADDMOD, RANGEN_SPLIT, RANGEN_JUMPAHEAD
  INTEGER(8) :: A0, C0, X0, AJ, CJ, XN

  ! IOR(...,1) force multiplicateur ET increment a etre IMPAIRS : les
  ! constantes de pcmatgen.f (MULT1=16838, IADD1=0) sont toutes deux
  ! PAIRES, ce qui rendrait A et C pairs et pourrait faire tomber
  ! X'=A*X+C sur un POINT FIXE (suite constante) pour certaines graines.
  ! A et C impairs excluent tout point fixe : (A-1) est alors pair, et
  ! (A-1)*X=-C (mod 2^32) n'a aucune solution quand C est impair.
  A0 = IOR(RANGEN_COMBINE(MULT), 1_8)
  C0 = IOR(RANGEN_COMBINE(IADD), 1_8)
  X0 = RANGEN_COMBINE(IRAN)

  CALL RANGEN_JUMPAHEAD(NJUMP, A0, C0, AJ, CJ)
  XN = RANGEN_ADDMOD(RANGEN_MULMOD(AJ,X0), CJ)

  CALL RANGEN_SPLIT(XN, IRANN)
  CALL RANGEN_SPLIT(AJ, IAJUMP)
  CALL RANGEN_SPLIT(CJ, ICJUMP)
END SUBROUTINE XJUMPM

! Une seule application de la recurrence (IA,IC) (potentiellement deja
! "sautee" via XJUMPM) a l'etat IRANN, produisant IRANM.
SUBROUTINE JUMPIT(IA, IC, IRANN, IRANM)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: IA(2), IC(2), IRANN(2)
  INTEGER, INTENT(OUT) :: IRANM(2)
  INTEGER(8) :: RANGEN_COMBINE, RANGEN_MULMOD, RANGEN_ADDMOD
  EXTERNAL   :: RANGEN_COMBINE, RANGEN_MULMOD, RANGEN_ADDMOD, RANGEN_SPLIT
  INTEGER(8) :: A, C, X, XM

  ! Meme correctif "impair" que XJUMPM (voir son commentaire).
  A = IOR(RANGEN_COMBINE(IA), 1_8)
  C = IOR(RANGEN_COMBINE(IC), 1_8)
  X = RANGEN_COMBINE(IRANN)
  XM = RANGEN_ADDMOD(RANGEN_MULMOD(A,X), C)
  CALL RANGEN_SPLIT(XM, IRANM)
END SUBROUTINE JUMPIT

! Etat partage entre SETRAN et PSRAND -- COMMON block Fortran77 (routines
! externes independantes, exactement comme pcmatgen.f les declare).
BLOCK DATA RANGEN_BLOCKDATA
  IMPLICIT NONE
  INTEGER :: GSEED(2), GMULT(2), GADD(2)
  COMMON /RANGEN_STATE/ GSEED, GMULT, GADD
  DATA GSEED /0, 0/, GMULT /0, 0/, GADD /0, 0/
END BLOCK DATA RANGEN_BLOCKDATA

! Fixe l'etat courant (ISEED) et memorise le "pas unitaire" (IA,IC)
! utilise par tous les appels suivants a PSRAND.
SUBROUTINE SETRAN(ISEED, IA, IC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ISEED(2), IA(2), IC(2)
  INTEGER :: GSEED(2), GMULT(2), GADD(2)
  COMMON /RANGEN_STATE/ GSEED, GMULT, GADD
  GSEED = ISEED
  GMULT = IA
  GADD  = IC
END SUBROUTINE SETRAN

! Avance l'etat courant d'un pas et renvoie un REEL uniforme dans [0,1).
! IDUMMY n'est jamais utilise par pcmatgen.f (toujours appele PSRAND(0))
! -- conserve pour matcher sa signature.
REAL FUNCTION PSRAND(IDUMMY)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IDUMMY
  INTEGER :: GSEED(2), GMULT(2), GADD(2)
  COMMON /RANGEN_STATE/ GSEED, GMULT, GADD
  INTEGER(8) :: RANGEN_COMBINE, RANGEN_MULMOD, RANGEN_ADDMOD
  EXTERNAL   :: RANGEN_COMBINE, RANGEN_MULMOD, RANGEN_ADDMOD, RANGEN_SPLIT
  INTEGER(8) :: A, C, X, XN

  A = IOR(RANGEN_COMBINE(GMULT), 1_8)
  C = IOR(RANGEN_COMBINE(GADD), 1_8)
  X = RANGEN_COMBINE(GSEED)
  XN = RANGEN_ADDMOD(RANGEN_MULMOD(A,X), C)
  CALL RANGEN_SPLIT(XN, GSEED)

  PSRAND = REAL(XN) / 4294967296.0
END FUNCTION PSRAND
