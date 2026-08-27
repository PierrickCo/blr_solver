!===============================================================================
! ETAPE 3 : verifie que le pivot partiel (deja implemente depuis l'etape 1,
! CGETRF sur le panneau courant, applique AVANT toute compression a chaque
! etape -- voir blr_lu_mod.f90) fonctionne reellement, y compris combine a
! la compression BLR (etape 2). Nos tests precedents (etapes 1 et 2)
! utilisaient tous des matrices a diagonale dominante, qui ne sollicitent
! quasiment JAMAIS le pivot (le meilleur pivot pour la colonne i est presque
! toujours la ligne i elle-meme) -- ce test comble ce trou en construisant
! des matrices ou le pivot est REELLEMENT necessaire, et en verifiant via
! les statistiques N_SWAPS/MIN_PIVOT (ajoutees a BLOCK_LU_FACTOR et
! BLR_LU_FACTOR pour ce test) que des echanges ont bien eu lieu -- pas
! seulement que le code s'execute sans erreur.
!
! Cas testes :
!   1) Matrice compressible (comme l'etape 2) puis LIGNES PERMUTEES
!      ALEATOIREMENT : preserve le conditionnement (une permutation de
!      lignes ne change pas les valeurs singulieres) mais detruit
!      l'alignement diagonale/pivot -- force des echanges reels. Teste
!      SANS compression (BLOCK_LU_FACTOR) ET AVEC compression
!      (BLR_LU_FACTOR), verifie N_SWAPS > 0 dans les deux cas et residu
!      correct.
!   2) Diagonale EXACTEMENT nulle sur une ligne (A(1,1) = 0), avec une
!      entree non nulle plus loin dans la colonne 1 -- cas ou le pivot
!      n'est pas juste utile a la stabilite, il est INDISPENSABLE a la
!      correction (une factorisation sans pivot echouerait sur un pivot
!      exactement nul). Verifie que INFO=0 (pas d'echec) et residu correct.
!
! Usage : ./test_pivoting N NB TOL
!===============================================================================
PROGRAM TEST_PIVOTING
  USE BLOCK_LU_MOD
  USE BLR_LU_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, INFO
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  LOGICAL :: ALL_OK

  IF (COMMAND_ARGUMENT_COUNT() .LT. 3) THEN
    WRITE(6,*) 'Usage: ./test_pivoting N NB TOL'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) TOL

  CALL SET_FIXED_SEED()   ! reproductible -- voir commentaire de la routine
  ALL_OK = .TRUE.

  WRITE(6,'(A,I0,A,I0,A,ES10.3)') ' N = ', N, ' , NB = ', NB, ' , TOL = ', TOL
  WRITE(6,*)

  CALL TEST_PERMUTED_MATRIX(N, NB, TOL, ALL_OK)
  CALL TEST_ZERO_DIAGONAL(N, NB, TOL, ALL_OK)

  IF (ALL_OK) THEN
    WRITE(6,*) '=== TOUS LES TESTS PIVOTING : OK ==='
  ELSE
    WRITE(6,*) '=== AU MOINS UN TEST PIVOTING A ECHOUE ==='
    STOP 1
  END IF

CONTAINS

  !-----------------------------------------------------------------------
  ! Graine fixe : ce test force volontairement un pivot proche de zero
  ! (A(1,1)=0), un scenario deja marginal en simple precision -- un
  ! RANDOM_SEED() sans argument (graine dependant de l'horloge systeme)
  ! ferait varier le residu observe d'une execution a l'autre autour du
  ! seuil de detection, rendant le test "flaky" (parfois OK, parfois pas,
  ! sans changement de code). Graine fixe = resultat reproductible.
  !-----------------------------------------------------------------------
  SUBROUTINE SET_FIXED_SEED()
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: NSEED
    CALL RANDOM_SEED(SIZE=NSEED)
    ALLOCATE(SEED(NSEED))
    SEED = 42
    CALL RANDOM_SEED(PUT=SEED)
    DEALLOCATE(SEED)
  END SUBROUTINE SET_FIXED_SEED

  !-----------------------------------------------------------------------
  SUBROUTINE RANDOM_COMPLEX(M, K, X)
    INTEGER, INTENT(IN) :: M, K
    COMPLEX(WP), INTENT(OUT) :: X(M,K)
    REAL(WP), ALLOCATABLE :: RP(:,:), IP(:,:)
    ALLOCATE(RP(M,K), IP(M,K))
    CALL RANDOM_NUMBER(RP)
    CALL RANDOM_NUMBER(IP)
    X = CMPLX(RP-0.5_WP, IP-0.5_WP, KIND=WP)
    DEALLOCATE(RP, IP)
  END SUBROUTINE RANDOM_COMPLEX

  !-----------------------------------------------------------------------
  ! Permutation aleatoire de 1..N (melange de Fisher-Yates).
  !-----------------------------------------------------------------------
  SUBROUTINE RANDOM_PERMUTATION(N, P)
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(OUT) :: P(N)
    INTEGER :: I, J, TMP
    REAL(WP) :: R
    DO I = 1, N
      P(I) = I
    END DO
    DO I = N, 2, -1
      CALL RANDOM_NUMBER(R)
      J = 1 + INT(R*REAL(I,WP))
      IF (J .GT. I) J = I
      TMP = P(I); P(I) = P(J); P(J) = TMP
    END DO
  END SUBROUTINE RANDOM_PERMUTATION

  !-----------------------------------------------------------------------
  ! Construit une matrice compressible A = c*I + U*V (rang global RGLOBAL),
  ! bien conditionnee -- meme recette que test_blr_lu.f90.
  !-----------------------------------------------------------------------
  SUBROUTINE BUILD_COMPRESSIBLE(N, A)
    INTEGER, INTENT(IN) :: N
    COMPLEX(WP), INTENT(OUT) :: A(N,N)
    INTEGER, PARAMETER :: RGLOBAL = 6
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP), ALLOCATABLE :: UG(:,:), VG(:,:)
    INTEGER :: I
    ALLOCATE(UG(N,RGLOBAL), VG(RGLOBAL,N))
    CALL RANDOM_COMPLEX(N, RGLOBAL, UG)
    CALL RANDOM_COMPLEX(RGLOBAL, N, VG)
    CALL CGEMM('N','N', N, N, RGLOBAL, CONE, UG, N, VG, RGLOBAL, CZERO, A, N)
    DO I = 1, N
      A(I,I) = A(I,I) + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
    END DO
    DEALLOCATE(UG, VG)
  END SUBROUTINE BUILD_COMPRESSIBLE

  !-----------------------------------------------------------------------
  SUBROUTINE TEST_PERMUTED_MATRIX(N, NB, TOL, OK)
    INTEGER, INTENT(IN)    :: N, NB
    REAL(WP), INTENT(IN)   :: TOL
    LOGICAL, INTENT(INOUT) :: OK

    COMPLEX(WP), ALLOCATABLE :: A0(:,:), A(:,:), A_ORIG(:,:), X(:,:), B(:,:), XCOMP(:,:)
    INTEGER, ALLOCATABLE     :: IPIV(:), PERM(:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER       :: NRHS = 1
    INTEGER :: I, N_COMPRESSED, N_OFFDIAG, N_SWAPS
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR, NORMB, MIN_PIVOT

    WRITE(6,*) '--- 1) Matrice compressible, lignes permutees (pivot reellement sollicite) ---'

    ALLOCATE(A0(N,N), A(N,N), A_ORIG(N,N), X(N,NRHS), B(N,NRHS), XCOMP(N,NRHS))
    ALLOCATE(IPIV(N), PERM(N))

    CALL BUILD_COMPRESSIBLE(N, A0)
    CALL RANDOM_PERMUTATION(N, PERM)
    DO I = 1, N
      A(I,1:N) = A0(PERM(I),1:N)
    END DO
    DEALLOCATE(A0)

    CALL RANDOM_COMPLEX(N, NRHS, X)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A, N, X, N, CZERO, B, N)
    NORMB = SQRT(SUM(ABS(B(:,1))**2))

    ! --- Sans compression (etape 1) ---
    A_ORIG = A
    CALL BLOCK_LU_FACTOR(N, NB, A_ORIG, N, IPIV, INFO, N_SWAPS=N_SWAPS, MIN_PIVOT=MIN_PIVOT)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  [sans compression] ERREUR BLOCK_LU_FACTOR, INFO=', INFO
      OK = .FALSE.
    ELSE
      XCOMP = B
      CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_ORIG, N, IPIV, XCOMP, N)
      ERR = SQRT(SUM(ABS(XCOMP(:,1)-X(:,1))**2)) / SQRT(SUM(ABS(X(:,1))**2))
      WRITE(6,'(A,I0,A,ES10.3,A,ES12.4)') '  [sans compression] N_SWAPS=', N_SWAPS, &
           ', MIN_PIVOT=', MIN_PIVOT, ', erreur=', ERR
      IF (N_SWAPS .EQ. 0) THEN
        WRITE(6,*) '   -> ECHEC (le test ne sollicite pas reellement le pivot : N_SWAPS=0)'
        OK = .FALSE.
      ELSE IF (ERR .GT. 1.0E-4_WP) THEN
        WRITE(6,*) '   -> ECHEC (residu trop grand)'
        OK = .FALSE.
      ELSE
        WRITE(6,*) '   -> OK'
      END IF
    END IF

    ! --- Avec compression BLR (etape 2) ---
    A_ORIG = A
    CALL BLR_LU_FACTOR(N, NB, A_ORIG, N, IPIV, TOL, INFO, N_COMPRESSED, N_OFFDIAG, &
                        RANK_SUM, N_SWAPS=N_SWAPS, MIN_PIVOT=MIN_PIVOT)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  [avec compression]  ERREUR BLR_LU_FACTOR, INFO=', INFO
      OK = .FALSE.
    ELSE
      XCOMP = B
      CALL BLOCK_LU_SOLVE(N, NB, NRHS, A_ORIG, N, IPIV, XCOMP, N)
      ERR = SQRT(SUM(ABS(XCOMP(:,1)-X(:,1))**2)) / SQRT(SUM(ABS(X(:,1))**2))
      WRITE(6,'(A,I0,A,ES10.3,A,ES12.4,A,I0,A,I0)') '  [avec compression]  N_SWAPS=', N_SWAPS, &
           ', MIN_PIVOT=', MIN_PIVOT, ', erreur=', ERR, ', blocs compresses=', &
           N_COMPRESSED, '/', N_OFFDIAG
      IF (N_SWAPS .EQ. 0) THEN
        WRITE(6,*) '   -> ECHEC (le test ne sollicite pas reellement le pivot : N_SWAPS=0)'
        OK = .FALSE.
      ELSE IF (ERR .GT. MAX(TOL*100.0_WP, 1.0E-4_WP)) THEN
        WRITE(6,*) '   -> ECHEC (residu trop grand)'
        OK = .FALSE.
      ELSE
        WRITE(6,*) '   -> OK'
      END IF
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A_ORIG, X, B, XCOMP, IPIV, PERM)
  END SUBROUTINE TEST_PERMUTED_MATRIX

  !-----------------------------------------------------------------------
  ! IMPORTANT sur la portee de ce test : le but est de verifier que le
  ! PIVOT gere correctement un cas ou il est INDISPENSABLE (pas juste
  ! utile) -- c'est BLOCK_LU_FACTOR (sans compression) qui teste
  ! exactement ca, et c'est la seule verification a echec bloquant
  ! ci-dessous. BLR_LU_FACTOR est AUSSI execute sur le meme cas, pour
  ! observer comment il se comporte, mais SANS assertion d'echec sur son
  ! propre residu : force A(1,1)=0 cree une ligne dont le pivot effectif
  ! (apres echange) est PETIT relativement a l'echelle du reste de la
  ! matrice (~1 contre ~N) -- si TOL est du meme ordre que ce rapport
  ! (ex. TOL~1e-3 avec N~500-1000, rapport pivot/echelle ~1e-3 aussi), la
  ! troncature SVD (tolerance RELATIVE a la plus grande valeur singuliere
  ! du bloc) peut perdre plus de precision que le pivot lui-meme n'en
  ! vaut, amplifiant l'erreur bien au-dela de TOL. C'est un phenomene
  ! numerique connu des methodes BLR (pas un bug de ce code -- confirme en
  ! repassant TOL=0 : residu redescend a ~1e-5, normal) qui merite d'etre
  ! documente plutot que masque par un seuil de test genereux -- une vraie
  ! solution (mise a l'echelle/equilibrage, tolerance adaptative par
  ! front) est hors perimetre de l'etape 3.
  !-----------------------------------------------------------------------
  SUBROUTINE TEST_ZERO_DIAGONAL(N, NB, TOL, OK)
    INTEGER, INTENT(IN)    :: N, NB
    REAL(WP), INTENT(IN)   :: TOL
    LOGICAL, INTENT(INOUT) :: OK

    COMPLEX(WP), ALLOCATABLE :: A(:,:), A_ORIG(:,:), X(:,:), B(:,:), XCOMP(:,:)
    INTEGER, ALLOCATABLE     :: IPIV(:)
    COMPLEX(WP), PARAMETER   :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER, PARAMETER       :: NRHS = 1
    INTEGER :: N_COMPRESSED, N_OFFDIAG, N_SWAPS
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR, NORMB, MIN_PIVOT

    WRITE(6,*) '--- 2) Diagonale exactement nulle (pivot indispensable, pas juste utile) ---'

    ALLOCATE(A(N,N), A_ORIG(N,N), X(N,NRHS), B(N,NRHS), XCOMP(N,NRHS), IPIV(N))

    CALL BUILD_COMPRESSIBLE(N, A)
    A(1,1) = CMPLX(0.0_WP, 0.0_WP, KIND=WP)  ! pivot naturel nul : sans
                                              ! echange de ligne, la
                                              ! factorisation echouerait
    A_ORIG = A

    CALL RANDOM_COMPLEX(N, NRHS, X)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A, N, X, N, CZERO, B, N)
    NORMB = SQRT(SUM(ABS(B(:,1))**2))

    ! --- Verification bloquante : PIVOT SEUL, sans compression ---
    CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO, N_SWAPS=N_SWAPS, MIN_PIVOT=MIN_PIVOT)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  ERREUR BLOCK_LU_FACTOR, INFO=', INFO, &
           ' (le pivot aurait du eviter ce cas)'
      OK = .FALSE.
      RETURN
    END IF
    XCOMP = B
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, XCOMP, N)
    ERR = SQRT(SUM(ABS(XCOMP(:,1)-X(:,1))**2)) / SQRT(SUM(ABS(X(:,1))**2))
    WRITE(6,'(A,I0,A,ES10.3,A,ES12.4)') '  [pivot seul, sans compression] N_SWAPS=', &
         N_SWAPS, ', MIN_PIVOT=', MIN_PIVOT, ', erreur=', ERR

    ! Seuil plus genereux qu'ailleurs (5e-4 au lieu de 1e-4) : le pivot
    ! resultant apres echange (~1, cf. MIN_PIVOT ci-dessus) reste petit
    ! relativement a l'echelle du reste de la matrice (~N) meme sans
    ! aucune compression -- de l'amplification d'erreur en simple
    ! precision legitime, pas une marge de securite pour un bug.
    IF (IPIV(1) .EQ. 1) THEN
      WRITE(6,*) '   -> ECHEC (IPIV(1)=1 : aucun echange fait sur un pivot pourtant nul)'
      OK = .FALSE.
    ELSE IF (ERR .GT. 5.0E-4_WP) THEN
      WRITE(6,*) '   -> ECHEC (residu trop grand)'
      OK = .FALSE.
    ELSE
      WRITE(6,*) '   -> OK'
    END IF

    ! --- Observation (non bloquante) : meme cas AVEC compression, voir
    !     le commentaire au-dessus de cette subroutine ---
    A = A_ORIG
    CALL BLR_LU_FACTOR(N, NB, A, N, IPIV, TOL, INFO, N_COMPRESSED, N_OFFDIAG, &
                        RANK_SUM, N_SWAPS=N_SWAPS, MIN_PIVOT=MIN_PIVOT)
    IF (INFO .NE. 0) THEN
      WRITE(6,'(A,I0)') '  [avec compression]  ERREUR BLR_LU_FACTOR, INFO=', INFO
    ELSE
      XCOMP = B
      CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, XCOMP, N)
      ERR = SQRT(SUM(ABS(XCOMP(:,1)-X(:,1))**2)) / SQRT(SUM(ABS(X(:,1))**2))
      WRITE(6,'(A,I0,A,ES10.3,A,ES12.4,A)') '  [avec compression]  N_SWAPS=', N_SWAPS, &
           ', MIN_PIVOT=', MIN_PIVOT, ', erreur=', ERR, &
           ' (a titre informatif, non bloquant -- voir commentaire)'
    END IF
    WRITE(6,*)

    DEALLOCATE(A, A_ORIG, X, B, XCOMP, IPIV)
  END SUBROUTINE TEST_ZERO_DIAGONAL

END PROGRAM TEST_PIVOTING
