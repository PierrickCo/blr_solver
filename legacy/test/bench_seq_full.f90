!===============================================================================
! BENCHMARK SEQUENTIEL : resolution d'un systeme N=10000, NRHS=100, sur
! trois solveurs (meme matrice, memes seconds membres, meme machine) --
! demande explicite de l'utilisateur pour voir concretement l'effet de la
! proposition 2 (stockage compresse) sur un cas realiste, meme methodologie
! que l'analyse de performance faite sur MUMPS.
!
!   1) Dense (BLOCK_LU_FACTOR/BLOCK_LU_SOLVE, etape 1) -- reference.
!   2) BLR, facteurs stockes DENSES (BLR_LU_FACTOR, etape 2/5 + BLOCK_LU_
!      SOLVE) -- la compression n'accelere QUE la factorisation, la
!      memoire du resultat stocke reste celle d'un tableau N x N dense.
!   3) BLR, facteurs stockes COMPRESSES (BLR_LU_FACTOR_COMPRESSED +
!      BLR_LU_SOLVE_COMPRESSED, proposition 2) -- c'est LA comparaison
!      demandee : memoire du stockage FINAL, une fois le tableau dense de
!      travail libere (voir DEALLOCATE(A) explicite ci-dessous, pour
!      mesurer ce qu'une application garderait reellement en memoire sur
!      la duree, pas le pic transitoire de la factorisation).
!
! Mesures par cas : temps d'ASSEMBLAGE (construction de la matrice/second
! membre), de FACTORISATION, de RESOLUTION (NRHS=100), erreur relative
! (||Xcalc-Xexact||/||Xexact||, norme de Frobenius sur les 100 colonnes),
! et memoire -- a la fois ANALYTIQUE (taille exacte des structures
! stockees, deja fiable depuis la proposition 2) et RSS REEL (/proc/self/
! status, bench_utils_mod.f90) pour se recouper.
!
! Matrice : meme recette "diagonale dominante + correction basse-rang"
! (rang global 6) que partout ailleurs dans ce projet, formule
! DETERMINISTE (pas de RANDOM_NUMBER pour A, comme les tests MPI) --
! permet une comparaison directe avec le benchmark MPI (meme matrice
! EXACTEMENT). Second membre B=A*Xexact avec Xexact aleatoire (a
! reproductibilite fixee).
!
! Usage : ./bench_seq [N] [NRHS] [NB] [TOL]  (defauts : 10000 100 100 1e-4)
!===============================================================================
PROGRAM BENCH_SEQ
  USE BLOCK_LU_MOD
  USE BLR_LU_MOD
  USE BLR_LU_STORAGE_MOD
  USE BENCH_UTILS_MOD
  IMPLICIT NONE

  INTEGER :: N, NRHS, NB, I
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  INTEGER, PARAMETER :: RGLOBAL = 6

  N = 10000; NRHS = 100; NB = 100; TOL = 1.0E-4_WP
  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) TOL
  END IF

  WRITE(6,'(A,I0,A,I0,A,I0,A,ES10.3)') ' === BENCHMARK SEQUENTIEL : N=', N, &
       ' , NRHS=', NRHS, ' , NB=', NB, ' , TOL=', TOL
  WRITE(6,*)

  CALL RUN_DENSE(N, NRHS, NB)
  CALL RUN_BLR_DENSE_STORAGE(N, NRHS, NB, TOL)
  CALL RUN_BLR_COMPRESSED_STORAGE(N, NRHS, NB, TOL)

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE BUILD_MATRIX(N, A, T_ASSEMBLY)
    INTEGER, INTENT(IN)        :: N
    COMPLEX(WP), INTENT(OUT)   :: A(N,N)
    DOUBLE PRECISION, INTENT(OUT) :: T_ASSEMBLY
    INTEGER :: GI, GJ, R
    COMPLEX(WP) :: AVAL
    DOUBLE PRECISION :: T0, T1

    T0 = OMP_GET_WTIME_L()
    DO GJ = 1, N
      DO GI = 1, N
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(GI,R)*LOWRANK_V(R,GJ)
        END DO
        IF (GI .EQ. GJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        A(GI,GJ) = AVAL
      END DO
    END DO
    T1 = OMP_GET_WTIME_L()
    T_ASSEMBLY = T1 - T0
  END SUBROUTINE BUILD_MATRIX

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

  DOUBLE PRECISION FUNCTION OMP_GET_WTIME_L()
    USE OMP_LIB, ONLY: OMP_GET_WTIME
    OMP_GET_WTIME_L = OMP_GET_WTIME()
  END FUNCTION OMP_GET_WTIME_L

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

  SUBROUTINE SET_FIXED_SEED()
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: NSEED
    CALL RANDOM_SEED(SIZE=NSEED)
    ALLOCATE(SEED(NSEED))
    SEED = 42
    CALL RANDOM_SEED(PUT=SEED)
    DEALLOCATE(SEED)
  END SUBROUTINE SET_FIXED_SEED

  SUBROUTINE REPORT(LABEL, T_ASM, T_COMPRESS, T_FACT_PURE, T_SOLVE, ERR, &
                     MEM_ANALYTIC_MB, RSS_PEAK_MB, RSS_FINAL_MB)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    DOUBLE PRECISION, INTENT(IN) :: T_ASM, T_COMPRESS, T_FACT_PURE, T_SOLVE
    REAL(WP), INTENT(IN) :: ERR
    DOUBLE PRECISION, INTENT(IN) :: MEM_ANALYTIC_MB, RSS_PEAK_MB, RSS_FINAL_MB
    WRITE(6,'(A)') ' --- ' // TRIM(LABEL) // ' ---'
    WRITE(6,'(A,F8.3,A)')   '   Temps assemblage              = ', T_ASM, ' s'
    WRITE(6,'(A,F8.3,A)')   '   Temps compression basse-rang  = ', T_COMPRESS, ' s'
    WRITE(6,'(A,F8.3,A)')   '   Temps factorisation (hors compression) = ', T_FACT_PURE, ' s'
    WRITE(6,'(A,F8.3,A)')   '   Temps resolution               = ', T_SOLVE, ' s'
    WRITE(6,'(A,ES12.4)')   '   Erreur relative                = ', ERR
    WRITE(6,'(A,F10.1,A)')  '   Memoire analytique (stockage final) = ', MEM_ANALYTIC_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)')  '   RSS pic (pendant factorisation)     = ', RSS_PEAK_MB, ' Mo'
    WRITE(6,'(A,F10.1,A)')  '   RSS final (apres liberation scratch)= ', RSS_FINAL_MB, ' Mo'
    WRITE(6,*)
  END SUBROUTINE REPORT

  FUNCTION KB_TO_MB(KB) RESULT(MB)
    INTEGER(8), INTENT(IN) :: KB
    DOUBLE PRECISION :: MB
    MB = REAL(KB,8) / 1024.0D0
  END FUNCTION KB_TO_MB

  FUNCTION FROB_ERR(N, NRHS, XCALC, XEX) RESULT(ERR)
    INTEGER, INTENT(IN) :: N, NRHS
    COMPLEX(WP), INTENT(IN) :: XCALC(N,NRHS), XEX(N,NRHS)
    REAL(WP) :: ERR
    ERR = SQRT(SUM(ABS(XCALC-XEX)**2)) / SQRT(SUM(ABS(XEX)**2))
  END FUNCTION FROB_ERR

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_DENSE(N, NRHS, NB)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    COMPLEX(WP), ALLOCATABLE :: A(:,:), XEX(:,:), B(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT, T_SOLVE, RSS_PEAK, RSS_FINAL
    INTEGER :: INFO
    REAL(WP) :: ERR

    CALL SET_FIXED_SEED()
    ALLOCATE(A(N,N))
    CALL BUILD_MATRIX(N, A, T_ASM)

    ALLOCATE(XEX(N,NRHS), B(N,NRHS), IPIV(N))
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A, N, XEX, N, CZERO, B, N)

    T0 = OMP_GET_WTIME_L()
    CALL BLOCK_LU_FACTOR(N, NB, A, N, IPIV, INFO)
    T1 = OMP_GET_WTIME_L()
    T_FACT = T1 - T0
    RSS_PEAK = KB_TO_MB(GET_RSS_KB())

    T0 = OMP_GET_WTIME_L()
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, B, N)
    T1 = OMP_GET_WTIME_L()
    T_SOLVE = T1 - T0

    ERR = FROB_ERR(N, NRHS, B, XEX)
    RSS_FINAL = KB_TO_MB(GET_RSS_KB())

    ! Pas de compression dans ce cas (etape 1 pure) : temps de compression
    ! nul par construction, pas mesure separement.
    CALL REPORT('1) Dense (etape 1)', T_ASM, 0.0D0, T_FACT, T_SOLVE, ERR, &
                8.0D0*REAL(N,8)*REAL(N,8)/1024.0D0/1024.0D0, RSS_PEAK, RSS_FINAL)

    DEALLOCATE(A, XEX, B, IPIV)
  END SUBROUTINE RUN_DENSE

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_BLR_DENSE_STORAGE(N, NRHS, NB, TOL)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    REAL(WP), INTENT(IN) :: TOL
    COMPLEX(WP), ALLOCATABLE :: A(:,:), XEX(:,:), B(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT_TOTAL, T_COMPRESS, T_SOLVE, RSS_PEAK, RSS_FINAL
    INTEGER :: INFO, N_COMPRESSED, N_OFFDIAG
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR

    CALL SET_FIXED_SEED()
    ALLOCATE(A(N,N))
    CALL BUILD_MATRIX(N, A, T_ASM)

    ALLOCATE(XEX(N,NRHS), B(N,NRHS), IPIV(N))
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A, N, XEX, N, CZERO, B, N)

    T0 = OMP_GET_WTIME_L()
    CALL BLR_LU_FACTOR(N, NB, A, N, IPIV, TOL, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                        T_COMPRESS=T_COMPRESS)
    T1 = OMP_GET_WTIME_L()
    T_FACT_TOTAL = T1 - T0
    RSS_PEAK = KB_TO_MB(GET_RSS_KB())

    T0 = OMP_GET_WTIME_L()
    CALL BLOCK_LU_SOLVE(N, NB, NRHS, A, N, IPIV, B, N)
    T1 = OMP_GET_WTIME_L()
    T_SOLVE = T1 - T0

    ERR = FROB_ERR(N, NRHS, B, XEX)
    RSS_FINAL = KB_TO_MB(GET_RSS_KB())

    CALL REPORT('2) BLR, facteurs stockes DENSES (etape 2/5)', T_ASM, T_COMPRESS, &
                T_FACT_TOTAL-T_COMPRESS, T_SOLVE, ERR, &
                8.0D0*REAL(N,8)*REAL(N,8)/1024.0D0/1024.0D0, RSS_PEAK, RSS_FINAL)
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   (compression : ', N_COMPRESSED, ' / ', &
         N_OFFDIAG, ' blocs, ', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F8.2,A,I0,A,F6.1,A)') '   (rang moyen = ', REAL(RANK_SUM,WP)/REAL(MAX(N_OFFDIAG,1),WP), &
         ' sur NB=', NB, ' -> memoire par bloc ~', &
         200.0_WP*REAL(RANK_SUM,WP)/(REAL(N_OFFDIAG,WP)*REAL(NB,WP)), ' % du bloc dense (UFAC+VFAC=2*rang*NB vs NB^2)'
    WRITE(6,*)

    DEALLOCATE(A, XEX, B, IPIV)
  END SUBROUTINE RUN_BLR_DENSE_STORAGE

  !-----------------------------------------------------------------------
  SUBROUTINE RUN_BLR_COMPRESSED_STORAGE(N, NRHS, NB, TOL)
    INTEGER, INTENT(IN) :: N, NRHS, NB
    REAL(WP), INTENT(IN) :: TOL
    COMPLEX(WP), ALLOCATABLE :: A(:,:), XEX(:,:), B(:,:)
    TYPE(BLR_FACTORED_LU_T) :: FACT
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    DOUBLE PRECISION :: T0, T1, T_ASM, T_FACT_PURE, T_COMPRESS, T_SOLVE, RSS_PEAK, RSS_FINAL
    INTEGER :: INFO, N_COMPRESSED, N_OFFDIAG
    INTEGER(8) :: RANK_SUM
    REAL(WP) :: ERR

    CALL SET_FIXED_SEED()
    ALLOCATE(A(N,N))
    CALL BUILD_MATRIX(N, A, T_ASM)

    ALLOCATE(XEX(N,NRHS), B(N,NRHS))
    CALL RANDOM_COMPLEX(N, NRHS, XEX)
    CALL CGEMM('N','N', N, NRHS, N, CONE, A, N, XEX, N, CZERO, B, N)

    T0 = OMP_GET_WTIME_L()
    CALL BLR_LU_FACTOR_COMPRESSED(N, NB, A, N, TOL, FACT, INFO, N_COMPRESSED, N_OFFDIAG, RANK_SUM, &
                                   T_FACTOR=T_FACT_PURE, T_COMPRESS=T_COMPRESS)
    T1 = OMP_GET_WTIME_L()
    RSS_PEAK = KB_TO_MB(GET_RSS_KB())

    ! liberer le tableau dense de travail -- c'est CE QUE GARDERAIT une
    ! application reelle sur la duree (voir note de tete de fichier)
    DEALLOCATE(A)
    RSS_FINAL = KB_TO_MB(GET_RSS_KB())

    T0 = OMP_GET_WTIME_L()
    CALL BLR_LU_SOLVE_COMPRESSED(FACT, NRHS, B, N)
    T1 = OMP_GET_WTIME_L()
    T_SOLVE = T1 - T0

    ERR = FROB_ERR(N, NRHS, B, XEX)

    CALL REPORT('3) BLR, facteurs stockes COMPRESSES (proposition 2)', T_ASM, T_COMPRESS, &
                T_FACT_PURE, T_SOLVE, ERR, &
                KB_TO_MB(BLR_FACTORED_LU_MEMORY_BYTES(FACT)/1024_8), RSS_PEAK, RSS_FINAL)
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') '   (compression : ', N_COMPRESSED, ' / ', &
         N_OFFDIAG, ' blocs, ', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F8.2,A,I0,A,F6.1,A)') '   (rang moyen = ', REAL(RANK_SUM,WP)/REAL(MAX(N_OFFDIAG,1),WP), &
         ' sur NB=', NB, ' -> memoire par bloc ~', &
         200.0_WP*REAL(RANK_SUM,WP)/(REAL(N_OFFDIAG,WP)*REAL(NB,WP)), ' % du bloc dense)'
    WRITE(6,'(A,F6.1,A)') '   (taux de compression GLOBAL du stockage final = ', &
         100.0_WP*REAL(BLR_FACTORED_LU_MEMORY_BYTES(FACT),WP)/(8.0_WP*REAL(N,WP)*REAL(N,WP)), &
         ' % du dense N x N)'

    DEALLOCATE(XEX, B)
  END SUBROUTINE RUN_BLR_COMPRESSED_STORAGE

END PROGRAM BENCH_SEQ
