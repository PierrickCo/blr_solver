!===============================================================================
! Micro-benchmark des noyaux du bloc diagonal (sans MPI) : compare, sur des
! blocs KW x KW complexes symetriques, la factorisation LU (CGETRF, celle du
! solveur LU) aux variantes LDL^T de LAPACK, et mesure les couts annexes du
! stockage compact et de la mise a jour de Schur diagonale.
!
!   getrf       CGETRF                         (LU, pivot partiel)
!   sytrf_rk    CSYTRF_RK                      (utilise aujourd'hui par le LDLt)
!   sytrf       CSYTRF                         (Bunch-Kaufman classique)
!   sytrf+conv  CSYTRF puis CSYCONVF('C')      (meme format de sortie que _RK)
!   rook        CSYTRF_ROOK                    (Bunch-Kaufman "rook")
!   aa          CSYTRF_AA                      (Aasen, L*T*L^T, T tridiagonale)
!   pack+unpk   CTPTTR puis CTRTTP             (cout d'un aller-retour du
!                                               stockage compact, par mise a jour)
!   gemm r      CGEMM KW x KW x r              (Schur diagonale, bloc entier)
!   gemmt r     CGEMMT KW x KW x r, triangle inferieur seul (MKL)
!
! Matrices : bloc diagonal de la matrice de demo (sym_testmat_mod.f90 :
! diagonale faible + anti-diagonale forte -> pivots 2x2, cas le plus
! defavorable pour Bunch-Kaufman) et bloc symetrique aleatoire. Temps = MINIMUM
! sur NREP repetitions (copie de la matrice hors chronometre). GF/s calcules
! avec les flops nominaux reels : 8/3 KW^3 (LU complexe), 4/3 KW^3 (LDLt
! complexe), 8 KW^2 r (GEMM), 4 KW^2 r (GEMMT).
!
! Usage : ./bench_diag_kernels [NREP] [RANK]     (defauts 5 8)
!         threads : OMP_NUM_THREADS / MKL_NUM_THREADS
!===============================================================================
PROGRAM BENCH_DIAG_KERNELS
  USE BLR_KINDS_MOD, ONLY: WP
  USE SYM_TESTMAT_MOD
  USE OMP_LIB
  IMPLICIT NONE
  INTEGER, PARAMETER :: NSIZES = 5
  INTEGER, PARAMETER :: SIZES(NSIZES) = (/ 200, 400, 600, 800, 1000 /)
  INTEGER :: NREP, RANKR, IS, KW, IM, IREP, INFO, LWORK, I, J
  COMPLEX(WP), ALLOCATABLE :: A0(:,:), A(:,:), E(:), WORK(:), AP(:), U(:,:), V(:,:)
  COMPLEX(WP) :: WQ(1)
  INTEGER, ALLOCATABLE :: IPIV(:)
  DOUBLE PRECISION :: T0, T, TMIN(9), FL
  CHARACTER(LEN=32) :: ARG
  CHARACTER(LEN=12) :: MATNAME
  REAL :: XR, XI
  CHARACTER(LEN=12), PARAMETER :: KNAMES(9) = (/ 'getrf       ', 'sytrf_rk    ', 'sytrf       ', &
       'sytrf+conv  ', 'rook        ', 'aa          ', 'pack+unpk   ', 'gemm r      ', 'gemmt r     ' /)

  NREP = 5; RANKR = 8
  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN; CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) NREP; END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN; CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) RANKR; END IF
  WRITE(6,'(A,I0,A,I0,A,I0)') ' Noyaux du bloc diagonal -- threads=', OMP_GET_MAX_THREADS(), ' , NREP=', NREP, &
       ' , rang Schur r=', RANKR

  DO IM = 1, 2
    IF (IM .EQ. 1) THEN
      MATNAME = 'demo (2x2)'
    ELSE
      MATNAME = 'aleatoire'
    END IF
    WRITE(6,'(/,A,A)') ' === Matrice : ', TRIM(MATNAME)
    WRITE(6,'(A6,9A13)') '   KW', (KNAMES(J), J=1,9)
    DO IS = 1, NSIZES
      KW = SIZES(IS)
      ALLOCATE(A0(KW,KW), A(KW,KW), E(KW), IPIV(KW), AP(KW*(KW+1)/2), U(KW,RANKR), V(KW,RANKR))
      IF (IM .EQ. 1) THEN
        CALL SYM_TESTMAT_SETUP(KW, KW)                 ! un seul bloc : anti-diagonale forte
        CALL SYM_GET_BLOCK(1, KW, 1, KW, A0)
      ELSE
        DO J = 1, KW
          DO I = J, KW
            CALL RANDOM_NUMBER(XR); CALL RANDOM_NUMBER(XI)
            A0(I,J) = CMPLX(XR-0.5, XI-0.5, KIND=WP); A0(J,I) = A0(I,J)
          END DO
        END DO
      END IF
      DO J = 1, RANKR
        DO I = 1, KW
          CALL RANDOM_NUMBER(XR); CALL RANDOM_NUMBER(XI)
          U(I,J) = CMPLX(XR, XI, KIND=WP); V(I,J) = CMPLX(XI, XR, KIND=WP)
        END DO
      END DO
      CALL CSYTRF_RK('L', KW, A, KW, E, IPIV, WQ, -1, INFO); LWORK = INT(REAL(WQ(1)))
      CALL CSYTRF('L', KW, A, KW, IPIV, WQ, -1, INFO);        LWORK = MAX(LWORK, INT(REAL(WQ(1))))
      CALL CSYTRF_ROOK('L', KW, A, KW, IPIV, WQ, -1, INFO);   LWORK = MAX(LWORK, INT(REAL(WQ(1))))
      CALL CSYTRF_AA('L', KW, A, KW, IPIV, WQ, -1, INFO);     LWORK = MAX(LWORK, INT(REAL(WQ(1))))
      ALLOCATE(WORK(MAX(LWORK,1)))
      TMIN = HUGE(1.0D0)
      DO IREP = 1, NREP
        A = A0; T0 = OMP_GET_WTIME(); CALL CGETRF(KW, KW, A, KW, IPIV, INFO)
        T = OMP_GET_WTIME()-T0; CALL CHK(INFO, 1); TMIN(1) = MIN(TMIN(1), T)
        A = A0; T0 = OMP_GET_WTIME(); CALL CSYTRF_RK('L', KW, A, KW, E, IPIV, WORK, LWORK, INFO)
        T = OMP_GET_WTIME()-T0; CALL CHK(INFO, 2); TMIN(2) = MIN(TMIN(2), T)
        A = A0; T0 = OMP_GET_WTIME(); CALL CSYTRF('L', KW, A, KW, IPIV, WORK, LWORK, INFO)
        T = OMP_GET_WTIME()-T0; CALL CHK(INFO, 3); TMIN(3) = MIN(TMIN(3), T)
        A = A0; T0 = OMP_GET_WTIME(); CALL CSYTRF('L', KW, A, KW, IPIV, WORK, LWORK, INFO)
        CALL CSYCONVF('L', 'C', KW, A, KW, E, IPIV, INFO)
        T = OMP_GET_WTIME()-T0; CALL CHK(INFO, 4); TMIN(4) = MIN(TMIN(4), T)
        A = A0; T0 = OMP_GET_WTIME(); CALL CSYTRF_ROOK('L', KW, A, KW, IPIV, WORK, LWORK, INFO)
        T = OMP_GET_WTIME()-T0; CALL CHK(INFO, 5); TMIN(5) = MIN(TMIN(5), T)
        A = A0; T0 = OMP_GET_WTIME(); CALL CSYTRF_AA('L', KW, A, KW, IPIV, WORK, LWORK, INFO)
        T = OMP_GET_WTIME()-T0; CALL CHK(INFO, 6); TMIN(6) = MIN(TMIN(6), T)
        A = A0; CALL CTRTTP('L', KW, A, KW, AP, INFO)
        T0 = OMP_GET_WTIME(); CALL CTPTTR('L', KW, AP, A, KW, INFO); CALL CTRTTP('L', KW, A, KW, AP, INFO)
        T = OMP_GET_WTIME()-T0; TMIN(7) = MIN(TMIN(7), T)
        A = A0; T0 = OMP_GET_WTIME()
        CALL CGEMM('N', 'T', KW, KW, RANKR, (-1.0_WP,0.0_WP), U, KW, V, KW, (1.0_WP,0.0_WP), A, KW)
        T = OMP_GET_WTIME()-T0; TMIN(8) = MIN(TMIN(8), T)
        A = A0; T0 = OMP_GET_WTIME()
        CALL CGEMMT('L', 'N', 'T', KW, RANKR, (-1.0_WP,0.0_WP), U, KW, V, KW, (1.0_WP,0.0_WP), A, KW)
        T = OMP_GET_WTIME()-T0; TMIN(9) = MIN(TMIN(9), T)
      END DO
      WRITE(6,'(I6,9F13.5)') KW, TMIN
      FL = REAL(KW,8)**3
      WRITE(6,'(A6,6F13.2,A13,2F13.2)') '  GF/s', 8.0D0/3.0D0*FL/TMIN(1)*1.0D-9, &
           (4.0D0/3.0D0*FL/TMIN(J)*1.0D-9, J=2,6), '-', &
           8.0D0*REAL(KW,8)**2*RANKR/TMIN(8)*1.0D-9, 4.0D0*REAL(KW,8)**2*RANKR/TMIN(9)*1.0D-9
      WRITE(6,'(A6,6F13.2)') ' t/LU', (TMIN(J)/TMIN(1), J=1,6)
      DEALLOCATE(A0, A, E, IPIV, AP, U, V, WORK)
    END DO
  END DO

CONTAINS

  SUBROUTINE CHK(INFO_, K)
    INTEGER, INTENT(IN) :: INFO_, K
    IF (INFO_ .NE. 0) WRITE(6,'(A,A,A,I0)') ' ATTENTION ', TRIM(KNAMES(K)), ' INFO=', INFO_
  END SUBROUTINE CHK

END PROGRAM BENCH_DIAG_KERNELS
