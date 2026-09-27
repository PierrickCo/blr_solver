!===============================================================================
! Interfaces explicites pour les routines BLAS/LAPACK/MKL appelees sans
! module dedie ailleurs dans ce dossier (CGEMM, CTRSM, CGETRF...) --
! elimine les avertissements ifx "explicit interface required" sans
! changer aucun site d'appel. Signatures copiees verbatim depuis les
! fichiers d'interface officiels de MKL (mkl_blas.fi, mkl_lapack.fi,
! mkl_service.fi, $MKLROOT/include/) plutot que d'inclure ces fichiers
! en entier (mkl_lapack.fi fait ~20000 lignes) dans chaque module.
!===============================================================================
MODULE MKL_INTERFACES_MOD
  IMPLICIT NONE

  ! CSWAP reste en simple declaration EXTERNAL (pas d'INTERFACE) : plusieurs
  ! sites d'appel lui passent un ELEMENT de tableau POINTER (BLK%UFAC(I,1),
  ! remappage arena) comme point de depart d'un parcours par pas INCX. C'est
  ! une association sequentielle interdite par la norme vers un argument
  ! muet assumed-size (*) des qu'une INTERFACE explicite impose la
  ! verification de forme (erreur ifx #8284) -- une interface implicite
  ! (EXTERNAL seul) ne fait pas cette verification et laisse passer ce motif
  ! F77 classique, exactement comme avant l'ajout de ce module. EXTERNAL
  ! seul suffit deja a satisfaire l'avertissement #8889.
  EXTERNAL :: CSWAP

  INTERFACE
    SUBROUTINE CGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      CHARACTER*1 :: TRANSA, TRANSB
      INTEGER     :: M, N, K, LDA, LDB, LDC
      COMPLEX     :: ALPHA, BETA
      COMPLEX     :: A(LDA,*), B(LDB,*), C(LDC,*)
    END SUBROUTINE CGEMM

    SUBROUTINE CTRSM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
      CHARACTER*1 :: SIDE, UPLO, TRANSA, DIAG
      INTEGER     :: M, N, LDA, LDB
      COMPLEX     :: ALPHA
      COMPLEX     :: A(LDA,*), B(LDB,*)
    END SUBROUTINE CTRSM

    SUBROUTINE CSCAL(N, CA, CX, INCX)
      INTEGER :: N, INCX
      COMPLEX :: CA, CX(*)
    END SUBROUTINE CSCAL

    SUBROUTINE CGERU(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
      INTEGER :: M, N, INCX, INCY, LDA
      COMPLEX :: ALPHA
      COMPLEX :: X(*), Y(*), A(LDA,*)
    END SUBROUTINE CGERU

    SUBROUTINE CGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO)
      CHARACTER*1 :: JOBU, JOBVT
      INTEGER     :: M, N, LDA, LDU, LDVT, LWORK, INFO
      REAL        :: S(*), RWORK(*)
      COMPLEX     :: A(LDA,*), U(LDU,*), VT(LDVT,*), WORK(*)
    END SUBROUTINE CGESVD

    SUBROUTINE CGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
      INTEGER :: M, N, LDA, LWORK, INFO
      COMPLEX :: A(LDA,*), TAU(*), WORK(*)
    END SUBROUTINE CGEQRF

    SUBROUTINE CUNGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
      INTEGER :: M, N, K, LDA, LWORK, INFO
      COMPLEX :: A(LDA,*), TAU(*), WORK(*)
    END SUBROUTINE CUNGQR

    SUBROUTINE CGETRF(M, N, A, LDA, IPIV, INFO)
      INTEGER :: M, N, LDA, INFO
      INTEGER :: IPIV(*)
      COMPLEX :: A(LDA,*)
    END SUBROUTINE CGETRF

    SUBROUTINE MKL_SET_NUM_THREADS(NTHRS)
      INTEGER :: NTHRS
    END SUBROUTINE MKL_SET_NUM_THREADS
  END INTERFACE

END MODULE MKL_INTERFACES_MOD
