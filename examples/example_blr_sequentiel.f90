!===============================================================================
! EXEMPLE CANONIQUE -- reference A/B (blr_matrix_mod.f90) : solveur LU
! sequentiel/OpenMP dont la matrice reste COMPRESSEE (Block Low-Rank) de
! l'assemblage jusqu'a la resolution, jamais dense en dehors des blocs
! diagonaux et des blocs hors-diagonale tombes en repli dense.
!
! But de ce fichier : montrer l'enchainement minimal des 3 appels qui
! forment ce solveur (assembler, factoriser, resoudre), sans toute
! l'instrumentation de mesure/validation croisee des drivers de test/bench
! (voir test/test_blr_matrix_progress.f90 pour la version qui verifie
! aussi la correction contre une reference CGESV independante).
!
! TOL controle la compression : TOL>0 donne la reference "A" (compresse) ;
! TOL=0 desactive la troncature basse-rang et donne la reference "B"
! (meme pipeline, aucun bloc ne se compresse) -- change juste la valeur
! passee ci-dessous pour comparer les deux.
!
! Compiler : voir INSTALL.md (fait partie de la ligne vivante, deja
! construit par `make`/`make -f Makefile.gcc`). Ce fichier n'est PAS un
! des 13 binaires de run_all_tests.sh -- a compiler a part si besoin :
!   gfortran -O2 -fopenmp -ffree-line-length-none -I../build-gcc \
!     examples/example_blr_sequentiel.f90 ../build-gcc/*.o -lblas -llapack -o example_blr
!
! Usage : ./example_blr [N] [NB] [NRHS] [TOL]   (defauts : 2000 200 4 1e-4)
!===============================================================================
PROGRAM EXAMPLE_BLR_SEQUENTIEL
  USE BLR_MATRIX_MOD
  IMPLICIT NONE

  INTEGER :: N, NB, NRHS, INFO, N_COMP, N_OFF, I, J
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B(:,:)
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACT

  N = 2000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP
  IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) TOL
  END IF

  WRITE(6,'(A,I0,A,I0,A,I0,A,ES9.2)') ' N=', N, ' NB=', NB, ' NRHS=', NRHS, ' TOL=', TOL

  ! --- 1) assemblage bloc par bloc + compression immediate --------------
  ! GET_BLOCK est un callback fourni par l'appelant (voir plus bas) : il
  ! ne calcule qu'UN bloc a la fois, jamais de tableau N x N dense complet
  ! en memoire, meme pour une tres grande matrice.
  CALL ASSEMBLE_BLR_MATRIX(N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, N_COMP, N_OFF, RANK_SUM)
  WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
  WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses : ', N_COMP, ' / ', N_OFF, ' (', &
       100.0_WP*REAL(N_COMP,WP)/REAL(MAX(N_OFF,1),WP), ' %)'

  ! --- 2) factorisation LU directement sur la matrice compressee --------
  CALL BLR_MATRIX_LU_FACTOR(A_BLR, TOL, INFO, N_COMP, N_OFF, RANK_SUM, T_FACT)
  IF (INFO .NE. 0) THEN
    WRITE(6,'(A,I0)') ' ERREUR BLR_MATRIX_LU_FACTOR, INFO=', INFO
    STOP 1
  END IF
  WRITE(6,'(A,F8.3,A)') ' Factorisation LU        : ', T_FACT, ' s'
  WRITE(6,'(A,F10.1,A)') ' Memoire (blocs BLR)     : ', &
       REAL(BLR_MATRIX_MEMORY_BYTES(A_BLR),8)/1024.0D0/1024.0D0, ' Mo'

  ! --- 3) resolution A*X=B, B ecrase par la solution ---------------------
  ! Second membre arbitraire (formule deterministe) -- dans un cas reel,
  ! B viendrait des donnees du probleme, pas d'une formule.
  ALLOCATE(B(N,NRHS))
  DO J = 1, NRHS
    DO I = 1, N
      B(I,J) = CMPLX(REAL(I,WP)*0.01_WP, -0.3_WP*REAL(J,WP), KIND=WP)
    END DO
  END DO
  CALL BLR_MATRIX_LU_SOLVE(A_BLR, NRHS, B, N)

  WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
  DO I = 1, MIN(5, N)
    WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(B(I,1)), ' + i*', AIMAG(B(I,1))
  END DO

  DEALLOCATE(B)

CONTAINS

  ! Matrice "diagonale + basse-rang" a diagonale dominante (ASSEMBLE_BLR_
  ! MATRIX suppose une matrice dense en entree ; SANS PIVOT, comme tous
  ! les solveurs de ce projet -- voir README.md, section "Limites").
  ! Remplacer ce callback par sa propre formule/assemblage pour un cas reel.
  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER, PARAMETER :: RGLOBAL = 6
    INTEGER :: II, JJ, R
    COMPLEX(WP) :: AVAL
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        AVAL = (0.0_WP,0.0_WP)
        DO R = 1, RGLOBAL
          AVAL = AVAL + LOWRANK_U(II,R)*LOWRANK_V(R,JJ)
        END DO
        IF (II .EQ. JJ) AVAL = AVAL + CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        BLK(II-ISTART+1, JJ-JSTART+1) = AVAL
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  COMPLEX(WP) FUNCTION LOWRANK_U(GI, R)
    INTEGER, INTENT(IN) :: GI, R
    LOWRANK_U = CMPLX(SIN(0.1_WP*REAL(GI*R,WP)), COS(0.07_WP*REAL(GI*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_U

  COMPLEX(WP) FUNCTION LOWRANK_V(R, GJ)
    INTEGER, INTENT(IN) :: R, GJ
    LOWRANK_V = CMPLX(COS(0.05_WP*REAL(GJ*R,WP)), SIN(0.03_WP*REAL(GJ*R,WP)), KIND=WP)
  END FUNCTION LOWRANK_V

END PROGRAM EXAMPLE_BLR_SEQUENTIEL
