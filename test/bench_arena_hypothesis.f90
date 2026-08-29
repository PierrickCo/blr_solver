!===============================================================================
! MICRO-BENCHMARK SYNTHETIQUE, PAS un test de correction (pas repris par
! run_all_tests.sh), PAS un changement au solveur reel -- ne touche AUCUN
! des 5 solveurs de reference ni blr_compress_mod.f90/mpi_blr_matrix_mod.f90.
!
! BUT : tester, a moindre risque, l'hypothese "un pool memoire (arena) par
! etape K reduirait la fragmentation observee sur example_blr_mpi_pivot"
! (voir la memoire du projet -- diagnostic memoire theorique-vs-RSS,
! mkl_free_buffers()/malloc_trim(0) tous deux ecartes comme cause
! principale) SANS toucher a BLR_BLOCK_T ni a la factorisation reelle --
! reproduit uniquement le PATTERN d'allocation suspect (des milliers de
! cycles allocate/deallocate de tailles VARIABLES, comme UFAC/VFAC d'un
! bloc BLR_BLOCK_T recompresse a chaque mise a jour de Schur), sous DEUX
! strategies comparables directement, MEME SEQUENCE de tailles (graine
! fixe) pour les deux :
!
!   BASELINE : un ALLOCATE/DEALLOCATE PAR bloc (ce que fait le solveur
!              reel aujourd'hui, via SCHUR_UPDATE_BLR_BLOCK/FAST_RECOMPRESS).
!   ARENA    : un SEUL grand ALLOCATE par etape K (taille = somme des
!              besoins de tous les blocs de l'etape, connue a l'avance ici
!              -- dans le vrai solveur il faudrait soit une passe de
!              dimensionnement prealable, soit une taille au pire cas),
!              chaque bloc "decoupe" une tranche de ce buffer plutot que
!              d'allouer le sien -- un seul DEALLOCATE en fin d'etape.
!
! Chaque etape imprime la RSS - la RSS de REFERENCE (avant toute boucle),
! meme methodologie que example_blr_mpi_pivot.f90 -- comparer directement
! les 2 courbes (deux executions, MODE=baseline puis MODE=arena) donne une
! reponse EMPIRIQUE avant d'investir dans un refactor du solveur reel.
!
! Usage : ./bench_arena_hypothesis MODE [NB] [NSTEPS] [BLOCKS_PER_STEP] [DENSE_PROB]
!   MODE            : "baseline" ou "arena"
!   NB              : taille de bloc simulee (defaut 500, comme les tests
!                     BLR+MPI a grande echelle de ce projet)
!   NSTEPS          : nombre d'etapes K simulees (defaut 40, ~NBLK a
!                     N=20000/NB=500)
!   BLOCKS_PER_STEP : nombre de blocs "recompresses" par etape (defaut 20)
!   DENSE_PROB      : probabilite qu'un bloc tombe en repli DENSE (rang=NB,
!                     gros bloc) plutot que bas-rang (rang petit, 5-50) --
!                     defaut 0.8, coherent avec le commentaire d'en-tete de
!                     example_blr_mpi_pivot.f90 ("la plupart des blocs
!                     hors-diagonale tombent en repli dense" pour cette
!                     matrice anti-diagonale/Cauchy peu compressible)
!===============================================================================
PROGRAM BENCH_ARENA_HYPOTHESIS
  USE BLR_KINDS_MOD, ONLY: WP
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE

  CHARACTER(LEN=16) :: MODE
  CHARACTER(LEN=32) :: ARG
  INTEGER :: NB, NSTEPS, BLOCKS_PER_STEP
  REAL :: DENSE_PROB
  INTEGER :: STEP, IBLK, TOTAL_ELEMS, OFFSET
  INTEGER, ALLOCATABLE :: RANKS(:)
  REAL :: U
  COMPLEX(WP), ALLOCATABLE :: UFAC(:,:), VFAC(:,:), ARENA(:)
  INTEGER(8) :: RSS_REF_KB, RSS_KB

  IF (COMMAND_ARGUMENT_COUNT() .LT. 1) THEN
    WRITE(6,*) 'Usage: ./bench_arena_hypothesis MODE [NB] [NSTEPS] [BLOCKS_PER_STEP] [DENSE_PROB]'
    WRITE(6,*) '  MODE = baseline ou arena'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, MODE)
  NB = 500; NSTEPS = 40; BLOCKS_PER_STEP = 20; DENSE_PROB = 0.8
  IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) NSTEPS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) BLOCKS_PER_STEP
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) DENSE_PROB
  END IF

  IF (MODE.NE.'baseline' .AND. MODE.NE.'arena') THEN
    WRITE(6,*) 'ERREUR : MODE doit etre "baseline" ou "arena"'
    STOP 1
  END IF

  ! Graine FIXE : les deux modes doivent voir EXACTEMENT la meme sequence
  ! de rangs, sinon la comparaison ne serait pas honnete.
  CALL SET_FIXED_SEED()

  WRITE(6,'(A,A,A,I0,A,I0,A,I0,A,F4.2)') ' MODE=', TRIM(MODE), ' NB=', NB, ' NSTEPS=', NSTEPS, &
       ' BLOCKS_PER_STEP=', BLOCKS_PER_STEP, ' DENSE_PROB=', DENSE_PROB

  ALLOCATE(RANKS(BLOCKS_PER_STEP))
  RSS_REF_KB = GET_RSS_KB()
  WRITE(6,'(A,F8.1,A)') '   RSS de reference (avant boucle) = ', REAL(RSS_REF_KB,8)/1024.0D0, ' Mo'

  DO STEP = 1, NSTEPS
    ! Meme tirage de rangs pour les 2 modes (RANDOM_SEED fixe + meme ordre
    ! d'appels RANDOM_NUMBER).
    DO IBLK = 1, BLOCKS_PER_STEP
      CALL RANDOM_NUMBER(U)
      IF (U .LT. DENSE_PROB) THEN
        RANKS(IBLK) = NB                    ! repli dense (gros bloc)
      ELSE
        CALL RANDOM_NUMBER(U)
        RANKS(IBLK) = 5 + INT(U * 45.0)      ! bas-rang, 5-50
      END IF
    END DO

    IF (MODE .EQ. 'baseline') THEN
      ! --- UN ALLOCATE/DEALLOCATE PAR BLOC (comportement actuel du solveur) ---
      DO IBLK = 1, BLOCKS_PER_STEP
        ALLOCATE(UFAC(NB, RANKS(IBLK)), VFAC(RANKS(IBLK), NB))
        CALL RANDOM_NUMBER_COMPLEX(UFAC)   ! force l'ecriture reelle des pages
        CALL RANDOM_NUMBER_COMPLEX(VFAC)
        DEALLOCATE(UFAC, VFAC)
      END DO
    ELSE
      ! --- UN SEUL ALLOCATE PAR ETAPE, blocs "decoupes" dans l'arena ---
      TOTAL_ELEMS = 0
      DO IBLK = 1, BLOCKS_PER_STEP
        TOTAL_ELEMS = TOTAL_ELEMS + 2 * NB * RANKS(IBLK)   ! UFAC + VFAC
      END DO
      ALLOCATE(ARENA(TOTAL_ELEMS))
      OFFSET = 0
      DO IBLK = 1, BLOCKS_PER_STEP
        CALL RANDOM_NUMBER_COMPLEX_FLAT(ARENA(OFFSET+1:OFFSET+NB*RANKS(IBLK)))
        OFFSET = OFFSET + NB*RANKS(IBLK)
        CALL RANDOM_NUMBER_COMPLEX_FLAT(ARENA(OFFSET+1:OFFSET+NB*RANKS(IBLK)))
        OFFSET = OFFSET + NB*RANKS(IBLK)
      END DO
      DEALLOCATE(ARENA)
    END IF

    RSS_KB = GET_RSS_KB()
    WRITE(6,'(A,I4,A,F8.1,A)') '   Etape ', STEP, ' : RSS - reference = ', &
         REAL(RSS_KB - RSS_REF_KB,8)/1024.0D0, ' Mo'
  END DO

CONTAINS

  SUBROUTINE SET_FIXED_SEED()
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: NSEED
    CALL RANDOM_SEED(SIZE=NSEED)
    ALLOCATE(SEED(NSEED))
    SEED = 42
    CALL RANDOM_SEED(PUT=SEED)
    DEALLOCATE(SEED)
  END SUBROUTINE SET_FIXED_SEED

  SUBROUTINE RANDOM_NUMBER_COMPLEX(A)
    COMPLEX(WP), INTENT(OUT) :: A(:,:)
    REAL(WP), ALLOCATABLE :: RP(:,:), IP(:,:)
    ALLOCATE(RP(SIZE(A,1),SIZE(A,2)), IP(SIZE(A,1),SIZE(A,2)))
    CALL RANDOM_NUMBER(RP)
    CALL RANDOM_NUMBER(IP)
    A = CMPLX(RP, IP, KIND=WP)
    DEALLOCATE(RP, IP)
  END SUBROUTINE RANDOM_NUMBER_COMPLEX

  SUBROUTINE RANDOM_NUMBER_COMPLEX_FLAT(A)
    COMPLEX(WP), INTENT(OUT) :: A(:)
    REAL(WP), ALLOCATABLE :: RP(:), IP(:)
    ALLOCATE(RP(SIZE(A)), IP(SIZE(A)))
    CALL RANDOM_NUMBER(RP)
    CALL RANDOM_NUMBER(IP)
    A = CMPLX(RP, IP, KIND=WP)
    DEALLOCATE(RP, IP)
  END SUBROUTINE RANDOM_NUMBER_COMPLEX_FLAT

END PROGRAM BENCH_ARENA_HYPOTHESIS
