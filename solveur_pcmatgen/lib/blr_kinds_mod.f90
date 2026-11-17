!===============================================================================
! Parametre de precision partage par l'ensemble du solveur. Doit rester le
! SEUL endroit ou WP est defini : deux modules redefinissant chacun leur
! propre WP, meme a la meme valeur, produisent une reference ambigue a la
! compilation des qu'ils sont USE-associes ensemble.
!
! Contient aussi les SEULS appels a l'API OpenMP de la bibliotheque
! (horloge, nombre de threads, numero de thread, coeurs visibles), en
! compilation conditionnelle OpenMP (lignes "!$ ", compilees seulement avec
! -qopenmp/-fopenmp) : la bibliotheque compile et tourne SANS OpenMP
! (Makefile OPENMP=0, MKL sequentiel, pas de libiomp5 -- voir CLUSTER.md :
! avec ifx, libiomp5 fait passer tous les ALLOCATE par l'allocateur du
! runtime OpenMP d'Intel, qui ne rend jamais la memoire liberee).
!===============================================================================
MODULE BLR_KINDS_MOD
!$ USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM, OMP_GET_NUM_PROCS
  IMPLICIT NONE
  INTEGER, PARAMETER :: WP = KIND(1.0)   ! simple precision (reel et complexe)

CONTAINS

  ! Horloge murale en secondes (remplace OMP_GET_WTIME).
  DOUBLE PRECISION FUNCTION BLR_WALLTIME()
    INTEGER(8) :: COUNT, RATE
    CALL SYSTEM_CLOCK(COUNT, RATE)
    BLR_WALLTIME = DBLE(COUNT) / DBLE(MAX(RATE, 1_8))
  END FUNCTION BLR_WALLTIME

  ! Nombre maximal de threads OpenMP (1 sans OpenMP).
  INTEGER FUNCTION BLR_MAX_THREADS()
    BLR_MAX_THREADS = 1
!$  BLR_MAX_THREADS = OMP_GET_MAX_THREADS()
  END FUNCTION BLR_MAX_THREADS

  ! Indice 1-base du thread courant (1 sans OpenMP).
  INTEGER FUNCTION BLR_THREAD_INDEX()
    BLR_THREAD_INDEX = 1
!$  BLR_THREAD_INDEX = OMP_GET_THREAD_NUM() + 1
  END FUNCTION BLR_THREAD_INDEX

  ! Coeurs visibles par le processus selon OpenMP (-1 sans OpenMP).
  INTEGER FUNCTION BLR_NUM_PROCS()
    BLR_NUM_PROCS = -1
!$  BLR_NUM_PROCS = OMP_GET_NUM_PROCS()
  END FUNCTION BLR_NUM_PROCS

END MODULE BLR_KINDS_MOD
