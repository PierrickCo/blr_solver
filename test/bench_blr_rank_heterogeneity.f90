!===============================================================================
! DIAGNOSTIC/EXPERIMENTATION -- matrice de test DELIBEREMENT HETEROGENE, pour
! valider (ou infirmer) l'idee d'une permutation des blocs BLR guidee par le
! rang mesure a l'assemblage (voir la discussion sur example_blr_mpi_pivot.f90/
! examples/analyze_block_ranks.py). Sur la matrice de demonstration canonique
! (noyau de Cauchy amorti), le diagnostic a montre que le rang d'un bloc ne
! depend QUE de la distance |I-J| au bloc diagonal -- deja parfaitement
! exploite par l'ordre naturel, rien a gagner d'une permutation (coefficient
! de variation du rang moyen entre blocs-lignes ~0.02-0.04, quasi nul).
!
! ICI, au contraire, la "durete" d'un bloc ne depend PAS de sa distance a la
! diagonale mais d'un groupe CACHE assigne a chaque indice de bloc par un hash
! deterministe SANS RAPPORT avec l'ordre naturel (CLUSTER_OF ci-dessous) --
! deux blocs-lignes voisins dans l'ordre naturel appartiennent typiquement a
! des groupes DIFFERENTS, et deux blocs du MEME groupe sont typiquement
! disperses loin l'un de l'autre. Formule :
!   - blocs (I,J) du MEME groupe cache -> "durs" : entrees pseudo-aleatoires
!     (HASH_NOISE) SANS structure de decroissance -- proches du rang plein,
!     quasi incompressibles.
!   - blocs (I,J) de groupes DIFFERENTS -> "faciles" : produit exterieur
!     U(I)*V(J) (rang EXACTEMENT 1 par construction), trivialement
!     compressibles quel que soit TOL.
! Un ordre naturel NE PEUT PAS regrouper les blocs du meme groupe (ils sont
! scattered par construction) -- c'est exactement le cas ou une permutation
! qui reordonne les blocs par groupe (mesure via le rang, sans connaitre le
! groupe a priori) devrait, si l'idee est valide, regrouper les blocs "durs"
! entre eux et reveler une forte heterogeneite du rang moyen PAR BLOC-LIGNE
! (contrairement au cas Cauchy amorti) -- verifiable directement avec
! examples/analyze_block_ranks.py sur le CSV produit ici.
!
! Ce fichier ne fait QUE l'assemblage + l'export CSV (DUMP_STRUCTURE) --
! PAS de factorisation/resolution, inutile pour ce diagnostic (voir le
! commentaire de tete d'example_blr_mpi_pivot.f90 pour le format CSV, repris
! ICI a l'identique : block_row,block_col,kind,rank).
!
! Usage : mpirun -np <PR*PC> ./bench_blr_rank_heterogeneity PR PC N NB TOL NCLUSTERS
!         (defauts : N=4000 NB=200 TOL=1e-4 NCLUSTERS=4)
!===============================================================================
PROGRAM BENCH_BLR_RANK_HETEROGENEITY
  USE MPI_BLR_MATRIX_MOD
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, IERR, MYID0, NCLUSTERS
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  INTEGER :: N_COMPRESSED, N_OFFDIAG
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM

  N = 4000; NB = 200; TOL = 1.0E-4_WP; NCLUSTERS = 4

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./bench_blr_rank_heterogeneity PR PC [N] [NB] [TOL] [NCLUSTERS]'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) PR
  CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) PC
  IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) N
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, ARG); READ(ARG,*) NB
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) NCLUSTERS
  END IF

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,I0)') &
         ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , TOL=', TOL, ' , NCLUSTERS=', NCLUSTERS
  END IF

  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE.)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF

  CALL DUMP_BLR_STRUCTURE('assembly')

  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Groupe cache assigne a un indice de bloc GLOBAL (1..NBLK), via un hash
  ! deterministe SANS RAPPORT avec l'ordre naturel -- constantes DIFFERENTES
  ! de HASH_NOISE ci-dessous pour eviter toute correlation involontaire.
  INTEGER FUNCTION CLUSTER_OF(IBLK)
    INTEGER, INTENT(IN) :: IBLK
    REAL(WP) :: S, FRAC
    S = SIN(REAL(IBLK,WP)*37.719_WP)*29143.671_WP
    FRAC = S - REAL(FLOOR(S),WP)
    CLUSTER_OF = MIN(INT(FRAC*REAL(NCLUSTERS,WP)), NCLUSTERS-1)
  END FUNCTION CLUSTER_OF

  ! Hash deterministe (I,J) -> [-1,1) -- meme recette qu'example_blr_mpi_
  ! pivot.f90/HASH_NOISE.
  REAL(WP) FUNCTION HASH_NOISE(GI, GJ)
    INTEGER, INTENT(IN) :: GI, GJ
    REAL(WP) :: S
    S = SIN(REAL(GI,WP)*12.9898_WP + REAL(GJ,WP)*78.233_WP) * 43758.5453_WP
    HASH_NOISE = 2.0_WP*(S - REAL(FLOOR(S),WP)) - 1.0_WP
  END FUNCTION HASH_NOISE

  ! Facteurs U/V du produit exterieur (rang EXACTEMENT 1) utilise entre
  ! blocs de groupes DIFFERENTS -- fonctions lisses, sans rapport avec
  ! HASH_NOISE/CLUSTER_OF.
  COMPLEX(WP) FUNCTION U_FUNC(GI)
    INTEGER, INTENT(IN) :: GI
    U_FUNC = CMPLX(SIN(REAL(GI,WP)*0.013_WP)+2.0_WP, COS(REAL(GI,WP)*0.007_WP), KIND=WP)
  END FUNCTION U_FUNC

  COMPLEX(WP) FUNCTION V_FUNC(GJ)
    INTEGER, INTENT(IN) :: GJ
    V_FUNC = CMPLX(COS(REAL(GJ,WP)*0.011_WP)+2.0_WP, SIN(REAL(GJ,WP)*0.005_WP), KIND=WP)
  END FUNCTION V_FUNC

  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ, IBLK, JBLK, CI, CJ
    DO JJ = JSTART, JEND
      JBLK = (JJ-1)/NB + 1
      CJ = CLUSTER_OF(JBLK)
      DO II = ISTART, IEND
        IBLK = (II-1)/NB + 1
        CI = CLUSTER_OF(IBLK)
        IF (CI .EQ. CJ) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(1.0_WP + 0.7_WP*HASH_NOISE(II,JJ), &
                                                 0.4_WP*HASH_NOISE(JJ,II), KIND=WP)
        ELSE
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(0.01_WP,0.0_WP,KIND=WP) * U_FUNC(II) * V_FUNC(JJ)
        END IF
        IF (II .EQ. JJ) BLK(II-ISTART+1, JJ-JSTART+1) = BLK(II-ISTART+1, JJ-JSTART+1) + CMPLX(50.0_WP, 0.0_WP, KIND=WP)
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  ! Meme format CSV qu'example_blr_mpi_pivot.f90/DUMP_BLR_STRUCTURE
  ! (block_row,block_col,kind,rank) -- reutilisable tel quel par
  ! examples/analyze_block_ranks.py ou plot_blr_structure.py.
  SUBROUTINE DUMP_BLR_STRUCTURE(STAGE_LABEL)
    CHARACTER(LEN=*), INTENT(IN) :: STAGE_LABEL
    INTEGER :: U, li, lj, I, J
    CHARACTER(LEN=96) :: FNAME
    WRITE(FNAME,'(A,A,A,I0,A)') 'blr_structure_', TRIM(STAGE_LABEL), '_rank', MYID0, '.csv'
    OPEN(NEWUNIT=U, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')
    WRITE(U,'(A)') 'block_row,block_col,kind,rank'
    DO I = 1, A_BLR%NBLK
      IF (ALLOCATED(A_BLR%DIAGBLK(I)%D)) THEN
        WRITE(U,'(I0,A,I0,A,I0)') I, ',', I, ',diagonal,', A_BLR%BEND(I)-A_BLR%BSTART(I)+1
      END IF
    END DO
    DO li = 1, A_BLR%NBLKR_LOC
      I = (li-1)*GRID%PR + GRID%MYROW + 1
      DO lj = 1, A_BLR%NBLKC_LOC
        J = (lj-1)*GRID%PC + GRID%MYCOL + 1
        IF (I .EQ. J) CYCLE
        IF (A_BLR%OFFBLK(li,lj)%IS_LOWRANK) THEN
          WRITE(U,'(I0,A,I0,A,I0)') I, ',', J, ',lowrank,', A_BLR%OFFBLK(li,lj)%RANK
        ELSE
          WRITE(U,'(I0,A,I0,A,I0)') I, ',', J, ',dense,', A_BLR%BEND(J)-A_BLR%BSTART(J)+1
        END IF
      END DO
    END DO
    CLOSE(U)
  END SUBROUTINE DUMP_BLR_STRUCTURE

END PROGRAM BENCH_BLR_RANK_HETEROGENEITY
