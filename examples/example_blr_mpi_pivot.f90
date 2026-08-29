!===============================================================================
! EXEMPLE CANONIQUE -- reference BLR+MPI, AVEC PIVOT PARTIEL CONFINE AU
! BLOC DIAGONAL (mpi_blr_matrix_mod.f90, MPI_BLR_MATRIX_LU_FACTOR_PIVOT) :
! meme philosophie que MUMPS/STRUMPACK (pivotage local au front, jamais
! global) -- permuter les lignes d'un bloc compresse ne coute que
! O(rang) (seul le facteur UFAC est touche), jamais de decompression.
!
! Contrairement a la variante sans pivot (example_blr_mpi_sans_pivot.f90),
! celle-ci tolere des matrices SANS diagonale dominante -- la matrice de
! demonstration ci-dessous est volontairement "anti-diagonale forte /
! diagonale naturelle faible" pour forcer un vrai pivotage.
!
! But de ce fichier : montrer l'enchainement minimal (assembler,
! factoriser AVEC pivot, resoudre AVEC pivot, rassembler) sans toute
! l'instrumentation de mesure/validation croisee des drivers de test (voir
! test/test_mpi_blr_matrix_pivot.f90 pour la version qui verifie aussi la
! correction contre une reference independante).
!
! Usage : mpirun -np <PR*PC> ./example_blr_mpi_pivot PR PC [N] [NB] [NRHS] [TOL] [USE_MOVE_ALLOC]
!         (defauts : N=4000 NB=200 NRHS=4 TOL=1e-4 USE_MOVE_ALLOC=0)
!         USE_MOVE_ALLOC=1 : SCHUR_UPDATE_BLR_BLOCK transfere la propriete
!         des tableaux via MOVE_ALLOC au lieu de "ABLK=NEWBLK" (qui
!         declenche une reallocation+copie implicite supplementaire) --
!         resultat numerique IDENTIQUE, bonne pratique defensive
!         independante du gain memoire (voir la memoire du projet :
!         mesure a N=20000, l'effet sur la RSS s'est avere negligeable,
!         mais le correctif reste correct et sans risque, garde par
!         defaut desactive pour ne rien changer au comportement existant).
!===============================================================================
PROGRAM EXAMPLE_BLR_MPI_PIVOT
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  USE BENCH_UTILS_MOD, ONLY: GET_RSS_KB
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  INTEGER :: USE_MOVE_ALLOC_INT
  LOGICAL :: USE_MOVE_ALLOC_FLAG

  N = 4000; NB = 200; NRHS = 4; TOL = 1.0E-4_WP; USE_MOVE_ALLOC_INT = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'Usage: mpirun -np PR*PC ./example_blr_mpi_pivot PR PC [N] [NB] [NRHS] [TOL]'
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
    CALL GET_COMMAND_ARGUMENT(5, ARG); READ(ARG,*) NRHS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, ARG); READ(ARG,*) TOL
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 7) THEN
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) USE_MOVE_ALLOC_INT
  END IF
  USE_MOVE_ALLOC_FLAG = (USE_MOVE_ALLOC_INT .NE. 0)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,L1)') &
       ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, &
       ' , USE_MOVE_ALLOC=', USE_MOVE_ALLOC_FLAG

  ! RSS DE REFERENCE : capturee AVANT toute allocation de donnees (A_BLR,
  ! B_LOCAL...) -- isole le cout FIXE du runtime deja charge en memoire a
  ! ce stade (Intel MPI, MKL, le Fortran runtime lui-meme) de la
  ! croissance imputable a NOTRE code. Comparer (RSS_apres - RSS_reference)
  ! a la valeur theorique est un test plus honnete que RSS_apres seule,
  ! qui melange les deux.
  CALL REPORT_MEMORY_MB(GRID, 'RSS de reference (avant donnees)', GET_RSS_KB()*1024_8)

  ! --- 1) assemblage BLR distribue (matrice "anti-diagonale forte / diagonale
  !     naturelle faible" -- force un vrai pivotage, chaotique/peu compressible :
  !     la plupart des blocs hors-diagonale tombent en repli dense) ---
  ! SHOW_PROGRESS=.TRUE. : barre de progression affichee par le processus 0
  ! (voir mpi_blr_matrix_mod.f90 -- approximative ici, l'assemblage n'est
  ! pas synchronise entre processus, contrairement aux 2 etapes suivantes).
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE.)
  ! N_COMPRESSED/N_OFFDIAG sont deja des sommes GLOBALES (MPI_ALLREDUCE
  ! interne a MPI_ASSEMBLE_BLR_MATRIX) -- identiques sur tous les
  ! processus, un seul affichage (rang 0) suffit, pas une approximation.
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Taux de compression (assemblage) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF
  ! Empreinte THEORIQUE (blocs BLR_BLOCK_T reellement stockes), PAS
  ! mesuree via l'OS.
  CALL REPORT_MEMORY_MB(GRID, 'Matrice A_BLR (compressee)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
  ! RSS REELLE (mesuree, /proc/self/status) -- forcement >= la valeur
  ! theorique ci-dessus : elle inclut aussi tout le reste du processus
  ! (workspaces SVD/QR transitoires par thread, runtime MPI/MKL charge en
  ! memoire, l'executable lui-meme, l'allocateur qui ne rend jamais toute
  ! la memoire liberee a l'OS...). Le RATIO RSS/theorique n'est PAS
  ! constant : a petite echelle (N~2000), l'overhead FIXE (runtime
  ! MPI/MKL, ~50-70 Mo/processus rien qu'au demarrage) domine largement
  ! et donne un ratio enorme (>60x, mesure) ; a plus grande echelle
  ! (N~20000), les donnees elles-memes dominent et le ratio se stabilise
  ! autour de 2.5-3x (mesure) -- c'est ce second regime qui sert de
  ! reference pour detecter une VRAIE anomalie (fuite, copie oubliee).
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres assemblage)', GET_RSS_KB()*1024_8)

  ! --- 2) factorisation LU, AVEC pivot confine au bloc diagonal ----------
  ! SHOW_PROGRESS=.TRUE. : ici K/NBLK est un pourcentage EXACT (boucle
  ! globalement synchronisee par les MPI_BCAST a chaque etape K).
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                       N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                       SHOW_PROGRESS=.FALSE., SHOW_RSS_TRACE=.TRUE., &
                                       USE_MOVE_ALLOC=USE_MOVE_ALLOC_FLAG)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
    WRITE(6,'(A,I0)') ' Nombre de permutations de lignes (pivot) : ', N_SWAPS
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Taux de compression (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
  END IF
  ! A_BLR est factorisee EN PLACE -- cette valeur peut differer de celle
  ! de l'assemblage (fill-in de la mise a jour de Schur).
  CALL REPORT_MEMORY_MB(GRID, 'Facteurs LU (compresses)', MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres factorisation)', GET_RSS_KB()*1024_8)

  ! --- 3) second membre B, distribue par blocs-lignes ET par colonnes sur
  !     NRHS -- formule deterministe, sans communication ---
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_B(A_BLR, GRID, NRHS_LOC, RSTART_MINE, B_LOCAL, MLOC)
  CALL REPORT_MEMORY_MB(GRID, 'Second membre B', 8_8*INT(MLOC,8)*INT(NRHS_LOC,8))
  CALL REPORT_MEMORY_MB(GRID, 'RSS reelle (apres allocation B)', GET_RSS_KB()*1024_8)

  ! --- 4) resolution 2D, AVEC pivot (applique la meme sequence de swaps
  !     que celle utilisee a la factorisation, deja diffusee a la grille) --
  ! SHOW_PROGRESS=.TRUE. : une seule barre pour les 2 sous-etapes
  ! (substitution avant + arriere), pourcentage exact (meme raison que la
  ! factorisation -- boucle globalement synchronisee).
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  T0 = MPI_WTIME()
  CALL MPI_BLR_MATRIX_LU_SOLVE_PIVOT_2D(GRID, A_BLR, NRHS, NRHS_LOC, B_LOCAL, MAX(MLOC,1), SHOW_PROGRESS=.TRUE.)
  T1 = MPI_WTIME()
  T_SOLVE = T1 - T0
  IF (MYID0.EQ.0) WRITE(6,'(A,F8.3,A)') ' Resolution : ', T_SOLVE, ' s'

  ALLOCATE(X_GLOBAL(N,NRHS))
  CALL GATHER_2D_RHS(GRID, N, NB, NRHS, B_LOCAL, MAX(MLOC,1), NRHS_LOC, 0, X_GLOBAL)

  IF (MYID0 .EQ. 0) THEN
    WRITE(6,'(A)') ' Resolution terminee. Premieres valeurs de X(:,1) :'
    DO I = 1, MIN(5, N)
      WRITE(6,'(A,I0,A,ES12.4,A,ES12.4)') '   X(', I, ',1) = ', REAL(X_GLOBAL(I,1)), ' + i*', AIMAG(X_GLOBAL(I,1))
    END DO
  END IF

  DEALLOCATE(B_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! A(GI,GJ) = 50 si GJ=N+1-GI (anti-diagonale forte), sinon un noyau de
  ! Cauchy amorti 1/(1+|N+1-GI-GJ|) -- diagonale naturelle petite face a
  ! l'anti-diagonale : force un vrai pivotage (contrairement a la matrice
  ! diagonale dominante de la variante sans pivot).
  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    REAL(WP) :: DF
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        IF (JJ .EQ. N+1-II) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(50.0_WP, 0.0_WP, KIND=WP)
        ELSE
          DF = 1.0_WP / (1.0_WP + ABS(REAL(N+1-II-JJ,WP)))
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(DF, 0.3_WP*DF, KIND=WP)
        END IF
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  SUBROUTINE FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: BL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, GRR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          GRR = RSTART + RR - 1
          BL(RS + (GII-AB%BSTART(II)), RR) = CMPLX(REAL(GII,WP)*0.01_WP, -0.3_WP*REAL(GRR,WP), KIND=WP)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_B

END PROGRAM EXAMPLE_BLR_MPI_PIVOT
