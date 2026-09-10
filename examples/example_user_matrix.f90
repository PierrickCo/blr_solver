!===============================================================================
! EXEMPLE MINIMAL -- point de depart pour brancher VOTRE PROPRE systeme
! lineaire A*X=B sur le solveur BLR+MPI (tout en RAM). Version out-of-core
! (RAM bornee) : voir example_user_matrix_ooc.f90, meme structure.
!
! CE QU'IL FAUT REMPLACER (marque ci-dessous) :
!   - GET_BLOCK : renvoie un bloc de VOTRE matrice A, appelee PLUSIEURS FOIS
!     par la bibliotheque pour des sous-blocs -- jamais pour A entiere d'un
!     coup (c'est ce qui permet a l'assemblage de compresser au fur et a
!     mesure, sans jamais materialiser A dense en RAM).
!   - FILL_B : renvoie VOTRE second membre b.
! RIEN D'AUTRE n'a besoin d'etre modifie -- assemblage/factorisation/
! resolution sont deja generiques par rapport a ce que fait GET_BLOCK.
!
! Le placeholder ci-dessous (noyau 1/(1+|i-j|), diagonale = N) est juste un
! exemple RUNNABLE pour verifier que tout compile et tourne AVANT de
! brancher vos vraies donnees -- il n'a aucune signification physique.
!
! COMPRESSIBILITE : la compression BLR ne fait gagner du temps/memoire QUE
! si les blocs hors-diagonale de VOTRE matrice ont un rang numerique
! faible (matrices issues d'un noyau a decroissance geometrique, methode
! des elements de frontiere, etc.). Si ce n'est pas le cas, les blocs
! tombent simplement en repli DENSE (comportement SUR, jamais de
! plantage), sans gain -- ajustez RATIO/TOL en consequence si besoin, voir
! example_blr_mpi_pivot.f90 pour ces parametres avances.
!
! VERIFICATION SANS SOLUTION CONNUE (CHECK_RESIDUAL) : le residu
! ||A*X-B||/||B|| ne suppose RIEN sur B -- fonctionne avec votre b reel,
! contrairement a une comparaison a une solution exacte (qui n'existe pas
! pour un vrai probleme). Desactive par defaut (cout CALCUL O(N^2*NRHS)
! reel, meme optimise BLAS+MPI -- voir COMPUTE_RESIDUAL) : activez-le pour
! verifier, desactivez-le une fois confiant pour ne pas payer ce cout a
! chaque run.
!
! RATIO (COMPRESS_WORTHWHILE_RATIO) et USE_RSVD (SVD complete vs
! randomisee adaptative) : memes parametres, memes valeurs par defaut, que
! example_blr_mpi_pivot.f90 -- voir son commentaire de tete pour le detail.
!
! VERIFIER QUE VOS BLOCS SE COMPRESSENT VRAIMENT : GET_BLOCK renvoie
! toujours un bloc DENSE (c'est la seule chose qu'elle sait faire) -- c'est
! COMPRESS_BLOCK, appelee automatiquement par l'assemblage, qui decide
! ENSUITE de le stocker compresse (SVD tronquee) ou de le garder dense
! (bloc juge "pas assez compressible", voir RATIO ci-dessus). Deux
! indicateurs affiches ci-dessous permettent de VERIFIER ce qui s'est
! reellement passe sur VOS blocs :
!   - "Blocs compresses" : COMBIEN de blocs ont ete stockes compresses
!     (IS_LOWRANK=.TRUE.) sur le total hors-diagonale -- une simple
!     proportion, qui peut rester a 100% meme si chaque bloc n'economise
!     presque rien (un rang proche du plafond RATIO*NB compresse encore
!     "techniquement" mais n'apporte quasi aucun gain memoire).
!   - "Taux de compression REEL" : la taille de stockage EFFECTIVE
!     (bas-rang + dense restant) rapportee a la taille dense equivalente,
!     basee sur le RANG MOYEN reellement trouve -- LA mesure qui repond
!     vraiment a "mes blocs entres en dense sont-ils bien compresses
!     ensuite, et de combien ?". Proche de 100% = quasi aucun gain reel
!     (meme si "Blocs compresses" affiche 100%) ; proche de 0% = tres
!     bonne compression.
!
! Usage : mpirun -np <PR*PC> ./example_user_matrix PR PC [N] [NB] [NRHS] [TOL] [RATIO] [USE_RSVD] [CHECK_RESIDUAL]
!         (defauts : N=4000 NB=200 NRHS=1 TOL=1e-4 RATIO=1/3 USE_RSVD=0 CHECK_RESIDUAL=0)
!===============================================================================
PROGRAM EXAMPLE_USER_MATRIX
  USE MPI_BLR_MATRIX_MOD
  USE MPI_LU_SOLVE_2D_MOD, ONLY: GATHER_2D_RHS
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(GRID_T) :: GRID
  INTEGER :: PR, PC, N, NB, NRHS, IERR, MYID0
  REAL(WP) :: TOL, RATIO
  CHARACTER(LEN=32) :: ARG
  TYPE(MPI_BLR_MATRIX_T) :: A_BLR
  COMPLEX(WP), ALLOCATABLE :: B_LOCAL(:,:), B_ORIG_LOCAL(:,:), X_GLOBAL(:,:)
  INTEGER :: MLOC, NRHS_LOC, RSTART_MINE, INFO, N_SWAPS
  INTEGER :: N_COMPRESSED, N_OFFDIAG, li, I
  INTEGER(8) :: RANK_SUM
  DOUBLE PRECISION :: T_ASM, T_FACTOR, T_SOLVE, T0, T1
  REAL(WP) :: RESID
  INTEGER :: USE_RSVD_INT, CHECK_RESIDUAL_INT
  LOGICAL :: USE_RSVD_FLAG, CHECK_RESIDUAL_FLAG

  N = 4000; NB = 200; NRHS = 1; TOL = 1.0E-4_WP
  RATIO = 1.0_WP/3.0_WP; USE_RSVD_INT = 0; CHECK_RESIDUAL_INT = 0

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID0, IERR)

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    IF (MYID0.EQ.0) WRITE(6,*) &
         'Usage: mpirun -np PR*PC ./example_user_matrix PR PC [N] [NB] [NRHS] [TOL] '// &
         '[RATIO] [USE_RSVD] [CHECK_RESIDUAL]'
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
    CALL GET_COMMAND_ARGUMENT(7, ARG); READ(ARG,*) RATIO
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 8) THEN
    CALL GET_COMMAND_ARGUMENT(8, ARG); READ(ARG,*) USE_RSVD_INT
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 9) THEN
    CALL GET_COMMAND_ARGUMENT(9, ARG); READ(ARG,*) CHECK_RESIDUAL_INT
  END IF
  USE_RSVD_FLAG = (USE_RSVD_INT .NE. 0)
  CHECK_RESIDUAL_FLAG = (CHECK_RESIDUAL_INT .NE. 0)

  ! Voir example_blr_mpi_pivot.f90/SET_COMPRESS_WORTHWHILE_RATIO pour le
  ! detail -- AVANT tout assemblage/factorisation.
  CALL SET_COMPRESS_WORTHWHILE_RATIO(RATIO)

  CALL GRID_INIT(GRID, MPI_COMM_WORLD, PR, PC, IERR)
  IF (IERR .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR : PR*PC doit egaler le nombre de processus lances'
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F6.3)') &
         ' Grille ', PR, 'x', PC, ' , N=', N, ' , NB=', NB, ' , NRHS=', NRHS, ' , TOL=', TOL, ' , RATIO=', RATIO
    IF (USE_RSVD_FLAG) THEN
      WRITE(6,'(A)') ' Compression : SVD randomisee adaptative (COMPRESS_BLOCK_RSVD_ADAPTIVE)'
    ELSE
      WRITE(6,'(A)') ' Compression : SVD complete (COMPRESS_BLOCK, CGESVD)'
    END IF
  END IF

  ! --- 1) assemblage BLR distribue -----------------------------------------
  CALL MPI_ASSEMBLE_BLR_MATRIX(GRID, N, NB, TOL, GET_BLOCK, A_BLR, T_ASM, &
                                N_COMPRESSED, N_OFFDIAG, RANK_SUM, SHOW_PROGRESS=.TRUE., USE_RSVD=USE_RSVD_FLAG)
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Assemblage+compression : ', T_ASM, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (bas-rang / dense equivalent) : ', &
         REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
  END IF

  ! --- 2) factorisation LU, pivot confine au bloc diagonal -----------------
  CALL MPI_BLR_MATRIX_LU_FACTOR_PIVOT(GRID, A_BLR, TOL, INFO, &
                                       N_COMPRESSED, N_OFFDIAG, RANK_SUM, N_SWAPS, T_FACTOR, &
                                       SHOW_PROGRESS=.TRUE.)
  IF (INFO .NE. 0) THEN
    IF (MYID0.EQ.0) WRITE(6,*) 'ERREUR MPI_BLR_MATRIX_LU_FACTOR_PIVOT, INFO=', INFO
    CALL MPI_FINALIZE(IERR)
    STOP 1
  END IF
  IF (MYID0.EQ.0) THEN
    WRITE(6,'(A,F8.3,A)') ' Factorisation LU : ', T_FACTOR, ' s'
    WRITE(6,'(A,I0,A,I0,A,F6.1,A)') ' Blocs compresses (fin factorisation) : ', N_COMPRESSED, ' / ', N_OFFDIAG, &
         ' (', 100.0_WP*REAL(N_COMPRESSED,WP)/REAL(MAX(N_OFFDIAG,1),WP), ' %)'
    WRITE(6,'(A,F6.2,A)') ' Taux de compression REEL (fin factorisation) : ', &
         REAL_COMPRESSION_PCT(N_COMPRESSED, N_OFFDIAG, RANK_SUM), ' %'
  END IF

  ! --- 3) second membre local, distribue comme A ---------------------------
  MLOC = 0
  DO li = 1, A_BLR%NBLKR_LOC
    I = (li-1)*GRID%PR + GRID%MYROW + 1
    MLOC = MLOC + (A_BLR%BEND(I)-A_BLR%BSTART(I)+1)
  END DO
  NRHS_LOC = RHS_NUMROC_MBM(NRHS, GRID%PC, GRID%MYCOL)
  CALL RHS_LOCAL_RANGE_MBM(NRHS, GRID%PC, GRID%MYCOL, RSTART_MINE, NRHS_LOC)

  ALLOCATE(B_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)), B_ORIG_LOCAL(MAX(MLOC,1),MAX(NRHS_LOC,1)))
  CALL FILL_B(A_BLR, GRID, NRHS_LOC, RSTART_MINE, B_LOCAL, MLOC)
  B_ORIG_LOCAL = B_LOCAL   ! B_LOCAL sera ecrase EN PLACE par la resolution -- garder une copie pour le residu

  ! --- 4) resolution 2D, avec pivot -----------------------------------------
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

  ! --- 5) verification ||A*X-B||/||B|| -- APPEL COLLECTIF (tous les rangs) -
  ! Optionnelle (CHECK_RESIDUAL=0 par defaut) : cout CALCUL reel meme
  ! optimise, voir le commentaire de tete.
  IF (CHECK_RESIDUAL_FLAG) THEN
    RESID = COMPUTE_RESIDUAL(X_GLOBAL, B_ORIG_LOCAL, GRID, N, NB, NRHS, NRHS_LOC, MLOC)
    IF (MYID0.EQ.0) WRITE(6,'(A,ES12.4)') ' Residu relatif ||A*X-B|| / ||B|| = ', RESID
  END IF

  DEALLOCATE(B_LOCAL, B_ORIG_LOCAL, X_GLOBAL)
  CALL MPI_FINALIZE(IERR)

CONTAINS

  ! Taux de compression REEL (stockage bas-rang effectif / stockage dense
  ! equivalent) -- PAS un simple pourcentage de blocs "compresses" (voir
  ! le commentaire de tete) : RANK_SUM mele le rang REEL des blocs
  ! compresses et un "rang virtuel"=NB pour chaque bloc en repli DENSE,
  ! d'ou le retrait de la contribution des blocs denses avant de calculer
  ! le ratio (meme correctif que example_blr_mpi_pivot.f90).
  REAL(WP) FUNCTION REAL_COMPRESSION_PCT(NCOMP, NOFF, RSUM)
    INTEGER, INTENT(IN)     :: NCOMP, NOFF
    INTEGER(8), INTENT(IN)  :: RSUM
    INTEGER :: NDENSE
    NDENSE = NOFF - NCOMP
    REAL_COMPRESSION_PCT = 100.0_WP * (2.0_WP*REAL(NB,WP)*REAL(RSUM,WP) - REAL(NB,WP)**2*REAL(NDENSE,WP)) &
         / (REAL(NB,WP)**2 * REAL(MAX(NOFF,1),WP))
  END FUNCTION REAL_COMPRESSION_PCT

  !===========================================================================
  ! REMPLACEZ CETTE FONCTION par votre propre matrice A.
  ! BLK(i,j) doit contenir A(ISTART+i-1, JSTART+j-1) -- appelee par blocs,
  ! jamais pour la matrice entiere. Le placeholder ci-dessous (noyau
  ! 1/(1+|i-j|), diagonale=N pour la dominance) n'est qu'un exemple
  ! RUNNABLE, sans signification physique.
  !===========================================================================
  SUBROUTINE GET_BLOCK(ISTART, IEND, JSTART, JEND, BLK)
    INTEGER, INTENT(IN)      :: ISTART, IEND, JSTART, JEND
    COMPLEX(WP), INTENT(OUT) :: BLK(IEND-ISTART+1, JEND-JSTART+1)
    INTEGER :: II, JJ
    DO JJ = JSTART, JEND
      DO II = ISTART, IEND
        IF (II .EQ. JJ) THEN
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(REAL(N,WP), 0.0_WP, KIND=WP)
        ELSE
          BLK(II-ISTART+1, JJ-JSTART+1) = CMPLX(1.0_WP/(1.0_WP+ABS(REAL(II-JJ,WP))), 0.0_WP, KIND=WP)
        END IF
      END DO
    END DO
  END SUBROUTINE GET_BLOCK

  !===========================================================================
  ! REMPLACEZ CETTE FONCTION par votre propre second membre b.
  ! Remplit les lignes GLOBALES possedees LOCALEMENT par ce processus (deja
  ! calculees via A_BLR%BSTART/BEND -- ne pas y toucher) pour les colonnes
  ! RSTART..RSTART+NRL-1. Placeholder : b = vecteur de 1.
  !===========================================================================
  SUBROUTINE FILL_B(AB, GD, NRL, RSTART, BL, ML)
    TYPE(MPI_BLR_MATRIX_T), INTENT(IN) :: AB
    TYPE(GRID_T), INTENT(IN)           :: GD
    INTEGER, INTENT(IN)                :: NRL, RSTART, ML
    COMPLEX(WP), INTENT(OUT)           :: BL(MAX(ML,1),MAX(NRL,1))
    INTEGER :: lli, II, GII, RR, RS
    RS = 1
    DO lli = 1, AB%NBLKR_LOC
      II = (lli-1)*GD%PR + GD%MYROW + 1
      DO GII = AB%BSTART(II), AB%BEND(II)
        DO RR = 1, NRL
          BL(RS + (GII-AB%BSTART(II)), RR) = CMPLX(1.0_WP, 0.0_WP, KIND=WP)
        END DO
      END DO
      RS = RS + (AB%BEND(II)-AB%BSTART(II)+1)
    END DO
  END SUBROUTINE FILL_B

  !===========================================================================
  ! ||A*X-B||/||B|| -- APPEL COLLECTIF (GATHER_2D_RHS/MPI_BCAST/
  ! MPI_ALLREDUCE a l'interieur, tous les rangs doivent l'appeler). A*X est
  ! calcule par blocs de lignes via GET_BLOCK + CGEMM (BLAS), lignes
  ! reparties entre TOUS les rangs. B_ORIG_LOCAL (VOTRE b reel, sauvegarde
  ! AVANT que la resolution ne l'ecrase en place) est rassemble en un
  ! tableau complet UNE SEULE FOIS ici -- jamais recalcule via une formule :
  ! ca marche donc que B vienne d'une formule, d'un fichier, ou de n'importe
  ! quelle autre source. Voir example_blr_mpi_pivot.f90/COMPUTE_RESIDUAL
  ! pour la version complete (avec barre de progression/RSS) dont celle-ci
  ! est la forme condensee.
  !===========================================================================
  REAL(WP) FUNCTION COMPUTE_RESIDUAL(XG, BORIG_LOCAL, GD, NN, NNB, NR, NRL, MLL)
    TYPE(GRID_T), INTENT(IN)   :: GD
    INTEGER, INTENT(IN)        :: NN, NNB, NR, NRL, MLL
    COMPLEX(WP), INTENT(INOUT) :: XG(NN,NR)
    COMPLEX(WP), INTENT(IN)    :: BORIG_LOCAL(MAX(MLL,1),MAX(NRL,1))
    COMPLEX(WP), ALLOCATABLE :: ABLOCK(:,:), AXBLOCK(:,:), BG(:,:)
    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    INTEGER :: ROWS_PER_PROC, RSTART_MINE, REND_MINE, RSTART_CHUNK, REND_CHUNK, NCHUNK
    INTEGER :: GI, R, IERR2
    REAL(WP) :: NUM2_LOCAL, DEN2_LOCAL, NUM2, DEN2

    ALLOCATE(BG(NN,NR))
    CALL GATHER_2D_RHS(GD, NN, NNB, NR, BORIG_LOCAL, MAX(MLL,1), NRL, 0, BG)
    CALL MPI_BCAST(XG, NN*NR, MPI_COMPLEX, 0, GD%COMM, IERR2)
    CALL MPI_BCAST(BG, NN*NR, MPI_COMPLEX, 0, GD%COMM, IERR2)

    ROWS_PER_PROC = (NN + GD%NPROCS - 1) / GD%NPROCS
    RSTART_MINE = GD%MYID*ROWS_PER_PROC + 1
    REND_MINE = MIN((GD%MYID+1)*ROWS_PER_PROC, NN)

    NUM2_LOCAL = 0.0_WP; DEN2_LOCAL = 0.0_WP
    DO RSTART_CHUNK = RSTART_MINE, REND_MINE, NNB
      REND_CHUNK = MIN(RSTART_CHUNK+NNB-1, REND_MINE)
      NCHUNK = REND_CHUNK - RSTART_CHUNK + 1
      ALLOCATE(ABLOCK(NCHUNK,NN), AXBLOCK(NCHUNK,NR))
      CALL GET_BLOCK(RSTART_CHUNK, REND_CHUNK, 1, NN, ABLOCK)
      CALL CGEMM('N','N', NCHUNK, NR, NN, CONE, ABLOCK, NCHUNK, XG, NN, CZERO, AXBLOCK, NCHUNK)
      DO R = 1, NR
        DO GI = 1, NCHUNK
          NUM2_LOCAL = NUM2_LOCAL + REAL(ABS(AXBLOCK(GI,R)-BG(RSTART_CHUNK+GI-1,R))**2,WP)
          DEN2_LOCAL = DEN2_LOCAL + REAL(ABS(BG(RSTART_CHUNK+GI-1,R))**2,WP)
        END DO
      END DO
      DEALLOCATE(ABLOCK, AXBLOCK)
    END DO
    DEALLOCATE(BG)

    CALL MPI_ALLREDUCE(NUM2_LOCAL, NUM2, 1, MPI_REAL, MPI_SUM, GD%COMM, IERR2)
    CALL MPI_ALLREDUCE(DEN2_LOCAL, DEN2, 1, MPI_REAL, MPI_SUM, GD%COMM, IERR2)
    COMPUTE_RESIDUAL = SQRT(NUM2) / SQRT(DEN2)
  END FUNCTION COMPUTE_RESIDUAL

END PROGRAM EXAMPLE_USER_MATRIX
