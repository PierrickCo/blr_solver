!===============================================================================
! Diffuse les panneaux L21/U12 deja compresses au lieu de diffuser dense
! puis compresser apres reception. Deux gains :
!   1) VOLUME DE COMMUNICATION : si le panneau compresse (rang R << KW),
!      on envoie UFAC (MxR) + VFAC (RxN) au lieu du bloc dense MxN -- un
!      gain propre au distribue, sans equivalent sequentiel.
!   2) CALCUL REDONDANT EVITE : diffuser un panneau dense fait recompresser
!      independamment (CGESVD) la meme donnee par chacun des P processus
!      recepteurs -- P appels CGESVD redondants sur des donnees
!      identiques. Compresser une fois a la source, avant diffusion,
!      elimine cette redondance completement, que le bloc soit
!      compressible ou non.
!
! MECANISME (BCAST_BLR_BLOCK, ci-dessous) : la source (qui a deja le bloc
! dense et vient de le compresser) diffuse d'abord un petit entete (2
! entiers : is_lowrank, rank) pour que les recepteurs sachent quoi
! allouer/recevoir (la taille du buffer n'est pas fixe d'avance). Puis
! diffusion des donnees elles-memes (UFAC+VFAC si basse-rang, DENSE
! sinon).
!
! PROPRIETE NUMERIQUE A CONNAITRE (pas un bug) : le seuil "rentable" de
! COMPRESS_BLOCK (RANG < 0.8*MIN(M,N)) est calibre pour un compromis en
! FLOPS (le GEMM basse-rang devient moins cher que le GEMM dense), pas
! pour un compromis en octets. Pour un bloc presque carre (M proche de N),
! le seuil ou la compression economise reellement de la memoire/bande
! passante est plus strict (RANG < M*N/(M+N), environ MIN(M,N)/2 pour un
! bloc carre) que le seuil "rentable pour les flops" -- un bloc peut donc
! etre marque "compresse" (et accelerer le GEMM) sans pour autant reduire
! le volume envoye ici, voire l'augmenter legerement (UFAC+VFAC > DENSE).
! Mesure sur N=60,NB=5 (blocs quasi carres en fin de matrice) : 101.2% du
! volume dense au lieu d'une reduction. Sans consequence pour un usage BLR
! realiste (NB << N, panneaux L21/U12 largement rectangulaires la plupart
! du temps) : N=120-200, NB=8-10 donne 66.9%-75.0% du volume dense, une
! reduction reelle et systematique. A garder en tete si NB descend
! en-dessous du rang typique des blocs.
!
! Structure generale identique aux autres variantes MPI+pivot de ce
! solveur (rassembler-factoriser-renvoyer le panneau pivote, propager les
! echanges de lignes), le calcul dense de U12 par CTRSM (petit bloc DKK,
! toujours diffuse dense -- KW<=NB, jamais interessant a compresser), et
! le dispatch SCHUR_UPDATE_BLOCK final. Seule la diffusion de LBUF/UBUF
! change de "dense puis compresser cote recepteur" a "compresser cote
! source puis diffuser compresse".
!===============================================================================
MODULE MPI_BLR_LU_PIVOT_BCAST_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  USE BLR_COMPRESS_MOD
  USE MPI_LU_PIVOT_MOD, ONLY: SWAP_ROWS_RANGE
  IMPLICIT NONE

  ! Compteurs des diffusions de blocs (BCAST_BLR_BLOCK) vues par ce
  ! processus, emetteur ou recepteur : appels MPI_BCAST et elements
  ! complexes transportes (diagnostic).
  INTEGER(8) :: BCAST_NCALL = 0_8, BCAST_ELEMS = 0_8

  ! Diffusion groupee des panneaux (BCAST_BLR_PANEL) : -1 = lire la
  ! variable d'environnement BLR_PANEL_PACK au premier appel, 1 = groupee
  ! (defaut), 0 = bloc par bloc.
  INTEGER :: PANEL_PACK_MODE = -1

CONTAINS

  !-----------------------------------------------------------------------
  ! Diffuse un BLR_BLOCK_T (M x N, deja pret sur le processus ROOT -- rang
  ! DANS COMM, pas rang global) depuis ROOT vers tous les membres de COMM.
  ! ROOT doit avoir BLK deja rempli (dense ou basse-rang) ; les AUTRES
  ! membres recoivent un BLK entierement neuf (alloue ici).
  !-----------------------------------------------------------------------
  ! CONCATENATION UFAC+VFAC (un appel MPI_BCAST de moins) : au lieu de
  ! diffuser UFAC (MxR) puis VFAC (RxN) separement (2 appels), les deux
  ! sont empiles dans un buffer 1D contigu BUF (taille M*R+R*N), diffuse en
  ! un seul appel -- 2 appels au lieu de 3 par panneau compresse (header +
  ! BUF, au lieu de header + UFAC + VFAC), donc environ un tiers de moins
  ! de rendez-vous collectifs sur le chemin dominant. RESHAPE cree une
  ! copie temporaire (empilage/depilage) -- cout negligeable (M*R+R*N,
  ! petit vu que R<<NB) compare a l'aller-retour MPI economise.
  ! T_ALLOC/T_MPICALL (optionnels, INTENT(INOUT), cumules par l'appelant --
  ! initialiser a 0 avant la premiere boucle) : mesurent separement (a) le
  ! cout du DEALLOCATE/ALLOCATE de UFAC/VFAC/DENSE cote recepteur (voir
  ! commentaire ci-dessous -- ce cout existe car BLK est INTENT(INOUT), pas
  ! INTENT(OUT)) et (b) le temps passe dans les appels MPI_BCAST eux-memes.
  ! Non fournis : comportement et cout inchanges.
  SUBROUTINE BCAST_BLR_BLOCK(COMM, ROOT, M, N, BLK, T_ALLOC, T_MPICALL)
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN)             :: COMM, ROOT, M, N
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLK
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: T_ALLOC, T_MPICALL

    INTEGER :: HEADER(2), IERR, MYID, NU, NV, NTOT
    COMPLEX(WP), ALLOCATABLE :: BUF(:)
    DOUBLE PRECISION :: TT0

    CALL MPI_COMM_RANK(COMM, MYID, IERR)

    IF (MYID .EQ. ROOT) THEN
      HEADER(1) = MERGE(1, 0, BLK%IS_LOWRANK)
      HEADER(2) = BLK%RANK
    END IF
    IF (PRESENT(T_MPICALL)) TT0 = MPI_WTIME()
    CALL MPI_BCAST(HEADER, 2, MPI_INTEGER, ROOT, COMM, IERR)
    BCAST_NCALL = BCAST_NCALL + 1
    IF (PRESENT(T_MPICALL)) T_MPICALL = T_MPICALL + (MPI_WTIME() - TT0)

    IF (MYID .NE. ROOT) THEN
      ! BLK est reutilise d'une etape de bloc K a l'autre par l'appelant
      ! (une seule variable au niveau de la boucle, pas une nouvelle a
      ! chaque K) -- contrairement a COMPRESS_BLOCK (BLK en INTENT(OUT),
      ! deallocation automatique des composantes allouables a l'entree),
      ! BLK est ici en INTENT(INOUT) : il faut liberer nous-memes une
      ! eventuelle allocation laissee par l'etape PRECEDENTE avant d'en
      ! creer une nouvelle, sous peine de plantage "already allocated".
      IF (PRESENT(T_ALLOC)) TT0 = MPI_WTIME()
      IF (ALLOCATED(BLK%DENSE)) DEALLOCATE(BLK%DENSE)
      IF (ALLOCATED(BLK%UFAC))  DEALLOCATE(BLK%UFAC)
      IF (ALLOCATED(BLK%VFAC))  DEALLOCATE(BLK%VFAC)
      IF (ALLOCATED(BLK%SV))    DEALLOCATE(BLK%SV)
      BLK%ORTHO = .FALSE.
      BLK%M = M
      BLK%N = N
      BLK%IS_LOWRANK = (HEADER(1) .EQ. 1)
      BLK%RANK = HEADER(2)
      IF (BLK%IS_LOWRANK) THEN
        ALLOCATE(BLK%UFAC(M,BLK%RANK), BLK%VFAC(BLK%RANK,N))
      ELSE
        ALLOCATE(BLK%DENSE(M,N))
      END IF
      IF (PRESENT(T_ALLOC)) T_ALLOC = T_ALLOC + (MPI_WTIME() - TT0)
    END IF

    IF (BLK%IS_LOWRANK) THEN
      NU = M*BLK%RANK
      NV = BLK%RANK*N
      NTOT = NU + NV
      IF (NTOT .GT. 0) THEN
        ALLOCATE(BUF(NTOT))
        IF (MYID .EQ. ROOT) THEN
          IF (NU .GT. 0) BUF(1:NU)      = RESHAPE(BLK%UFAC, [NU])
          IF (NV .GT. 0) BUF(NU+1:NTOT) = RESHAPE(BLK%VFAC, [NV])
        END IF
        IF (PRESENT(T_MPICALL)) TT0 = MPI_WTIME()
        CALL MPI_BCAST(BUF, NTOT, MPI_COMPLEX, ROOT, COMM, IERR)
        BCAST_NCALL = BCAST_NCALL + 1; BCAST_ELEMS = BCAST_ELEMS + NTOT
        IF (PRESENT(T_MPICALL)) T_MPICALL = T_MPICALL + (MPI_WTIME() - TT0)
        IF (MYID .NE. ROOT) THEN
          IF (NU .GT. 0) BLK%UFAC = RESHAPE(BUF(1:NU),      [M,BLK%RANK])
          IF (NV .GT. 0) BLK%VFAC = RESHAPE(BUF(NU+1:NTOT), [BLK%RANK,N])
        END IF
        DEALLOCATE(BUF)
      END IF
    ELSE
      IF (M*N .GT. 0) THEN
        IF (PRESENT(T_MPICALL)) TT0 = MPI_WTIME()
        CALL MPI_BCAST(BLK%DENSE, M*N, MPI_COMPLEX, ROOT, COMM, IERR)
        BCAST_NCALL = BCAST_NCALL + 1; BCAST_ELEMS = BCAST_ELEMS + INT(M,8)*INT(N,8)
        IF (PRESENT(T_MPICALL)) T_MPICALL = T_MPICALL + (MPI_WTIME() - TT0)
      END IF
    END IF
  END SUBROUTINE BCAST_BLR_BLOCK

  !-----------------------------------------------------------------------
  ! Diffusion GROUPEE d'un panneau : les NBK blocs BLKS(IDX(1:NBK)) partent
  ! en deux appels collectifs (entetes, puis donnees empilees) au lieu de
  ! deux par bloc. Le nombre de rendez-vous collectifs par etape ne depend
  ! plus du nombre de blocs du panneau : c'est ce qui compte quand la
  ! grille de processus grandit (latence par appel x nombre d'appels).
  !
  ! ROOTS(B) = rang DANS COMM du processus qui possede le bloc B, MS/NS ses
  ! dimensions. Tous les membres de COMM appellent avec la meme liste.
  !   - une seule racine  : MPI_BCAST des entetes, MPI_BCAST des donnees ;
  !   - plusieurs racines : MPI_ALLREDUCE des entetes (chaque processus
  !     remplit ceux de ses blocs), MPI_ALLGATHERV des donnees.
  ! Les blocs recus sont alloues ici (comme BCAST_BLR_BLOCK) ; ceux de la
  ! racine ne sont pas modifies. Octets transportes identiques a la
  ! diffusion bloc par bloc : resultat identique au bit pres.
  !
  ! Repli sur BCAST_BLR_BLOCK bloc par bloc si le panneau empile depasse
  ! 2^31-1 elements (compteurs MPI 32 bits), ou si BLR_PANEL_PACK=0
  ! (comparaison avant/apres).
  !-----------------------------------------------------------------------
  SUBROUTINE BCAST_BLR_PANEL(COMM, NBK, IDX, ROOTS, MS, NS, BLKS, T_ALLOC, T_MPICALL)
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN)              :: COMM, NBK
    INTEGER, INTENT(IN)              :: IDX(*), ROOTS(*), MS(*), NS(*)
    TYPE(BLR_BLOCK_T), INTENT(INOUT) :: BLKS(:)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: T_ALLOC, T_MPICALL

    INTEGER, ALLOCATABLE :: HDR(:,:), HLOC(:,:), CNT(:), DSP(:)
    INTEGER(8), ALLOCATABLE :: LEN8(:), OFF8(:), RUN8(:)
    COMPLEX(WP), ALLOCATABLE :: BUF(:)
    INTEGER(8) :: TOT8, O8
    INTEGER :: IERR, MYID, NP, B, P, R, M, N, NU, NV
    LOGICAL :: SINGLE
    DOUBLE PRECISION :: TT0
    CHARACTER(LEN=8) :: ENVT

    IF (NBK .LE. 0) RETURN
    IF (PANEL_PACK_MODE .LT. 0) THEN
      CALL GET_ENVIRONMENT_VARIABLE('BLR_PANEL_PACK', ENVT)
      PANEL_PACK_MODE = MERGE(0, 1, TRIM(ENVT) .EQ. '0')
    END IF
    IF (PANEL_PACK_MODE .EQ. 0) THEN
      CALL PER_BLOCK()
      RETURN
    END IF

    CALL MPI_COMM_RANK(COMM, MYID, IERR)
    CALL MPI_COMM_SIZE(COMM, NP, IERR)
    SINGLE = ALL(ROOTS(1:NBK) .EQ. ROOTS(1))

    ! --- entetes : (basse-rang ?, rang) de chaque bloc ---
    ALLOCATE(HDR(2,NBK))
    HDR = 0
    DO B = 1, NBK
      IF (ROOTS(B) .EQ. MYID) THEN
        HDR(1,B) = MERGE(1, 0, BLKS(IDX(B))%IS_LOWRANK)
        HDR(2,B) = BLKS(IDX(B))%RANK
      END IF
    END DO
    IF (PRESENT(T_MPICALL)) TT0 = MPI_WTIME()
    IF (SINGLE) THEN
      CALL MPI_BCAST(HDR, 2*NBK, MPI_INTEGER, ROOTS(1), COMM, IERR)
    ELSE
      ALLOCATE(HLOC(2,NBK))
      HLOC = HDR
      CALL MPI_ALLREDUCE(HLOC, HDR, 2*NBK, MPI_INTEGER, MPI_SUM, COMM, IERR)
      DEALLOCATE(HLOC)
    END IF
    BCAST_NCALL = BCAST_NCALL + 1
    IF (PRESENT(T_MPICALL)) T_MPICALL = T_MPICALL + (MPI_WTIME() - TT0)

    ! --- longueur de chaque bloc et position dans le tampon empile :
    !     les blocs d'un meme processus se suivent (ordre de MPI_ALLGATHERV) ---
    ALLOCATE(LEN8(NBK), OFF8(NBK), RUN8(0:NP-1), CNT(0:NP-1), DSP(0:NP-1))
    RUN8 = 0_8
    DO B = 1, NBK
      IF (HDR(1,B) .EQ. 1) THEN
        LEN8(B) = (INT(MS(B),8) + INT(NS(B),8))*INT(HDR(2,B),8)
      ELSE
        LEN8(B) = INT(MS(B),8)*INT(NS(B),8)
      END IF
      RUN8(ROOTS(B)) = RUN8(ROOTS(B)) + LEN8(B)
    END DO
    TOT8 = SUM(RUN8)
    IF (TOT8 .GT. INT(HUGE(1),8)) THEN
      DEALLOCATE(HDR, LEN8, OFF8, RUN8, CNT, DSP)
      CALL PER_BLOCK()
      RETURN
    END IF
    O8 = 0_8
    DO P = 0, NP-1
      CNT(P) = INT(RUN8(P)); DSP(P) = INT(O8)
      O8 = O8 + RUN8(P)
      RUN8(P) = INT(DSP(P),8)
    END DO
    DO B = 1, NBK
      OFF8(B) = RUN8(ROOTS(B))
      RUN8(ROOTS(B)) = RUN8(ROOTS(B)) + LEN8(B)
    END DO

    ! --- blocs recus : allocation a la taille annoncee ---
    IF (PRESENT(T_ALLOC)) TT0 = MPI_WTIME()
    DO B = 1, NBK
      IF (ROOTS(B) .EQ. MYID) CYCLE
      P = IDX(B)
      IF (ALLOCATED(BLKS(P)%DENSE)) DEALLOCATE(BLKS(P)%DENSE)
      IF (ALLOCATED(BLKS(P)%UFAC))  DEALLOCATE(BLKS(P)%UFAC)
      IF (ALLOCATED(BLKS(P)%VFAC))  DEALLOCATE(BLKS(P)%VFAC)
      IF (ALLOCATED(BLKS(P)%SV))    DEALLOCATE(BLKS(P)%SV)
      BLKS(P)%ORTHO = .FALSE.
      BLKS(P)%M = MS(B)
      BLKS(P)%N = NS(B)
      BLKS(P)%IS_LOWRANK = (HDR(1,B) .EQ. 1)
      BLKS(P)%RANK = HDR(2,B)
      IF (BLKS(P)%IS_LOWRANK) THEN
        ALLOCATE(BLKS(P)%UFAC(MS(B),HDR(2,B)), BLKS(P)%VFAC(HDR(2,B),NS(B)))
      ELSE
        ALLOCATE(BLKS(P)%DENSE(MS(B),NS(B)))
      END IF
    END DO
    IF (PRESENT(T_ALLOC)) T_ALLOC = T_ALLOC + (MPI_WTIME() - TT0)

    IF (TOT8 .GT. 0_8) THEN
      ALLOCATE(BUF(TOT8))
      ! --- la racine de chaque bloc l'empile a sa place ---
      DO B = 1, NBK
        IF (ROOTS(B) .NE. MYID .OR. LEN8(B) .EQ. 0_8) CYCLE
        P = IDX(B); M = MS(B); N = NS(B)
        IF (HDR(1,B) .EQ. 1) THEN
          R = HDR(2,B); NU = M*R; NV = R*N
          IF (NU .GT. 0) CALL COPYV(NU, BLKS(P)%UFAC, BUF(OFF8(B)+1))
          IF (NV .GT. 0) CALL COPYV(NV, BLKS(P)%VFAC, BUF(OFF8(B)+NU+1))
        ELSE
          CALL COPYV(M*N, BLKS(P)%DENSE, BUF(OFF8(B)+1))
        END IF
      END DO
      IF (PRESENT(T_MPICALL)) TT0 = MPI_WTIME()
      IF (SINGLE) THEN
        CALL MPI_BCAST(BUF, INT(TOT8), MPI_COMPLEX, ROOTS(1), COMM, IERR)
      ELSE
        CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_COMPLEX, BUF, CNT, DSP, MPI_COMPLEX, COMM, IERR)
      END IF
      BCAST_NCALL = BCAST_NCALL + 1; BCAST_ELEMS = BCAST_ELEMS + TOT8
      IF (PRESENT(T_MPICALL)) T_MPICALL = T_MPICALL + (MPI_WTIME() - TT0)
      ! --- les autres le depilent ---
      DO B = 1, NBK
        IF (ROOTS(B) .EQ. MYID .OR. LEN8(B) .EQ. 0_8) CYCLE
        P = IDX(B); M = MS(B); N = NS(B)
        IF (HDR(1,B) .EQ. 1) THEN
          R = HDR(2,B); NU = M*R; NV = R*N
          IF (NU .GT. 0) CALL COPYV(NU, BUF(OFF8(B)+1), BLKS(P)%UFAC)
          IF (NV .GT. 0) CALL COPYV(NV, BUF(OFF8(B)+NU+1), BLKS(P)%VFAC)
        ELSE
          CALL COPYV(M*N, BUF(OFF8(B)+1), BLKS(P)%DENSE)
        END IF
      END DO
      DEALLOCATE(BUF)
    END IF
    DEALLOCATE(HDR, LEN8, OFF8, RUN8, CNT, DSP)

  CONTAINS

    SUBROUTINE PER_BLOCK()
      INTEGER :: BB
      DO BB = 1, NBK
        CALL BCAST_BLR_BLOCK(COMM, ROOTS(BB), MS(BB), NS(BB), BLKS(IDX(BB)), T_ALLOC, T_MPICALL)
      END DO
    END SUBROUTINE PER_BLOCK

  END SUBROUTINE BCAST_BLR_PANEL

  ! Copie de L elements (tableaux de forme quelconque, vus a plat).
  SUBROUTINE COPYV(L, SRC, DST)
    INTEGER, INTENT(IN)      :: L
    COMPLEX(WP), INTENT(IN)  :: SRC(*)
    COMPLEX(WP), INTENT(OUT) :: DST(*)
    DST(1:L) = SRC(1:L)
  END SUBROUTINE COPYV

  !-----------------------------------------------------------------------
  ! Meme structure generale que MPI_LU_PIVOT_FACTOR (mpi_lu_pivot_mod.f90)
  ! sauf la toute fin de chaque etape de bloc : compresse LBUF/UBUF cote
  ! source avant diffusion (au lieu de diffuser dense puis compresser cote
  ! recepteur).
  !-----------------------------------------------------------------------
  ! T_BCAST (optionnel) : temps cumule (max/processus) passe dans les 2
  ! appels BCAST_BLR_BLOCK. Un rang non-source qui atteint BCAST_BLR_BLOCK
  ! est deja bloque en attente pendant que la source (et elle seule)
  ! termine sa compression avant de pouvoir emettre quoi que ce soit -- ce
  ! temps d'attente est comptabilise dans T_BCAST (pas dans T_COMPRESS, qui
  ! ne mesure que le calcul de la source) : cette serialisation source-puis-
  ! diffusion est le cout a mettre en balance avec le gain de volume de
  ! communication documente en tete de fichier.
  SUBROUTINE MPI_BLR_LU_PIVOT_BCAST_FACTOR(GRID, N, NB, A_LOCAL, MLOC, NLOC, TOL, &
                                            IPIV_GLOBAL, INFO, N_COMPRESSED, N_TOTAL_PANELS, &
                                            RANK_SUM, BYTES_SENT, BYTES_DENSE_EQUIV, N_SWAPS, MIN_PIVOT, &
                                            T_COMPRESS, T_BCAST)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)        :: GRID
    INTEGER, INTENT(IN)             :: N, NB, MLOC, NLOC
    COMPLEX(WP), INTENT(INOUT)      :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    REAL(WP), INTENT(IN)            :: TOL
    INTEGER, INTENT(OUT)            :: IPIV_GLOBAL(N)
    INTEGER, INTENT(OUT)            :: INFO
    INTEGER, INTENT(OUT)            :: N_COMPRESSED, N_TOTAL_PANELS
    INTEGER(8), INTENT(OUT)         :: RANK_SUM
    INTEGER(8), INTENT(OUT)         :: BYTES_SENT        ! diagnostic : volume REEL envoye (compresse)
    INTEGER(8), INTENT(OUT)         :: BYTES_DENSE_EQUIV ! diagnostic : ce que la diffusion dense aurait envoye
    INTEGER, INTENT(OUT), OPTIONAL  :: N_SWAPS
    REAL(WP), INTENT(OUT), OPTIONAL :: MIN_PIVOT
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: T_COMPRESS, T_BCAST

    COMPLEX(WP), ALLOCATABLE :: PANEL(:,:), RBUF(:,:), SBUF(:,:), DKK(:,:), LBUF(:,:), UBUF(:,:)
    TYPE(BLR_BLOCK_T) :: LBLK, UBLK
    INTEGER :: IPIV_LOCAL(NB), IPIV_PANEL_GLOBAL(NB)
    INTEGER :: NBLK, K, KSTART, KEND, KW, MPANEL
    INTEGER :: PROW_K, PCOL_K, OWNER_RANK, RS, CS, RSTART0, NROWS_MINE
    INTEGER :: R, RMLOC_R, RSTART0_R, RNROWS_R, I, J, GI
    INTEGER :: LEFT_COUNT, RIGHT_COUNT, RTRAIL0, CTRAIL0, NROWTRAIL, NCOLTRAIL
    INTEGER :: R1, R2, PROW2, IINFO, INFO_STEP, IERR, NS
    INTEGER :: NC_LOCAL, NP_LOCAL
    INTEGER(8) :: RSUM_LOCAL, BYTES_LOCAL, BYTES_DENSE_LOCAL
    INTEGER :: STATUS(MPI_STATUS_SIZE)
    REAL(WP) :: MP
    DOUBLE PRECISION :: TC_LOCAL, TB_LOCAL, TC0

    INFO = 0
    NS = 0
    MP = HUGE(1.0_WP)
    NC_LOCAL = 0
    NP_LOCAL = 0
    RSUM_LOCAL = 0_8
    BYTES_LOCAL = 0_8
    BYTES_DENSE_LOCAL = 0_8
    TC_LOCAL = 0.0D0
    TB_LOCAL = 0.0D0
    NBLK = (N + NB - 1) / NB

    DO K = 1, NBLK
      KSTART = (K-1)*NB + 1
      KEND   = MIN(K*NB, N)
      KW     = KEND - KSTART + 1
      MPANEL = N - KSTART + 1

      PROW_K = INDXG2P(KSTART, NB, GRID%PR)
      PCOL_K = INDXG2P(KSTART, NB, GRID%PC)
      OWNER_RANK = PROW_K*GRID%PC + PCOL_K

      INFO_STEP = 0

      ! === A+B+C : rassembler le panneau, factoriser (CGETRF), renvoyer
      !     (rassembler-factoriser-renvoyer le panneau) ===
      IF (GRID%MYCOL .EQ. PCOL_K) THEN
        RSTART0    = NUMROC(KSTART-1, NB, GRID%MYROW, GRID%PR)
        NROWS_MINE = MLOC - RSTART0
        CS         = INDXG2L(KSTART, NB, GRID%PC)

        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(PANEL(MPANEL,KW))
          IF (NROWS_MINE .GT. 0) THEN
            DO J = 1, KW
              DO I = 1, NROWS_MINE
                GI = INDXL2G(RSTART0+I, NB, GRID%MYROW, GRID%PR)
                PANEL(GI-KSTART+1,J) = A_LOCAL(RSTART0+I, CS+J-1)
              END DO
            END DO
          END IF
          DO R = 0, GRID%PR-1
            IF (R .EQ. PROW_K) CYCLE
            RMLOC_R   = NUMROC(N, NB, R, GRID%PR)
            RSTART0_R = NUMROC(KSTART-1, NB, R, GRID%PR)
            RNROWS_R  = RMLOC_R - RSTART0_R
            IF (RNROWS_R .LE. 0) CYCLE
            ALLOCATE(RBUF(RNROWS_R,KW))
            CALL MPI_RECV(RBUF, RNROWS_R*KW, MPI_COMPLEX, R, 10, GRID%COL_COMM, STATUS, IERR)
            DO J = 1, KW
              DO I = 1, RNROWS_R
                GI = INDXL2G(RSTART0_R+I, NB, R, GRID%PR)
                PANEL(GI-KSTART+1,J) = RBUF(I,J)
              END DO
            END DO
            DEALLOCATE(RBUF)
          END DO

          CALL CGETRF(MPANEL, KW, PANEL, MPANEL, IPIV_LOCAL, IINFO)
          IF (IINFO .NE. 0) INFO_STEP = KSTART - 1 + IINFO
          DO J = 1, KW
            IPIV_PANEL_GLOBAL(J) = IPIV_LOCAL(J)
            MP = MIN(MP, ABS(PANEL(J,J)))
          END DO

          IF (NROWS_MINE .GT. 0) THEN
            DO J = 1, KW
              DO I = 1, NROWS_MINE
                GI = INDXL2G(RSTART0+I, NB, GRID%MYROW, GRID%PR)
                A_LOCAL(RSTART0+I, CS+J-1) = PANEL(GI-KSTART+1,J)
              END DO
            END DO
          END IF
          DO R = 0, GRID%PR-1
            IF (R .EQ. PROW_K) CYCLE
            RMLOC_R   = NUMROC(N, NB, R, GRID%PR)
            RSTART0_R = NUMROC(KSTART-1, NB, R, GRID%PR)
            RNROWS_R  = RMLOC_R - RSTART0_R
            IF (RNROWS_R .LE. 0) CYCLE
            ALLOCATE(SBUF(RNROWS_R,KW))
            DO J = 1, KW
              DO I = 1, RNROWS_R
                GI = INDXL2G(RSTART0_R+I, NB, R, GRID%PR)
                SBUF(I,J) = PANEL(GI-KSTART+1,J)
              END DO
            END DO
            CALL MPI_SEND(SBUF, RNROWS_R*KW, MPI_COMPLEX, R, 11, GRID%COL_COMM, IERR)
            DEALLOCATE(SBUF)
          END DO
          DEALLOCATE(PANEL)
        ELSE
          IF (NROWS_MINE .GT. 0) THEN
            ALLOCATE(SBUF(NROWS_MINE,KW))
            SBUF = A_LOCAL(RSTART0+1:MLOC, CS:CS+KW-1)
            CALL MPI_SEND(SBUF, NROWS_MINE*KW, MPI_COMPLEX, PROW_K, 10, GRID%COL_COMM, IERR)
            DEALLOCATE(SBUF)

            ALLOCATE(RBUF(NROWS_MINE,KW))
            CALL MPI_RECV(RBUF, NROWS_MINE*KW, MPI_COMPLEX, PROW_K, 11, GRID%COL_COMM, STATUS, IERR)
            A_LOCAL(RSTART0+1:MLOC, CS:CS+KW-1) = RBUF
            DEALLOCATE(RBUF)
          END IF
        END IF
      END IF

      CALL MPI_BCAST(INFO_STEP, 1, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IF (INFO_STEP .NE. 0) THEN
        INFO = INFO_STEP
        RETURN
      END IF
      CALL MPI_BCAST(IPIV_PANEL_GLOBAL, KW, MPI_INTEGER, OWNER_RANK, GRID%COMM, IERR)
      IPIV_GLOBAL(KSTART:KEND) = IPIV_PANEL_GLOBAL(1:KW)

      ! === E : propager les echanges de lignes (rassembler-factoriser-renvoyer le panneau) ===
      LEFT_COUNT  = NUMROC(KSTART-1, NB, GRID%MYCOL, GRID%PC)
      RIGHT_COUNT = NLOC - NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      DO I = 1, KW
        R1 = KSTART - 1 + I
        R2 = KSTART - 1 + IPIV_PANEL_GLOBAL(I)
        IF (R1 .NE. R2) THEN
          NS = NS + 1
          PROW2 = INDXG2P(R2, NB, GRID%PR)
          IF (LEFT_COUNT .GT. 0) THEN
            CALL SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW_K, PROW2, &
                                  1, LEFT_COUNT)
          END IF
          IF (RIGHT_COUNT .GT. 0) THEN
            CALL SWAP_ROWS_RANGE(GRID, NB, A_LOCAL, MLOC, NLOC, R1, R2, PROW_K, PROW2, &
                                  NLOC-RIGHT_COUNT+1, NLOC)
          END IF
        END IF
      END DO

      ! === F : U12 par CTRSM, puis compression A LA
      !     SOURCE + diffusion COMPRESSEE (nouveaute de cette proposition) ===
      RTRAIL0   = NUMROC(KEND, NB, GRID%MYROW, GRID%PR)
      CTRAIL0   = NUMROC(KEND, NB, GRID%MYCOL, GRID%PC)
      NROWTRAIL = MLOC - RTRAIL0
      NCOLTRAIL = NLOC - CTRAIL0

      IF (GRID%MYROW .EQ. PROW_K) THEN
        ALLOCATE(DKK(KW,KW))
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CS = INDXG2L(KSTART, NB, GRID%PC)
          DKK = A_LOCAL(RS:RS+KW-1, CS:CS+KW-1)
        END IF
        CALL MPI_BCAST(DKK, KW*KW, MPI_COMPLEX, PCOL_K, GRID%ROW_COMM, IERR)
        IF (NCOLTRAIL .GT. 0) THEN
          RS = INDXG2L(KSTART, NB, GRID%PR)
          CALL CTRSM('Left', 'Lower', 'No transpose', 'Unit', KW, NCOLTRAIL, &
                      (1.0_WP,0.0_WP), DKK, KW, A_LOCAL(RS,CTRAIL0+1), MLOC)
        END IF
        DEALLOCATE(DKK)
      END IF

      ! --- L : compresser cote source (MYCOL==PCOL_K) AVANT diffusion ---
      IF (NROWTRAIL .GT. 0) THEN
        IF (GRID%MYCOL .EQ. PCOL_K) THEN
          ALLOCATE(LBUF(NROWTRAIL,KW))
          CS = INDXG2L(KSTART, NB, GRID%PC)
          LBUF = A_LOCAL(RTRAIL0+1:MLOC, CS:CS+KW-1)
          TC0 = MPI_WTIME()
          CALL COMPRESS_BLOCK(NROWTRAIL, KW, LBUF, NROWTRAIL, TOL, LBLK)
          TC_LOCAL = TC_LOCAL + (MPI_WTIME() - TC0)
          DEALLOCATE(LBUF)
          NP_LOCAL = NP_LOCAL + 1
          BYTES_DENSE_LOCAL = BYTES_DENSE_LOCAL + 8_8*INT(NROWTRAIL,8)*KW
          IF (LBLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RSUM_LOCAL = RSUM_LOCAL + LBLK%RANK
            BYTES_LOCAL = BYTES_LOCAL + 8_8*(INT(NROWTRAIL,8)*LBLK%RANK + INT(LBLK%RANK,8)*KW)
          ELSE
            RSUM_LOCAL = RSUM_LOCAL + MIN(NROWTRAIL,KW)
            BYTES_LOCAL = BYTES_LOCAL + 8_8*INT(NROWTRAIL,8)*KW
          END IF
        END IF
        TC0 = MPI_WTIME()
        CALL BCAST_BLR_BLOCK(GRID%ROW_COMM, PCOL_K, NROWTRAIL, KW, LBLK)
        TB_LOCAL = TB_LOCAL + (MPI_WTIME() - TC0)
      END IF

      ! --- U : compresser cote source (MYROW==PROW_K) AVANT diffusion ---
      IF (NCOLTRAIL .GT. 0) THEN
        IF (GRID%MYROW .EQ. PROW_K) THEN
          ALLOCATE(UBUF(KW,NCOLTRAIL))
          RS = INDXG2L(KSTART, NB, GRID%PR)
          UBUF = A_LOCAL(RS:RS+KW-1, CTRAIL0+1:NLOC)
          TC0 = MPI_WTIME()
          CALL COMPRESS_BLOCK(KW, NCOLTRAIL, UBUF, KW, TOL, UBLK)
          TC_LOCAL = TC_LOCAL + (MPI_WTIME() - TC0)
          DEALLOCATE(UBUF)
          NP_LOCAL = NP_LOCAL + 1
          BYTES_DENSE_LOCAL = BYTES_DENSE_LOCAL + 8_8*INT(KW,8)*NCOLTRAIL
          IF (UBLK%IS_LOWRANK) THEN
            NC_LOCAL = NC_LOCAL + 1
            RSUM_LOCAL = RSUM_LOCAL + UBLK%RANK
            BYTES_LOCAL = BYTES_LOCAL + 8_8*(INT(KW,8)*UBLK%RANK + INT(UBLK%RANK,8)*NCOLTRAIL)
          ELSE
            RSUM_LOCAL = RSUM_LOCAL + MIN(KW,NCOLTRAIL)
            BYTES_LOCAL = BYTES_LOCAL + 8_8*INT(KW,8)*NCOLTRAIL
          END IF
        END IF
        TC0 = MPI_WTIME()
        CALL BCAST_BLR_BLOCK(GRID%COL_COMM, PROW_K, KW, NCOLTRAIL, UBLK)
        TB_LOCAL = TB_LOCAL + (MPI_WTIME() - TC0)
      END IF

      IF (NROWTRAIL.GT.0 .AND. NCOLTRAIL.GT.0) THEN
        CALL SCHUR_UPDATE_BLOCK(NROWTRAIL, NCOLTRAIL, KW, LBLK, UBLK, &
                                 A_LOCAL(RTRAIL0+1,CTRAIL0+1), MLOC)
      END IF
    END DO

    CALL MPI_ALLREDUCE(NC_LOCAL, N_COMPRESSED, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(NP_LOCAL, N_TOTAL_PANELS, 1, MPI_INTEGER, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(RSUM_LOCAL, RANK_SUM, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(BYTES_LOCAL, BYTES_SENT, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(BYTES_DENSE_LOCAL, BYTES_DENSE_EQUIV, 1, MPI_INTEGER8, MPI_SUM, GRID%COMM, IERR)
    IF (PRESENT(N_SWAPS)) N_SWAPS = NS
    IF (PRESENT(MIN_PIVOT)) THEN
      CALL MPI_ALLREDUCE(MP, MIN_PIVOT, 1, MPI_REAL, MPI_MIN, GRID%COMM, IERR)
    END IF
    IF (PRESENT(T_COMPRESS)) THEN
      CALL MPI_ALLREDUCE(TC_LOCAL, T_COMPRESS, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
    IF (PRESENT(T_BCAST)) THEN
      CALL MPI_ALLREDUCE(TB_LOCAL, T_BCAST, 1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    END IF
  END SUBROUTINE MPI_BLR_LU_PIVOT_BCAST_FACTOR

END MODULE MPI_BLR_LU_PIVOT_BCAST_MOD
