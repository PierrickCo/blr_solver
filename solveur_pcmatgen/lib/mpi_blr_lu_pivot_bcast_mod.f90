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
END MODULE MPI_BLR_LU_PIVOT_BCAST_MOD
