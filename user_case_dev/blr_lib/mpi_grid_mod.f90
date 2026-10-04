!===============================================================================
! Infrastructure de distribution bloc-cyclique 2D (a la ScaLAPACK/BLACS),
! base commune du chemin MPI du solveur -- aucune factorisation ici,
! uniquement : grille de processus, formules de correspondance d'indices
! local<->global, et une routine de rassemblement (GATHER) pour valider la
! distribution avant de construire quoi que ce soit dessus.
!
! POURQUOI 2D BLOC-CYCLIQUE : c'est le schema utilise par ScaLAPACK
! (PxGETRF), avec un meilleur equilibrage de charge et un volume de
! communication moindre qu'une distribution 1D par blocs de lignes, au
! prix d'une complexite plus grande dans les modules qui s'appuient dessus.
!
! GRILLE DE PROCESSUS : PR x PC processus (PR lignes, PC colonnes de la
! grille -- PR*PC = nombre total de processus MPI). Convention retenue
! partout dans ce solveur : numerotation "ligne-major" -- rang MPI =
! MYROW*PC + MYCOL, MYROW = rang/PC, MYCOL = MOD(rang,PC).
!
! DEUX SOUS-COMMUNICATEURS PAR PROCESSUS (utilises par les modules de
! factorisation : la diffusion du panneau L se fait le long d'une ligne de
! la grille, celle du panneau U le long d'une colonne) :
!   ROW_COMM : tous les processus qui partagent le MEME MYROW (varie MYCOL)
!   COL_COMM : tous les processus qui partagent le MEME MYCOL (varie MYROW)
!
! DECOUPAGE BLOC-CYCLIQUE 1D (applique separement aux lignes ET aux
! colonnes) : le bloc global d'indice IBLK (0-based) est possede par le
! processus IBLK mod NPROCS (NPROCS = PR pour les lignes, PC pour les
! colonnes) -- "cyclique" car apres avoir distribue une premiere serie de
! blocs a tous les processus, la distribution recommence au premier. Les
! fonctions NUMROC/INDXL2G/INDXG2L/INDXG2P ci-dessous sont les formules
! standard de ScaLAPACK pour ce decoupage.
!===============================================================================
MODULE MPI_GRID_MOD
  USE BLR_KINDS_MOD, ONLY: WP, BLR_MAX_THREADS, BLR_NUM_PROCS
  USE MPI
  IMPLICIT NONE

  TYPE :: GRID_T
    INTEGER :: COMM        ! communicateur MPI de la grille (souvent MPI_COMM_WORLD)
    INTEGER :: NPROCS      ! PR*PC
    INTEGER :: PR, PC      ! dimensions de la grille
    INTEGER :: MYID        ! rang dans COMM
    INTEGER :: MYROW, MYCOL
    INTEGER :: ROW_COMM    ! communicateur des processus de MEME MYROW
    INTEGER :: COL_COMM    ! communicateur des processus de MEME MYCOL
  END TYPE GRID_T

CONTAINS

  !-----------------------------------------------------------------------
  ! Initialise la grille PR x PC a partir de COMM (typiquement
  ! MPI_COMM_WORLD). Echoue proprement (IERR renvoye) si PR*PC ne
  ! correspond pas au nombre de processus reellement lances.
  !-----------------------------------------------------------------------
  SUBROUTINE GRID_INIT(GRID, COMM, PR, PC, IERR)
    TYPE(GRID_T), INTENT(OUT) :: GRID
    INTEGER, INTENT(IN)       :: COMM, PR, PC
    INTEGER, INTENT(OUT)      :: IERR

    GRID%COMM = COMM
    GRID%PR = PR
    GRID%PC = PC
    CALL MPI_COMM_SIZE(COMM, GRID%NPROCS, IERR)
    IF (GRID%NPROCS .NE. PR*PC) THEN
      IERR = -1
      RETURN
    END IF
    CALL MPI_COMM_RANK(COMM, GRID%MYID, IERR)

    GRID%MYROW = GRID%MYID / PC
    GRID%MYCOL = MOD(GRID%MYID, PC)

    CALL MPI_COMM_SPLIT(COMM, GRID%MYROW, GRID%MYCOL, GRID%ROW_COMM, IERR)
    CALL MPI_COMM_SPLIT(COMM, GRID%MYCOL, GRID%MYROW, GRID%COL_COMM, IERR)
  END SUBROUTINE GRID_INIT

  !-----------------------------------------------------------------------
  ! Diagnostic OMP_NUM_THREADS / affinite CPU -- a appeler une fois, juste
  ! apres MPI_INIT, par tout programme MPI+OpenMP destine a tourner sur un
  ! cluster (pas necessaire pour les petits tests de correction). Rend
  ! visible, des le lancement, la configuration reelle (processus MPI par
  ! noeud x threads OpenMP par processus, comparee aux coeurs disponibles)
  ! plutot que de la laisser se manifester indirectement par une
  ! performance degradee.
  !
  ! COMM : communicateur global (typiquement MPI_COMM_WORLD), utilise pour
  ! detecter quels rangs partagent la meme memoire physique (meme noeud)
  ! via MPI_COMM_SPLIT_TYPE(MPI_COMM_TYPE_SHARED) -- une seule ligne de
  ! diagnostic est imprimee par noeud (par le rang 0 de chaque groupe
  ! partage), pas par processus, pour rester lisible avec beaucoup de
  ! rangs.
  !
  ! DETECTION DE LA SUR-SOUSCRIPTION : un simple "processus_sur_noeud x
  ! OMP_NUM_THREADS vs coeurs vus par CE processus" donne un faux positif
  ! des qu'un pinning CPU est actif (ex : Intel MPI epingle par defaut
  ! chaque rang sur un sous-ensemble disjoint de coeurs -- avec 4 rangs sur
  ! une machine a 8 coeurs, chacun epingle sur 2 coeurs disjoints, "coeurs
  ! vus par ce processus" vaut 2, et 4x2=8 declencherait a tort une alerte
  ! alors que les 8 coeurs physiques sont exactement et exclusivement
  ! couverts).
  !
  ! Au lieu de multiplier un nombre de coeurs "par processus" (qui ne dit
  ! pas si les processus se partagent les memes coeurs ou en occupent des
  ! disjoints), le masque d'affinite CPU reel de chaque processus
  ! (/proc/self/status, "Cpus_allowed") est lu, puis son UNION sur tous les
  ! rangs du noeud calculee (MPI_ALLREDUCE MAX bit-a-bit) -- le nombre de
  ! coeurs distincts reellement occupes sur le noeud, pinning ou pas -- et
  ! comparee a la somme des threads demandes par tous les rangs du noeud.
  ! Ce test est correct dans les deux regimes :
  !   - pinning actif (masques disjoints) : UNION = somme des tailles de
  !     masques -> reflete le vrai total de coeurs occupes.
  !   - pas de pinning (tous les rangs voient le noeud entier, masques
  !     identiques) : UNION = taille d'UN masque (le noeud entier) -> la
  !     somme des threads demandes (elle, pas dedupliquee) depasse vite ce
  !     total des que plusieurs processus tournent sans restriction.
  !-----------------------------------------------------------------------
  SUBROUTINE REPORT_THREAD_AFFINITY(COMM)
    INTEGER, INTENT(IN) :: COMM
    INTEGER, PARAMETER :: MAXBITS = 2048
    INTEGER :: IERR, SHM_COMM, RANKS_ON_NODE, SHM_RANK, NAME_LEN
    INTEGER :: NTHREADS, NPROCS_VISIBLE, TOTAL_THREADS_NODE, TOTAL_CORES_USED
    INTEGER :: MYMASK(0:MAXBITS-1), UNIONMASK(0:MAXBITS-1)
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: NODE_NAME

    CALL MPI_GET_PROCESSOR_NAME(NODE_NAME, NAME_LEN, IERR)

    CALL MPI_COMM_SPLIT_TYPE(COMM, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, SHM_COMM, IERR)
    CALL MPI_COMM_SIZE(SHM_COMM, RANKS_ON_NODE, IERR)
    CALL MPI_COMM_RANK(SHM_COMM, SHM_RANK, IERR)

    NTHREADS       = BLR_MAX_THREADS()
    NPROCS_VISIBLE = BLR_NUM_PROCS()

    CALL GET_CPU_AFFINITY_MASK(MYMASK, MAXBITS)
    ! sans OpenMP : coeurs du masque d'affinite de ce processus
    IF (NPROCS_VISIBLE .LT. 1) NPROCS_VISIBLE = COUNT(MYMASK .NE. 0)
    CALL MPI_ALLREDUCE(MYMASK, UNIONMASK, MAXBITS, MPI_INTEGER, MPI_MAX, SHM_COMM, IERR)
    TOTAL_CORES_USED = SUM(UNIONMASK)
    CALL MPI_ALLREDUCE(NTHREADS, TOTAL_THREADS_NODE, 1, MPI_INTEGER, MPI_SUM, SHM_COMM, IERR)

    IF (SHM_RANK .EQ. 0) THEN
      WRITE(6,'(A)') '--------------------------------------------------------------'
      WRITE(6,'(A,A)')  ' [AFFINITE] Noeud                           = ', TRIM(NODE_NAME)
      WRITE(6,'(A,I0)') ' [AFFINITE] Processus MPI sur ce noeud      = ', RANKS_ON_NODE
      WRITE(6,'(A,I0)') ' [AFFINITE] OMP_NUM_THREADS (ce processus)  = ', NTHREADS
      WRITE(6,'(A,I0)') ' [AFFINITE] Coeurs vus (ce processus)       = ', NPROCS_VISIBLE
      WRITE(6,'(A,I0)') ' [AFFINITE] Coeurs distincts occupes/noeud  = ', TOTAL_CORES_USED
      WRITE(6,'(A,I0)') ' [AFFINITE] Total threads demandes / noeud  = ', TOTAL_THREADS_NODE
      IF (TOTAL_THREADS_NODE .GT. TOTAL_CORES_USED) THEN
        WRITE(6,'(A)') ' [AFFINITE] ATTENTION : sur-souscription reelle sur ce noeud'
        WRITE(6,'(A)') '            (total threads demandes > coeurs distincts occupes).'
        WRITE(6,'(A)') '            Fixez OMP_NUM_THREADS = coeurs_par_noeud / processus_MPI_par_noeud,'
        WRITE(6,'(A)') '            ou activez/verifiez le pinning CPU (ex: --cpu-bind sous Slurm).'
      ELSE
        WRITE(6,'(A)') ' [AFFINITE] OK : pas de sur-souscription (masques d''affinite pris en compte).'
      END IF
      WRITE(6,'(A)') '--------------------------------------------------------------'
      CALL FLUSH(6)
    END IF

    CALL MPI_COMM_FREE(SHM_COMM, IERR)
  END SUBROUTINE REPORT_THREAD_AFFINITY

  !-----------------------------------------------------------------------
  ! Affiche l'empreinte memoire THEORIQUE (calculee a partir des tailles
  ! de tableaux, PAS mesuree via l'OS -- voir GET_RSS_KB/bench_utils_mod.f90
  ! pour la RSS reelle, qui inclut l'overhead allocateur et les copies
  ! transitoires) d'une quantite donnee (matrice, facteurs LU, second
  ! membre...), en Mo, sous TROIS formes : le TOTAL agrege sur toute la
  ! grille (utile pour la capacite globale du cluster), le MAXIMUM
  ! par-processus (utile pour verifier que CHAQUE noeud/processus tient
  ! dans sa RAM disponible -- c'est la contrainte qui bloque reellement un
  ! job, pas le total), et le MINIMUM par-processus (compare au max, permet
  ! de juger si la charge est bien REPARTIE entre processus -- max proche
  ! du min = repartition equilibree ; max >> min = desequilibre, un ou
  ! quelques processus concentrent la memoire pendant que d'autres restent
  ! sous-charges). BYTES_LOCAL : empreinte de CE processus SEUL, calculee
  ! par l'appelant (ex: MLOC*NLOC*8 pour un tableau dense COMPLEX(WP), ou
  ! MPI_BLR_MATRIX_MEMORY_BYTES(A_BLR) pour du BLR compresse) -- reduite
  ! ici via MPI_ALLREDUCE (SUM, MAX et MIN).
  !-----------------------------------------------------------------------
  SUBROUTINE REPORT_MEMORY_MB(GRID, LABEL, BYTES_LOCAL)
    TYPE(GRID_T), INTENT(IN)        :: GRID
    CHARACTER(LEN=*), INTENT(IN)    :: LABEL
    INTEGER(8), INTENT(IN)          :: BYTES_LOCAL
    DOUBLE PRECISION :: MB_LOCAL, MB_TOTAL, MB_MAX, MB_MIN
    INTEGER :: IERR

    MB_LOCAL = REAL(BYTES_LOCAL,8) / 1024.0D0 / 1024.0D0
    CALL MPI_ALLREDUCE(MB_LOCAL, MB_TOTAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(MB_LOCAL, MB_MAX,   1, MPI_DOUBLE_PRECISION, MPI_MAX, GRID%COMM, IERR)
    CALL MPI_ALLREDUCE(MB_LOCAL, MB_MIN,   1, MPI_DOUBLE_PRECISION, MPI_MIN, GRID%COMM, IERR)
    IF (GRID%MYID .EQ. 0) THEN
      WRITE(6,'(A,A,A,F10.1,A,F10.1,A,F10.1,A)') '   ', LABEL, ' : total = ', MB_TOTAL, &
           ' Mo , min/processus = ', MB_MIN, ' Mo , max/processus = ', MB_MAX, ' Mo'
    END IF
  END SUBROUTINE REPORT_MEMORY_MB

  !-----------------------------------------------------------------------
  ! Lit le masque d'affinite CPU du processus courant (/proc/self/status,
  ! champ "Cpus_allowed", format hexadecimal eventuellement groupe par
  ! virgules pour les masques > 32 bits) et le decode en tableau 0/1
  ! MASK(0:MAXBITS-1) (bit i = 1 si le coeur logique i est autorise pour ce
  ! processus). Specifique Linux (format /proc). Si le fichier ou le champ
  ! est introuvable, renvoie un masque nul (le diagnostic degrade
  ! proprement, il n'affichera aucune alerte plutot que de planter).
  !-----------------------------------------------------------------------
  SUBROUTINE GET_CPU_AFFINITY_MASK(MASK, MAXBITS)
    INTEGER, INTENT(IN)  :: MAXBITS
    INTEGER, INTENT(OUT) :: MASK(0:MAXBITS-1)
    INTEGER :: UNIT, IOS, POS, NC, IC, BITBASE, HV
    CHARACTER(LEN=512) :: LINE, HEXSTR, CLEANHEX
    CHARACTER :: CH
    LOGICAL :: FOUND

    MASK = 0
    HEXSTR = ' '
    FOUND = .FALSE.
    OPEN(NEWUNIT=UNIT, FILE='/proc/self/status', STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) RETURN
    DO
      READ(UNIT,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) EXIT
      IF (LINE(1:13) .EQ. 'Cpus_allowed:') THEN
        HEXSTR = ADJUSTL(LINE(14:))
        FOUND = .TRUE.
        EXIT
      END IF
    END DO
    CLOSE(UNIT)
    IF (.NOT. FOUND) RETURN

    NC = 0
    CLEANHEX = ' '
    DO POS = 1, LEN_TRIM(HEXSTR)
      CH = HEXSTR(POS:POS)
      IF (CH .NE. ',' .AND. CH .NE. ' ') THEN
        NC = NC + 1
        CLEANHEX(NC:NC) = CH
      END IF
    END DO

    BITBASE = 0
    DO POS = NC, 1, -1
      CH = CLEANHEX(POS:POS)
      HV = HEXDIGIT_VALUE(CH)
      DO IC = 0, 3
        IF (BITBASE+IC .LE. MAXBITS-1 .AND. BTEST(HV, IC)) MASK(BITBASE+IC) = 1
      END DO
      BITBASE = BITBASE + 4
    END DO
  END SUBROUTINE GET_CPU_AFFINITY_MASK

  PURE FUNCTION HEXDIGIT_VALUE(CH) RESULT(V)
    CHARACTER, INTENT(IN) :: CH
    INTEGER :: V
    SELECT CASE (CH)
    CASE ('0':'9')
      V = ICHAR(CH) - ICHAR('0')
    CASE ('a':'f')
      V = ICHAR(CH) - ICHAR('a') + 10
    CASE ('A':'F')
      V = ICHAR(CH) - ICHAR('A') + 10
    CASE DEFAULT
      V = 0
    END SELECT
  END FUNCTION HEXDIGIT_VALUE

  !-----------------------------------------------------------------------
  ! Nombre d'elements (lignes OU colonnes, meme formule dans les 2 cas)
  ! possedes localement par le processus IPROC (0-based) sur NPROCS, pour
  ! une dimension globale N decoupee en blocs de taille NB. Formule
  ! standard ScaLAPACK NUMROC (ISRCPROC=0 -- le premier bloc appartient
  ! toujours au processus 0).
  !-----------------------------------------------------------------------
  PURE FUNCTION NUMROC(N, NB, IPROC, NPROCS) RESULT(LOCAL_N)
    INTEGER, INTENT(IN) :: N, NB, IPROC, NPROCS
    INTEGER :: LOCAL_N, NBLOCKS, EXTRABLOCKS

    NBLOCKS = N / NB
    LOCAL_N = (NBLOCKS/NPROCS) * NB
    EXTRABLOCKS = MOD(NBLOCKS, NPROCS)
    IF (IPROC .LT. EXTRABLOCKS) THEN
      LOCAL_N = LOCAL_N + NB
    ELSE IF (IPROC .EQ. EXTRABLOCKS) THEN
      LOCAL_N = LOCAL_N + MOD(N, NB)
    END IF
  END FUNCTION NUMROC

  !-----------------------------------------------------------------------
  ! Indice LOCAL (1-based) -> indice GLOBAL (1-based), pour un processus
  ! IPROC (0-based) sur NPROCS, bloc de taille NB.
  !-----------------------------------------------------------------------
  PURE FUNCTION INDXL2G(LOCAL_IDX, NB, IPROC, NPROCS) RESULT(GLOBAL_IDX)
    INTEGER, INTENT(IN) :: LOCAL_IDX, NB, IPROC, NPROCS
    INTEGER :: GLOBAL_IDX
    GLOBAL_IDX = NPROCS*NB*((LOCAL_IDX-1)/NB) + MOD(LOCAL_IDX-1,NB) + IPROC*NB + 1
  END FUNCTION INDXL2G

  !-----------------------------------------------------------------------
  ! Indice GLOBAL (1-based) -> indice LOCAL (1-based) SUR LE PROCESSUS QUI
  ! LE POSSEDE (appelant responsable de verifier via INDXG2P que c'est
  ! bien lui le proprietaire -- sinon resultat sans sens).
  !-----------------------------------------------------------------------
  PURE FUNCTION INDXG2L(GLOBAL_IDX, NB, NPROCS) RESULT(LOCAL_IDX)
    INTEGER, INTENT(IN) :: GLOBAL_IDX, NB, NPROCS
    INTEGER :: LOCAL_IDX
    LOCAL_IDX = NB*((GLOBAL_IDX-1)/(NB*NPROCS)) + MOD(GLOBAL_IDX-1,NB) + 1
  END FUNCTION INDXG2L

  !-----------------------------------------------------------------------
  ! Indice GLOBAL (1-based) -> processus proprietaire (0-based).
  !-----------------------------------------------------------------------
  PURE FUNCTION INDXG2P(GLOBAL_IDX, NB, NPROCS) RESULT(IPROC)
    INTEGER, INTENT(IN) :: GLOBAL_IDX, NB, NPROCS
    INTEGER :: IPROC
    IPROC = MOD((GLOBAL_IDX-1)/NB, NPROCS)
  END FUNCTION INDXG2P

  !-----------------------------------------------------------------------
  ! Rassemble sur ROOT (rang dans GRID%COMM) une matrice locale
  ! bloc-cyclique 2D (A_LOCAL, LLDxLLDC utiles, MLOC x NLOC elements) en la
  ! matrice GLOBALE A_GLOBAL (N x N, seule valide sur ROOT en sortie).
  ! Usage : validation de la distribution et diagnostics -- pas destine a
  ! etre utilise a chaque iteration d'un vrai calcul (couteux, communication
  ! en etoile vers ROOT).
  !-----------------------------------------------------------------------
  SUBROUTINE GATHER_2D(GRID, N, NB, A_LOCAL, MLOC, NLOC, ROOT, A_GLOBAL)
    TYPE(GRID_T), INTENT(IN)   :: GRID
    INTEGER, INTENT(IN)        :: N, NB, MLOC, NLOC, ROOT
    COMPLEX(WP), INTENT(IN)    :: A_LOCAL(MLOC,NLOC)
    COMPLEX(WP), INTENT(OUT)   :: A_GLOBAL(N,N)

    COMPLEX(WP), ALLOCATABLE :: RECVBUF(:,:)
    INTEGER :: P, RMLOC, RNLOC, RMYROW, RMYCOL, I, J, GI, GJ, IERR
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    IF (GRID%MYID .EQ. ROOT) THEN
      ! Copie directe de la propre part de ROOT
      DO J = 1, NLOC
        GJ = INDXL2G(J, NB, GRID%MYCOL, GRID%PC)
        DO I = 1, MLOC
          GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
          A_GLOBAL(GI,GJ) = A_LOCAL(I,J)
        END DO
      END DO
      ! Reception de la part de chaque autre processus
      DO P = 0, GRID%NPROCS-1
        IF (P .EQ. ROOT) CYCLE
        RMYROW = P / GRID%PC
        RMYCOL = MOD(P, GRID%PC)
        RMLOC = NUMROC(N, NB, RMYROW, GRID%PR)
        RNLOC = NUMROC(N, NB, RMYCOL, GRID%PC)
        IF (RMLOC.EQ.0 .OR. RNLOC.EQ.0) CYCLE
        ALLOCATE(RECVBUF(RMLOC,RNLOC))
        CALL MPI_RECV(RECVBUF, RMLOC*RNLOC, MPI_COMPLEX, P, 0, GRID%COMM, &
                       STATUS, IERR)
        DO J = 1, RNLOC
          GJ = INDXL2G(J, NB, RMYCOL, GRID%PC)
          DO I = 1, RMLOC
            GI = INDXL2G(I, NB, RMYROW, GRID%PR)
            A_GLOBAL(GI,GJ) = RECVBUF(I,J)
          END DO
        END DO
        DEALLOCATE(RECVBUF)
      END DO
    ELSE
      IF (MLOC.GT.0 .AND. NLOC.GT.0) THEN
        CALL MPI_SEND(A_LOCAL, MLOC*NLOC, MPI_COMPLEX, ROOT, 0, GRID%COMM, IERR)
      END IF
    END IF
  END SUBROUTINE GATHER_2D

  !-----------------------------------------------------------------------
  ! Rassemble sur ROOT un vecteur distribue PAR LIGNES uniquement (NUMROC
  ! sur PR, un element par MYROW) et REPLIQUE sur tous les MYCOL. Utilisee
  ! par legacy/src/mpi_lu_solve_mod.f90. Seuls les processus MYCOL==0
  ! participent a l'envoi (evite d'envoyer PC fois la meme donnee) ; ROOT
  ! doit lui-meme etre un processus MYCOL==0 (toujours vrai pour ROOT=0
  ! avec la convention de numerotation ligne-major).
  !-----------------------------------------------------------------------
  SUBROUTINE GATHER_ROW_VECTOR(GRID, N, NB, Y_LOCAL, MLOC, ROOT, Y_GLOBAL)
    TYPE(GRID_T), INTENT(IN) :: GRID
    INTEGER, INTENT(IN)      :: N, NB, MLOC, ROOT
    COMPLEX(WP), INTENT(IN)  :: Y_LOCAL(MLOC)
    COMPLEX(WP), INTENT(OUT) :: Y_GLOBAL(N)

    COMPLEX(WP), ALLOCATABLE :: RECVBUF(:)
    INTEGER :: R, RMLOC, I, GI, IERR, SRC_RANK
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    IF (GRID%MYCOL .NE. 0) RETURN

    IF (GRID%MYID .EQ. ROOT) THEN
      DO I = 1, MLOC
        GI = INDXL2G(I, NB, GRID%MYROW, GRID%PR)
        Y_GLOBAL(GI) = Y_LOCAL(I)
      END DO
      DO R = 0, GRID%PR-1
        SRC_RANK = R*GRID%PC   ! rang du processus (R,0) -- numerotation ligne-major
        IF (SRC_RANK .EQ. ROOT) CYCLE
        RMLOC = NUMROC(N, NB, R, GRID%PR)
        IF (RMLOC .EQ. 0) CYCLE
        ALLOCATE(RECVBUF(RMLOC))
        CALL MPI_RECV(RECVBUF, RMLOC, MPI_COMPLEX, SRC_RANK, 1, GRID%COMM, &
                       STATUS, IERR)
        DO I = 1, RMLOC
          GI = INDXL2G(I, NB, R, GRID%PR)
          Y_GLOBAL(GI) = RECVBUF(I)
        END DO
        DEALLOCATE(RECVBUF)
      END DO
    ELSE
      IF (MLOC .GT. 0) THEN
        CALL MPI_SEND(Y_LOCAL, MLOC, MPI_COMPLEX, ROOT, 1, GRID%COMM, IERR)
      END IF
    END IF
  END SUBROUTINE GATHER_ROW_VECTOR

END MODULE MPI_GRID_MOD
