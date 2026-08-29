!===============================================================================
! ARENA (pool memoire) POUR LE STOCKAGE PERSISTANT DES BLOCS BLR HORS-
! DIAGONALE PENDANT UNE FACTORISATION MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA
! (voir mpi_blr_matrix_mod.f90) -- correctif "lourd" identifie et documente
! dans la memoire du projet apres une investigation methodique de l'ecart
! RSS-vs-theorique (~39 Mo/processus non expliques a N=20000, 4 hypotheses
! plus legeres testees et ecartees : cache MKL, malloc_trim, MOVE_ALLOC,
! seuil mmap -- voir examples/README.md).
!
! PRINCIPE : au lieu de laisser CHAQUE bloc BLR_BLOCK_T faire son propre
! ALLOCATE/DEALLOCATE a chaque recompression (des milliers de cycles de
! tailles VARIABLES pendant une factorisation, cause de la fragmentation
! de tas observee), on reserve de GROS buffers ("chunks") a l'avance et on
! y DECOUPE des vues (pointeurs) pour chaque bloc -- glibc ne voit alors
! qu'une poignee de GRANDES allocations (une par chunk), jamais les
! milliers de petites.
!
! POURQUOI UN ARENA COUVRANT TOUTE LA FACTORISATION (pas un arena par
! etape K, comme envisage initialement) : le stockage d'un bloc doit
! SURVIVRE au-dela d'une seule etape K (il est relu/remis a jour a des
! etapes ulterieures) -- un arena vide a chaque fin d'etape detruirait a
! tort des blocs encore necessaires.
!
! COMPROMIS ASSUME (mesure, pas suppose -- voir examples/README.md pour
! les chiffres) : cet arena ne RECUPERE JAMAIS l'espace d'un bloc
! remplace par une nouvelle version (pas de compaction, cf. discussion
! sur les limites fondamentales des allocateurs C/Fortran -- un pointeur
! deja distribue ne peut pas etre deplace). L'arena peut donc consommer
! PLUS de memoire au pic qu'une gestion bloc-par-bloc, en echange de
! l'elimination de la fragmentation de tas -- c'est exactement le
! compromis a mesurer en comparant aux executables existants (jamais
! modifies, voir mpi_blr_matrix_mod.f90).
!
! POURQUOI CHUNKS(:) EST UN TABLEAU DE TAILLE FIXE (ARENA_MAX_CHUNKS), PAS
! ALLOCATABLE : agrandir un tableau ALLOCATABLE de types derives a
! composante ALLOCATABLE (via une reaffectation) declenche la MEME
! reallocation-et-copie implicite que "ABLK = NEWBLK" (voir
! blr_compress_mod.f90/SCHUR_UPDATE_BLR_BLOCK) -- ca DEPLACERAIT le
! buffer BUF de chaque chunk deja actif, invalidant tous les pointeurs
! deja distribues vers ce chunk. Un tableau de taille FIXE (jamais
! redimensionne) evite ce piege : seul BUF de chaque chunk INDIVIDUEL est
! alloue (une fois, jamais redeplace ensuite), le tableau CHUNKS(:) qui
! les contient ne bouge jamais.
!===============================================================================
MODULE BLR_ARENA_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  IMPLICIT NONE

  ! Taille standard d'un chunk : 2 M elements COMPLEX(WP) = 16 Mo. Choisi
  ! assez grand pour contenir de nombreux blocs BLR typiques (NB~500,
  ! rang bas-rang ~5-50 -> quelques dizaines de Ko/bloc ; un bloc dense de
  ! repli NBxNB ~2 Mo tient aussi confortablement dedans).
  INTEGER(8), PARAMETER :: ARENA_CHUNK_ELEMS = 2097152_8

  ! Borne haute sur le nombre de chunks -- genereuse et JAMAIS
  ! redimensionnee (voir commentaire de tete). Un chunk non active ne
  ! coute quasiment rien (descripteur ALLOCATABLE non alloue).
  INTEGER, PARAMETER :: ARENA_MAX_CHUNKS = 4000

  ! NOTE FORTRAN : TARGET n'est PAS un attribut valide sur une composante
  ! de type derive (seulement sur une variable/argument declare(e)) --
  ! BUF herite du statut "target" de la variable BLR_ARENA_T qui le
  ! contient (voir ARENA_CARVE ci-dessous, dummy INTENT(INOUT), TARGET) :
  ! c'est a l'APPELANT de declarer sa propre variable BLR_ARENA_T (ou le
  ! type derive qui la contient) avec l'attribut TARGET pour que les
  ! pointeurs renvoyes par ARENA_CARVE restent valides apres l'appel.
  TYPE :: ARENA_CHUNK_T
    COMPLEX(WP), ALLOCATABLE :: BUF(:)
    INTEGER(8) :: USED = 0_8
  END TYPE ARENA_CHUNK_T

  TYPE :: BLR_ARENA_T
    TYPE(ARENA_CHUNK_T) :: CHUNKS(ARENA_MAX_CHUNKS)
    INTEGER :: NCHUNKS = 0
  END TYPE BLR_ARENA_T

CONTAINS

  !-----------------------------------------------------------------------
  ! Decoupe une vue de NELEMS elements COMPLEX(WP) dans l'arena -- BUMP-
  ! POINTER dans le chunk courant s'il reste de la place, sinon active un
  ! nouveau chunk (ou un chunk dedie si NELEMS depasse la taille standard
  ! d'un chunk, ex. un tres gros bloc de repli dense). Ne libere JAMAIS
  ! rien individuellement (voir compromis en tete de fichier) -- seul
  ! ARENA_FREE (fin de factorisation) libere tout d'un coup.
  !-----------------------------------------------------------------------
  ! TARGET sur ARENA (dummy) : necessaire pour que les pointeurs renvoyes
  ! (PTR, et transitivement toute vue BLK%UFAC/VFAC/DENSE construite a
  ! partir de lui) soient valides -- l'appelant DOIT, de son cote,
  ! declarer sa propre variable BLR_ARENA_T (ou le type qui la contient,
  ! ex. MPI_BLR_MATRIX_ARENA_T) avec l'attribut TARGET, sans quoi le
  ! pointeur ne resterait valide que le temps de CET appel (voir la norme
  ! Fortran sur l'association de cible via un argument).
  FUNCTION ARENA_CARVE(ARENA, NELEMS) RESULT(PTR)
    TYPE(BLR_ARENA_T), INTENT(INOUT), TARGET :: ARENA
    INTEGER(8), INTENT(IN)           :: NELEMS
    COMPLEX(WP), POINTER             :: PTR(:)
    INTEGER :: C

    NULLIFY(PTR)
    IF (NELEMS .LE. 0_8) RETURN

    IF (NELEMS .GT. ARENA_CHUNK_ELEMS) THEN
      ! Plus gros qu'un chunk standard -- chunk dedie, taille exacte.
      IF (ARENA%NCHUNKS .GE. ARENA_MAX_CHUNKS) THEN
        WRITE(6,'(A)') 'ERREUR ARENA_CARVE : ARENA_MAX_CHUNKS depasse (augmenter la constante).'
        STOP 1
      END IF
      ARENA%NCHUNKS = ARENA%NCHUNKS + 1
      C = ARENA%NCHUNKS
      ALLOCATE(ARENA%CHUNKS(C)%BUF(NELEMS))
      ARENA%CHUNKS(C)%USED = NELEMS
      PTR => ARENA%CHUNKS(C)%BUF
      RETURN
    END IF

    IF (ARENA%NCHUNKS .EQ. 0) THEN
      ARENA%NCHUNKS = 1
      ALLOCATE(ARENA%CHUNKS(1)%BUF(ARENA_CHUNK_ELEMS))
    END IF

    C = ARENA%NCHUNKS
    IF (ARENA%CHUNKS(C)%USED + NELEMS .GT. ARENA_CHUNK_ELEMS) THEN
      IF (ARENA%NCHUNKS .GE. ARENA_MAX_CHUNKS) THEN
        WRITE(6,'(A)') 'ERREUR ARENA_CARVE : ARENA_MAX_CHUNKS depasse (augmenter la constante).'
        STOP 1
      END IF
      ARENA%NCHUNKS = ARENA%NCHUNKS + 1
      C = ARENA%NCHUNKS
      ALLOCATE(ARENA%CHUNKS(C)%BUF(ARENA_CHUNK_ELEMS))
    END IF

    PTR => ARENA%CHUNKS(C)%BUF(ARENA%CHUNKS(C)%USED+1 : ARENA%CHUNKS(C)%USED+NELEMS)
    ARENA%CHUNKS(C)%USED = ARENA%CHUNKS(C)%USED + NELEMS
  END FUNCTION ARENA_CARVE

  !-----------------------------------------------------------------------
  ! CARVE-OU-REUTILISE : comme ARENA_CARVE, mais reutilise FLAT (un
  ! pointeur PERSISTANT que l'appelant conserve d'un appel a l'autre, ex.
  ! BLR_BLOCK_ARENA_T%UFAC_FLAT) s'il est deja assez grand -- ne carve un
  ! NOUVEAU chunk QUE si FLAT est encore non-associe ou trop petit.
  !
  ! Correctif du gaspillage x28 mesure avec l'arena "pleine factorisation"
  ! en carve-systematique (voir la memoire du projet) : un bloc dont la
  ! taille de stockage est STABLE ou BORNEE d'une etape K a l'autre (ex. le
  ! repli dense, toujours de taille MI*NJ FIXE pour un slot de bloc donne
  ! -- MI/NJ ne dependent que de la partition, jamais de K ; ou le rang
  ! d'un bloc bas-rang, fluctuant mais borne par RMAX) finit par stabiliser
  ! sa capacite apres quelques etapes et REUTILISE ensuite le MEME buffer
  ! indefiniment, au lieu d'en carver un nouveau a CHAQUE etape --
  ! strategie "high-water-mark" : la capacite ne fait QUE croitre (jamais
  ! re-carvee plus petite), simple et sans risque, meme si pas parfaitement
  ! optimale (un bloc dont le rang grimpe une fois puis redescend garde la
  ! capacite maximale). Ne re-fragmente PAS le tas glibc (on reste dans de
  ! gros chunks arena) -- seul le GASPILLAGE (memoire arena non reutilisee)
  ! est reduit.
  !
  ! SECURITE (aliasing) : l'appelant reste responsable de ne JAMAIS
  ! reutiliser FLAT en place quand l'ancien contenu de FLAT doit encore
  ! etre LU par un calcul qui ECRIT simultanement dans le resultat (ex. une
  ! expansion bas-rang -> dense ne doit PAS ecrire dans le buffer dense
  ! reutilise avant d'avoir fini de lire les anciens UFAC/VFAC -- ceux-ci
  ! sont de toute facon des buffers DIFFERENTS de DENSE_FLAT, donc sans
  ! risque ; voir SCHUR_UPDATE_BLR_BLOCK_ARENA/blr_compress_mod.f90 pour le
  ! cas "deja dense -> dense", qui evite le probleme en mettant a jour EN
  ! PLACE via CGEMM, sans jamais relire-puis-ecrire le meme buffer via deux
  ! arguments (dummy) DIFFERENTS qui aliaseraient la meme memoire).
  !
  ! CAP_MAX (optionnel) : PLAFOND de croissance en elements -- quand un
  ! agrandissement est necessaire, la NOUVELLE capacite carvee est
  ! MAX(NELEMS, 2*ancienne_capacite) PLAFONNE a CAP_MAX (croissance
  ! amortie classique, facteur 2, comme un tableau dynamique) plutot que
  ! carvee exactement a NELEMS -- borne le gaspillage cumule (somme des
  ! buffers successivement abandonnes lors d'agrandissements repetes) a
  ! ~2x la taille finale, MEME dans le pire cas.
  !
  ! POURQUOI PAS "carver directement CAP_MAX des le premier appel" (option
  ! plus simple envisagee, ecartee) : un bloc qui tente le bas-rang puis
  ! finit par basculer en dense (rang qui depasse le plafond de
  ! compression, ex. RMAX ~ 80% de la dimension du bloc -- voir
  ! COMPRESS_WORTHWHILE_RATIO, blr_compress_mod.f90) paierait ALORS la
  ! capacite MAXIMALE en pure perte (orpheline des la bascule vers dense,
  ! voir SCHUR_UPDATE_BLR_BLOCK_ARENA), potentiellement PLUS cher que le
  ! bloc dense final lui-meme. La croissance doublante, elle, ne paie
  ! JAMAIS plus que ~2x le besoin REEL a n'importe quel instant, y compris
  ! pour un bloc qui echoue tot.
  !-----------------------------------------------------------------------
  FUNCTION ARENA_CARVE_REUSE(ARENA, FLAT, NELEMS, CAP_MAX) RESULT(PTR)
    TYPE(BLR_ARENA_T), INTENT(INOUT), TARGET  :: ARENA
    COMPLEX(WP), POINTER, INTENT(INOUT)       :: FLAT(:)
    INTEGER(8), INTENT(IN)                    :: NELEMS
    INTEGER(8), INTENT(IN), OPTIONAL          :: CAP_MAX
    COMPLEX(WP), POINTER                      :: PTR(:)
    INTEGER(8) :: NEWCAP, CEILING

    IF (NELEMS .LE. 0_8) THEN
      NULLIFY(PTR)
      RETURN
    END IF
    IF ((.NOT. ASSOCIATED(FLAT)) .OR. (SIZE(FLAT,KIND=8) .LT. NELEMS)) THEN
      IF (.NOT. ASSOCIATED(FLAT)) THEN
        NEWCAP = NELEMS   ! premiere carve : taille exacte, pas de sur-allocation
      ELSE
        NEWCAP = MAX(NELEMS, 2_8*SIZE(FLAT,KIND=8))   ! agrandissement : double (amorti)
      END IF
      CEILING = NELEMS
      IF (PRESENT(CAP_MAX)) CEILING = MAX(CAP_MAX, NELEMS)   ! garde-fou : jamais < besoin immediat
      NEWCAP = MIN(NEWCAP, CEILING)
      FLAT => ARENA_CARVE(ARENA, NEWCAP)
    END IF
    PTR => FLAT(1:NELEMS)
  END FUNCTION ARENA_CARVE_REUSE

  !-----------------------------------------------------------------------
  ! Empreinte memoire ANALYTIQUE de l'arena (octets) -- somme des chunks
  ! REELLEMENT actives (les slots de CHUNKS(:) au-dela de NCHUNKS restent
  ! non alloues, cout nul). INCLUT l'espace gaspille par les anciennes
  ! versions de blocs jamais recuperees (voir compromis en tete de
  ! fichier) -- c'est un chiffre HONNETE, pas juste la donnee utile.
  !-----------------------------------------------------------------------
  FUNCTION ARENA_BYTES(ARENA) RESULT(NBYTES)
    TYPE(BLR_ARENA_T), INTENT(IN) :: ARENA
    INTEGER(8) :: NBYTES
    INTEGER :: C
    NBYTES = 0_8
    DO C = 1, ARENA%NCHUNKS
      IF (ALLOCATED(ARENA%CHUNKS(C)%BUF)) NBYTES = NBYTES + 8_8*INT(SIZE(ARENA%CHUNKS(C)%BUF),8)
    END DO
  END FUNCTION ARENA_BYTES

  !-----------------------------------------------------------------------
  ! DIAGNOSTIC TEMPORAIRE (voir MPI_BLR_MATRIX_ARENA_MEMORY_DIAG,
  ! mpi_blr_matrix_mod.f90) -- somme CHUNKS(:)%USED (octets REELLEMENT
  ! carves, qu'ils soient encore atteignables depuis un bloc ou deja
  ! orphelins suite a un agrandissement -- voir ARENA_CARVE_REUSE), PAR
  ! OPPOSITION a ARENA_BYTES qui compte la taille COMPLETE de chaque chunk
  ! ALLOUE (y compris la partie jamais carvee en fin de chunk). La
  ! difference ARENA_BYTES-ARENA_USED_BYTES isole la fragmentation INTERNE
  ! aux chunks (fin de chunk non comblee) ; la difference
  ! ARENA_USED_BYTES-(USED+SLACK du diagnostic bloc-par-bloc) isole les
  ! buffers orphelins suite a un AGRANDISSEMENT (ARENA_CARVE_REUSE carve un
  ! nouveau buffer PLUS GRAND et abandonne l'ancien, plus petit, sans
  ! aucune trace conservee -- hypothese a confirmer/infirmer avant de
  ! coder un correctif).
  !-----------------------------------------------------------------------
  FUNCTION ARENA_USED_BYTES(ARENA) RESULT(NBYTES)
    TYPE(BLR_ARENA_T), INTENT(IN) :: ARENA
    INTEGER(8) :: NBYTES
    INTEGER :: C
    NBYTES = 0_8
    DO C = 1, ARENA%NCHUNKS
      IF (ALLOCATED(ARENA%CHUNKS(C)%BUF)) NBYTES = NBYTES + 8_8*ARENA%CHUNKS(C)%USED
    END DO
  END FUNCTION ARENA_USED_BYTES

  SUBROUTINE ARENA_FREE(ARENA)
    TYPE(BLR_ARENA_T), INTENT(INOUT) :: ARENA
    INTEGER :: C
    DO C = 1, ARENA%NCHUNKS
      IF (ALLOCATED(ARENA%CHUNKS(C)%BUF)) DEALLOCATE(ARENA%CHUNKS(C)%BUF)
    END DO
    ARENA%NCHUNKS = 0
  END SUBROUTINE ARENA_FREE

END MODULE BLR_ARENA_MOD
