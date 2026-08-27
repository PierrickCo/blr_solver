!===============================================================================
! ETAPE 6b : produit matrice-vecteur y = A*x avec A distribuee en
! bloc-cyclique 2D (module MPI_GRID_MOD, etape 6a). Objectif : valider les
! DEUX communicateurs ROW_COMM/COL_COMM sur une operation simple (un seul
! aller-retour de communication) avant de les reutiliser dans la
! factorisation distribuee, bien plus complexe (etapes 6c-6d).
!
! ALGORITHME (schema standard bloc-cyclique 2D, le meme principe -- diffuser
! le long des colonnes, reduire le long des lignes -- que ScaLAPACK PxGEMV) :
!
!   1) x est conceptuellement decoupe PAR COLONNES (meme bloc NB que les
!      colonnes de A), et fourni en entree UNIQUEMENT par les processus
!      MYROW=0 (X_LOCAL de taille NLOC = NUMROC(N,NB,MYCOL,PC) sur ces
!      processus -- valeur ignoree ailleurs).
!   2) DIFFUSION le long de chaque colonne de la grille (COL_COMM, racine
!      = le processus MYROW=0 de cette colonne) : apres cet appel, TOUS
!      les processus (quel que soit MYROW) possedent le morceau de x
!      correspondant a LEUR PROPRE MYCOL.
!   3) PRODUIT LOCAL : chaque processus calcule Y_PARTIAL = A_LOCAL*X_LOCAL
!      (CGEMV), une contribution PARTIELLE aux lignes qu'il possede (ne
!      voit qu'UN SEUL bloc de colonnes parmi les PC qui contribuent a ces
!      memes lignes).
!   4) REDUCTION le long de chaque ligne de la grille (ROW_COMM, somme) :
!      les PC contributions partielles pour les memes lignes globales
!      (une par MYCOL, 0..PC-1) sont sommees -- ALLREDUCE (pas REDUCE) :
!      tous les processus de cette ligne finissent avec le meme resultat,
!      pas besoin d'une diffusion separee ensuite.
!
! Y_LOCAL en sortie est donc distribue PAR LIGNES SEULEMENT (NUMROC sur
! PR), REPLIQUE sur tous les MYCOL -- voir GATHER_ROW_VECTOR dans
! MPI_GRID_MOD pour le rassembler en vue de la verification.
!===============================================================================
MODULE MPI_MATVEC_MOD
  USE BLR_KINDS_MOD, ONLY: WP
  USE MPI_GRID_MOD
  IMPLICIT NONE

CONTAINS

  SUBROUTINE MPI_MATVEC_2D(GRID, MLOC, NLOC, A_LOCAL, X_LOCAL, Y_LOCAL)
    INCLUDE 'mpif.h'
    TYPE(GRID_T), INTENT(IN)   :: GRID
    INTEGER, INTENT(IN)        :: MLOC, NLOC
    COMPLEX(WP), INTENT(IN)    :: A_LOCAL(MAX(MLOC,1),MAX(NLOC,1))
    COMPLEX(WP), INTENT(INOUT) :: X_LOCAL(MAX(NLOC,1))  ! valide en entree
                                                         ! seulement si MYROW=0
    COMPLEX(WP), INTENT(OUT)   :: Y_LOCAL(MAX(MLOC,1))

    COMPLEX(WP), PARAMETER :: CONE = (1.0_WP,0.0_WP), CZERO = (0.0_WP,0.0_WP)
    COMPLEX(WP) :: Y_PARTIAL(MAX(MLOC,1))
    INTEGER :: IERR

    ! 2) Diffusion de x le long de la colonne de grille (COL_COMM).
    ! MYROW=0 est toujours le rang 0 DANS COL_COMM : COL_COMM a ete cree
    ! par MPI_COMM_SPLIT(...,COLOR=MYCOL,KEY=MYROW) (voir GRID_INIT), et
    ! la norme MPI garantit un ordre de rangs croissant avec KEY -- donc
    ! KEY=MYROW=0 => rang 0 dans le communicateur resultant.
    IF (NLOC .GT. 0) THEN
      CALL MPI_BCAST(X_LOCAL, NLOC, MPI_COMPLEX, 0, GRID%COL_COMM, IERR)
    END IF

    ! 3) Produit local (BLAS quick-return correct si MLOC=0 ou NLOC=0).
    Y_PARTIAL = CZERO
    IF (MLOC.GT.0 .AND. NLOC.GT.0) THEN
      CALL CGEMV('No transpose', MLOC, NLOC, CONE, A_LOCAL, MLOC, &
                  X_LOCAL, 1, CZERO, Y_PARTIAL, 1)
    END IF

    ! 4) Reduction le long de la ligne de grille (ROW_COMM).
    IF (MLOC .GT. 0) THEN
      CALL MPI_ALLREDUCE(Y_PARTIAL, Y_LOCAL, MLOC, MPI_COMPLEX, MPI_SUM, &
                          GRID%ROW_COMM, IERR)
    END IF
  END SUBROUTINE MPI_MATVEC_2D

END MODULE MPI_MATVEC_MOD
