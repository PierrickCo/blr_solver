!===============================================================================
! Analyse de la ligne de commande -- SPLIT_KV (parsing --cle=valeur) et
! PRINT_USAGE (message d'aide), utilisees par example_blr_mpi_pivot_user.f90.
! Aucun etat, aucune dependance au solveur -- reutilisable telle quelle.
!===============================================================================
MODULE CLI_OPTIONS_MOD
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SPLIT_KV, PRINT_USAGE

CONTAINS

  ! Coupe "--cle=valeur" en KEYOUT="cle"/VALOUT="valeur" (le "--" initial
  ! est deja verifie par l'appelant). Sans "=" (option booleenne nue),
  ! VALOUT="1" -- equivalent a "--cle=1".
  SUBROUTINE SPLIT_KV(ARGSTR, KEYOUT, VALOUT)
    CHARACTER(LEN=*), INTENT(IN)  :: ARGSTR
    CHARACTER(LEN=*), INTENT(OUT) :: KEYOUT, VALOUT
    INTEGER :: EQPOS
    EQPOS = INDEX(ARGSTR, '=')
    IF (EQPOS .EQ. 0) THEN
      KEYOUT = ARGSTR(3:)
      VALOUT = '1'
    ELSE
      KEYOUT = ARGSTR(3:EQPOS-1)
      VALOUT = ARGSTR(EQPOS+1:)
    END IF
  END SUBROUTINE SPLIT_KV

  ! Message d'aide -- affiche sur --help/-h, ligne de commande invalide,
  ! ou option/argument non reconnu.
  SUBROUTINE PRINT_USAGE()
    WRITE(6,'(A)') 'Usage : mpirun -np PR*PC ./example_blr_mpi_pivot_user PR PC [N] [NB] [NRHS] [TOL] [OPTIONS...]'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  PR PC          grille de processus MPI (obligatoire, PR*PC = nb de processus lances)'
    WRITE(6,'(A)') '  N NB NRHS TOL  positionnels optionnels, DANS CET ORDRE (defauts 4000 200 4 1e-4)'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Options nommees (apres les positionnels ci-dessus, dans N''IMPORTE QUEL ORDRE ;'
    WRITE(6,'(A)') '--cle=valeur, ou --cle seul pour une option booleenne, equivalent a --cle=1) :'
    WRITE(6,'(A)') '  --n=ENTIER              equivalent au positionnel N (utile pour ne changer que NRHS/TOL...)'
    WRITE(6,'(A)') '  --nb=ENTIER             equivalent au positionnel NB'
    WRITE(6,'(A)') '  --nrhs=ENTIER           equivalent au positionnel NRHS'
    WRITE(6,'(A)') '  --tol=REEL              equivalent au positionnel TOL'
    WRITE(6,'(A)') '  --ratio=REEL            (defaut 0.333) COMPRESS_WORTHWHILE_RATIO'
    WRITE(6,'(A)') '  --use-rsvd=0|1          (defaut 0) SVD randomisee adaptative au lieu de CGESVD'
    WRITE(6,'(A)') '  --check-residual=0|1    (defaut 0) verification ||A*X-B||/||B||'
    WRITE(6,'(A)') '  --profile-usage=0|1     (defaut 0) trace RSS/CPU par etape (usage_profile_rank<N>.csv)'
    WRITE(6,'(A)') '  --dump-structure=0|1    (defaut 0) export structure BLR (blr_structure_<etape>_rank<N>.csv)'
    WRITE(6,'(A)') '  --save-matrix=0|1       (defaut 0) sauvegarde la matrice assemblee'
    WRITE(6,'(A)') '  --save-rhs=0|1          (defaut 0) sauvegarde le second membre'
    WRITE(6,'(A)') '  --save-factored=0|1     (defaut 0) sauvegarde la matrice factorisee'
    WRITE(6,'(A)') '  --save-solution=0|1     (defaut 0) sauvegarde la solution'
    WRITE(6,'(A)') '  --load-matrix=0|1       (defaut 0) relit la matrice sauvegardee (saute l''assemblage)'
    WRITE(6,'(A)') '  --load-factored=0|1     (defaut 0) relit la matrice factorisee (saute assemblage+factorisation,'
    WRITE(6,'(A)') '                          prioritaire sur --load-matrix)'
    WRITE(6,'(A)') '  --load-rhs=0|1          (defaut 0) relit le second membre sauvegarde'
    WRITE(6,'(A)') '  --checkpoint-dir=CHEMIN (defaut .) repertoire des fichiers de checkpoint (doit deja exister)'
    WRITE(6,'(A)') '  --nrhs-batch=ENTIER     (defaut 0 = pas de decoupage) colonnes RHS traitees par lot'
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'Exemple :'
    WRITE(6,'(A)') '  mpirun -np 4 ./example_blr_mpi_pivot_user 2 2 4000 200 8 1e-4 \'
    WRITE(6,'(A)') '      --check-residual --save-factored --checkpoint-dir=ckpt'
  END SUBROUTINE PRINT_USAGE

END MODULE CLI_OPTIONS_MOD
