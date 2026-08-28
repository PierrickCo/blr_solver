# Exemples canoniques

Cinq programmes courts, un par paradigme du projet (voir le tableau de
[`README.md`](../README.md)), montrant l'enchainement minimal d'appels
pour chaque solveur -- sans l'instrumentation de mesure ni la validation
croisee contre une reference independante des drivers de `test/`/`bench/`
(qui restent la source de verite pour la correction, voir
`../run_all_tests.sh`).

| Fichier | Reference | Paradigme |
|---|---|---|
| `example_blr_sequentiel.f90` | A/B | Sequentiel/OpenMP, matrice compressee (BLR) de bout en bout |
| `example_dense_sequentiel.f90` | C | Sequentiel/OpenMP, LU dense par blocs, pivot complet |
| `example_mpi_dense.f90` | MPI | MPI, LU dense distribue, pivot, resolution 2D |
| `example_blr_mpi_sans_pivot.f90` | BLR+MPI | MPI, stockage BLR distribue, sans pivot (diagonale dominante) |
| `example_blr_mpi_pivot.f90` | BLR+MPI | MPI, stockage BLR distribue, pivot confine au bloc diagonal |

## Compiler

```bash
cd ..                    # depuis examples/, remonter a blr_solver/
make examples             # chaine Intel (ifx+MKL+Intel MPI)
# ou
make -f Makefile.gcc examples   # chaine gcc+BLAS/LAPACK+OpenMPI
```

Chaque exemple se construit comme binaire a la racine de `blr_solver/`
(ou `build-gcc/` pour la chaine gcc), en reutilisant directement les
`.o`/`.mod` deja compiles par `make`/`make -f Makefile.gcc` -- compiler
le projet principal d'abord.

## Lancer

```bash
./example_blr_sequentiel [N] [NB] [NRHS] [TOL]
./example_dense_sequentiel [N] [NB]
mpirun -np 4 ./example_mpi_dense 2 2 [N] [NB] [NRHS]
mpirun -np 4 ./example_blr_mpi_sans_pivot 2 2 [N] [NB] [NRHS] [TOL]
mpirun -np 4 ./example_blr_mpi_pivot 2 2 [N] [NB] [NRHS] [TOL]
```

Chaque exemple assemble une matrice de demonstration (formule analytique,
pas de fichier d'entree), factorise, resout un second membre arbitraire,
et affiche les premieres valeurs de la solution -- point de depart a
adapter (callback `GET_BLOCK`, second membre) pour un probleme reel.

## Barres de progression et taux de compression (BLR+MPI avec pivot)

`example_blr_mpi_pivot.f90` active `SHOW_PROGRESS=.TRUE.` sur les 3 etapes
couteuses (assemblage, factorisation, resolution) et affiche le taux de
compression (`N_COMPRESSED/N_OFFDIAG`) apres l'assemblage et a la fin de
la factorisation. Ces deux compteurs sont deja des sommes GLOBALES
(agregees par `MPI_ALLREDUCE` a l'interieur de `mpi_blr_matrix_mod.f90`) --
identiques sur tous les processus, affiches une seule fois (rang 0).

La barre de progression, elle, n'a la meme garantie que pour la
factorisation et la resolution : ces deux boucles sont globalement
synchronisees (un `MPI_BCAST` a chaque etape K), donc `K/NBLK` est un
pourcentage d'avancement EXACT. L'assemblage, en revanche, n'a aucune
synchronisation entre processus (chacun assemble ses blocs locaux de son
cote) -- la barre affichee par le rang 0 n'est donc qu'une approximation
(mais fidele, le decoupage bloc-cyclique equilibrant la charge entre
processus). `MPI_ASSEMBLE_BLR_MATRIX` (etape d'assemblage) est PARTAGEE
par les 2 variantes -- `example_blr_mpi_sans_pivot.f90` beneficie donc
deja de sa barre de progression. En revanche, `SHOW_PROGRESS` n'a ete
ajoute qu'aux routines de factorisation/resolution PIVOT
(`MPI_BLR_MATRIX_LU_FACTOR_PIVOT`/`_SOLVE_PIVOT_2D`) -- les equivalents
NOPIV (`MPI_BLR_MATRIX_LU_FACTOR_NOPIV`/`_SOLVE_2D`) n'ont pas encore cet
argument optionnel.
