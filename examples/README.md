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

## Barres de progression et taux de compression

Les 5 exemples activent `SHOW_PROGRESS=.TRUE.` sur chacune de leurs etapes
couteuses (assemblage, factorisation, resolution) -- couverture complete
des 5 solveurs de reference, pas seulement BLR+MPI. Seule exception :
`BLOCK_LU_SOLVE` (reference C, dense sequentiel) n'a pas de barre --
elle resout en un seul appel `CTRSM` par triangle, sans boucle par bloc a
instrumenter (deja tres rapide, `O(N^2)`, pas `O(N^3)` comme la
factorisation).

Les deux exemples BLR+MPI (`example_blr_mpi_sans_pivot.f90`,
`example_blr_mpi_pivot.f90`) affichent en plus le taux de compression
(`N_COMPRESSED/N_OFFDIAG`) apres l'assemblage et a la fin de la
factorisation. Ces deux compteurs sont deja des sommes GLOBALES (agregees
par `MPI_ALLREDUCE` a l'interieur de `mpi_blr_matrix_mod.f90`) --
identiques sur tous les processus, affiches une seule fois (rang 0).

**Deux mecanismes d'affichage distincts** (voir `src/progress_bar_mod.f90`) :
- **Sequentiel/OpenMP** (`example_blr_sequentiel.f90`,
  `example_dense_sequentiel.f90`) : `PRINT_PROGRESS_BAR`, une barre a
  rafraichissement en place (retour chariot `CHAR(13)`) -- fonctionne
  correctement, attache a un vrai terminal.
- **MPI** (les 3 autres exemples, et les modules `mpi_lu_pivot_mod.f90`/
  `mpi_lu_solve_2d_mod.f90`/`mpi_blr_matrix_mod.f90` sous-jacents) :
  `PRINT_PROGRESS_MILESTONE`, un affichage par PALIERS (une ligne complete
  tous les 5%) -- `mpirun` relaie le stdout de chaque processus via un
  pipe, pas un pseudo-terminal, donc le retour chariot n'y est PAS
  interprete comme "ecraser la ligne precedente" : la barre en place
  finissait par polluer la sortie d'une ligne par mise a jour (bug
  constate et corrige) -- d'ou l'affichage par paliers, robuste quel que
  soit le lanceur.

Pour les 3 boucles de factorisation/resolution MPI (reference MPI et
BLR+MPI, avec ou sans pivot), le pourcentage affiche est EXACT : chaque
etape K est globalement synchronisee (un `MPI_BCAST`/diffusion a chaque
K), donc `K/NBLK` (ou l'equivalent pour la resolution, 2 sous-etapes)
represente le meme etat sur tous les processus. Seul l'assemblage BLR+MPI
fait exception : aucune synchronisation entre processus a cette etape
(chacun assemble ses blocs locaux de son cote) -- l'affichage du rang 0
n'est donc qu'une approximation (mais fidele, le decoupage bloc-cyclique
equilibrant la charge entre processus).
