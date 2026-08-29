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

## Empreinte memoire theorique (solveurs MPI)

Les 3 exemples MPI (`example_mpi_dense.f90`, `example_blr_mpi_sans_pivot.f90`,
`example_blr_mpi_pivot.f90`) affichent, a 3 moments cles, l'empreinte
memoire THEORIQUE (calculee a partir des tailles de tableaux -- PAS
mesuree via l'OS, voir `GET_RSS_KB`/`CLUSTER.md` pour la RSS reelle) :

1. **Matrice** (dense ou compressee BLR) juste apres l'assemblage.
2. **Facteurs LU** juste apres la factorisation -- pour la reference MPI
   (dense), c'est EXACTEMENT la meme valeur qu'a l'assemblage (facteurs
   ecrasent la matrice en place, aucune compression ne peut changer la
   taille) ; pour BLR+MPI, cette valeur peut DIFFERER de celle de
   l'assemblage a cause du remplissage ("fill-in") introduit par la mise
   a jour de Schur -- observe concretement sur `example_blr_mpi_pivot`
   (matrice 1.9 Mo -> facteurs 2.0 Mo a N=1000/NB=100/grille 2x2).
3. **Second membre B** juste apres son allocation/remplissage.

Chaque ligne affiche DEUX nombres (`REPORT_MEMORY_MB`, nouvelle routine
partagee de `mpi_grid_mod.f90`) : le **total** agrege sur toute la grille
(`MPI_ALLREDUCE`/`MPI_SUM`, utile pour la capacite globale du cluster) et
le **maximum par-processus** (`MPI_ALLREDUCE`/`MPI_MAX`, c'est LUI qui
determine si le job tient dans la RAM de chaque noeud -- la contrainte
qui bloque reellement un job, pas le total). Pour BLR+MPI, la memoire de
la matrice/des facteurs reutilise directement `MPI_BLR_MATRIX_MEMORY_
BYTES` (deja existante) ; pour la reference MPI dense et le second
membre (les deux BLR+MPI inclus), c'est une simple formule
`lignes_locales * colonnes_locales * 8 octets` (`COMPLEX(WP)` simple
precision).

**Comparaison a la RSS reelle** (`example_blr_mpi_pivot.f90` uniquement,
via `GET_RSS_KB()`/`bench_utils_mod.f90`, meme mecanisme `REPORT_MEMORY_
MB`) : la RSS mesuree est TOUJOURS superieure a la valeur theorique
(elle inclut le runtime MPI/MKL charge en memoire, les workspaces
SVD/QR transitoires par thread, l'executable lui-meme...), mais le
RATIO n'est PAS constant -- mesure sur ce projet : ~65x a N=2000 (l'
overhead fixe du runtime domine a cette echelle), ~2.5-3x a N=20000
(les donnees elles-memes dominent). C'est ce second regime (grande
echelle) qui sert de reference pour detecter une vraie anomalie (fuite,
copie transitoire oubliee) -- pas un ratio absolu independant de N.

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

## Investigation de l'ecart RSS-vs-theorique (`example_blr_mpi_pivot.f90`)

A N=20000/NB=500/grille 2x2, la factorisation ajoute ~49 Mo/processus de
RSS reelle alors que le stockage BLR (theorique) ne grossit que de ~10 Mo
(remplissage/fill-in) -- un ecart de ~39 Mo non explique par les donnees
elles-memes. Quatre hypotheses ciblees ont ete testees, TOUTES ecartees
ou negligeables :

1. **Cache interne MKL** (`mkl_free_buffers()`, appele apres la
   factorisation) -- ~15% de l'ecart recupere, pas la cause principale.
2. **Fragmentation en haut de tas** (`malloc_trim(0)`, glibc) -- ~0%
   recupere (il ne peut agir que sur la memoire libre CONTIGUE en haut du
   tas, pas les trous au milieu).
3. **Reallocation implicite Fortran** (`ABLK = NEWBLK` declenche une
   reallocation+copie supplementaire meme si `NEWBLK` est deja a la bonne
   taille) -- corrige via `MOVE_ALLOC` (voir `SCHUR_UPDATE_BLR_BLOCK`,
   `blr_compress_mod.f90`, argument optionnel `USE_MOVE_ALLOC`,
   propage par `MPI_BLR_MATRIX_LU_FACTOR_PIVOT`) -- resultat numerique
   IDENTIQUE (verifie bit-a-bit), mais ~0% de gain memoire mesure.
4. **Seuil `mmap` glibc** (`mallopt(M_MMAP_THRESHOLD, ...)`, testé puis
   retire -- voir l'historique git) -- l'appel reussit (code retour 1)
   mais ~0% d'effet, tres probablement parce que le runtime Fortran
   d'Intel (`ifx`) ne route pas les `ALLOCATE`/`DEALLOCATE` par le
   `malloc()` de glibc que ce reglage controle (non confirme -- non
   teste sous gfortran, qui route typiquement `ALLOCATE` directement par
   `malloc()`).

**Conclusion** : un tracage RSS par etape K (`SHOW_RSS_TRACE`, argument
optionnel de `MPI_BLR_MATRIX_LU_FACTOR_PIVOT`, actif dans cet exemple)
montre une croissance CONTINUE qui RALENTIT et se stabilise vers la fin
de la factorisation -- une signature de fragmentation proportionnelle au
VOLUME de travail (le nombre de blocs mis a jour a chaque etape diminue
mecaniquement quand la sous-matrice traINante retrecit), pas d'une fuite
qui grossirait sans lien avec le travail restant. Aucun allocateur
standard (glibc, jemalloc, tcmalloc) ne peut compacter (deplacer) de la
memoire encore utilisee -- seule une reecriture manuelle (arena/pool
dedie sur toute la duree de la factorisation, pas juste une etape --
le stockage des blocs doit survivre au-dela d'une seule etape K)
offrirait une garantie reelle, avec son propre compromis (memoire
gaspillee par les anciennes versions de blocs jamais recuperees) --
chantier volontairement laisse de cote pour l'instant, le rapport
effort/gain n'etant pas demontre.

`USE_MOVE_ALLOC` est conserve (defaut `.FALSE.`, aucun changement de
comportement) car c'est une bonne pratique defensive independante du
gain memoire observe ici. `SHOW_RSS_TRACE` est conserve comme outil de
diagnostic permanent (defaut `.FALSE.`, cout nul).
