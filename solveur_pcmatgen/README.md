# solveur_pcmatgen

Solveur BLR (block low-rank) + MPI, **LU avec pivot** ou **LDLᵀ** au choix à
l'exécution. La matrice A et le second membre B sont générés par `PCMATGEN`,
le générateur de test de ScaLAPACK. Dossier autonome : il contient tout le
code nécessaire et rien d'autre.

## Compilation

```bash
source /opt/intel/oneapi/setvars.sh
make                       # -> ./solveur_blr
make OPENMP=0 CHECK=0      # -> ./solveur_blr_noomp_nocheck (production MPI pur)
```

| Option     | Effet                                                                |
|------------|----------------------------------------------------------------------|
| `OPENMP=0` | sans OpenMP, MKL séquentiel (recommandé en MPI pur)                  |
| `CHECK=0`  | sans vérification des bornes de tableaux (plus rapide)               |
| `FAST=1`   | `-O3` + jeu d'instructions de la machine (`FAST_ARCH`, défaut `-xHost`) |

Chaque combinaison a son propre exécutable et son dossier `build*`.
`make clean` supprime tous les produits de compilation.

## Exécution

```bash
mpirun -np PR*PC ./solveur_blr PR PC [N] [NB] [NRHS] [TOL] [--option=valeur ...]
mpirun -np 4 ./solveur_blr 2 2 4000 200 8 1e-4 --solver=ldlt --check-residual
./solveur_blr --help       # liste complète des options
```

Options principales (booléennes : `--option` seul vaut `--option=1`) :

| Option | Défaut | Rôle |
|---|---|---|
| `--solver=lu\|ldlt` | `lu` | LU avec pivot (A générale) ou LDLᵀ Bunch-Kaufman (A symétrique) |
| `--use-rsvd` | 1 | SVD randomisée adaptative (0 : SVD complète) |
| `--ortho-recomp` | 1 | recompression incrémentale du complément de Schur |
| `--diag-spread` | 1 | distribution décalée des blocs-colonnes |
| `--lookahead` | 1 | anticipation d'une étape de factorisation |
| `--panel-pack` | 1 | diffusion groupée des panneaux |
| `--descent-mode=0\|1\|2` | 0 | résolution : automatique, colonne, distribuée |
| `--nrhs-batch=K` | 0 | seconds membres traités par lots de K colonnes |
| `--check-residual` | 0 | affiche ‖A·X − B‖ / ‖B‖ |
| `--profile-usage`, `--trace-fact`, `--dump-work`, `--dump-structure` | 0 | diagnostics (fichiers par processus) |
| `--malloc-trim` | 1 | rend la mémoire libérée au système à chaque étape |
| `--save-matrix/-factored/-rhs/-solution`, `--load-matrix/-factored/-rhs`, `--checkpoint-dir`, `--reload-solution` | 0 / `.` | sauvegarde et reprise |

Le ratio de compression (`COMPRESS_WORTHWHILE_RATIO` = 1/3) et les graines de
génération sont fixés dans le code.

## Génération de A et B

`app/matrix_gen_mod.f90` :

- **A** : `PCMATGEN`, diagonale dominante, graine 1 ; générale (`'N'`) pour le
  LU, symétrique (`'S'`) pour le LDLᵀ.
- **B** : `PCMATGEN`, matrice générale N × NRHS, graine 999999, générée
  **directement** (pas de B = A·X) : la solution exacte n'est pas connue, on
  vérifie avec `--check-residual`.

Chaque processus appelle `PCMATGEN` sur une grille virtuelle 1×1 et ne génère
que ses propres blocs : le résultat est indépendant de la grille MPI.

Une matrice aléatoire de `PCMATGEN` n'a pas de structure bas-rang : tous les
blocs restent denses (0 % compressés). Le code sert à tester le solveur sur
une matrice non structurée, pas à mesurer le gain de la compression.

## Contenu

```
Makefile
app/        solveur_blr.f90 (programme), matrix_gen_mod (A et B), options_mod,
            checkpoint_io_mod, profiling_mod, residual_mod
lib/        bibliothèque BLR+MPI (copie allégée de ../src et ../ldlt/src :
            sans les variantes arena, fixedcap, out-of-core ni LU sans pivot)
scalapack/  pcmatgen.f (copie de scalapack_ooc_ref) et scalapack_randgen.f90
            (générateur aléatoire interne, absent de MKL)
```

`scalapack_randgen.f90` réimplémente le générateur aléatoire interne de
ScaLAPACK : les valeurs sont reproductibles mais pas identiques au bit près à
celles d'un vrai ScaLAPACK avec la même graine.

Dépendances externes : ifx, Intel MPI et MKL (dont ScaLAPACK/BLACS, qui
fournissent `NUMROC`, `ICEIL` et `PXERBLA` utilisés par `PCMATGEN`).
