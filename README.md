# blr_solver

Solveur direct **Block Low-Rank (BLR)** écrit *from scratch*, en Fortran, pour
des systèmes linéaires denses complexes en simple précision — LU avec
pivotage, séquentiel, OpenMP et MPI. Construit pas à pas dans un but
d'apprentissage : comprendre en profondeur comment fonctionne un solveur BLR
(à la MUMPS/STRUMPACK) plutôt que d'utiliser ces bibliothèques comme des
boîtes noires.

Ce n'est **pas** un concurrent de MUMPS ou STRUMPACK — c'est un projet
pédagogique, sans support de la sparsité générale ni de nombreuses
fonctionnalités de production. Il reste néanmoins un solveur réel, validé
numériquement à chaque étape, avec des résultats de performance mesurés
(pas supposés) sur chaque optimisation.

## Ce qui existe

Cinq solveurs de référence, tous validés contre une solution exacte
indépendante et comparables entre eux (même recette de matrice de test) :

| Référence | Description | Fichiers |
|---|---|---|
| **A** | BLR compressé dès l'assemblage, jamais dense (« pipeline en 5 étapes ») | `blr_matrix_mod.f90`, `bench_blr_pipeline.f90` (`TOL>0`) |
| **B** | Même pipeline, sans compression (`TOL=0`) — isole le gain dû à la compression seule | `bench_blr_pipeline.f90` (`TOL=0`) |
| **MPI** | LU dense distribué avec pivot, résolution 2D complète (aucune réplication résiduelle) | `mpi_lu_pivot_mod.f90`, `mpi_lu_solve_2d_mod.f90`, `bench_mpi.f90` |
| **C** | LU dense monolithique, sans découpage en blocs (référence MKL pure) | `block_lu_mod.f90`, `bench_seq.f90` |
| **BLR+MPI** | Solveur MPI dont le stockage local par processus reste compressé de bout en bout (jamais dense hors blocs diagonaux), panneaux diffusés déjà compressés | `mpi_blr_matrix_mod.f90`, `test_mpi_blr_matrix.f90` (sans pivot), `test_mpi_blr_matrix_pivot.f90` (avec pivot confiné au bloc diagonal) |

**BLR+MPI** est la pièce la plus aboutie du projet : c'est le seul module où
le stockage *distribué* lui-même reste compressé (tous les modules MPI plus
anciens gardent une matrice locale dense, la compression n'y accélère que le
calcul). Il existe en deux variantes :
- **sans pivot** (`MPI_BLR_MATRIX_LU_FACTOR_NOPIV`) — matrices à diagonale
  dominante uniquement.
- **avec pivot partiel confiné au bloc diagonal**
  (`MPI_BLR_MATRIX_LU_FACTOR_PIVOT`) — même philosophie que MUMPS/STRUMPACK
  (pivotage local au front, pas global) ; permuter les lignes d'un bloc
  compressé ne coûte que `O(rang)` (seul `UFAC` est touché), jamais de
  décompression.

## Limites connues, assumées

- **Complexité asymptotique** : le découpage en blocs est **plat** (taille
  fixe, pas de clustering hiérarchique type HODLR/H-matrix) — la
  factorisation reste `O(N³)` en ordre, avec une constante bien plus petite
  grâce à la compression (mesuré : BLR+MPI ×2-2.4 plus rapide et ×2.7-4.5
  moins gourmand en mémoire que la référence MPI dense, à N=10000-20000).
  Voir la mémoire du projet pour la localisation précise du coût
  (`FAST_RECOMPRESS`, 83-86% du temps de factorisation).
- **Pivotage confiné, pas global** : le pivot ne cherche que dans les lignes
  du bloc diagonal courant — plus faible qu'un pivotage global, mais c'est
  exactement le compromis utilisé en production par MUMPS/STRUMPACK
  (pivotage confiné au front).
- **Pas de sparsité générale** : matrices denses uniquement (la compression
  BLR exploite la structure basse-rang, pas la sparsité).
- **Code mis de côté** : une tentative LDLt+BLR a été abandonnée (voir
  `legacy/README.md`) — `CLASYF` utilise un schéma de mise à jour différée
  incompatible avec l'approche « factoriser dense, compresser après » qui a
  fonctionné pour LU.

## Démarrer

Voir [`INSTALL.md`](INSTALL.md) pour compiler, et [`CLUSTER.md`](CLUSTER.md)
pour un déploiement sur un vrai cluster (affinité CPU/OpenMP, scripts de
balayage de performance, exemple Slurm). Pour une prise en main rapide de
l'API de chaque solveur, voir [`examples/`](examples/README.md) (un
programme court par paradigme, sans l'instrumentation des drivers de test).

```bash
source /opt/intel/oneapi/setvars.sh
make
./run_all_tests.sh
```

Le dépôt est vérifié en continu par une CI GitHub Actions (chaîne
gcc+BLAS/LAPACK+OpenMPI, voir [`.github/workflows/ci.yml`](.github/workflows/ci.yml))
à chaque push/pull request sur `main`.

## Structure du dépôt

```
blr_solver/
├── src/               modules vivants (les 5 references + BLR+MPI)
├── test/              drivers de validation et de benchmark, un par module
├── examples/          programmes courts, un par paradigme -- prise en main
│                      rapide de l'API (voir examples/README.md)
├── bench/             scripts de balayage reutilisables (taille de bloc,
│                      tolerance, scalabilite) -- voir CLUSTER.md
├── legacy/            code mis de cote (LDLt, ancienne lignee MPI a
│                      stockage dense) -- voir legacy/README.md
├── .github/workflows/ CI (chaine gcc, build + run_all_tests.sh a chaque push)
├── run_all_tests.sh   compile + lance le jeu de tests cure, resume PASS/FAIL
├── Makefile           chaine Intel (ifx+MKL+Intel MPI), reference
├── Makefile.gcc       chaine alternative (gfortran+BLAS/LAPACK+OpenMPI)
├── INSTALL.md
├── CLUSTER.md
└── LICENSE            MIT
```

## Licence

MIT, voir [`LICENSE`](LICENSE). Auteur : Dr Pierrick Cordel.
