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
| **BLR+MPI** | Solveur MPI dont le stockage local par processus reste compressé de bout en bout (jamais dense hors blocs diagonaux), panneaux diffusés déjà compressés | `mpi_blr_matrix_mod.f90`, `test_mpi_blr_matrix.f90` (sans pivot), `test_mpi_blr_matrix_pivot.f90` (avec pivot confiné au bloc diagonal), `test_mpi_blr_matrix_pivot_arena.f90` (même algorithme, allocateur arena pour le stockage bas-rang), `test_mpi_blr_matrix_pivot_ooc.f90` (même algorithme, out-of-core) |

**Distribution des blocs** (BLR+MPI, LU et LDLᵀ) : blocs-lignes
cycliques, blocs-colonnes en cyclique **décalé** par défaut (`DIAG_SPREAD`
des assemblages, `BLR_COL_OWNER` dans `mpi_blr_matrix_mod.f90`) : les blocs
diagonaux passent par tous les processus, au lieu de PR·PC/PGCD(PR,PC) en
cyclique pur (2 sur 4 en 2×2). Mesuré (N=16000, 1 thread) : mémoire max
par processus −37 à −47 % (LU 2×2, NB=800 : 52,1 → 27,8 Mo), mémoire
totale et temps inchangés ; identique au cyclique pur si PGCD(PR,PC)=1.
Les checkpoints enregistrent la distribution et les anciens restent
relisibles (voir [`examples/README.md`](examples/README.md)).

**Recompression incrémentale du Schur hors diagonale** (LU et LDLᵀ,
depuis 2026-11-03). La mise à jour d'un bloc bas-rang par un produit
bas-rang, suivie de sa recompression, est le coût dominant de la
factorisation.

- **Invariant gardé :** un bloc issu d'une recompression garde `UFAC` à
  colonnes orthonormées et `VFAC = diag(SV)·Wᵀ`, avec W orthonormée (champs
  `ORTHO` et `SV` de `BLR_BLOCK_T`).
- **Mise à jour suivante** (`FAST_RECOMPRESS_ORTHO`) : la mise à jour est
  projetée sur ces bases ; seul le résidu est orthonormalisé, puis vient
  la SVD du petit noyau. La matrice et la troncature sont les mêmes que
  la recompression complète ; seuls les arrondis diffèrent.
- **Mesuré** sur une matrice BEM à N = 40 230 (LDLᵀ, NB = 702,
  1 processus) : factorisation 90 → 67-69 s (−24 %), mémoire inchangée.
- **Comparaison :** `ORTHO_RECOMP = .FALSE.` rétablit la recompression
  complète.
- **Diagnostic :** `SCHUR_PROF = .TRUE.` active le profil par cas, étape de
  recompression et rangs (variables `SPROF_*` de `blr_compress_mod.f90`,
  séquentiel). `T_BREAKDOWN_MIN` de `MPI_BLR_SYM_LDLT_FACTOR` donne le
  minimum par processus de chaque poste, à comparer au maximum
  `T_BREAKDOWN`.

**Mémoire de la factorisation** (LU et LDLᵀ, depuis 2026-11-05) :
- **Tampons de panneau libérés à chaque étape.** Les tampons des
  blocs-lignes et blocs-colonnes déjà traités ne sont plus gardés jusqu'à
  la fin de la factorisation.
- **Espace de travail dimensionné pour la factorisation.**
  `BLR_WORKSPACE_ALLOC(WS, NB, FACT_WORKSPACE_RANK(NB))` alloue des tableaux
  NB × rang maximal au lieu de NB × NB.
- **Mesuré** sur la matrice BEM à N = 40 230 (2×2) : pic par processus
  487-498 → 342 Mo (−30 %), résultat identique au bit près.
- **Diagnostic :** `FACT_MEM_PROF = .TRUE.` imprime le profil mémoire par
  étape du LDLᵀ (lignes `[MEMFACT]`).

**BLR+MPI** est la pièce la plus aboutie du projet : c'est le seul module où
le stockage *distribué* lui-même reste compressé (tous les modules MPI plus
anciens gardent une matrice locale dense, la compression n'y accélère que le
calcul). Il existe en quatre variantes, toutes validées :
- **sans pivot** (`MPI_BLR_MATRIX_LU_FACTOR_NOPIV`) — matrices à diagonale
  dominante uniquement.
- **avec pivot partiel confiné au bloc diagonal**
  (`MPI_BLR_MATRIX_LU_FACTOR_PIVOT`) — même philosophie que MUMPS/STRUMPACK
  (pivotage local au front, pas global) ; permuter les lignes d'un bloc
  compressé ne coûte que `O(rang)` (seul `UFAC` est touché), jamais de
  décompression.
- **avec pivot, stockage arena** (`MPI_BLR_MATRIX_LU_FACTOR_PIVOT_ARENA`,
  type `MPI_BLR_MATRIX_ARENA_T`) — même algorithme numérique que la variante
  ci-dessus (résultats identiques), mais les blocs bas-rang hors-diagonale
  sont carvés dans un pool mémoire dédié (`blr_arena_mod.f90`) au lieu d'un
  `ALLOCATE`/`DEALLOCATE` individuel à chaque recompression, pour contenir
  la fragmentation de tas observée à grande échelle (voir la mémoire du
  projet pour l'investigation complète). Compromis MESURÉ, pas supposé : à
  N=50000/NB=500, RSS réelle inférieure de 11.5 % à la référence classique
  et factorisation 30 % plus rapide, contre un assemblage 19 % plus lent et
  une empreinte théorique plus élevée (le prix du compromis mémoire de
  l'arena, jamais récupérée pour les anciennes versions de blocs) — l'écart
  se creuse avec N, l'arena n'a pas encore d'avantage net en RSS absolue à
  plus petite échelle (N=20000).
- **avec pivot, out-of-core** (`MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC`, type
  `MPI_BLR_MATRIX_OOC_T`) — même algorithme numérique, mais la matrice, les
  facteurs LU et le second membre vivent sur un fichier scratch par
  processus (jamais en RAM en permanence), la mémoire résidente étant
  bornée par un paramètre utilisateur (`MAX_MB_PER_PROCESS`) via deux
  caches LRU (voir `examples/README.md` pour le détail). Compromis MESURÉ :
  le budget mémoire doit rester au-dessus de l'ensemble de travail réel du
  processus, sous peine de faire thrasher le cache — à N=20000, un budget
  de 200 Mo perdait 68 % du temps de factorisation en accès disque, contre
  seulement quelques % à 800 Mo (RSS réelle quasi identique dans les deux
  cas, ~130 Mo) ; le compromis n'est donc pas « moins de RAM contre plus de
  temps » de façon linéaire, mais un vrai seuil à ne pas franchir.

## Limites connues, assumées

- **Complexité asymptotique** : le découpage en blocs est **plat** (taille
  fixe, pas de clustering hiérarchique type HODLR/H-matrix) — la
  factorisation reste `O(N³)` en ordre, avec une constante bien plus petite
  grâce à la compression (mesuré : BLR+MPI ×2-2.4 plus rapide et ×2.7-4.5
  moins gourmand en mémoire que la référence MPI dense, à N=10000-20000).
  Le coût dominant est la recompression du Schur hors diagonale : 83-86 %
  du temps de factorisation, et 94 % du Schur sur une matrice BEM à
  N = 40 230, avant la recompression incrémentale qui en retire environ un
  quart.
- **Pivotage confiné, pas global** : le pivot ne cherche que dans les lignes
  du bloc diagonal courant — plus faible qu'un pivotage global, mais c'est
  exactement le compromis utilisé en production par MUMPS/STRUMPACK
  (pivotage confiné au front).
- **Compression à l'assemblage par SVD randomisée adaptative, par défaut**
  (depuis 2026-10-03, `COMPRESS_BLOCK_RSVD_ADAPTIVE`) : mesurée ×9-12 plus
  rapide que `CGESVD` sur l'assemblage complet à N=20000, même rang et même
  erreur. Garantie heuristique (pas de preuve formelle), avec une zone
  morte étroite près du plafond `RATIO` où elle replie en dense là où
  `CGESVD` aurait compressé (jamais de corruption). `--use-rsvd=0` dans les
  exemples, ou `USE_RSVD=.FALSE.` dans l'API, rétablit `CGESVD`.
- **Pas de sparsité générale** : matrices denses uniquement (la compression
  BLR exploite la structure basse-rang, pas la sparsité).
- **Code mis de côté** : une première tentative LDLt+BLR a été abandonnée
  (voir `legacy/README.md`) — `CLASYF` utilise un schéma de mise à jour
  différée incompatible avec l'approche « factoriser dense, compresser
  après ». Un LDLᵀ BLR+MPI fonctionnel existe désormais dans
  [`ldlt/`](ldlt/README.md) (matrices complexes symétriques, Bunch-Kaufman
  borné confiné au bloc diagonal, qui contourne ce problème).

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
à chaque push/pull request sur `master` (ou `main`).

## Structure du dépôt

```
blr_solver/
├── src/               modules vivants (les 5 references + BLR+MPI)
├── test/              drivers de validation et de benchmark, un par module
├── examples/          programmes courts, un par paradigme -- prise en main
│                      rapide de l'API (voir examples/README.md)
├── tools/             sweep_incore.sh (recherche des meilleurs parametres,
│                      voir CLUSTER.md), analyze_deps.py
├── ldlt/              LDLt BLR+MPI pour matrices complexes symetriques
│                      (Bunch-Kaufman confine) -- voir ldlt/README.md
├── legacy/            code mis de cote (LDLt, ancienne lignee MPI a
│                      stockage dense) -- voir legacy/README.md
├── .github/workflows/ CI (chaine gcc, build + run_all_tests.sh a chaque push)
├── run_all_tests.sh   compile + lance le jeu de tests cure, resume PASS/FAIL
├── Makefile           chaine Intel oneAPI (ifx+MKL+Intel MPI), reference
├── Makefile.gcc       chaine alternative (gfortran+BLAS/LAPACK+OpenMPI)
├── Makefile.ifort2019 chaine alternative (ifort/Parallel Studio XE 2019 --
│                      non validee en conditions reelles, voir INSTALL.md)
├── INSTALL.md
├── CLUSTER.md
└── LICENSE            MIT
```

## Licence

MIT, voir [`LICENSE`](LICENSE). Auteur : Dr Pierrick Cordel.
