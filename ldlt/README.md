# ldlt/ — LDLᵀ BLR+MPI pour matrices complexes symétriques

Pendant symétrique du LU BLR+MPI avec pivot confiné (`src/mpi_blr_matrix_mod.f90`,
`MPI_BLR_MATRIX_LU_FACTOR_PIVOT`) : matrices **complexes symétriques**
(Aᵀ = A, *pas* hermitiennes — transposée simple partout, jamais de
conjugaison), simple précision, MPI 2D bloc-cyclique + compression BLR.
Trois variantes, comme pour le LU : **classique** (en mémoire), **arena**
(pool mémoire contre la fragmentation de tas) et **out-of-core**.

Le code de `src/` n'est pas modifié : ce répertoire réutilise ses `.o`/`.mod`
déjà construits (même principe que `legacy/`).

## Algorithme

- **Stockage** : triangle inférieur uniquement. Blocs hors-diagonale I>J
  (dense ou bas-rang) ; blocs I<J jamais assemblés ni alloués (`GET_BLOCK`
  n'est appelé que pour I≥J). Blocs **diagonaux en stockage compact**
  (triangle inférieur seul, KW·(KW+1)/2 éléments, format LAPACK `'L'`) :
  dépliés dans un unique espace de travail NB×NB le temps d'une mise à jour
  de Schur ou de `CSYTRF_RK`, et diffusés compacts (moitié du volume).
- **Pivot** : Bunch-Kaufman borné (`CSYTRF_RK`, pivots 1×1 et 2×2) **confiné
  au bloc diagonal**, même philosophie que le pivot partiel confiné du LU (et
  que MUMPS/STRUMPACK). Le bloc diagonal K est factorisé *en entier*, une
  fois toutes ses mises à jour de Schur reçues : aucune mise à jour différée
  ne traverse la compression — c'est ce qui avait fait échouer la première
  tentative LDLt+BLR (`legacy/src/blr_ldlt_mod.f90`, `CLASYF`).
- **Propagation de la permutation symétrique P_K sans décompression** :
  lignes des L(K,J<K) déjà calculés → échange de lignes de `UFAC` ; colonnes
  des A(I>K,K) → échange de colonnes de `VFAC`. Factorisation obtenue :
  Q·A·Qᵀ = L·D·Lᵀ, Q = diag(P_Kᵀ).
- **Étape K** : `CSYTRF_RK` sur (K,K) → IPIV et D_K diffusés à toute la grille ;
  W(I,K) = A(I,K)·P_K·L_KK⁻ᵀ sur la colonne de grille de K, diffusé (déjà
  compressé) le long des lignes ; L(I,K) = W·D_K⁻¹ ; L(J,K)ᵀ diffusé le long
  des colonnes depuis le processus (ligne(J), colonne(J)) — schéma du
  Cholesky ScaLAPACK ; Schur A(I,J) −= W(I,K)·L(J,K)ᵀ sur I≥J seulement.
- **Résolution** : Q·B → descente L (algorithme du LU) → D (local, aucune
  communication) → remontée Lᵀ « left-looking » (une réduction le long d'une
  colonne de grille + une diffusion par étape) → Qᵀ.

La convention de `CSYTRF_RK` (Pᵀ·A·P = L·D·Lᵀ, échanges (k,|IPIV(k)|)
appliqués dans l'ordre k=1..n, bloc 2×2 signalé par IPIV(k)<0, terme
hors-diagonal dans E(k)) a été **vérifiée numériquement** avant d'écrire la
propagation, pas supposée d'après la documentation.

Compression à l'assemblage : **SVD randomisée adaptative par défaut**
(`USE_RSVD=.TRUE.`, comme tout le projet depuis 2026-10-03) ; `USE_RSVD=.FALSE.`
/ `--use-rsvd=0` rétablit `CGESVD`. Tous les tests tournent avec le RSVD.

## Variantes

| Variante | Module | Spécificité |
|---|---|---|
| classique | `mpi_blr_ldlt_mod.f90` | type `MPI_BLR_SYM_T` ; aussi `MPI_BLR_SYM_SAVE`/`_LOAD` (checkpoint) |
| arena | `mpi_blr_ldlt_arena_mod.f90` | blocs hors-diagonale carvés dans `BLR_ARENA_T` ; Schur via `SCHUR_UPDATE_BLR_BLOCK_ARENA` ; opérations de panneau sur une copie ordinaire puis recopie dans le **même** buffer (rang inchangé) |
| out-of-core | `mpi_blr_ldlt_ooc_mod.f90` | matrice, facteur et second membre sur un fichier scratch par processus ; caches LRU `CACHE_A`/`CACHE_D` de `blr_ooc_mod.f90` ; bloc diagonal compact stocké comme un bloc dense NP×1 (aucune modification de `blr_ooc_mod.f90`) ; second membre par paquets (`OOC_RHS_T`, `CONSUME_X_BATCH`) |

## Fichiers

| Fichier | Rôle |
|---|---|
| `src/mpi_blr_ldlt_mod.f90` | version classique + primitives de bloc (échange de colonnes, TRSM transposé, D⁻¹, transposée, produit transposé, stockage compact) + checkpoint |
| `src/mpi_blr_ldlt_arena_mod.f90` | variante arena |
| `src/mpi_blr_ldlt_ooc_mod.f90` | variante out-of-core |
| `src/sym_testmat_mod.f90` | matrice de test complexe symétrique (noyau de Helmholtz sur le cercle + diagonale/anti-diagonale réglables) |
| `test/test_mpi_blr_ldlt.f90` | validation : erreur vs Xexact, référence dense `CGESV` indépendante, reconstruction directe ‖Q·A·Qᵀ − L·D·Lᵀ‖ (grille 1×1) |
| `test/test_mpi_blr_ldlt_arena.f90` | arena vs classique (erreur, écart, mémoire) |
| `test/test_mpi_blr_ldlt_ooc.f90` | OOC vs classique, second membre exact entrée par entrée, cache résident observé |
| `test/bench_ldlt_vs_lu.f90` | LDLᵀ vs LU BLR+MPI sur la même matrice (temps, mémoire, erreur, phases de résolution) |
| `examples/example_blr_mpi_ldlt.f90` | exemple avec les options du LU, y compris sauvegarde/reprise |

## Compiler et lancer

```bash
source /opt/intel/oneapi/setvars.sh
make                     # projet principal d'abord (depuis blr_solver/)
make -C ldlt             # puis le LDLt -> ldlt/bin/
mpirun -np 4 ldlt/bin/example_blr_mpi_ldlt 2 2 4000 200 8 1e-4 --check-error
```

Chaîne gcc : `make -f Makefile.gcc && make -C ldlt -f Makefile.gcc`
(→ `ldlt/build-gcc/`, LAPACK ≥ 3.7 requis pour `CSYTRF_RK`).
`./run_all_tests.sh [intel|gcc]` compile et teste aussi ce répertoire
(9 tests LDLᵀ : classique sur 3 régimes de pivots, exemple, benchmark,
arena, OOC avec évictions réelles, cycle sauvegarde → reprise).

Options de l'exemple : `--use-move-alloc --eps --ratio --use-rsvd
--check-error --check-residual --profile-usage --dump-structure --nrhs-batch`,
sauvegarde/reprise `--save-matrix --save-factored --save-rhs --save-solution
--load-matrix --load-factored --load-rhs --checkpoint-dir --ntailbloc` (même
sens que pour le LU, voir `--help`), plus les paramètres de la matrice de
démo `--wavenumber --strong --diag`. Les CSV de `--profile-usage` /
`--dump-structure` ont le même format que ceux du LU
(`examples/plot_usage_profile.py`, `examples/plot_blr_structure.py`).

## Validation (2026-10-03)

- Grilles 1×1, 1×2, 2×1, 2×2, 1×4, 4×1, 2×3, 3×2 ; trois régimes de pivots
  (tout 2×2, mélange 1×1/2×2, tout 1×1 sans échange) ; avec et sans
  compression.
- Sans compression : reconstruction de la factorisation ~1e-7 (N=200) à
  9e-7 (N=1000), erreur sur la solution au niveau de `CGESV` dense ou
  meilleure.
- Avec TOL=1e-4 : écart à `CGESV` ~1e-5 à 8e-5 selon N, identique (à
  l'arrondi près) d'une grille à l'autre.
- Arena et OOC : même erreur que la version classique, écart entre variantes
  ~1.5e-5 (les tirages aléatoires du RSVD diffèrent d'un run à l'autre), même
  nombre d'échanges de pivot. OOC : résultat **identique** quel que soit le
  budget (1, 3, 64 Mo) — l'éviction ne change pas les calculs.
- Reprise : `--load-factored` / `--load-matrix` redonnent bit pour bit
  l'erreur et le résidu du run d'origine ; checkpoint incompatible détecté sur
  chaque rang, arrêt propre de tous les rangs.
- Chaînes Intel (ifx+MKL+Intel MPI) et gcc (gfortran+LAPACK de
  référence+OpenMPI) : mêmes résultats.

## Mesures LDLᵀ vs LU (`bench_ldlt_vs_lu`, grille 2×2, OMP_NUM_THREADS=2, TOL=1e-4 sauf mention, RSVD)

Ratios LDLᵀ/LU (< 1 = LDLᵀ meilleur), stockage diagonal compact :

| Cas | Assemblage | Factorisation | Résolution | Mémoire facteurs | Erreur LU / LDLᵀ |
|---|---|---|---|---|---|
| N=8000, NB=400 | ×0.48 | ×0.81 | ×0.79 | ×0.50 | 1.9e-4 / 1.9e-4 |
| N=8000, NB=400, TOL=0 (dense) | ×0.55 | ×0.64 | ×1.29 | ×0.50 | 1.9e-5 / 1.5e-5 |
| N=16000, NB=200 | ×0.51 | ×0.46 | ×0.87 | ×0.50 | 1.0e-4 / 1.2e-4 |
| N=16000, NB=400 | ×0.52 | ×0.82 | ×0.84 | ×0.50 | 2.7e-4 / 2.7e-4 |
| N=16000, NB=800 | ×0.53 | ×1.09 | ×1.22 | ×0.50 | 3.3e-4 / 3.2e-4 |
| N=32000, NB=400 | ×0.57 | ×0.45 | ×0.73 | ×0.50 | 2.1e-4 / 2.5e-4 |
| N=32000, NB=800 | ×0.54 | ×0.93 | ×1.21 | ×0.50 | 5.3e-4 / 5.2e-4 |
| N=48000, NB=800 | ×0.54 | ×0.90 | ×0.86 | ×0.50 | 4.6e-4 / 4.6e-4 |

Temps absolus de factorisation (s), pour comparer **chaque solveur à son
meilleur NB** :

| N | LU NB=400 | LU NB=800 | LDLᵀ NB=200 | LDLᵀ NB=400 | LDLᵀ NB=800 |
|---|---|---|---|---|---|
| 16000 | 0.479 | **0.403** | 1.020 | **0.392** | 0.440 |
| 32000 | 5.917 | **1.728** | — | 2.674 | **1.613** |

Lecture honnête :
- **Mémoire ×0.50 partout** (stockage compact des blocs diagonaux ; avant
  compaction, seulement ×0.89 à TOL=1e-4 : les blocs diagonaux pleins
  dominaient).
- **Assemblage ×0.5** (moitié des blocs à évaluer et compresser).
- **Factorisation : ×0.45 à NB fixe petit/moyen, mais seulement ×0.90–0.97
  au meilleur NB de chaque solveur.** À grand NB, le coût est dominé par les
  blocs diagonaux ; `CSYTRF_RK` (Bunch-Kaufman, moins parallèle que
  `CGETRF`) et le dépliage/repliage du stockage compact mangent l'essentiel
  du gain de flops.
- **Résolution** : proche du LU (×0.73 à ×1.29) ; un premier écart ×2.2
  mesuré sur un seul run s'est révélé être du bruit (re-mesuré ×1.02–1.07).

## Limites connues, assumées

- **Gain de factorisation modeste au NB optimal** (voir ci-dessus) : le
  levier serait un `CSYTRF_RK` plus efficace sur gros blocs (variante
  blocked/parallèle) ou une mise à jour de Schur diagonale limitée au
  triangle inférieur — non fait.
- **Remontée Lᵀ** : copie pleine largeur `XFULL` (MLOC × NRHS) des inconnues
  déjà calculées — mémoire RHS ×PC pendant cette phase, bornée par
  `--nrhs-batch` (et par paquet en OOC).
- **OOC** : le budget `MAX_MB_PER_PROCESS` s'applique **à chacun des deux
  caches** (diagonal et hors-diagonale), comme dans le LU OOC — le résident
  peut atteindre 2× le budget. Coût fixe d'environ 2.5–3 ms par étape K
  (écritures disque des blocs de panneau évincés), indépendant du budget :
  ×1.6–2.2 sur la factorisation à N=8000 sur ce poste, du même ordre que le
  LU OOC (0.236 s contre 0.217 s à NB=200). Mêmes mises en garde de
  dimensionnement que le LU OOC (thrashing sous l'ensemble de travail).
- **Arena** : `ARENA_BYTES` compte la capacité des chunks réservés (~16 Mo
  au minimum par processus dans les tests), pas seulement les données
  utiles ; le gain attendu porte sur la RSS réelle à grande échelle (comme
  pour le LU arena), non mesuré ici.
- **Déséquilibre de charge sur grille carrée** (hérité de la distribution
  bloc-cyclique, aussi présent dans le LU) : les blocs diagonaux et leur
  `CSYTRF_RK` tombent tous sur les processus (r,r) — mémoire max/processus
  ×0.55 au lieu de ×0.50 en dense.
- **Matrices mal conditionnées** : comme pour le LU, l'erreur de compression
  (~TOL) est amplifiée par le conditionnement. Exemple : la variante « pivots
  mixtes » du test (`STRONG=1 DIAG=0.1`, `CGESV` dense seul ~6e-4) dépasse
  1e-3 à TOL=1e-4 alors que la factorisation est exacte au niveau de TOL —
  il faut resserrer TOL (1e-6 dans `run_all_tests.sh`).
- **Pas de checkpoint pour les variantes arena/OOC** (seulement la version
  classique, via l'exemple).
