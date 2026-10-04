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
  mis à jour **directement dans le stockage compact** (une `CGEMV` par
  colonne quand l'un des facteurs de Schur est bas-rang), dépliés dans un
  unique espace de travail NB×NB seulement pour la factorisation du bloc ou
  une mise à jour dense, et diffusés compacts (moitié du volume).
- **Pivot** : Bunch-Kaufman (`CSYTRF` puis `CSYCONVF`, qui convertit la
  sortie au format de `CSYTRF_RK` ; pivots 1×1 et 2×2) **confiné
  au bloc diagonal**, même philosophie que le pivot partiel confiné du LU (et
  que MUMPS/STRUMPACK). Le bloc diagonal K est factorisé *en entier*, une
  fois toutes ses mises à jour de Schur reçues : aucune mise à jour différée
  ne traverse la compression — c'est ce qui avait fait échouer la première
  tentative LDLt+BLR (`legacy/src/blr_ldlt_mod.f90`, `CLASYF`).
- **Propagation de la permutation symétrique P_K sans décompression** :
  lignes des L(K,J<K) déjà calculés → échange de lignes de `UFAC` ; colonnes
  des A(I>K,K) → échange de colonnes de `VFAC`. Factorisation obtenue :
  Q·A·Qᵀ = L·D·Lᵀ, Q = diag(P_Kᵀ).
- **Étape K** : Bunch-Kaufman sur (K,K) → IPIV et D_K diffusés à toute la grille ;
  W(I,K) = A(I,K)·P_K·L_KK⁻ᵀ sur la colonne de grille de K, diffusé (déjà
  compressé) le long des lignes ; L(I,K) = W·D_K⁻¹ ; L(J,K)ᵀ diffusé le long
  des colonnes depuis le processus (ligne(J), colonne(J)) — schéma du
  Cholesky ScaLAPACK ; Schur A(I,J) −= W(I,K)·L(J,K)ᵀ sur I≥J seulement.
- **Distribution** : blocs-lignes cycliques (ligne de grille MOD(I−1,PR)) ;
  blocs-colonnes en cyclique **décalé** par défaut (`DIAG_SPREAD=.TRUE.` à
  l'assemblage, les trois variantes) : J → MOD(J−1 + ⌊(J−1)/L⌋, PC), avec
  L = PPCM(PR,PC). Avec le cyclique pur, le bloc diagonal K tombe sur
  (MOD(K−1,PR), MOD(K−1,PC)) : si PGCD(PR,PC) = g > 1, seuls PR·PC/g
  processus possèdent des blocs diagonaux (2 sur 4 en 2×2, 4 sur 8 en
  2×4). Le décalage d'une colonne à chaque période les fait passer par les
  PR·PC processus (vérifié jusqu'à 16×18), sans changer le nombre de
  blocs-colonnes par colonne de grille. Identique au cyclique pur si g = 1
  (dont PR=1 ou PC=1). Le checkpoint enregistre la distribution (numéro
  magique 7311 + CMAP).
- **Anticipation** (`LOOKAHEAD=.TRUE.` par défaut, version classique) : à
  l'étape K, la colonne K+1 reçoit sa mise à jour de Schur en premier ; le
  propriétaire de (K+1,K+1) le factorise aussitôt et lance des diffusions
  **non bloquantes** (`MPI_IBCAST`) de IPIV/D à la grille et du bloc compact
  à sa colonne de grille, qui traite tout de suite son panneau (échanges de
  colonnes + TRSM) ; le reste du Schur de l'étape K recouvre le tout, et
  l'étape K+1 commence par un `MPI_WAITALL`. Résultat numérique identique
  (mêmes opérations, seul l'ordre change).
- **Résolution** : Q·B → descente L → D (local, aucune communication) →
  remontée Lᵀ → Qᵀ. Deux schémas pour la descente et la remontée
  (`DESCENT_MODE` de `MPI_BLR_SYM_LDLT_SOLVE_2D`, choix automatique par
  défaut) :
  - **mode 1, « colonne »** : à l'étape K, X_K est rassemblé sur toute la
    largeur des seconds membres sur la colonne de grille propriétaire de K,
    qui calcule seule les contributions puis les disperse (descente) ; la
    remontée garde une copie pleine largeur du second membre (XFULL) et
    réduit vers le propriétaire du bloc diagonal. Communications
    proportionnelles à NRHS, une seule colonne de grille calcule par étape ;
  - **mode 2, « distribué »** : X_K (ma tranche de seconds membres) est
    diffusé le long de la colonne de grille, la colonne L(·,K) de mes lignes
    est diffusée **compressée** le long de la ligne de grille (un paquet :
    en-têtes + données), et chaque processus applique L (descente) ou Lᵀ
    (remontée, suivie d'une réduction le long de la colonne de grille) à
    **ses** seconds membres. Toutes les colonnes de grille calculent, le
    volume de L ne dépend pas de NRHS, et la remontée n'a plus de copie
    pleine largeur ;
  - **automatique** : mode 2 si NRHS > 1,5 × W, où W est la largeur
    effective de L (stockage des blocs I>J divisé par la somme de leurs
    lignes, 4,5 à 7,8 sur la matrice de démo).

Conventions **vérifiées numériquement**, pas supposées d'après la
documentation : Pᵀ·A·P = L·D·Lᵀ, échanges (k,|IPIV(k)|) appliqués dans
l'ordre k=1..n ; bloc 2×2 en (k,k+1) ⇔ terme hors-diagonal E(k) ≠ 0. Ce
dernier critère est le seul commun aux deux formats : `CSYTRF_RK` marque un
2×2 par IPIV(k), IPIV(k+1) < 0, alors que `CSYTRF` + `CSYCONVF` donne
IPIV(k) = k, IPIV(k+1) < 0 (un premier essai fondé sur le signe de IPIV(k)
faussait toute la factorisation, détecté par le test de reconstruction).

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
| `test/bench_ldlt_vs_lu.f90` | LDLᵀ vs LU BLR+MPI sur la même matrice (temps, mémoire, erreur, découpage de la factorisation et phases de résolution) |
| `test/bench_diag_kernels.f90` | micro-benchmark des noyaux du bloc diagonal (`CGETRF`, variantes Bunch-Kaufman et Aasen, stockage compact, `CGEMM`/`CGEMMT`) — Intel seulement |
| `test/sweep_driver.f90` | pilote de `tools/sweep_incore.sh` |
| `examples/example_blr_mpi_ldlt.f90` | exemple avec les options du LU, y compris sauvegarde/reprise |

## Compiler et lancer

```bash
source /opt/intel/oneapi/setvars.sh
make                     # projet principal d'abord (depuis blr_solver/)
make -C ldlt             # puis le LDLt -> ldlt/bin/
mpirun -np 4 ldlt/bin/example_blr_mpi_ldlt 2 2 4000 200 8 1e-4 --check-error
```

Chaîne gcc : `make -f Makefile.gcc && make -C ldlt -f Makefile.gcc`
(→ `ldlt/build-gcc/`, LAPACK ≥ 3.7 requis pour `CSYCONVF`).
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

## Optimisations du 2026-10-03 (mesurées)

Micro-benchmark des noyaux du bloc diagonal (1 thread, bloc 800, rang 8) :
- `CSYTRF_RK` prend 0.65× le temps de `CGETRF` (deux fois moins de calculs,
  mais 94 contre 121 GF/s) ; `CSYTRF` + `CSYCONVF` 0.61× (0.82× à 2 threads,
  contre 0.95× pour `_RK`) → adopté. Aasen (`CSYTRF_AA`) n'est pas meilleur.
- Mise à jour de Schur d'un bloc diagonal compact : dépliage + `CGEMM` +
  repliage 0.77 ms ; `CGEMV` par colonne directement dans le stockage compact
  0.21 ms, résultat identique → adopté.

Découpage de la factorisation (N=16000, grille 2×4 × 1 thread, s, MAX sur les
processus), LDLᵀ avant → après :

| Poste | NB=400 | NB=800 |
|---|---|---|
| factorisation des blocs diagonaux | 0.019 → 0.018 | 0.067 → 0.055 |
| Schur des blocs diagonaux | 0.072 → 0.022 | 0.150 → 0.027 |
| attentes MPI | 0.232 → 0.176 | 0.508 → 0.364 |
| **factorisation totale** (LU : 0.39 / 0.56) | **0.310 → 0.251** | **0.552 → 0.400** |

L'attente MPI reste le premier poste (60 à 90 % de la boucle, LU compris) :
les blocs diagonaux sont traités par une partie seulement des processus
pendant que les autres attendent — c'est la cible de l'anticipation
(« lookahead »). La remontée Lᵀ en pleine largeur ramène la résolution
LDLᵀ au niveau du LU ou en dessous : ×0.93–0.99 (2×4, NB=800, contre ×1.41
avant), ×0.80–0.86 (1×8, NB=600, contre ×1.9) ; la remontée ne coûte plus
que 0.016–0.045 s, la descente L domine désormais.

## Anticipation (2026-10-04, mesurée)

Factorisation LDLᵀ (s, N=16000, TOL=1e-4, RSVD, 1 thread, 3 runs
alternés ; erreurs **identiques au chiffre près** avec et sans) :

| Grille, NB | sans | avec |
|---|---|---|
| 2×2, NB=400 | 0.351–0.402 | 0.336–0.364 |
| 2×2, NB=800 | 0.372–0.462 | 0.375–0.466 |
| 2×4, NB=400 | 0.409–0.430 | 0.401–0.437 |
| 2×4, NB=800 | 0.490–0.516 | 0.496–0.506 |

Gain de 0 à 7 % seulement sur cette machine (4 cœurs physiques, 8 avec
hyperthreading : la grille 2×4 y est sursouscrite et les rangs en attente
active volent des cycles aux rangs qui calculent). Un découpage fin des
attentes, par processus (instrumentation temporaire, grille 2×2), montre
pourquoi : les diffusions elles-mêmes (W le long des lignes, Lᵀ le long des
colonnes) ne coûtent que 0.05–0.09 s ; l'essentiel de l'attente vient du
**déséquilibre de charge** — les processus (0,0) et (1,1), qui possèdent
tous les blocs diagonaux (factorisation de Bunch-Kaufman + Schur dense),
calculent 0.20–0.25 s contre 0.06–0.17 s pour (0,1) et (1,0) — combiné à la
synchronisation de chaque étape par les diffusions bloquantes du panneau.
L'anticipation de profondeur 1 retire la factorisation du bloc diagonal du
chemin critique mais ne rééquilibre rien. Sur un vrai réseau (latence
inter-nœuds), le recouvrement des diffusions devrait davantage compter ;
l'option est conservée par défaut (jamais mesurablement plus lente). La
section suivante corrige ce diagnostic : une fois la charge rééquilibrée, le
temps ne baisse toujours pas.

## Mesures LDLᵀ vs LU avant ces optimisations (`bench_ldlt_vs_lu`, grille 2×2, OMP_NUM_THREADS=2, TOL=1e-4 sauf mention, RSVD)

Ratios LDLᵀ/LU (< 1 = LDLᵀ meilleur), stockage diagonal compact (dépliage à
chaque mise à jour, `CSYTRF_RK`, remontée left-looking) :

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
  au meilleur NB de chaque solveur** avant optimisation. L'explication
  initiale (« `CSYTRF_RK` moins efficace que `CGETRF` ») était fausse en
  temps absolu : le découpage a montré que le surcoût venait surtout du
  dépliage du stockage compact à chaque mise à jour de Schur, et des attentes
  MPI qu'il allongeait (voir la section précédente).
- **Résolution** : de ×0.73 à ×1.29 avant optimisation, et jusqu'à ×1.7–1.9
  sur les grilles larges (2×4, 1×8) — corrigé par la remontée pleine largeur.

## Répartition des blocs diagonaux (2026-10-04, mesurée)

N=16000, TOL=1e-4, RSVD, 1 thread, anticipation active ; facteurs, max par
processus (Mo), cyclique pur → décalé (identique au chiffre près d'un run
à l'autre) :

| Grille, NB | mémoire max/processus | factorisation (s, 3 runs) |
|---|---|---|
| 2×2, NB=400 | 14.8 → **8.9** | 0.263–0.268 → 0.263–0.278 |
| 2×2, NB=800 | 26.1 → **14.0** | 0.353–0.419 → 0.382–0.430 |
| 2×4, NB=400 | 7.5 → **4.4** | 0.269–0.351 → 0.296–0.366 |
| 2×4, NB=800 | 13.1 → **8.2** | 0.471–0.504 → 0.482–0.517 |
| 2×2, N=40000, NB=1000 | 82.9 → **45.1** | 1.49–1.71 → 1.40–1.61 |

La mémoire maximale par processus baisse de 37 à 46 %, et le calcul est
rééquilibré (instrumentation temporaire, 2×2, NB=800 : 0.20 / 0.06 / 0.06 /
0.20 s → 0.11 / 0.11 / 0.12 / 0.12 s). Le temps de factorisation, lui, ne
change pas au bruit près : il est borné par la **chaîne séquentielle des
factorisations de blocs diagonaux** (chaque bloc K+1 attend le Schur de
l'étape K, qui attend la factorisation de K) — leur somme vaut à elle seule
≈ 0.23 s sur une boucle de 0.35 s, et le Schur hors-diagonale, compressé,
est trop léger (≈ 0.05 s en tout) pour la recouvrir. Ce plancher croît comme
N·NB² : c'est ce qui fait préférer NB=400–600 au balayage de paramètres.

## Résolution distribuée (2026-10-04, mesurée)

Résolution LDLᵀ (s), N=16000, NB=400, TOL=1e-4 (W=4,5), 1 thread, 2 runs :

| Grille | NRHS | mode 1 (colonne) | mode 2 (distribué) | gain |
|---|---|---|---|---|
| 2×2 | 1 | 0,019–0,022 | 0,031–0,032 | ×0,65 |
| 2×2 | 8 | 0,044–0,045 | 0,038–0,043 | ≈ |
| 2×2 | 32 | 0,102–0,105 | 0,061–0,064 | ×1,7 |
| 2×2 | 128 | 0,315–0,325 | 0,140–0,144 | ×2,3 |
| 2×2 | 512 | 1,38–1,41 | 0,52–0,54 | ×2,6 |
| 2×2 | 1024 | 2,65–2,79 | 0,88–1,04 | ×2,9 |
| 2×4 | 1 | 0,037–0,038 | 0,065–0,068 | ×0,55 |
| 2×4 | 128 | 0,34–0,35 | 0,16–0,19 | ×2,0 |
| 2×4 | 512 | 1,35 | 0,41–0,46 | ×3,1 |
| 2×4 | 1024 | 2,70–2,81 | 0,85–0,89 | ×3,2 |

Descente seule à NRHS=1024 : 1,86–1,96 → 0,46–0,53 s ; remontée :
0,77–0,97 → 0,34–0,48 s. Le choix automatique prend le mode 1 en dessous
d'environ 7 seconds membres sur cette matrice. La répartition des blocs
diagonaux (`DIAG_SPREAD`) ne change pas le temps de résolution, quel que
soit NRHS (mesuré jusqu'à 1024, écarts dans le bruit).

## Limites connues, assumées

- **Attentes MPI dominantes** (60 à 90 % de la factorisation, LU compris) :
  ni l'anticipation de profondeur 1 (0 à 7 %, version classique seulement)
  ni la répartition des blocs diagonaux (calcul rééquilibré, temps inchangé)
  ne les font baisser, car elles reflètent la chaîne séquentielle des
  factorisations denses de blocs diagonaux (voir « Répartition des blocs
  diagonaux »). Leviers restants : factorisation diagonale plus rapide
  (multithreadée chez le propriétaire, ou NB plus petit), ou découpage du
  bloc diagonal en sous-blocs pour raccourcir la chaîne.
- **Remontée Lᵀ en mode 1** : copie pleine largeur `XFULL` (MLOC × NRHS) du
  second membre — mémoire RHS ×PC pendant cette phase, bornée par
  `--nrhs-batch` (et par paquet en OOC). Le mode 2 n'en a pas besoin.
- **Résolution distribuée (mode 2)** : version classique seulement ; les
  variantes arena et OOC gardent le mode 1.
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
- **Déséquilibre de charge sur grille carrée** : corrigé par la
  distribution décalée des colonnes (voir « Répartition des blocs
  diagonaux »), pour le LDLᵀ comme pour le LU.
- **Matrices mal conditionnées** : comme pour le LU, l'erreur de compression
  (~TOL) est amplifiée par le conditionnement. Exemple : la variante « pivots
  mixtes » du test (`STRONG=1 DIAG=0.1`, `CGESV` dense seul ~6e-4) dépasse
  1e-3 à TOL=1e-4 alors que la factorisation est exacte au niveau de TOL —
  il faut resserrer TOL (1e-6 dans `run_all_tests.sh`).
- **Pas de checkpoint pour les variantes arena/OOC** (seulement la version
  classique, via l'exemple).
