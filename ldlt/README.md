# ldlt/ — LDLᵀ BLR+MPI pour matrices complexes symétriques

Pendant symétrique du LU BLR+MPI avec pivot confiné (`src/mpi_blr_matrix_mod.f90`,
`MPI_BLR_MATRIX_LU_FACTOR_PIVOT`) : matrices **complexes symétriques**
(Aᵀ = A, *pas* hermitiennes — transposée simple partout, jamais de
conjugaison), simple précision, MPI 2D bloc-cyclique + compression BLR.

Le code de `src/` n'est pas modifié : ce répertoire réutilise ses `.o`/`.mod`
déjà construits (même principe que `legacy/`).

## Algorithme

- **Stockage** : triangle inférieur uniquement — blocs diagonaux denses +
  blocs I>J (dense ou bas-rang). Les blocs I<J ne sont jamais assemblés ni
  alloués ; `GET_BLOCK` n'est appelé que pour I≥J.
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

## Fichiers

| Fichier | Rôle |
|---|---|
| `src/mpi_blr_ldlt_mod.f90` | type `MPI_BLR_SYM_T`, `MPI_ASSEMBLE_BLR_SYM`, `MPI_BLR_SYM_LDLT_FACTOR`, `MPI_BLR_SYM_LDLT_SOLVE_2D`, primitives de bloc (échange de colonnes, TRSM transposé, D⁻¹, transposée, produit transposé) |
| `src/sym_testmat_mod.f90` | matrice de test complexe symétrique (noyau de Helmholtz sur le cercle + diagonale/anti-diagonale réglables) |
| `test/test_mpi_blr_ldlt.f90` | validation : erreur vs Xexact, référence dense `CGESV` indépendante, reconstruction directe ‖Q·A·Qᵀ − L·D·Lᵀ‖ (grille 1×1) |
| `test/bench_ldlt_vs_lu.f90` | LDLᵀ vs LU BLR+MPI sur la même matrice (temps, mémoire, erreur) |
| `examples/example_blr_mpi_ldlt.f90` | exemple avec les options du LU (hors sauvegarde/reprise) |

## Compiler et lancer

```bash
source /opt/intel/oneapi/setvars.sh
make                     # projet principal d'abord (depuis blr_solver/)
make -C ldlt             # puis le LDLt -> ldlt/bin/
mpirun -np 4 ldlt/bin/example_blr_mpi_ldlt 2 2 4000 200 8 1e-4 --check-error
```

Chaîne gcc : `make -f Makefile.gcc && make -C ldlt -f Makefile.gcc`
(→ `ldlt/build-gcc/`, LAPACK ≥ 3.7 requis pour `CSYTRF_RK`).
`./run_all_tests.sh [intel|gcc]` compile et teste aussi ce répertoire.

Options de l'exemple : `--use-move-alloc --eps --ratio --use-rsvd
--check-error --check-residual --profile-usage --dump-structure --nrhs-batch`
(même sens que pour le LU, voir `--help`), plus les paramètres de la matrice
de démo `--wavenumber --strong --diag`. Les CSV de `--profile-usage` /
`--dump-structure` ont le même format que ceux du LU
(`examples/plot_usage_profile.py`, `examples/plot_blr_structure.py`).

## Validation (2026-10-03)

- Grilles 1×1, 1×2, 2×1, 2×2, 1×4, 4×1, 2×3, 3×2 ; N jusqu'à 2000 dans le
  test ; avec et sans compression ; trois régimes de pivots (tout 2×2,
  mélange 1×1/2×2, tout 1×1 sans échange).
- Sans compression : reconstruction de la factorisation ~1e-7 (N=200) à
  4e-7 (N=600), erreur sur la solution au niveau de `CGESV` dense ou
  meilleure.
- Avec TOL=1e-4 : écart à `CGESV` ~1e-5 à 8e-5 selon N, identique (à
  l'arrondi près) d'une grille à l'autre.
- Chaînes Intel (ifx+MKL+Intel MPI) et gcc (gfortran+LAPACK de
  référence+OpenMPI) : mêmes résultats.

## Mesures LDLᵀ vs LU (`bench_ldlt_vs_lu`, grille 2×2, OMP_NUM_THREADS=2, matrice par défaut, `--use-rsvd`)

| Cas | Assemblage | Factorisation | Résolution | Mémoire facteurs (total) | Erreur LU / LDLᵀ |
|---|---|---|---|---|---|
| N=8000, NB=400, TOL=0 (dense) | ×0.53 | ×0.61 | ×1.20 | ×0.525 | 1.9e-5 / 1.5e-5 |
| N=8000, NB=400, TOL=1e-6 | ×0.52 | ×0.71 | ×0.69 | ×0.87 | 2.1e-5 / 2.1e-5 |
| N=8000, NB=400, TOL=1e-4 | ×0.53 | ×0.82 | ×0.95 | ×0.89 | 1.9e-4 / 1.9e-4 |
| N=16000, NB=400, TOL=1e-4 | ×0.54 | ×0.77 | ×1.01 | ×0.85 | 2.7e-4 / 2.7e-4 |
| N=16000, NB=800, TOL=1e-4 | ×0.53 | ×0.92 | ×2.18 | ×0.94 | 3.3e-4 / 3.2e-4 |

(ratios LDLᵀ/LU, < 1 = LDLᵀ meilleur)

## Limites connues, assumées

- **Gain mémoire ×0.5 seulement en régime peu compressible.** Quand les blocs
  hors-diagonale sont très compressés, la mémoire est dominée par les blocs
  **diagonaux denses**, que le LDLᵀ stocke aussi en entier (`CSYTRF_RK`
  travaille sur un tableau plein, pas en stockage compact) : à N=8000/NB=400,
  TOL=1e-4, les 20 blocs diagonaux pèsent ~25 Mo sur ~30. Le gain tend vers
  ×0.5 quand les rangs augmentent (TOL=0 : ×0.525).
- **Résolution non optimisée** : la remontée Lᵀ fait une réduction + une
  diffusion par bloc et garde une copie pleine largeur `XFULL`
  (MLOC × NRHS) des inconnues déjà calculées — plus lente que le LU à grand
  NB (×2.2 à NB=800), mémoire RHS ×PC pendant cette phase (bornée par
  `--nrhs-batch`). Temps absolus petits devant la factorisation ici.
- **Déséquilibre de charge sur grille carrée** (hérité de la distribution
  bloc-cyclique, aussi présent dans le LU) : les blocs diagonaux denses et
  leur `CSYTRF_RK` tombent tous sur les processus (r,r).
- **Matrices mal conditionnées** : comme pour le LU, l'erreur de compression
  (~TOL) est amplifiée par le conditionnement. Exemple : la variante « pivots
  mixtes » du test (`STRONG=1 DIAG=0.1`, `CGESV` dense seul ~6e-4) dépasse
  1e-3 à TOL=1e-4 alors que la factorisation est exacte au niveau de TOL —
  il faut resserrer TOL (1e-6 dans `run_all_tests.sh`).
- **Pas de variante OOC / arena / checkpoint** pour l'instant (une seule
  version, en mémoire, demandée pour ce premier jet).
