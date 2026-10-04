# Déploiement sur un vrai cluster

Ce document règle une question précise, testée et validée durant le
développement : **comment éviter la sur-souscription de threads sur un
cluster réel**, et comment retrouver les réglages optimaux (taille de bloc,
nombre de processus) sans repartir de zéro.

## 1. Affinité CPU / `OMP_NUM_THREADS`

Chaque binaire de ce projet destiné à un usage réel (`bench_mpi`,
`test_mpi_blr_matrix`, `test_mpi_blr_matrix_pivot`) appelle
`REPORT_THREAD_AFFINITY` juste après l'initialisation MPI. Il affiche, une
fois par nœud physique :

```
--------------------------------------------------------------
 [AFFINITE] Noeud                           = <hostname>
 [AFFINITE] Processus MPI sur ce noeud      = <N>
 [AFFINITE] OMP_NUM_THREADS (ce processus)  = <T>
 [AFFINITE] Coeurs vus (ce processus)       = <P>
 [AFFINITE] Coeurs distincts occupes/noeud  = <union>
 [AFFINITE] Total threads demandes / noeud  = <somme>
 [AFFINITE] OK : pas de sur-souscription (masques d'affinite pris en compte).
--------------------------------------------------------------
```

**Toujours vérifier cette sortie au premier lancement sur un nouveau
cluster.** Le diagnostic lit le masque d'affinité CPU réel de chaque
processus (`/proc/self/status`) et compare le nombre de **cœurs
distincts réellement occupés sur le nœud** (union des masques) au
**nombre total de threads demandés** — ce test reste correct que le
lanceur MPI applique ou non un pinning par défaut (Intel MPI pine par
défaut, `I_MPI_PIN=1` ; d'autres lanceurs peuvent ne rien faire).

**Formule de réglage manuel**, si vous lancez sans `mpirun`/`srun` géré
automatiquement par un script :

```
OMP_NUM_THREADS = coeurs_par_noeud / processus_MPI_par_noeud
```

Fixez-la explicitement (`export OMP_NUM_THREADS=...`) plutôt que de
compter sur la détection automatique de MKL, qui ne connaît pas le nombre
de processus MPI colocalisés sur le nœud.

## 2. Exemple de script Slurm

```bash
#!/bin/bash
#SBATCH --job-name=blr_mpi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8        # ajuster : coeurs_du_noeud / ntasks-per-node
#SBATCH --time=00:30:00

source /opt/intel/oneapi/setvars.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --cpu-bind=cores ./bin/test_mpi_blr_matrix 2 2 20000 500 100 1e-4
```

`--cpu-bind=cores` (ou l'option équivalente du lanceur utilisé) est un
complément plus robuste que le seul réglage de `OMP_NUM_THREADS`,
particulièrement sur un nœud partagé entre plusieurs jobs.

## 3. Recherche des meilleurs paramètres (`tools/sweep_incore.sh`)

Les anciens scripts `bench/` ont été retirés (commit `29ee771`). Le script
`tools/sweep_incore.sh` les remplace : il cherche automatiquement les
meilleurs réglages des solveurs **en mémoire** (LU et LDLᵀ, variantes
classique et arena, pas d'out-of-core) sur la matrice de démo complexe
symétrique, via le pilote `ldlt/bin/sweep_driver`, en phases successives :

1. **TOL** : le plus grand TOL dont l'erreur reste ≤ `ERR_MAX` ;
2. **NB** : la taille de bloc minimisant l'objectif ;
3. **grille / threads** : PR × PC × `OMP_NUM_THREADS` = `CORES`, plusieurs formes ;
4. **RATIO** (`COMPRESS_WORTHWHILE_RATIO`) ;
5. **arena** au meilleur réglage (conseillée seulement si la RSS baisse).

Chaque configuration est lancée `REPEATS` fois (médiane) ; en phase 2, parmi
les NB à moins de `NEAR` (5 %) du meilleur temps on retient le moins gourmand
en mémoire ; en phases 3-4 on ne quitte le réglage par défaut que pour un gain
> `MIN_GAIN` (3 %) — garde-fous contre le bruit de mesure. Sorties dans
`OUTDIR` : `runs.csv`, `summary.csv`, `decisions.txt`, `best_params.txt`
(meilleurs paramètres, commandes d'exemple conseillées, et meilleures
configurations mesurées par critère : factorisation, total, mémoire), `logs/`.
`./tools/sweep_incore.sh --report OUTDIR` régénère le rapport sans relancer.

```bash
make && make -C ldlt
./tools/sweep_incore.sh                                # poste local : N=20000, CORES=$(nproc)
N=100000 CORES=300 LAUNCH="srun -n" LAUNCH_THREADS_OPT="--cpus-per-task=" ./tools/sweep_incore.sh
```

Toutes les options sont des variables d'environnement documentées en tête du
script (`NB_LIST`, `TOL_LIST`, `THREADS_LIST`, `OBJECTIVE`, `ERR_MAX`,
`TIMEOUT`...). Exemple de soumission Slurm : `tools/sweep_incore.slurm`.

**Coût** : ~20 configurations × `REPEATS` × solveurs (≈120 runs par défaut).
Chaque run calcule aussi B = A·Xexact par produit brut (non chronométré,
O(N²/PR) évaluations du noyau par processus). À N=100000, prévoir plusieurs
heures ; réduire `REPEATS`, `NB_LIST` ou `SOLVERS` si besoin.

## 4. Taille de bloc (`NB`) : ce qui la détermine

Pour `mpi_blr_matrix_mod.f90` (BLR+MPI), la taille de bloc optimale résulte
d'un compromis entre deux coûts qui évoluent en sens opposé avec `NB` :
- la **recompression** (dominante, 83-86% du temps de factorisation)
  décroît en `O(N³/NB²)` — favorise un `NB` grand ;
- le **LU du bloc diagonal** croît en `O(N·NB²)` — favorise un `NB` petit.

Ni l'un ni l'autre ne dépend du nombre de threads OpenMP (aucune boucle
`!$OMP` dans ce module — le calcul est distribué par MPI, pas threadé en
interne). Le point de croisement mesuré sur la machine de développement
était `NB≈500` à N=10000-20000 ; **relancez `tools/sweep_incore.sh` sur le
cluster cible** plutôt que de réutiliser cette valeur telle quelle — le
point de croisement dépend du débit mono-cœur de la machine, qui diffère
d'un cluster à l'autre.

## 5. Compilation pour les vrais calculs : OpenMP, mémoire et vérifications

**Runs MPI purs (1 thread par processus), calculs de production :
compiler avec `OPENMP=0 CHECK=0`, plus `FAST=1` (sur un nœud de calcul)
quand l'assemblage pèse lourd dans le temps total — voir « Optimisation »
ci-dessous.**

```bash
O="OPENMP=0 CHECK=0 FAST=1"
make $O && make -C ldlt $O              # -> bin_noomp_nocheck_fast/, ldlt/bin_noomp_nocheck_fast/
make -C user_case_dev $O                # -> example_blr_mpi_pivot_user_noomp_nocheck_fast
CHECK=0 FAST=1 ./run_all_tests.sh intel-noomp   # meme suite de tests, meme compilation
```

Les trois options sont indépendantes et combinables ; chaque combinaison a
ses propres répertoires de sortie (suffixes `_noomp`, `_nocheck`, `_fast`,
dans cet ordre), et les Makefiles `ldlt/` et `user_case_dev/` prennent les
mêmes options que le Makefile racine.

Pourquoi : avec ifx, dès que `libiomp5` (le runtime OpenMP d'Intel, amené
par MKL threadée) est lié, **tous les `ALLOCATE` passent par l'allocateur
interne de ce runtime, qui ne rend jamais la mémoire libérée au système**.
Or la factorisation réalloue sans cesse les facteurs U/V des blocs
bas-rang à des rangs croissants. `OPENMP=0` compile sans OpenMP et lie MKL
séquentielle : les `ALLOCATE` passent par glibc, et chaque étape de
factorisation rend la mémoire libérée par `malloc_trim(0)`
(`BLR_TRIM_HEAP`, désactivable par `BLR_MALLOC_TRIM=0`).

Mesuré (N=40000, NB=500, grille 2×2, matrice peu compressible K=5000 dont
le rang moyen passe de 8 à 53 pendant la factorisation ; RSS max par
processus en fin de factorisation) :

| | défaut (OpenMP, MKL threadée) | `OPENMP=0` | `OPENMP=0`, `BLR_MALLOC_TRIM=0` |
|---|---|---|---|
| LDLᵀ | 816 Mo | **449 Mo** (−45 %) | 516 Mo |
| LU | 1332 Mo | **721 Mo** (−46 %) | 832 Mo |
| facteurs seuls (théorique) + RSS de départ | LDLᵀ 344 Mo, LU 618 Mo | | |

Temps de factorisation LU : 178 / 177 / 177 s (identiques). LDLᵀ : 67 à
91 s selon les runs dans les trois configurations (variabilité de la
machine de test, pas d'écart systématique mesuré).

Non concernés : la chaîne gcc (`Makefile.gcc`, libgomp n'intercepte pas
les allocations ; `BLR_TRIM_HEAP` y agit directement : 460 → 300 Mo en fin
de factorisation sur le même cas) et `Makefile.ifort2019` (non modifié).
**Mode hybride MPI + threads** : `OPENMP=0` ne s'applique pas tel quel ;
garder la compilation par défaut (les couches de threads GNU ou TBB de MKL,
qui éviteraient `libiomp5`, n'ont pas été testées).

**Vérifications de bornes (`CHECK=0`).** Par défaut, tous les Makefiles
compilent avec `-check bounds` (`-fcheck=bounds` en gcc, option aussi
disponible dans `Makefile.gcc` et `ldlt/Makefile.gcc`), utile au
développement mais coûteux : mesuré sur le même cas (K=8), assemblage ×2,4
(LDLᵀ) à ×3,5 (LU), factorisation LU ×1,7, résolution ×1,1 à ×1,4.
`CHECK=0` les retire. Garder la compilation par défaut pour développer et
déboguer : un dépassement de tableau y est signalé au lieu de corrompre
silencieusement la mémoire.

**Optimisation (`FAST=1`).** Par défaut `-O2` générique ; `FAST=1` passe en
`-O3` avec le jeu d'instructions de la machine qui compile (`FAST_ARCH`,
défaut `-xHost`, `-march=native` en gcc). Mesuré (N=40000, NB=500, sans
OpenMP ni vérification des bornes, une mesure par cas) :

| `-O2` → `FAST=1` | LU K=8 | LDLᵀ K=8 | LU K=5000 | LDLᵀ K=5000 |
|---|---|---|---|---|
| assemblage | 20,8 → 11,7 s | 15,0 → 5,9 s | 32,6 → 14,0 s | 15,2 → 8,5 s |
| factorisation | 2,1 → 2,9 s | 1,6 → 1,7 s | 172,7 → 159,5 s | 80,9 → 97,8 s |
| résolution (500 seconds membres) | 6,3 → 7,8 s | 1,9 → 2,0 s | 14,4 → 15,6 s | 6,8 → 8,6 s |

L'assemblage gagne nettement (l'évaluation du noyau de la matrice est
vectorisée) ; factorisation et résolution varient de −8 % à +30 % selon le
cas (sur cette machine, des runs identiques de la factorisation LDLᵀ K=5000
ont varié de 65 à 91 s : les écarts de 10 à 20 % ne sont pas
significatifs). Le gain net dépend donc du poids de l'assemblage : à
mesurer sur sa propre matrice (`bench_steps`, en compilant les deux
variantes). **`-xHost` vise le processeur de
la machine de compilation** : compiler sur un nœud de calcul, ou imposer
l'architecture des nœuds, par exemple `make $O FAST_ARCH=-xCORE-AVX512`
(un binaire compilé pour une architecture plus récente que celle des nœuds
s'arrête au démarrage).

**Mesurer soi-même** : `ldlt/bin*/bench_steps` (temps, RSS de fin et pic
par étape, taux de compression, décomposition de la RSS) :

```bash
mpirun -np 4 ldlt/bin_noomp_nocheck_fast/bench_steps LDLT 2 2 40000 500 500 3 1e-4 5000 4
```

