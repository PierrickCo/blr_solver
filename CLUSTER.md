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
- la **recompression** (dominante, 83-86% du temps de factorisation,
  environ −25 % depuis la recompression incrémentale de 2026-11-03)
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

**Coût en temps de `OPENMP=0` : aucun.** Mesure isolée (N=40000, NB=500,
K=8, 3 paquets de 500 seconds membres, `CHECK=0 FAST=1` partout, 4
compilations **entrelacées** et répétées deux fois) :

| | LU : fact. / résolution | LU : RSS fin / pic | LDLᵀ : fact. / résolution | LDLᵀ : RSS fin / pic |
|---|---|---|---|---|
| OpenMP, MKL threadée | 2,7–3,1 / 22,7–25,4 s | 149 / 363 Mo | 1,95 / 6,0–6,1 s | 119 / 178 Mo |
| OpenMP, MKL séquentielle | 3,1–3,7 / 24,0–25,6 s | 147 / 362 Mo | 1,8–1,9 / 5,6–6,3 s | 117 / 176 Mo |
| `OPENMP=0` | 3,0–3,2 / 22,6–22,9 s | 133 / 316 Mo | 1,3–1,9 / 4,9–6,1 s | 106 / 162 Mo |
| `OPENMP=0`, `BLR_MALLOC_TRIM=0` | 2,9–3,1 / 22,7–23,4 s | 136 / 316 Mo | 1,3–1,8 / 5,1–6,2 s | 108 / 163 Mo |

Ni MKL séquentielle, ni l'allocateur glibc, ni `malloc_trim` ne coûtent de
temps au-delà du bruit (10 à 15 % entre deux runs identiques) ; la
mémoire baisse dès cette matrice compressible (−10 % en fin de
factorisation, −13 % au pic de résolution) et de 45 % sur la matrice peu
compressible. **Pour comparer des temps, entrelacer les variantes et
répéter** : sur la machine de test, une même compilation a varié de 30 %
entre le matin et le soir (LU, résolution : 19,5 puis 22,7 à 25,4 s), ce
qui suffit à fabriquer de faux écarts entre runs lancés à des heures
différentes. Testé et écarté : relever `MALLOC_MMAP_THRESHOLD_` et
`MALLOC_TRIM_THRESHOLD_` (64 Mo ou 1 Go) ralentit la résolution de 17 %
(LDLᵀ) à 30 % (LU) pour un gain mémoire négligeable.

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

## 6. Préparer un passage sur cluster sans cluster (trace et compteurs)

`BLR_TRACE_FACT=1` (variable d'environnement, tout pilote du LDLᵀ :
`bench_steps`, exemples, `bem_solve`) fait imprimer par la factorisation
LDLᵀ un résumé **indépendant des temps**, donc valable même avec plus de
processus que de cœurs :

- **Travail du Schur par processus**, estimé en opérations
  (`WORK_COUNT`, `SCHUR_WORK`) : moyenne, maximum, et somme des maximums
  par étape.
- **Trois indicateurs** :
  - le déséquilibre global, maximum / moyenne ;
  - la perte de synchronisation par étape, somme des maximums par étape /
    maximum des totaux ;
  - une borne de l'accélération, P × moyenne / somme des maximums.
- **Diffusions de blocs par processus** : nombre d'appels et éléments
  transportés (`BCAST_NCALL`, `BCAST_ELEMS`).

Jusqu'à 8 processus s'ajoute le détail des temps : communications
ventilées par processus, puis Schur et communications à chaque étape.

Émulation sur un portable à 4 cœurs (matrice de test, N = 12 000, NB = 150,
80 blocs ; `I_MPI_WAIT_MODE=1 I_MPI_PIN=off`, les temps n'ont alors pas de
sens) :

| Grille | Déséquilibre global | Synchronisation par étape | Borne d'accélération | Diffusions par processus | Travail par diffusion |
|---|---|---|---|---|---|
| 2 × 2 | 2,3 % | 2,4 % | 3,8 / 4 | 6 320 | 5,2e6 |
| 4 × 4 | 7,7 % | 7,4 % | 13,8 / 16 | 3 160 | 2,6e6 |
| 6 × 6 | 14,9 % | 10,0 % | 28,5 / 36 | 2 107 | 1,7e6 |
| 8 × 8 | 25,9 % | 12,3 % | 45,3 / 64 | 1 580 | 1,3e6 |

**Lecture :**
- **Déséquilibre.** À 64 processus sur 80 blocs, le déséquilibre global
  est la première perte (26 %), devant la synchronisation (12 %).
- **Messages.** Le nombre de diffusions par processus ne décroît que comme
  1/√P : le travail utile par message est divisé par 4 de 2×2 à 8×8.
- **Limite.** Ces chiffres dépendent du nombre de blocs par ligne de
  grille (10 ici en 8×8). À refaire à la taille visée.

**Entiers 32 bits.** Les résolutions LU et LDLᵀ s'arrêtent avec un message
(`BLR_CHECK_COUNT32`) si un compteur MPI (lignes locales × seconds membres,
taille de bloc × seconds membres, colonne de L diffusée) dépasse 2³¹−1 :
réduire alors le lot de seconds membres. La factorisation elle-même n'a pas
de compteur à risque (blocs au plus NB × NB).

## 7. Simulateur de répartition (`tools/blr_sim/simulate.py`)

`BLR_DUMP_WORK=1` fait écrire par chaque processus un relevé
`blr_work_rank<id>.bin` : étape K, bloc (I, J) et travail estimé de chaque
mise à jour de Schur. Ce relevé ne dépend pas de la grille du calcul. Le
simulateur rejoue ensuite la charge de chaque processus pour n'importe
quelle grille et n'importe quelle règle de répartition, en quelques
secondes et sans MPI :

```bash
BLR_DUMP_WORK=1 mpirun -np 4 ldlt/bin_noomp_nocheck_fast/bench_steps LDLT 2 2 28600 100 4 1 1e-4 5000 4 2
python3 tools/blr_sim/simulate.py . --grilles 8x8,16x16,32x32,44x45 --repartitions actuel,cyclique,pondere
```

**Validation.** Avec un relevé pris en 2×2 (N = 12 000, 80 blocs), le
simulateur retrouve exactement les indicateurs mesurés directement :
- en 4×4 : 7,3 % de déséquilibre, 7,1 % de synchronisation, borne 13,9 ;
- en 8×8 : 21,1 %, 13,0 %, borne 46,8.

**Étude à 286 blocs.** C'est le nombre de blocs attendu à N = 1 000 000
avec NB = 3 500. La matrice de test utilisée (N = 28 600, NB = 100) a la
bonne structure, mais pas les rangs d'un vrai problème à un million.

| Grille | Répartition | Déséquilibre | Synchronisation | Borne d'accélération | Diffusions par processus | Si regroupées |
|---|---|---|---|---|---|---|
| 8 × 8 | actuelle | 5,2 % | 6,0 % | 57 / 64 | 20 378 | 1 124 |
| 10 × 10 | actuelle | 7,8 % | 7,3 % | 87 / 100 | 16 302 | 1 121 |
| 16 × 16 | actuelle | 17,4 % | 8,7 % | 201 / 256 | 10 189 | 1 108 |
| 16 × 16 | pondérée | 8,6 % | 46,7 % | 161 / 256 | | |
| 32 × 32 | actuelle | 46,7 % | 7,5 % | 649 / 1 024 | 5 094 | 1 076 |
| 32 × 32 | pondérée | 20,2 % | 78,1 % | 478 / 1 024 | | |
| 44 × 45 | actuelle | 69,9 % | 6,3 % | 1 096 / 1 980 | 3 664 | 1 053 |
| 44 × 45 | pondérée | 46,9 % | 120 % | 612 / 1 980 | | |

**Lecture :**
- **La répartition actuelle** (cyclique décalée) reste la meilleure. La
  cyclique simple donne la même borne à 1 % près.
- **Une répartition pondérée par le coût total fait pire.** Elle réduit le
  déséquilibre global mais regroupe le travail dans les mêmes étapes, et la
  perte de synchronisation explose. Il faut équilibrer étape par étape, ce
  que le cyclique fait naturellement.
- **L'efficacité de la charge baisse avec la grille** : 90 % à 64
  processus, 78 % à 256, 63 % à 1 024, 55 % à 1 980. La cause est la
  granularité : 6 à 9 blocs par ligne de grille. Les remèdes sont plus de
  blocs (NB plus petit) ou moins de processus avec des threads.
- **Regrouper les diffusions par panneau** divise le nombre d'appels par 18
  en 8×8 et par 3,5 en 44×45.

## 8. Diffusion groupée des panneaux (`BCAST_BLR_PANEL`)

**Avant.** Chaque bloc d'un panneau partait en deux appels collectifs :
un en-tête, puis ses données. Le nombre d'appels par étape croissait avec
le nombre de blocs du panneau.

**Maintenant.** Tous les blocs d'un panneau, vus d'un communicateur de
ligne ou de colonne, partent en deux appels : les en-têtes, puis les
données empilées.
- Panneau à une seule racine (W le long des lignes en LDLᵀ, L et U en
  LU) : deux `MPI_BCAST`.
- Panneau à plusieurs racines (Lᵀ le long des colonnes en LDLᵀ) : un
  `MPI_ALLREDUCE` pour les en-têtes, un `MPI_ALLGATHERV` pour les données.
- Les envois restent bloquants. Les octets transportés sont les mêmes :
  le résultat est identique au bit près.

**Routines concernées :** `MPI_BLR_SYM_LDLT_FACTOR` (classique),
`MPI_BLR_MATRIX_LU_FACTOR_PIVOT` et `MPI_BLR_MATRIX_LU_FACTOR_NOPIV`. Les
variantes arène, hors mémoire et capacité fixe diffusent encore bloc par
bloc.

**Garde-fous.**
- Si un panneau empilé dépasse 2³¹−1 éléments, la routine revient à la
  diffusion bloc par bloc.
- `BLR_PANEL_PACK=0` force la diffusion bloc par bloc, pour comparer.

**Mesures** (`bem_solver`, LDLᵀ, grille 2×2, appels de diffusion par
processus, `BLR_TRACE_FACT=1`) :

| Cas | Bloc par bloc | Groupée | Solution |
|---|---|---|---|
| 6 816 inconnues, 2×2 | 552 | 90 | identique au bit près |
| 6 816 inconnues, 4×1 et 1×4 | 690 | 89 | identique au bit près |
| 40 230 inconnues, 2×2 | 3 306 | 226 | identique au bit près |

- **Éléments transportés inchangés** (9,345e7 par processus à 40k).
- **Temps de factorisation à 40k inchangé sur le portable** : 28,1 / 28,7 /
  31,8 s bloc par bloc, 28,4 / 28,7 / 33,9 s groupée, passes alternées.
  Avec 4 processus en mémoire partagée, un appel collectif ne coûte presque
  rien ; le gain est attendu sur un réseau, quand la grille grandit.
- **Pic mémoire à 40k inchangé** : 366-368 Mo par processus groupée,
  369-373 bloc par bloc. Le tampon empilé d'un panneau est transitoire.
- **LU avec pivot** (3 cas de 6 816 à 12 291 inconnues) : solutions
  identiques au bit près.
- `run_all_tests.sh intel-noomp` : tous les tests passent, bornes
  vérifiées.

**Suite possible :** rendre ces deux appels non bloquants, pour que le
panneau de l'étape suivante parte pendant les mises à jour de Schur.
