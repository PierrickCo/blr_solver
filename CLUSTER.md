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
