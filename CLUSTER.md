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

## 3. Scripts de balayage (`bench/`)

Trois scripts réutilisables, avec sortie CSV et timeout de sécurité par
run (`RUN_TIMEOUT`, défaut 300s) :

```bash
# Taille de bloc optimale (voir aussi la section 4 ci-dessous)
./bench/sweep_blocksize.sh [N] [NRHS] [TOL] [PR] [PC] [NB_LIST]
#   defauts : 10000 100 1e-4 2 2 "50 100 200 300 500 750 1000 1500 2000"

# Reaction au taux de compression (A et BLR+MPI uniquement -- seuls
# solveurs ou TOL a un sens)
./bench/sweep_tolerance.sh [N] [NB] [NRHS] [PR] [PC] [TOL_LIST]
#   defauts : 10000 500 100 2 2 "1e-1 1e-2 1e-3 1e-4 1e-6 1e-8 1e-10"

# Scalabilite forte : threads OpenMP (A/B/C) + grille MPI (MPI/BLR+MPI)
./bench/sweep_scaling.sh [N] [NRHS] [NB] [TOL]
#   defauts : 10000 100 500 1e-4
```

**Adaptation cluster** : surchargez `CORES_PER_NODE` si le script tourne
sur un nœud de connexion différent des nœuds de calcul, ou si un
scheduler restreint les cœurs réellement alloués :

```bash
CORES_PER_NODE=32 ./bench/sweep_blocksize.sh 20000 100 1e-4 4 4
```

Si le cluster utilise Slurm plutôt que `mpirun` directement, remplacer
`mpirun -np N` par `srun -n N` dans les trois scripts (un seul type
d'appel à changer, voir `bench/common.sh`).

## 4. Taille de bloc (`NB`) : ce qui la détermine

Pour `mpi_blr_matrix_mod.f90` (BLR+MPI), la taille de bloc optimale résulte
d'un compromis entre deux coûts qui évoluent en sens opposé avec `NB` :
- la **recompression** (dominante, 83-86% du temps de factorisation)
  décroît en `O(N³/NB²)` — favorise un `NB` grand ;
- le **LU du bloc diagonal** croît en `O(N·NB²)` — favorise un `NB` petit.

Ni l'un ni l'autre ne dépend du nombre de threads OpenMP (aucune boucle
`!$OMP` dans ce module — le calcul est distribué par MPI, pas threadé en
interne). Le point de croisement mesuré sur la machine de développement
était `NB≈500` à N=10000-20000 ; **relancez `sweep_blocksize.sh` sur le
cluster cible** plutôt que de réutiliser cette valeur telle quelle — le
point de croisement dépend du débit mono-cœur de la machine, qui diffère
d'un cluster à l'autre.
