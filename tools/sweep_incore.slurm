#!/bin/bash
#===============================================================================
# Exemple de soumission Slurm de tools/sweep_incore.sh (N ~ 100000, ~300
# coeurs). A ADAPTER : partition, nombre de noeuds / coeurs par noeud,
# modules, duree. Lancer depuis blr_solver/ :  sbatch tools/sweep_incore.slurm
#
# L'allocation reserve CORES coeurs ; chaque run du balayage utilise
#   srun -n <PR*PC> --cpus-per-task=<OMP_NUM_THREADS>
# a l'interieur de cette allocation (PR*PC*threads = CORES).
#===============================================================================
#SBATCH --job-name=blr_sweep
#SBATCH --nodes=5                  # ex. 5 noeuds x 60 coeurs = 300 coeurs
#SBATCH --ntasks=300
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=blr_sweep_%j.out

set +u
source /opt/intel/oneapi/setvars.sh      # ou : module load intel-oneapi ...
set -u

export N=100000
export CORES=${SLURM_NTASKS:-300}
export LAUNCH="srun -n"
export LAUNCH_THREADS_OPT="--cpus-per-task="
export THREADS_LIST="1 2 4"
export REPEATS=2                          # 3 si le budget temps le permet
export TIMEOUT=7200
export OUTDIR=sweep_results_N${N}_${SLURM_JOB_ID:-local}

./tools/sweep_incore.sh
