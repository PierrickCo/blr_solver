#!/bin/bash
#===============================================================================
# Script de lancement pour examples/example_blr_mpi_pivot.f90 -- le solveur
# le plus abouti du projet : BLR+MPI avec pivot partiel confine au bloc
# diagonal (voir README.md, section "Ce qui existe"). Compile si besoin,
# fixe OMP_NUM_THREADS selon la convention validee du projet
# (voir CLUSTER.md), puis lance mpirun.
#
# Usage : ./examples/run_example_blr_mpi_pivot.sh [PR] [PC] [N] [NB] [NRHS] [TOL]
#   (a lancer depuis blr_solver/, ou directement -- le script se replace
#   tout seul dans blr_solver/)
#   defauts : PR=2 PC=2 N=4000 NB=200 NRHS=4 TOL=1e-4
#
# Variables d'environnement optionnelles :
#   TOOLCHAIN=intel|gcc   chaine de compilation (defaut intel, voir INSTALL.md)
#   CORES_PER_NODE=N      coeurs REELLEMENT alloues a ce job (defaut: nproc) --
#                         a surcharger sous un scheduler (Slurm --cpus-per-task
#                         etc.), ne jamais se fier a nproc seul sur un noeud
#                         partage (voir bench/common.sh, meme convention).
#   MPIRUN="srun -n"      remplacer le lanceur si mpirun n'est pas celui du
#                         cluster (defaut: "mpirun -np").
#===============================================================================
set -u
BLR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BLR_ROOT"

PR="${1:-2}"
PC="${2:-2}"
N="${3:-4000}"
NB="${4:-200}"
NRHS="${5:-4}"
TOL="${6:-1e-4}"

TOOLCHAIN="${TOOLCHAIN:-intel}"
if [ "$TOOLCHAIN" = "gcc" ]; then
  MAKEFILE="Makefile.gcc"
  BINARY="build-gcc/example_blr_mpi_pivot"
else
  MAKEFILE="Makefile"
  BINARY="example_blr_mpi_pivot"
fi

# Source l'environnement Intel oneAPI si besoin (sans effet si deja fait,
# ou si TOOLCHAIN=gcc -- idempotent, meme garde `set +u`/`set -u` que
# bench/common.sh : setvars.sh reference des variables non definies).
if [ "$TOOLCHAIN" = "intel" ] && [ -f /opt/intel/oneapi/setvars.sh ]; then
  set +u
  source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1
  set -u
fi

if [ ! -x "$BINARY" ]; then
  echo "=== $BINARY absent, compilation (chaine: $TOOLCHAIN) ==="
  make -f "$MAKEFILE" examples || { echo "ECHEC DE COMPILATION"; exit 1; }
fi

NPROCS=$((PR * PC))
CORES_PER_NODE="${CORES_PER_NODE:-$(nproc)}"
export OMP_NUM_THREADS=$((CORES_PER_NODE / NPROCS))
[ "$OMP_NUM_THREADS" -lt 1 ] && export OMP_NUM_THREADS=1

MPIRUN="${MPIRUN:-mpirun -np}"

echo "Grille ${PR}x${PC} (${NPROCS} processus), N=$N NB=$NB NRHS=$NRHS TOL=$TOL, OMP_NUM_THREADS=$OMP_NUM_THREADS"
$MPIRUN "$NPROCS" "./$BINARY" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL"
