#!/bin/bash
#===============================================================================
# Balayage de la tolerance de compression TOL, a N/NB/NRHS/grille fixes, sur
# les 2 solveurs pour lesquels TOL a un sens (A -- OpenMP compresse -- et
# BLR+MPI sans pivot) : montre comment temps, memoire, taux de compression
# et rang moyen reagissent quand on exige plus ou moins de precision par
# bloc. B/MPI-dense/C sont hors-sujet ici (TOL fixe ou inexistante).
#
# Usage : ./bench/sweep_tolerance.sh [N] [NB] [NRHS] [PR] [PC] [TOL_LIST]
#   defauts : N=10000 NB=500 NRHS=100 PR=2 PC=2
#   TOL_LIST (entre guillemets) : defaut "1e-1 1e-2 1e-3 1e-4 1e-6 1e-8 1e-10"
#     (TOL=0 volontairement exclue : degenere en "B", deja couvert par
#     sweep_blocksize.sh)
#
# Sortie : bench/results_tolerance.csv (solver,N,NB,NRHS,TOL,compression_pct,
#   avg_rank,t_assemble,t_factor,t_solve,mem_mb,error)
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-10000}"
NB="${2:-500}"
NRHS="${3:-100}"
PR="${4:-2}"
PC="${5:-2}"
TOL_LIST_STR="${6:-1e-1 1e-2 1e-3 1e-4 1e-6 1e-8 1e-10}"
read -r -a TOL_LIST <<< "$TOL_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

CSV="bench/results_tolerance.csv"
echo "solver,N,NB,NRHS,TOL,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_mb,error" > "$CSV"
echo "Balayage TOL : N=$N NB=$NB NRHS=$NRHS grille=${PR}x${PC} OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "TOL testees : ${TOL_LIST[*]}"

for TOL in "${TOL_LIST[@]}"; do
  echo "=== TOL=$TOL ==="

  # --- A : OpenMP compresse ---
  log=$(run_and_log "tolerance_A_TOL${TOL}" ./bench_blr_pipeline "$N" "$NRHS" "$NB" "$TOL")
  comp_pct=$(floats_on_line "Blocs compresses" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  t_asm=$(grep_num "Temps assemblage+compression" "$log")
  t_fac=$(grep_num "Temps factorisation =" "$log")
  t_sol=$(grep_num "Temps resolution =" "$log")
  mem=$(grep_num "Memoire des facteurs L/U" "$log")
  err=$(grep_num "Erreur relative ||Xcalc-Xexact" "$log")
  echo "A,$N,$NB,$NRHS,$TOL,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem,$err" >> "$CSV"

  # --- BLR+MPI sans pivot ---
  log=$(run_and_log "tolerance_BLRMPI_TOL${TOL}" mpirun -np "$NPROCS" ./test_mpi_blr_matrix "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL")
  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  err=$(grep_num "Erreur relative ||X_distribue" "$log")
  echo "BLR+MPI,$N,$NB,$NRHS,$TOL,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
