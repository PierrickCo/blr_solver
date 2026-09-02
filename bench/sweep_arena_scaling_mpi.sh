#!/bin/bash
#===============================================================================
# BENCHMARK ARENA 5/5 -- balayage du NOMBRE DE TACHES MPI (grille PRxPC), a
# N/NB/EPS/RATIO FIXES -- montre comment le temps par etape, la memoire
# (theorique et RSS, max/processus) et leur ratio evoluent avec le nombre de
# processus pour le solveur BLR+MPI arena. OMP_NUM_THREADS est ajuste a
# chaque grille pour eviter la sur-souscription (voir bench/common.sh).
#
# Usage : ./bench/sweep_arena_scaling_mpi.sh [N] [NB] [NRHS] [TOL] [EPS] [RATIO] [GRID_LIST]
#   defauts : N=20000 NB=200 NRHS=8 TOL=1e-4 EPS=0 RATIO=0.333333
#   GRID_LIST (entre guillemets, format "PRxPC") : defaut "1x1 1x2 2x1 2x2"
#
# Sortie : bench/results_arena_scaling_mpi.csv (NPROCS,PR,PC,N,NB,NRHS,TOL,
#   EPS,RATIO,compression_pct,avg_rank,t_assemble,t_factor,t_solve,
#   mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur)
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-20000}"
NB="${2:-200}"
NRHS="${3:-8}"
TOL="${4:-1e-4}"
EPS="${5:-0}"
RATIO="${6:-0.333333}"
GRID_LIST_STR="${7:-1x1 1x2 2x1 2x2}"
read -r -a GRID_LIST <<< "$GRID_LIST_STR"

BIN="./bin/bench_blr_mpi_pivot_arena"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_arena)..."
  make -f Makefile bin/bench_blr_mpi_pivot_arena
fi

CSV="bench/results_arena_scaling_mpi.csv"
echo "NPROCS,PR,PC,N,NB,NRHS,TOL,EPS,RATIO,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "[ARENA 5/5] Balayage nombre de taches MPI : N=$N NB=$NB NRHS=$NRHS TOL=$TOL EPS=$EPS RATIO=$RATIO"
echo "Grilles testees : ${GRID_LIST[*]}"

for GRID in "${GRID_LIST[@]}"; do
  PR="${GRID%x*}"
  PC="${GRID#*x}"
  NPROCS=$((PR*PC))
  OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")
  export OMP_NUM_THREADS
  echo "=== Grille ${PR}x${PC} (NPROCS=$NPROCS, OMP_NUM_THREADS=$OMP_NUM_THREADS) ==="
  log=$(run_and_log "arena_scaling_mpi_${PR}x${PC}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS" "$RATIO")

  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem_theo=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  rss=$(floats_on_line "Memoire ANALYTIQUE" "$log" | tail -1)
  ratio_rss=$(grep_num "Ratio RSS" "$log")
  err=$(grep_num "Erreur relative" "$log")

  echo "$NPROCS,$PR,$PC,$N,$NB,$NRHS,$TOL,$EPS,$RATIO,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem_theo,$rss,$ratio_rss,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
