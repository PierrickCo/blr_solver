#!/bin/bash
#===============================================================================
# BENCHMARK ARENA 4/5 -- balayage de la TAILLE DU SYSTEME (N), a NB/EPS/
# RATIO/grille FIXES -- montre comment le temps par etape, la memoire
# (theorique et RSS) et leur ratio evoluent avec N pour le solveur BLR+MPI
# arena. EPS=0 par defaut (matrice pleinement compressible) pour isoler
# l'effet de N seul, sans confondre avec l'effet de compressibilite deja
# caracterise par sweep_arena_compressibility.sh.
#
# Usage : ./bench/sweep_arena_scaling_n.sh [NB] [NRHS] [TOL] [EPS] [RATIO] [PR] [PC] [N_LIST]
#   defauts : NB=200 NRHS=8 TOL=1e-4 EPS=0 RATIO=0.333333 PR=2 PC=2
#   N_LIST (entre guillemets) : defaut "2000 5000 10000 20000 40000"
#
# Sortie : bench/results_arena_scaling_n.csv (N,NB,NRHS,TOL,EPS,RATIO,
#   compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,
#   rss_mb,ratio_rss_theorique,erreur)
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

NB="${1:-200}"
NRHS="${2:-8}"
TOL="${3:-1e-4}"
EPS="${4:-0}"
RATIO="${5:-0.333333}"
PR="${6:-2}"
PC="${7:-2}"
N_LIST_STR="${8:-2000 5000 10000 20000 40000}"
read -r -a N_LIST <<< "$N_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

BIN="./bin/bench_blr_mpi_pivot_arena"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_arena)..."
  make -f Makefile bin/bench_blr_mpi_pivot_arena
fi

CSV="bench/results_arena_scaling_n.csv"
echo "N,NB,NRHS,TOL,EPS,RATIO,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "[ARENA 4/5] Balayage taille du systeme : NB=$NB NRHS=$NRHS TOL=$TOL EPS=$EPS RATIO=$RATIO grille=${PR}x${PC}"
echo "N testees : ${N_LIST[*]}"

for N in "${N_LIST[@]}"; do
  echo "=== N=$N ==="
  log=$(run_and_log "arena_scaling_n_N${N}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS" "$RATIO")

  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem_theo=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  rss=$(floats_on_line "Memoire ANALYTIQUE" "$log" | tail -1)
  ratio_rss=$(grep_num "Ratio RSS" "$log")
  err=$(grep_num "Erreur relative" "$log")

  echo "$N,$NB,$NRHS,$TOL,$EPS,$RATIO,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem_theo,$rss,$ratio_rss,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
