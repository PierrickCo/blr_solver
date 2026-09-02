#!/bin/bash
#===============================================================================
# BENCHMARK ARENA 2/5 -- balayage de la TAILLE DE BLOC (NB), a N/EPS/RATIO/
# grille FIXES -- montre comment NB affecte le temps par etape, la memoire
# (theorique et RSS) et leur ratio pour le solveur BLR+MPI arena.
#
# Usage : ./bench/sweep_arena_blocksize.sh [N] [NRHS] [TOL] [EPS] [RATIO] [PR] [PC] [NB_LIST]
#   defauts : N=20000 NRHS=8 TOL=1e-4 EPS=0 RATIO=0.333333 PR=2 PC=2
#   NB_LIST (entre guillemets) : defaut "100 200 400 800 2000"
#
# Sortie : bench/results_arena_blocksize.csv (NB,N,NRHS,TOL,EPS,RATIO,
#   compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,
#   rss_mb,ratio_rss_theorique,erreur)
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-20000}"
NRHS="${2:-8}"
TOL="${3:-1e-4}"
EPS="${4:-0}"
RATIO="${5:-0.333333}"
PR="${6:-2}"
PC="${7:-2}"
NB_LIST_STR="${8:-100 200 400 800 2000}"
read -r -a NB_LIST <<< "$NB_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

BIN="./bin/bench_blr_mpi_pivot_arena"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_arena)..."
  make -f Makefile bin/bench_blr_mpi_pivot_arena
fi

CSV="bench/results_arena_blocksize.csv"
echo "NB,N,NRHS,TOL,EPS,RATIO,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "[ARENA 2/5] Balayage taille de bloc : N=$N NRHS=$NRHS TOL=$TOL EPS=$EPS RATIO=$RATIO grille=${PR}x${PC}"
echo "NB testees : ${NB_LIST[*]}"

for NB in "${NB_LIST[@]}"; do
  echo "=== NB=$NB ==="
  log=$(run_and_log "arena_blocksize_NB${NB}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS" "$RATIO")

  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem_theo=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  rss=$(floats_on_line "Memoire ANALYTIQUE" "$log" | tail -1)
  ratio_rss=$(grep_num "Ratio RSS" "$log")
  err=$(grep_num "Erreur relative" "$log")

  echo "$NB,$N,$NRHS,$TOL,$EPS,$RATIO,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem_theo,$rss,$ratio_rss,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
