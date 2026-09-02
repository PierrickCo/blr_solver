#!/bin/bash
#===============================================================================
# BENCHMARK ARENA 3/5 -- balayage de COMPRESS_WORTHWHILE_RATIO, a N/NB/EPS/
# grille FIXES, sur une matrice MOYENNEMENT compressible (ni 100% compresse,
# ni totalement dense) -- EPS est un PARAMETRE explicite (defaut=1e-6, situe
# dans la zone de transition mesuree par sweep_arena_compressibility.sh sous
# RATIO=1/3, ~13% de compression -- voir memoire du projet), a ajuster selon
# ce que l'utilisateur considere "moyennement compressible" pour son cas.
#
# Usage : ./bench/sweep_arena_ratio.sh [N] [NB] [NRHS] [TOL] [EPS] [PR] [PC] [RATIO_LIST]
#   defauts : N=8000 NB=200 NRHS=8 TOL=1e-4 EPS=1e-6 PR=2 PC=2
#   RATIO_LIST (entre guillemets) : defaut "0.1 0.2 0.3 0.333333 0.4 0.5 0.6 0.7 0.8 0.9"
#
# Sortie : bench/results_arena_ratio.csv (RATIO,RMAX,EPS,N,NB,NRHS,TOL,
#   compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,
#   rss_mb,ratio_rss_theorique,erreur)
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-8000}"
NB="${2:-200}"
NRHS="${3:-8}"
TOL="${4:-1e-4}"
EPS="${5:-1e-6}"
PR="${6:-2}"
PC="${7:-2}"
RATIO_LIST_STR="${8:-0.1 0.2 0.3 0.333333 0.4 0.5 0.6 0.7 0.8 0.9}"
read -r -a RATIO_LIST <<< "$RATIO_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

BIN="./bin/bench_blr_mpi_pivot_arena"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_arena)..."
  make -f Makefile bin/bench_blr_mpi_pivot_arena
fi

CSV="bench/results_arena_ratio.csv"
echo "RATIO,RMAX,EPS,N,NB,NRHS,TOL,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "[ARENA 3/5] Balayage COMPRESS_WORTHWHILE_RATIO (matrice moyennement compressible, EPS=$EPS) : N=$N NB=$NB NRHS=$NRHS TOL=$TOL grille=${PR}x${PC}"
echo "RATIO testes : ${RATIO_LIST[*]}"

for RATIO in "${RATIO_LIST[@]}"; do
  echo "=== RATIO=$RATIO ==="
  log=$(run_and_log "arena_ratio_R${RATIO}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS" "$RATIO")

  rmax=$(grep_num "RMAX =" "$log")
  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem_theo=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  rss=$(floats_on_line "Memoire ANALYTIQUE" "$log" | tail -1)
  ratio_rss=$(grep_num "Ratio RSS" "$log")
  err=$(grep_num "Erreur relative" "$log")

  echo "$RATIO,$rmax,$EPS,$N,$NB,$NRHS,$TOL,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem_theo,$rss,$ratio_rss,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
