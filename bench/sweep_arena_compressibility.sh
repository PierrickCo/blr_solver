#!/bin/bash
#===============================================================================
# BENCHMARK ARENA 1/5 -- balayage de la COMPRESSIBILITE (EPS), UNIQUEMENT
# autour de la zone MARGINALE (transition compressible/incompressible) --
# reconnaissance prealable (voir memoire du projet) : sous RATIO=1/3 par
# defaut, N=8000, NB=200, cette matrice passe de 100% compresse a EPS=1e-7
# a ~0% a EPS=1e-5, la transition la plus raide se situant entre 1e-7 et
# 1e-6 -- le balayage par defaut se concentre sur cette fenetre.
#
# Usage : ./bench/sweep_arena_compressibility.sh [N] [NB] [NRHS] [TOL] [RATIO] [PR] [PC] [EPS_LIST]
#   defauts : N=8000 NB=200 NRHS=8 TOL=1e-4 RATIO=0.333333 PR=2 PC=2
#   EPS_LIST (entre guillemets) : defaut "1e-7 3e-7 6e-7 1e-6 2e-6 4e-6 6e-6 1e-5"
#
# Sortie : bench/results_arena_compressibility.csv (EPS,N,NB,NRHS,TOL,RATIO,
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
RATIO="${5:-0.333333}"
PR="${6:-2}"
PC="${7:-2}"
EPS_LIST_STR="${8:-1e-7 3e-7 6e-7 1e-6 2e-6 4e-6 6e-6 1e-5}"
read -r -a EPS_LIST <<< "$EPS_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

BIN="./bin/bench_blr_mpi_pivot_arena"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_arena)..."
  make -f Makefile bin/bench_blr_mpi_pivot_arena
fi

CSV="bench/results_arena_compressibility.csv"
echo "EPS,N,NB,NRHS,TOL,RATIO,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "[ARENA 1/5] Balayage compressibilite (zone marginale) : N=$N NB=$NB NRHS=$NRHS TOL=$TOL RATIO=$RATIO grille=${PR}x${PC}"
echo "EPS testees : ${EPS_LIST[*]}"

for EPS in "${EPS_LIST[@]}"; do
  echo "=== EPS=$EPS ==="
  log=$(run_and_log "arena_compressibility_EPS${EPS}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS" "$RATIO")

  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem_theo=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  rss=$(floats_on_line "Memoire ANALYTIQUE" "$log" | tail -1)
  ratio_rss=$(grep_num "Ratio RSS" "$log")
  err=$(grep_num "Erreur relative" "$log")

  echo "$EPS,$N,$NB,$NRHS,$TOL,$RATIO,$comp_pct,$avg_rank,$t_asm,$t_fac,$t_sol,$mem_theo,$rss,$ratio_rss,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
