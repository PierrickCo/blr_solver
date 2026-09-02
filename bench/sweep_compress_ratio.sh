#!/bin/bash
#===============================================================================
# Balayage de COMPRESS_WORTHWHILE_RATIO (blr_compress_mod.f90, voir
# SET_COMPRESS_WORTHWHILE_RATIO) a EPS/N/NB/NRHS/TOL/grille FIXES, DANS LE
# CAS MARGINAL identifie par sweep_compressibility.sh : EPS proche du point
# de transition compressible/incompressible, ou une grande partie des blocs
# atteignent un rang final proche du plafond RMAX -- c'est LA ou RATIO
# influe le plus (voir test/bench_blr_mpi_pivot_ratio.f90 pour la
# justification complete).
#
# Usage : ./bench/sweep_compress_ratio.sh [N] [NB] [NRHS] [TOL] [EPS] [PR] [PC] [RATIO_LIST]
#   defauts : N=8000 NB=200 NRHS=8 TOL=1e-4 EPS=1e-6 PR=2 PC=2
#   RATIO_LIST (entre guillemets) : defaut "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
#
# Sortie : bench/results_compress_ratio.csv (RATIO,RMAX,EPS,N,NB,NRHS,TOL,
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
RATIO_LIST_STR="${8:-0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9}"
read -r -a RATIO_LIST <<< "$RATIO_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

BIN="./bin/bench_blr_mpi_pivot_ratio"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_ratio)..."
  make -f Makefile bin/bench_blr_mpi_pivot_ratio
fi

CSV="bench/results_compress_ratio.csv"
echo "RATIO,RMAX,EPS,N,NB,NRHS,TOL,compression_pct,avg_rank,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "Balayage COMPRESS_WORTHWHILE_RATIO (cas marginal EPS=$EPS) : N=$N NB=$NB NRHS=$NRHS TOL=$TOL grille=${PR}x${PC} OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "RATIO testes : ${RATIO_LIST[*]}"

for RATIO in "${RATIO_LIST[@]}"; do
  echo "=== RATIO=$RATIO ==="
  log=$(run_and_log "compress_ratio_R${RATIO}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS" "$RATIO")

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
