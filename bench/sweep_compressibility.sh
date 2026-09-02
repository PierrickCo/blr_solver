#!/bin/bash
#===============================================================================
# Balayage de la COMPRESSIBILITE de la matrice (parametre EPS -- bruit
# plein-rang additif, voir test/bench_blr_mpi_pivot_compressibility.f90), a
# N/NB/NRHS/TOL/grille FIXES, sur le solveur BLR+MPI de reference AVEC PIVOT
# (MPI_BLR_MATRIX_LU_FACTOR_PIVOT/_SOLVE_PIVOT_2D -- JAMAIS modifie).
#
# EPS=0 : noyau de Cauchy amorti pur, tres compressible. EPS croissant :
# degrade progressivement la compressibilite (de plus en plus de blocs
# tombent en repli dense a TOL fixe) -- montre l'impact sur RSS reelle,
# memoire theorique (+ratio des deux), precision numerique et temps par
# etape (assemblage/factorisation/resolution).
#
# Usage : ./bench/sweep_compressibility.sh [N] [NB] [NRHS] [TOL] [PR] [PC] [EPS_LIST]
#   defauts : N=8000 NB=200 NRHS=8 TOL=1e-4 PR=2 PC=2
#   EPS_LIST (entre guillemets) : defaut "0 1e-4 1e-3 1e-2 1e-1 1 10"
#
# Sortie : bench/results_compressibility.csv (EPS,N,NB,NRHS,TOL,
#   compression_pct,avg_rank,n_swaps,t_assemble,t_factor,t_solve,
#   mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur)
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-8000}"
NB="${2:-200}"
NRHS="${3:-8}"
TOL="${4:-1e-4}"
PR="${5:-2}"
PC="${6:-2}"
EPS_LIST_STR="${7:-0 1e-4 1e-3 1e-2 1e-1 1 10}"
read -r -a EPS_LIST <<< "$EPS_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

BIN="./bin/bench_blr_mpi_pivot_compressibility"
if [ ! -x "$BIN" ]; then
  echo "$BIN absent -- compilation (make bin/bench_blr_mpi_pivot_compressibility)..."
  make -f Makefile bin/bench_blr_mpi_pivot_compressibility
fi

CSV="bench/results_compressibility.csv"
echo "EPS,N,NB,NRHS,TOL,compression_pct,avg_rank,n_swaps,t_assemble,t_factor,t_solve,mem_theorique_mb,rss_mb,ratio_rss_theorique,erreur" > "$CSV"
echo "Balayage compressibilite (EPS) : N=$N NB=$NB NRHS=$NRHS TOL=$TOL grille=${PR}x${PC} OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "EPS testees : ${EPS_LIST[*]}"

for EPS in "${EPS_LIST[@]}"; do
  echo "=== EPS=$EPS ==="
  log=$(run_and_log "compressibility_EPS${EPS}" mpirun -np "$NPROCS" "$BIN" "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL" "$EPS")

  comp_pct=$(floats_on_line "Compression :" "$log" | tail -1)
  avg_rank=$(grep_num "Rang moyen =" "$log")
  n_swaps=$(grep_num "Nombre de swaps" "$log")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem_theo=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  rss=$(floats_on_line "Memoire ANALYTIQUE" "$log" | tail -1)
  ratio=$(grep_num "Ratio RSS" "$log")
  err=$(grep_num "Erreur relative" "$log")

  echo "$EPS,$N,$NB,$NRHS,$TOL,$comp_pct,$avg_rank,$n_swaps,$t_asm,$t_fac,$t_sol,$mem_theo,$rss,$ratio,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
