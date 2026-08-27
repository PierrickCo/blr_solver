#!/bin/bash
#===============================================================================
# Scalabilite FORTE (taille de probleme fixe, ressources croissantes) des 4
# solveurs de reference :
#   - A/B/C (OpenMP pur, un seul processus) : on fait varier OMP_NUM_THREADS.
#   - MPI dense+2D et BLR+MPI (distribues) : on fait varier la grille PRxPC
#     (donc le nombre de processus), avec OMP_NUM_THREADS RECALCULE a
#     chaque config = coeurs_par_noeud / processus (cf. REPORT_THREAD_
#     AFFINITY, mpi_grid_mod.f90 -- chaque run affiche son propre diagnostic
#     d'affinite, a verifier surtout la premiere fois sur un vrai cluster).
#
# Usage : ./bench/sweep_scaling.sh [N] [NRHS] [NB] [TOL]
#   defauts : N=10000 NRHS=100 NB=500 TOL=1e-4 (NB=500 = optimum trouve en
#   session pour le module BLR+MPI a cette echelle -- ajuster si N change
#   beaucoup, cf. sweep_blocksize.sh)
#
# Grilles MPI testees (filtrees a <= CORES_PER_NODE) : 1x1 1x2 2x1 2x2 1x4
#   4x1 2x3 3x2 1x8 8x1 4x2 2x4 -- adapter GRID_LIST ci-dessous pour un
#   vrai cluster multi-noeuds (grilles plus grandes, PR*PC > coeurs d'UN
#   noeud -- alors OMP_NUM_THREADS doit rester borne par les coeurs D'UN
#   SEUL noeud, pas le total cluster : voir la memoire projet sur
#   OMP_NUM_THREADS/affinite avant de pousser au-dela d'un noeud).
#
# Sorties : bench/results_scaling_openmp.csv, bench/results_scaling_mpi.csv
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-10000}"
NRHS="${2:-100}"
NB="${3:-500}"
TOL="${4:-1e-4}"

# --- (a) balayage OpenMP : threads = 1,2,4,8,... jusqu'a CORES_PER_NODE ---
THREADS_LIST=()
t=1
while [ "$t" -le "$CORES_PER_NODE" ]; do
  THREADS_LIST+=("$t")
  t=$((t*2))
done
[ "${THREADS_LIST[-1]}" -ne "$CORES_PER_NODE" ] && THREADS_LIST+=("$CORES_PER_NODE")

CSV_OMP="bench/results_scaling_openmp.csv"
echo "solver,threads,N,NRHS,NB,TOL,t_factor,t_solve,mem_mb,error,speedup_vs_1thread" > "$CSV_OMP"
echo "Scalabilite OpenMP : N=$N NRHS=$NRHS NB=$NB TOL=$TOL, threads testes : ${THREADS_LIST[*]}"

BASE_A=""; BASE_B=""; BASE_C=""

for th in "${THREADS_LIST[@]}"; do
  export OMP_NUM_THREADS=$th
  echo "=== OMP_NUM_THREADS=$th ==="

  log=$(run_and_log "scaling_A_T${th}" ./bench_blr_pipeline "$N" "$NRHS" "$NB" "$TOL")
  t_fac=$(grep_num "Temps factorisation =" "$log"); t_sol=$(grep_num "Temps resolution =" "$log")
  mem=$(grep_num "Memoire des facteurs L/U" "$log"); err=$(grep_num "Erreur relative ||Xcalc-Xexact" "$log")
  [ "$th" -eq 1 ] && BASE_A=$t_fac
  sp=$(awk -v b="${BASE_A:-$t_fac}" -v x="$t_fac" 'BEGIN{if(x>0) printf "%.2f", b/x; else print ""}')
  echo "A,$th,$N,$NRHS,$NB,$TOL,$t_fac,$t_sol,$mem,$err,$sp" >> "$CSV_OMP"

  log=$(run_and_log "scaling_B_T${th}" ./bench_blr_pipeline "$N" "$NRHS" "$NB" 0)
  t_fac=$(grep_num "Temps factorisation =" "$log"); t_sol=$(grep_num "Temps resolution =" "$log")
  mem=$(grep_num "Memoire des facteurs L/U" "$log"); err=$(grep_num "Erreur relative ||Xcalc-Xexact" "$log")
  [ "$th" -eq 1 ] && BASE_B=$t_fac
  sp=$(awk -v b="${BASE_B:-$t_fac}" -v x="$t_fac" 'BEGIN{if(x>0) printf "%.2f", b/x; else print ""}')
  echo "B,$th,$N,$NRHS,$NB,0,$t_fac,$t_sol,$mem,$err,$sp" >> "$CSV_OMP"

  log=$(run_and_log "scaling_C_T${th}" ./bench_seq "$N" "$NRHS" "$N" "$TOL")
  t_fac=$(grep_num_after "Direct dense, sans decoupage" "Temps factorisation (hors compression)" "$log")
  t_sol=$(grep_num_after "Direct dense, sans decoupage" "Temps resolution" "$log")
  mem=$(grep_num_after "Direct dense, sans decoupage" "Memoire analytique" "$log")
  err=$(grep_num_after "Direct dense, sans decoupage" "Erreur relative" "$log")
  [ "$th" -eq 1 ] && BASE_C=$t_fac
  sp=$(awk -v b="${BASE_C:-$t_fac}" -v x="$t_fac" 'BEGIN{if(x>0) printf "%.2f", b/x; else print ""}')
  echo "C,$th,$N,$NRHS,$N,$TOL,$t_fac,$t_sol,$mem,$err,$sp" >> "$CSV_OMP"
done
echo "Resultats OpenMP : $CSV_OMP"
echo

# --- (b) balayage grille MPI : PRxPC, filtree a PR*PC <= CORES_PER_NODE ---
GRID_LIST="1x1 1x2 2x1 2x2 1x4 4x1 2x3 3x2 1x8 8x1 4x2 2x4"
CSV_MPI="bench/results_scaling_mpi.csv"
echo "solver,grid,nprocs,N,NRHS,NB,TOL,t_factor,t_solve,mem_mb,error,speedup_vs_min_procs" > "$CSV_MPI"
echo "Scalabilite MPI : N=$N NRHS=$NRHS NB=$NB TOL=$TOL, grilles candidates : $GRID_LIST (filtrees a <= $CORES_PER_NODE coeurs)"

BASEMPI_A=""; BASEMPI_B=""; FIRST_MPI=1; FIRST_BLRMPI=1

for grid in $GRID_LIST; do
  PR="${grid%x*}"; PC="${grid#*x}"
  NPROCS=$((PR*PC))
  [ "$NPROCS" -gt "$CORES_PER_NODE" ] && continue
  export OMP_NUM_THREADS
  OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")
  echo "=== grille ${PR}x${PC} ($NPROCS proc, OMP_NUM_THREADS=$OMP_NUM_THREADS) ==="

  log=$(run_and_log "scaling_MPI_${PR}x${PC}" mpirun -np "$NPROCS" ./bench_mpi "$PR" "$PC" "$N" "$NRHS" "$NB" "$TOL")
  t_fac=$(grep_num "Temps factorisation (hors compression, distribuee)" "$log")
  t_sol=$(grep_num "Temps resolution (2D)" "$log")
  mem=$(grep_num "Memoire PAR PROCESSUS pendant factor." "$log")
  err=$(grep_num_after "resolution DISTRIBUEE" "Erreur relative" "$log")
  if [ "$FIRST_MPI" -eq 1 ] && [ -n "$t_fac" ]; then BASEMPI_A="$t_fac"; FIRST_MPI=0; fi
  sp=$(awk -v b="${BASEMPI_A:-$t_fac}" -v x="$t_fac" 'BEGIN{if(x>0) printf "%.2f", b/x; else print ""}')
  echo "MPI,${PR}x${PC},$NPROCS,$N,$NRHS,$NB,$TOL,$t_fac,$t_sol,$mem,$err,$sp" >> "$CSV_MPI"

  log=$(run_and_log "scaling_BLRMPI_${PR}x${PC}" mpirun -np "$NPROCS" ./test_mpi_blr_matrix "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL")
  read -r _ t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  err=$(grep_num "Erreur relative ||X_distribue" "$log")
  if [ "$FIRST_BLRMPI" -eq 1 ] && [ -n "$t_fac" ]; then BASEMPI_B="$t_fac"; FIRST_BLRMPI=0; fi
  sp=$(awk -v b="${BASEMPI_B:-$t_fac}" -v x="$t_fac" 'BEGIN{if(x>0) printf "%.2f", b/x; else print ""}')
  echo "BLR+MPI,${PR}x${PC},$NPROCS,$N,$NRHS,$NB,$TOL,$t_fac,$t_sol,$mem,$err,$sp" >> "$CSV_MPI"
done
echo "Resultats MPI : $CSV_MPI"
