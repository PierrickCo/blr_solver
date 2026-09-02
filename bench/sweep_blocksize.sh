#!/bin/bash
#===============================================================================
# Balayage de la taille de bloc NB, a N/NRHS/TOL/grille fixes, sur les 3
# solveurs qui exposent NB comme parametre : A et B (OpenMP, meme binaire
# bench_blr_pipeline, A=TOL>0 compresse / B=TOL=0 plein bloc), MPI dense+2D
# (bench_mpi, cas 1 SEULEMENT -- cas 2/3 ignores, cf. commentaire plus bas),
# et BLR+MPI sans pivot (test_mpi_blr_matrix). Generalise le balayage
# manuel fait en session (NB=100/300/500/750/1000 a N=10000/20000, qui a
# trouve NB~500 optimal pour ce dernier -- voir memoire projet) pour le
# rejouer facilement a d'autres N/machines/clusters.
#
# Usage : ./bench/sweep_blocksize.sh [N] [NRHS] [TOL] [PR] [PC] [NB_LIST]
#   N, NRHS, TOL, PR, PC : comme d'habitude (defauts 10000 100 1e-4 2 2)
#   NB_LIST : liste de tailles de bloc separees par des espaces, ENTRE
#     GUILLEMETS (defaut "50 100 200 300 500 750 1000 1500 2000")
#   Exemple sur un vrai cluster, noeud a 32 coeurs, grille 4x4 :
#     CORES_PER_NODE=32 ./bench/sweep_blocksize.sh 20000 100 1e-4 4 4
#
# ATTENTION cout : bench_mpi.f90 lance TOUJOURS ses 3 cas (dense+pivot,
# BLR+pivot dense-puis-compression, BLR+pivot compresse) -- seul le cas 1
# nous interesse ici (reference "MPI"), mais les cas 2/3 s'executent quand
# meme et peuvent etre lents a grand N (deja observe : cas 2 peut depasser
# 10 min a N=20000). RUN_TIMEOUT (defaut 300s, cf. common.sh) coupe le
# calcul avant que ca ne tourne dans le vide -- le cas 1 a TOUJOURS fini et
# ete imprime avant que le cas 2 ne demarre, donc ses chiffres restent
# fiables meme si le process entier est tue ensuite.
#
# Sortie : bench/results_blocksize.csv
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source bench/common.sh

N="${1:-10000}"
NRHS="${2:-100}"
TOL="${3:-1e-4}"
PR="${4:-2}"
PC="${5:-2}"
NB_LIST_STR="${6:-50 100 200 300 500 750 1000 1500 2000}"
read -r -a NB_LIST <<< "$NB_LIST_STR"

NPROCS=$((PR*PC))
export OMP_NUM_THREADS
OMP_NUM_THREADS=$(omp_for_nprocs "$NPROCS")

CSV="bench/results_blocksize.csv"
echo "solver,N,NB,NRHS,TOL,t_assemble,t_factor,t_solve,mem_mb,error" > "$CSV"
echo "Balayage NB : N=$N NRHS=$NRHS TOL=$TOL grille=${PR}x${PC} OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "NB testes : ${NB_LIST[*]}"

for NB in "${NB_LIST[@]}"; do
  echo "=== NB=$NB ==="

  # --- A : OpenMP, compresse (TOL>0) ---
  log=$(run_and_log "blocksize_A_NB${NB}" ./bin/bench_blr_pipeline "$N" "$NRHS" "$NB" "$TOL")
  t_asm=$(grep_num "Temps assemblage+compression" "$log")
  t_fac=$(grep_num "Temps factorisation =" "$log")
  t_sol=$(grep_num "Temps resolution =" "$log")
  mem=$(grep_num "Memoire des facteurs L/U" "$log")
  err=$(grep_num "Erreur relative ||Xcalc-Xexact" "$log")
  echo "A,$N,$NB,$NRHS,$TOL,$t_asm,$t_fac,$t_sol,$mem,$err" >> "$CSV"

  # --- B : OpenMP, plein bloc (TOL=0) ---
  log=$(run_and_log "blocksize_B_NB${NB}" ./bin/bench_blr_pipeline "$N" "$NRHS" "$NB" 0)
  t_asm=$(grep_num "Temps assemblage+compression" "$log")
  t_fac=$(grep_num "Temps factorisation =" "$log")
  t_sol=$(grep_num "Temps resolution =" "$log")
  mem=$(grep_num "Memoire des facteurs L/U" "$log")
  err=$(grep_num "Erreur relative ||Xcalc-Xexact" "$log")
  echo "B,$N,$NB,$NRHS,0,$t_asm,$t_fac,$t_sol,$mem,$err" >> "$CSV"

  # --- MPI dense+2D (cas 1 de bench_mpi) ---
  log=$(run_and_log "blocksize_MPI_NB${NB}" mpirun -np "$NPROCS" ./bin/bench_mpi "$PR" "$PC" "$N" "$NRHS" "$NB" "$TOL")
  t_asm=$(grep_num "Temps assemblage (max/processus)" "$log")
  t_fac=$(grep_num "Temps factorisation (hors compression, distribuee)" "$log")
  t_sol=$(grep_num "Temps resolution (2D)" "$log")
  mem=$(grep_num "Memoire PAR PROCESSUS pendant factor." "$log")
  err=$(grep_num_after "resolution DISTRIBUEE" "Erreur relative" "$log")
  echo "MPI,$N,$NB,$NRHS,$TOL,$t_asm,$t_fac,$t_sol,$mem,$err" >> "$CSV"

  # --- BLR+MPI sans pivot (test_mpi_blr_matrix) ---
  log=$(run_and_log "blocksize_BLRMPI_NB${NB}" mpirun -np "$NPROCS" ./bin/test_mpi_blr_matrix "$PR" "$PC" "$N" "$NB" "$NRHS" "$TOL")
  read -r t_asm t_fac t_sol <<< "$(floats_on_line "Temps assemblage" "$log" | tr '\n' ' ')"
  mem=$(floats_on_line "Memoire ANALYTIQUE" "$log" | head -1)
  err=$(grep_num "Erreur relative ||X_distribue" "$log")
  echo "BLR+MPI,$N,$NB,$NRHS,$TOL,$t_asm,$t_fac,$t_sol,$mem,$err" >> "$CSV"
done

echo
echo "Resultats : $CSV"
