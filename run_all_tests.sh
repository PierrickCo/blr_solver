#!/bin/bash
#===============================================================================
# Compile puis lance le jeu de tests de correction CURE de la ligne vivante
# (13 binaires -- voir INSTALL.md pour le detail individuel de chacun),
# avec un resume PASS/FAIL clair a la fin -- pour ne pas avoir a connaitre
# chaque nom de binaire/syntaxe pour verifier qu'une installation marche.
#
# Ne couvre PAS le code de legacy/ (voir legacy/README.md, meme principe
# manuellement si besoin) ni les scripts de performance de bench/ (voir
# CLUSTER.md -- ce script teste la CORRECTION, pas la performance).
#
# Usage : ./run_all_tests.sh [intel|gcc]   (depuis blr_solver/, apres avoir
#         source l'environnement -- voir INSTALL.md)
#   intel (defaut) : Makefile (ifx+MKL+Intel MPI), binaires dans bin/.
#   gcc            : Makefile.gcc (gfortran+BLAS/LAPACK+OpenMPI), binaires
#                     dans build-gcc/ -- chemin utilise par la CI (voir
#                     .github/workflows/ci.yml), aucune licence Intel requise.
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"

TOOLCHAIN="${1:-intel}"
if [ "$TOOLCHAIN" = "gcc" ]; then
  MAKEFILE="Makefile.gcc"
  BINDIR="build-gcc"
else
  MAKEFILE="Makefile"
  BINDIR="bin"
fi

echo "=== Compilation (chaine: $TOOLCHAIN) ==="
if ! make -f "$MAKEFILE" -j4; then
  echo "ECHEC DE COMPILATION -- voir la sortie ci-dessus."
  exit 1
fi
echo

LOGDIR="run_all_tests_logs"
mkdir -p "$LOGDIR"
RESULTS=()
TIMEOUT_S=120

# run_test NAME CMD... -- detecte PASS/FAIL/WARN via les marqueurs "OK"/
# "ECHEC"/"SUSPECT" que chaque binaire de test imprime lui-meme.
run_test() {
  local name="$1"; shift
  local log="$LOGDIR/$name.log"
  timeout "$TIMEOUT_S" "$@" > "$log" 2>&1
  local rc=$?
  if grep -q "ECHEC" "$log"; then
    RESULTS+=("FAIL  $name")
  elif grep -q "SUSPECT" "$log"; then
    RESULTS+=("WARN  $name (voir $log)")
  elif grep -q " OK" "$log"; then
    RESULTS+=("PASS  $name")
  elif [ $rc -ne 0 ]; then
    RESULTS+=("FAIL  $name (exit $rc, voir $log)")
  else
    RESULTS+=("?     $name (aucun verdict detecte, voir $log)")
  fi
}

# run_bench NAME CMD... -- pour bench_seq/bench_mpi, qui n'impriment pas de
# marqueur OK/ECHEC (ce sont des outils de mesure, pas des tests) : verifie
# juste que ca tourne sans planter ET que l'erreur relative rapportee reste
# raisonnable (<1e-3).
run_bench() {
  local name="$1"; shift
  local log="$LOGDIR/$name.log"
  timeout "$TIMEOUT_S" "$@" > "$log" 2>&1
  local rc=$?
  local err
  err=$(grep -m1 "Erreur relative" "$log" | grep -oE '[0-9]+\.[0-9]+[eE][-+]?[0-9]+' | tail -1)
  if [ $rc -ne 0 ]; then
    RESULTS+=("FAIL  $name (exit $rc, voir $log)")
  elif [ -n "$err" ] && awk -v e="$err" 'BEGIN{exit !(e+0<1e-3)}'; then
    RESULTS+=("PASS  $name (erreur=$err)")
  else
    RESULTS+=("?     $name (erreur=${err:-non trouvee}, voir $log)")
  fi
}

echo "=== Tests sequentiels / OpenMP ==="
run_test test_block_lu            "$BINDIR/test_block_lu" 500 32
run_test test_blr_compress        "$BINDIR/test_blr_compress"
run_test test_blr_workspace       "$BINDIR/test_blr_workspace"
run_test test_blr_matrix_progress "$BINDIR/test_blr_matrix_progress"
run_test test_openmp_dense        "$BINDIR/test_openmp_dense" 1000 100

echo "=== Tests MPI (grille 2x2) ==="
run_test test_mpi_grid             mpirun -np 4 "$BINDIR/test_mpi_grid" 2 2 500 50
run_test test_mpi_lu_pivot         mpirun -np 4 "$BINDIR/test_mpi_lu_pivot" 2 2 500 50
run_test test_mpi_lu_solve_2d      mpirun -np 4 "$BINDIR/test_mpi_lu_solve_2d" 2 2 500 50 20
run_test test_mpi_blr_matrix       mpirun -np 4 "$BINDIR/test_mpi_blr_matrix" 2 2 1000 100 20 1e-4
run_test test_mpi_blr_matrix_pivot mpirun -np 4 "$BINDIR/test_mpi_blr_matrix_pivot" 2 2 1000 100 20 1e-4
run_test test_mpi_blr_matrix_pivot_arena mpirun -np 4 "$BINDIR/test_mpi_blr_matrix_pivot_arena" 2 2 1000 100 20 1e-4
run_test test_mpi_blr_matrix_pivot_ooc mpirun -np 4 "$BINDIR/test_mpi_blr_matrix_pivot_ooc" 2 2 1000 100 20 1e-4 8 7 "$BINDIR/blr_ooc_scratch_ci"

echo "=== Benchmarks (fumee -- petite echelle) ==="
run_bench bench_seq  "$BINDIR/bench_seq" 1000 20 100 1e-4
run_bench bench_mpi  mpirun -np 4 "$BINDIR/bench_mpi" 2 2 1000 20 100 1e-4
run_test  bench_blr_pipeline "$BINDIR/bench_blr_pipeline" 1000 20 100 1e-4

echo
echo "=== Resume ==="
NFAIL=0
NWARN=0
for r in "${RESULTS[@]}"; do
  echo "  $r"
  case "$r" in
    FAIL*) NFAIL=$((NFAIL+1)) ;;
    WARN*|\?*) NWARN=$((NWARN+1)) ;;
  esac
done
echo

if [ "$NFAIL" -eq 0 ] && [ "$NWARN" -eq 0 ]; then
  echo "Tous les tests sont passes."
  exit 0
elif [ "$NFAIL" -eq 0 ]; then
  echo "Aucun echec, mais $NWARN resultat(s) a verifier manuellement (voir $LOGDIR/)."
  exit 0
else
  echo "$NFAIL test(s) en echec, $NWARN a verifier -- voir $LOGDIR/ pour le detail."
  exit 1
fi
