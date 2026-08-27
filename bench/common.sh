#!/bin/bash
#===============================================================================
# Fonctions partagees par les scripts bench/sweep_*.sh -- a SOURCER (jamais
# executer directement) depuis un script lance avec cwd = blr_solver/.
#
# ADAPTATION CLUSTER : quand ces scripts tournent enfin sur un vrai cluster
# multi-noeuds (et pas seulement cette machine WSL2 de developpement), deux
# choses a verifier/adapter ICI :
#   1) CORES_PER_NODE ci-dessous doit refleter les coeurs REELLEMENT
#      alloues a CE job (respecter un scheduler comme Slurm --cpus-per-task,
#      pas juste `nproc` qui peut voir le noeud entier meme si le job n'en a
#      qu'une partie) -- surcharger la variable d'environnement du meme nom
#      avant d'appeler un script (`CORES_PER_NODE=16 ./sweep_blocksize.sh`).
#   2) `mpirun -np N` ci-dessous peut devoir devenir `srun -n N` (ou
#      `mpirun` avec les options d'un hostfile) selon le lanceur du cluster
#      -- un seul point de changement si besoin : la fonction mpi_run().
#
# Chaque run individuel reste protege par un timeout (RUN_TIMEOUT) --
# discipline du projet : ne jamais laisser un calcul tourner dans le vide.
# Les scripts appelants peuvent l'augmenter pour un vrai balayage large
# echelle (`RUN_TIMEOUT=1800 ./sweep_scaling.sh`).
#===============================================================================

BLR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOGDIR="$BLR_ROOT/bench/logs"
mkdir -p "$LOGDIR"

# Source l'environnement Intel oneAPI (mpiifx/mpiicx/MKL) une seule fois --
# sans effet si deja fait par l'appelant (source est idempotent ici).
# `set +u`/`set -u` autour : setvars.sh reference des variables non
# definies en interne -- sous `set -u` (active par les scripts sweep_*.sh
# AVANT de sourcer ce fichier), ca fait avorter silencieusement TOUT le
# script appelant (aucune sortie, exit 0) -- piege classique des scripts
# d'environnement vendeur non prevus pour tourner sous `set -u`.
if [ -f /opt/intel/oneapi/setvars.sh ]; then
  set +u
  source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1
  set -u
fi

CORES_PER_NODE="${CORES_PER_NODE:-$(nproc)}"
RUN_TIMEOUT="${RUN_TIMEOUT:-300}"

# NOTE ADAPTATION CLUSTER : tous les scripts sweep_*.sh appellent
# `mpirun -np N ...` directement (pas de fonction wrapper -- `timeout`,
# utilise par run_and_log ci-dessous, execute ses arguments comme un
# PROCESSUS, il ne peut pas voir une fonction shell). Sur un cluster
# Slurm ou `mpirun` n'est pas le lanceur voulu, remplacer par un simple
# chercher/remplacer "mpirun -np" -> "srun -n" dans les 3 scripts.

# run_and_log LABEL CMD... -- execute CMD, journalise stdout+stderr dans
# bench/logs/LABEL.log, protege par RUN_TIMEOUT. Affiche le chemin du log
# (a capturer par l'appelant) ; signale en stderr tout timeout/echec SANS
# arreter le script (un point de balayage rate ne doit pas perdre tous les
# autres -- les valeurs correspondantes resteront simplement vides dans le
# CSV, visibles et non silencieuses).
run_and_log() {
  local label="$1"; shift
  local logfile="$LOGDIR/${label}.log"
  timeout "$RUN_TIMEOUT" "$@" > "$logfile" 2>&1
  local rc=$?
  if [ $rc -eq 124 ]; then
    echo "  [TIMEOUT apres ${RUN_TIMEOUT}s] $label" >&2
  elif [ $rc -ne 0 ]; then
    echo "  [ECHEC rc=$rc] $label (voir $logfile)" >&2
  fi
  echo "$logfile"
}

# _extract_num -- lit une ligne sur stdin, renvoie le premier nombre
# (entier/flottant/notation E) APRES le DERNIER '=' de la ligne. Utilise
# par grep_num/grep_num_after ci-dessous.
_extract_num() {
  sed -E 's/^.*=[[:space:]]*//' | grep -oE '[-+]?[0-9]+\.?[0-9]*([EeDd][-+]?[0-9]+)?' | head -1
}

# grep_num PATTERN FILE -- trouve la PREMIERE ligne de FILE contenant
# PATTERN (chaine litterale), en extrait le nombre suivant le dernier '='.
# Vide (pas d'erreur) si PATTERN absent -- ex: run tue par timeout avant
# d'atteindre cette ligne, a gerer cote appelant (CSV avec case vide).
grep_num() {
  local pattern="$1" file="$2"
  grep -m1 -F -- "$pattern" "$file" 2>/dev/null | _extract_num
}

# grep_num_after MARKER PATTERN FILE -- comme grep_num, mais ne cherche
# PATTERN qu'APRES la premiere ligne contenant MARKER (utile quand le meme
# libelle PATTERN apparait plusieurs fois dans un log, ex: bench_mpi.f90
# imprime "Erreur relative" a la fois pour l'ancien chemin [ANCIEN] et le
# nouveau [NOUVEAU] -- MARKER = le titre de bloc juste au-dessus permet de
# viser la bonne occurrence sans ambiguite).
grep_num_after() {
  local marker="$1" pattern="$2" file="$3"
  local ln
  ln=$(grep -n -F -- "$marker" "$file" 2>/dev/null | head -1 | cut -d: -f1)
  [ -z "$ln" ] && return
  tail -n "+$ln" "$file" | grep -m1 -F -- "$pattern" | _extract_num
}

# floats_on_line PATTERN FILE -- renvoie TOUS les nombres flottants
# presents sur la PREMIERE ligne contenant PATTERN, un par ligne (utile
# pour les lignes qui rapportent plusieurs valeurs en une fois, ex:
# "Temps assemblage = X s , factorisation = Y s , resolution = Z s").
floats_on_line() {
  local pattern="$1" file="$2"
  grep -m1 -F -- "$pattern" "$file" 2>/dev/null | grep -oE '[-+]?[0-9]+\.[0-9]+([EeDd][-+]?[0-9]+)?'
}

# omp_for_nprocs NPROCS -- calcule OMP_NUM_THREADS = coeurs_par_noeud /
# processus_MPI_par_noeud (convention etablie et validee via
# REPORT_THREAD_AFFINITY, cf. mpi_grid_mod.f90) ; jamais < 1.
omp_for_nprocs() {
  local nprocs="$1"
  local t=$((CORES_PER_NODE / nprocs))
  [ "$t" -lt 1 ] && t=1
  echo "$t"
}
