#!/bin/bash
#===============================================================================
# Recherche des meilleurs parametres des solveurs BLR+MPI EN MEMOIRE (LU et
# LDLt, variantes classique et arena -- pas d'out-of-core) sur la matrice de
# demo complexe symetrique (ldlt/src/sym_testmat_mod.f90), via le pilote
# ldlt/bin/sweep_driver. Concu pour etre lance tel quel sur un cluster.
#
# Balayage en phases successives, chaque phase fixant le meilleur reglage
# trouve avant de passer a la suivante (un balayage complet serait
# combinatoire) :
#   1. TOL           : le plus GRAND TOL dont l'erreur (mediane) <= ERR_MAX
#   2. NB            : taille de bloc minimisant l'objectif (OBJECTIVE)
#   3. grille/threads: PR x PC x OMP_NUM_THREADS = CORES, plusieurs formes
#   4. RATIO         : COMPRESS_WORTHWHILE_RATIO (seuil de compression)
#   5. arena         : variante arena au meilleur reglage (si TEST_ARENA=1)
# Chaque configuration est lancee REPEATS fois ; on retient la MEDIANE (un
# run isole peut etre tres bruite). Compression : SVD randomisee (defaut).
#
# Sorties (dans OUTDIR) : runs.csv (chaque run), summary.csv (medianes par
# configuration), best_params.txt (meilleurs parametres + commandes
# conseillees), logs/ (sortie brute de chaque run).
#
# Usage (depuis blr_solver/, apres `make && make -C ldlt`) :
#   ./tools/sweep_incore.sh                         # defauts ci-dessous
#   N=100000 CORES=300 LAUNCH="srun -n" LAUNCH_THREADS_OPT="--cpus-per-task=" ./tools/sweep_incore.sh
#   ./tools/sweep_incore.sh --report OUTDIR         # regenere best_params.txt
#
# Variables (toutes optionnelles) :
#   N=20000            taille de la matrice
#   NRHS=8             colonnes de second membre
#   CORES=$(nproc)     coeurs a utiliser (= MPI x OMP_NUM_THREADS)
#   SOLVERS="lu ldlt"  solveurs a regler (ldlt : matrices Aᵀ=A seulement)
#   TOL_LIST="1e-3 1e-4 1e-5 1e-6"
#   NB_LIST=auto       auto : NB tels que 20 <= N/NB <= 100 parmi
#                      200 300 400 600 800 1000 1500 2000 3000 4000 5000
#   THREADS_LIST="1 2 4"   OMP_NUM_THREADS essayes (diviseurs de CORES)
#   BASE_THREADS=2     threads utilises pendant les phases 1-2
#   RATIO_LIST="0.2 0.333 0.5"
#   ERR_MAX=1e-3       erreur relative maximale acceptee ||X-Xexact||/||Xexact||
#   OBJECTIVE=total    total (assemblage+factorisation+resolution) | factor
#                      | asm_fact
#   REPEATS=3          repetitions par configuration (mediane)
#   MIN_GAIN=0.03      phases 3-4 : on ne quitte le reglage par defaut (grille
#                      carree, RATIO=1/3) que pour un gain > MIN_GAIN (anti-bruit)
#   TEST_ARENA=1       phase 5 (variante arena)
#   LAUNCH="mpirun -np"  lanceur MPI ("srun -n" sous Slurm)
#   LAUNCH_THREADS_OPT=""  option du lanceur suivie du nombre de threads, ex.
#                      "--cpus-per-task=" sous Slurm (vide pour mpirun)
#   NEAR=0.05          phase NB : parmi les NB a moins de NEAR du meilleur
#                      temps, on retient celui qui consomme le moins de memoire
#   MPI_EXTRA=""       options supplementaires du lanceur (ex. --oversubscribe)
#   TIMEOUT=3600       limite par run (s)
#   OUTDIR=sweep_results_<date>
#   DRIVER=ldlt/bin/sweep_driver
#===============================================================================
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."

N=${N:-20000}
NRHS=${NRHS:-8}
CORES=${CORES:-$(nproc)}
SOLVERS=${SOLVERS:-"lu ldlt"}
TOL_LIST=${TOL_LIST:-"1e-3 1e-4 1e-5 1e-6"}
NB_LIST=${NB_LIST:-auto}
THREADS_LIST=${THREADS_LIST:-"1 2 4"}
BASE_THREADS=${BASE_THREADS:-2}
RATIO_LIST=${RATIO_LIST:-"0.2 0.333 0.5"}
ERR_MAX=${ERR_MAX:-1e-3}
OBJECTIVE=${OBJECTIVE:-total}
REPEATS=${REPEATS:-3}
MIN_GAIN=${MIN_GAIN:-0.03}
TEST_ARENA=${TEST_ARENA:-1}
LAUNCH=${LAUNCH:-"mpirun -np"}
MPI_EXTRA=${MPI_EXTRA:-""}
TIMEOUT=${TIMEOUT:-3600}
LAUNCH_THREADS_OPT=${LAUNCH_THREADS_OPT:-""}
NEAR=${NEAR:-0.05}
REPORT_ONLY=0
if [ "${1:-}" = "--report" ]; then
  REPORT_ONLY=1
  OUTDIR=${2:?usage : $0 --report OUTDIR}
  [ -f "$OUTDIR/decisions.txt" ] || { echo "$OUTDIR/decisions.txt absent"; exit 1; }
  # parametres du balayage d'origine (ecrits en tete de decisions.txt)
  eval "$(grep '^#PARAMS ' "$OUTDIR/decisions.txt" | sed 's/^#PARAMS //')"
fi
OUTDIR=${OUTDIR:-sweep_results_$(date +%Y%m%d_%H%M%S)}
DRIVER=${DRIVER:-ldlt/bin/sweep_driver}

if [ "$REPORT_ONLY" = "0" ] && [ ! -x "$DRIVER" ]; then
  echo "Pilote $DRIVER absent -- compilation (make && make -C ldlt)..."
  make -j4 >/dev/null && make -C ldlt -j4 >/dev/null || { echo "ECHEC de compilation"; exit 1; }
fi

if [ "$NB_LIST" = "auto" ]; then
  NB_LIST=""
  for nb in 200 300 400 600 800 1000 1500 2000 3000 4000 5000; do
    r=$(( N / nb ))
    if [ $r -ge 20 ] && [ $r -le 100 ]; then NB_LIST="$NB_LIST $nb"; fi
  done
  NB_LIST=${NB_LIST# }
  [ -z "$NB_LIST" ] && NB_LIST="200 400 800"
fi

# Threads : seulement les diviseurs de CORES
TL=""
for t in $THREADS_LIST; do [ $(( CORES % t )) -eq 0 ] && TL="$TL $t"; done
THREADS_LIST=${TL# }
[ -z "$THREADS_LIST" ] && THREADS_LIST=1
[ $(( CORES % BASE_THREADS )) -ne 0 ] && BASE_THREADS=1

RUNS="$OUTDIR/runs.csv"
SUMMARY="$OUTDIR/summary.csv"
BEST="$OUTDIR/best_params.txt"
DECISIONS="$OUTDIR/decisions.txt"
if [ "$REPORT_ONLY" = "0" ]; then
mkdir -p "$OUTDIR/logs"
echo "phase,solver,pr,pc,threads,n,nb,nrhs,tol,ratio,rep,status,info,t_asm,t_fact,t_solve,t_total,mem_max_mb,mem_tot_mb,rss_max_mb,err" > "$RUNS"
echo "phase,solver,pr,pc,threads,n,nb,nrhs,tol,ratio,nok,t_asm,t_fact,t_solve,t_total,objective,mem_max_mb,mem_tot_mb,rss_max_mb,err" > "$SUMMARY"
fi

# Echec = valeur numerique BIG ("inf" vaudrait 0 pour awk).
BIG=1e300
median() { sort -g | awk '{a[NR]=$1} END{ if (NR==0) {print "nan"} else if (NR%2) {print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2} }'; }

# Grille "la plus carree" pour P processus, et deux variantes (plus large /
# plus haute) -- lignes "PR PC".
grid_shapes() {
  local P=$1 pr0=1 d
  for ((d=1; d*d<=P; d++)); do [ $((P % d)) -eq 0 ] && pr0=$d; done
  local lo=1 hi=$P
  for ((d=1; d<=pr0/2; d++)); do [ $((P % d)) -eq 0 ] && lo=$d; done
  for ((d=P; d>=2*pr0; d--)); do [ $((P % d)) -eq 0 ] && hi=$d; done
  { echo "$pr0 $((P/pr0))"; echo "$lo $((P/lo))"; [ $hi -le $P ] && echo "$hi $((P/hi))"; } | awk '!seen[$0]++'
}
square_grid() { grid_shapes "$1" | head -1; }

# run_cfg PHASE SOLVER PR PC THREADS NB TOL RATIO
# Lance REPEATS fois, ajoute chaque run a runs.csv, la mediane a summary.csv,
# et imprime "objectif err" (objectif = inf si un run a echoue).
run_cfg() {
  local phase=$1 solver=$2 pr=$3 pc=$4 th=$5 nb=$6 tol=$7 ratio=$8
  local np=$((pr*pc)) tmp; tmp=$(mktemp)
  local rep line status
  for ((rep=1; rep<=REPEATS; rep++)); do
    local log="$OUTDIR/logs/${phase}_${solver}_${pr}x${pc}_t${th}_nb${nb}_tol${tol}_r${ratio}_rep${rep}.log"
    OMP_NUM_THREADS=$th timeout "$TIMEOUT" $LAUNCH $np ${LAUNCH_THREADS_OPT:+${LAUNCH_THREADS_OPT}$th} $MPI_EXTRA "$DRIVER" "$solver" $pr $pc $N $nb $NRHS $tol $ratio 1 < /dev/null > "$log" 2>&1
    line=$(grep -m1 '^SWEEP_RESULT' "$log" | sed 's/= */=/g')
    if [ -n "$line" ]; then
      echo "$line" | awk -v phase="$phase" -v rep="$rep" '
        { for (i=2;i<=NF;i++) { split($i,kv,"="); v[kv[1]]=kv[2] }
          st = (v["info"]==0 && v["err"]==v["err"]+0) ? "ok" : "fail"
          tot = v["t_asm"]+v["t_fact"]+v["t_solve"]
          printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%s,%.4f,%s,%s,%s,%s\n", phase, v["solver"], v["pr"], v["pc"],
                 v["threads"], v["n"], v["nb"], v["nrhs"], v["tol"], v["ratio"], rep, st, v["info"], v["t_asm"],
                 v["t_fact"], v["t_solve"], tot, v["mem_max_mb"], v["mem_tot_mb"], v["rss_max_mb"], v["err"] }' >> "$tmp"
    else
      echo "$phase,$solver,$pr,$pc,$th,$N,$nb,$NRHS,$tol,$ratio,$rep,fail,,,,,,,,," >> "$tmp"
    fi
  done
  cat "$tmp" >> "$RUNS"
  local nok; nok=$(awk -F, '$12=="ok"' "$tmp" | wc -l)
  if [ "$nok" -lt "$REPEATS" ]; then
    echo "$phase,$solver,$pr,$pc,$th,$N,$nb,$NRHS,$tol,$ratio,$nok,,,,,$BIG,,,," >> "$SUMMARY"
    echo "$BIG $BIG $BIG"; rm -f "$tmp"; return
  fi
  local ta tf ts tt mm mt rss er obj
  ta=$(cut -d, -f14 "$tmp" | median); tf=$(cut -d, -f15 "$tmp" | median); ts=$(cut -d, -f16 "$tmp" | median)
  tt=$(cut -d, -f17 "$tmp" | median); mm=$(cut -d, -f18 "$tmp" | median); mt=$(cut -d, -f19 "$tmp" | median)
  rss=$(cut -d, -f20 "$tmp" | median); er=$(cut -d, -f21 "$tmp" | median)
  case "$OBJECTIVE" in
    factor)   obj=$tf ;;
    asm_fact) obj=$(awk -v a="$ta" -v b="$tf" 'BEGIN{printf "%.4f", a+b}') ;;
    *)        obj=$tt ;;
  esac
  echo "$phase,$solver,$pr,$pc,$th,$N,$nb,$NRHS,$tol,$ratio,$nok,$ta,$tf,$ts,$tt,$obj,$mm,$mt,$rss,$er" >> "$SUMMARY"
  printf "    %-10s %3sx%-3s t=%-2s nb=%-5s tol=%-6s ratio=%-5s | asm %8.3f fact %8.3f solve %7.3f | mem %8.1f Mo/proc | err %s\n" \
         "$solver" $pr $pc $th $nb $tol $ratio $ta $tf $ts $mm $er >&2
  echo "$obj $er $mm"; rm -f "$tmp"
}

lt() { awk -v a="$1" -v b="$2" 'BEGIN{ exit !((a+0) < (b+0)) }'; }
# a < b*(1-MIN_GAIN) : amelioration significative
lt_gain() { awk -v a="$1" -v b="$2" -v g="$MIN_GAIN" 'BEGIN{ exit !((a+0) < (b+0)*(1-g)) }'; }
le_errmax() { awk -v e="$1" -v m="$ERR_MAX" 'BEGIN{ exit !((e+0) <= (m+0)) }'; }

if [ "$REPORT_ONLY" = "0" ]; then
echo "=== Balayage en memoire : N=$N NRHS=$NRHS CORES=$CORES solveurs=[$SOLVERS] objectif=$OBJECTIVE ERR_MAX=$ERR_MAX ==="
echo "    TOL=[$TOL_LIST] NB=[$NB_LIST] threads=[$THREADS_LIST] RATIO=[$RATIO_LIST] repetitions=$REPEATS MIN_GAIN=$MIN_GAIN"
echo "    resultats -> $OUTDIR"
{ echo "#PARAMS N=$N NRHS=$NRHS CORES=$CORES OBJECTIVE=$OBJECTIVE ERR_MAX=$ERR_MAX REPEATS=$REPEATS LAUNCH='$LAUNCH' LAUNCH_THREADS_OPT='$LAUNCH_THREADS_OPT'"; } > "$DECISIONS"

for SOLVER in $SOLVERS; do
  echo; echo "##### $SOLVER"
  P0=$((CORES / BASE_THREADS)); read PR0 PC0 <<< "$(square_grid $P0)"
  NBARR=($NB_LIST); NB0=${NBARR[$(( ${#NBARR[@]} / 2 ))]}

  echo "  [1/5] TOL (grille ${PR0}x${PC0}, $BASE_THREADS threads, NB=$NB0)"
  BEST_TOL=""; SMALLEST_TOL=""
  for tol in $TOL_LIST; do
    read obj er mm <<< "$(run_cfg tol $SOLVER $PR0 $PC0 $BASE_THREADS $NB0 $tol 0.333)"
    SMALLEST_TOL=$tol
    if le_errmax "$er"; then
      if [ -z "$BEST_TOL" ] || awk -v a="$tol" -v b="$BEST_TOL" 'BEGIN{exit !(a+0>b+0)}'; then BEST_TOL=$tol; fi
    fi
  done
  if [ -z "$BEST_TOL" ]; then
    BEST_TOL=$SMALLEST_TOL
    echo "  ATTENTION : aucun TOL ne respecte ERR_MAX=$ERR_MAX -- on garde le dernier essaye ($BEST_TOL)"
  fi
  echo "  -> TOL retenu : $BEST_TOL"

  echo "  [2/5] NB (TOL=$BEST_TOL)"
  BEST_NB=""; BEST_OBJ=$BIG; NBRES=""
  for nb in $NB_LIST; do
    read obj er mm <<< "$(run_cfg nb $SOLVER $PR0 $PC0 $BASE_THREADS $nb $BEST_TOL 0.333)"
    if le_errmax "$er"; then
      NBRES="$NBRES$nb $obj $mm\n"
      if lt "$obj" "$BEST_OBJ"; then BEST_OBJ=$obj; BEST_NB=$nb; fi
    fi
  done
  if [ -n "$NBRES" ]; then
    # parmi les NB a moins de NEAR du meilleur temps : le moins gourmand en memoire
    read BEST_NB BEST_OBJ <<< "$(printf "$NBRES" | awk -v b="$BEST_OBJ" -v t="$NEAR"         '$2 <= b*(1+t) { if (best=="" || $3 < bm) { best=$1; bo=$2; bm=$3 } } END { print best, bo }')"
  fi
  [ -z "$BEST_NB" ] && BEST_NB=$NB0
  echo "  -> NB retenu : $BEST_NB (objectif $BEST_OBJ s ; a moins de $NEAR du meilleur temps, memoire minimale)"

  echo "  [3/5] grille / threads (NB=$BEST_NB, TOL=$BEST_TOL)"
  BEST_PR=$PR0; BEST_PC=$PC0; BEST_TH=$BASE_THREADS; BEST_MM=$BIG
  for th in $THREADS_LIST; do
    P=$((CORES / th))
    while read pr pc; do
      read obj er mm <<< "$(run_cfg grid $SOLVER $pr $pc $th $BEST_NB $BEST_TOL 0.333)"
      # memoire de reference = celle de la config de depart, tant qu'elle reste la meilleure
      if [ "$pr" = "$PR0" ] && [ "$pc" = "$PC0" ] && [ "$th" = "$BASE_THREADS" ] && \
         [ "$BEST_PR" = "$PR0" ] && [ "$BEST_PC" = "$PC0" ] && [ "$BEST_TH" = "$BASE_THREADS" ]; then BEST_MM=$mm; fi
      if le_errmax "$er" && lt_gain "$obj" "$BEST_OBJ"; then BEST_OBJ=$obj; BEST_PR=$pr; BEST_PC=$pc; BEST_TH=$th; BEST_MM=$mm; fi
    done <<< "$(grid_shapes $P)"
  done
  echo "  -> grille retenue : ${BEST_PR}x${BEST_PC}, OMP_NUM_THREADS=$BEST_TH (objectif $BEST_OBJ s)"

  echo "  [4/5] RATIO"
  BEST_RATIO=0.333
  for ratio in $RATIO_LIST; do
    [ "$ratio" = "0.333" ] && continue          # deja mesure en phase 3
    read obj er mm <<< "$(run_cfg ratio $SOLVER $BEST_PR $BEST_PC $BEST_TH $BEST_NB $BEST_TOL $ratio)"
    # RATIO ne change le calcul que s'il fait basculer des blocs dense <-> bas-rang :
    # memoire identique (a 1 % pres) = aucun effet reel, l'ecart de temps est du bruit.
    if le_errmax "$er" && lt_gain "$obj" "$BEST_OBJ" && \
       awk -v a="$mm" -v b="$BEST_MM" 'BEGIN{ d=a-b; if (d<0) d=-d; exit !(d > 0.01*b) }'; then
      BEST_OBJ=$obj; BEST_RATIO=$ratio
    fi
  done
  echo "  -> RATIO retenu : $BEST_RATIO (objectif $BEST_OBJ s ; change seulement si la structure compressee change)"

  # Ligne de reference (meilleure config) dans summary.csv
  REF=$(awk -F, -v s=$SOLVER -v pr=$BEST_PR -v pc=$BEST_PC -v th=$BEST_TH -v nb=$BEST_NB -v r=$BEST_RATIO '
        $2==s && $3==pr && $4==pc && $5==th && $7==nb && ($10+0)==(r+0) && ($16+0)<1e299 {l=$0} END{print l}' "$SUMMARY")
  ARENA_LINE=""
  if [ "$TEST_ARENA" = "1" ]; then
    echo "  [5/5] variante arena (meme reglage)"
    run_cfg arena ${SOLVER}_arena $BEST_PR $BEST_PC $BEST_TH $BEST_NB $BEST_TOL $BEST_RATIO > /dev/null
    ARENA_LINE=$(awk -F, -v s=${SOLVER}_arena '$1=="arena" && $2==s {l=$0} END{print l}' "$SUMMARY")
  fi
  echo "$SOLVER|$BEST_TOL|$BEST_NB|$BEST_PR|$BEST_PC|$BEST_TH|$BEST_RATIO|$REF|$ARENA_LINE" >> "$DECISIONS"
done
fi

#------------------------------------------------------------------ rapport
{
  echo "Meilleurs parametres -- balayage en memoire ($(date))"
  echo "N=$N NRHS=$NRHS CORES=$CORES objectif=$OBJECTIVE ERR_MAX=$ERR_MAX repetitions=$REPEATS (medianes)"
  echo "Matrice de demo complexe symetrique (sym_testmat_mod.f90), compression SVD randomisee."
  echo
  while IFS='|' read s tol nb pr pc th ratio ref arena; do
    IFS=, read -r _ _ _ _ _ _ _ _ _ _ _ ta tf ts tt obj mm mt rss er <<< "$ref"
    echo "== $s =="
    echo "  TOL=$tol  NB=$nb  grille=${pr}x${pc}  OMP_NUM_THREADS=$th  RATIO=$ratio"
    printf "  assemblage %.3f s | factorisation %.3f s | resolution %.3f s | total %.3f s\n" "$ta" "$tf" "$ts" "$tt"
    printf "  facteurs %.1f Mo/proc (%.1f Mo au total) | RSS max %.1f Mo/proc | erreur %s\n" "$mm" "$mt" "$rss" "$er"
    if [ -n "$arena" ]; then
      IFS=, read -r _ _ _ _ _ _ _ _ _ _ _ aa af as at aobj amm amt arss aer <<< "$arena"
      if [ -n "$at" ]; then
        printf "  arena : total %.3f s (x%.2f) | RSS max %.1f Mo/proc (x%.2f) | erreur %s\n" "$at" \
               "$(awk -v a=$at -v b=$tt 'BEGIN{print a/b}')" "$arss" "$(awk -v a=$arss -v b=$rss 'BEGIN{print a/b}')" "$aer"
        if awk -v a=$arss -v b=$rss -v x=$at -v y=$tt 'BEGIN{exit !(a<0.95*b && x<1.10*y)}'; then
          echo "  -> arena CONSEILLEE (RSS reduite d'au moins 5 %, temps degrade de moins de 10 %)"
        else
          echo "  -> arena non conseillee a cette taille (pas de gain de RSS suffisant)"
        fi
      else
        echo "  arena : echec (voir $OUTDIR/logs/arena_*)"
      fi
    fi
    if [ "${s#ldlt}" != "$s" ]; then EX=ldlt/bin/example_blr_mpi_ldlt; else EX=bin/example_blr_mpi_pivot; fi
    echo "  commande conseillee :"
    echo "    OMP_NUM_THREADS=$th $LAUNCH $((pr*pc)) $EX $pr $pc $N $nb $NRHS $tol --ratio=$ratio"
    echo
  done < <(grep -v '^#' "$DECISIONS")
  echo "Comparaison des solveurs (meilleur reglage de chacun) :"
  while IFS='|' read s tol nb pr pc th ratio ref arena; do
    IFS=, read -r _ _ _ _ _ _ _ _ _ _ _ ta tf ts tt obj mm mt rss er <<< "$ref"
    printf "  %-6s objectif %10.3f s | facteurs %9.1f Mo/proc | erreur %s\n" "$s" "$obj" "$mm" "$er"
  done < <(grep -v '^#' "$DECISIONS")
  echo
  echo "Meilleures configurations MESUREES par critere (toutes phases, erreur <= ERR_MAX, medianes) :"
  for s in $(grep -v '^#' "$DECISIONS" | cut -d'|' -f1); do
    for crit in "t_fact:13:factorisation" "t_total:15:assemblage+factorisation+resolution" "mem_max_mb:17:memoire des facteurs par processus"; do
      IFS=: read cname col label <<< "$crit"
      awk -F, -v s="$s" -v c="$col" -v m="$ERR_MAX" -v lab="$label" 'NR>1 && $2==s && $11>0 && ($20+0)<=(m+0) && $16+0<1e299 {
            if (best=="" || ($c+0) < bv) { best=$0; bv=$c+0 } }
          END { if (best!="") { split(best,f,",");
            printf "  %-6s min %-36s : %10.3f | grille %sx%s, %s threads, NB=%s, TOL=%s, RATIO=%s (fact %.3f s, total %.3f s, %.1f Mo/proc, err %s)\n",
                   s, lab, bv, f[3], f[4], f[5], f[7], f[9], f[10], f[13], f[15], f[17], f[20] } }' "$SUMMARY"
    done
  done
  echo
  echo "Rappels : LDLt seulement si la matrice est symetrique (Aᵀ=A) ; refaire ce balayage avec la"
  echo "vraie matrice et sur la machine cible (l'optimum de NB depend de N et du materiel)."
} > "$BEST"
echo; cat "$BEST"
