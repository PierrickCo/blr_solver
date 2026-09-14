#!/usr/bin/env python3
"""
analyze_block_ranks.py -- DIAGNOSTIC SEUL (aucune permutation implementee) :
quantifie la variabilite du rang des blocs hors-diagonale d'une matrice BLR,
a partir des CSV deja produits par example_blr_mpi_pivot(--dump-structure)
(blr_structure_assembly_rank*.csv), pour juger si une reorganisation des
blocs (permutation guidee par le rang) a une chance d'apporter un gain
avant d'investir dans son implementation.

Usage:
    python3 examples/analyze_block_ranks.py [directory] [stage]
    (defauts : directory='.', stage='assembly')
"""
import csv
import glob
import os
import sys
import statistics


def load_stage(directory, stage):
    grid = {}
    nblk = 0
    for path in glob.glob(os.path.join(directory, f'blr_structure_{stage}_rank*.csv')):
        with open(path, 'r', newline='') as f:
            for row in csv.DictReader(f):
                r, c = int(row['block_row']), int(row['block_col'])
                grid[(r, c)] = (row['kind'], int(row['rank']))
                nblk = max(nblk, r, c)
    return nblk, grid


def effective_rank(kind, rank, nb_diag):
    # Pour un bloc DENSE, le "rang effectif" au sens stockage est min(M,N)
    # (repli plein, meme convention que RANK_SUM dans l'exemple Fortran).
    if kind == 'dense':
        return nb_diag
    return rank


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    stage = sys.argv[2] if len(sys.argv) > 2 else 'assembly'

    nblk, grid = load_stage(directory, stage)
    if nblk == 0:
        print(f"Aucun CSV blr_structure_{stage}_rank*.csv trouve dans {directory!r}.")
        print("Lancer d'abord : mpirun -np <P> ./bin/example_blr_mpi_pivot PR PC N NB NRHS TOL --dump-structure")
        return 1

    # Rang diagonal de reference (taille de bloc dense typique) pour
    # convertir un repli "dense" en rang effectif -- prend le NB du
    # premier bloc diagonal trouve (tous identiques sauf le dernier,
    # potentiellement plus petit si N non multiple de NB).
    nb_diag = max(rank for (r, c), (kind, rank) in grid.items() if kind == 'diagonal' and r == c)

    offdiag = [(r, c, kind, rank) for (r, c), (kind, rank) in grid.items() if r != c]
    ranks = [effective_rank(k, rk, nb_diag) for (_, _, k, rk) in offdiag]
    n_dense = sum(1 for (_, _, k, _) in offdiag if k == 'dense')
    n_lowrank = len(offdiag) - n_dense

    print(f"=== Diagnostic rang des blocs -- stage={stage!r}, NBLK={nblk}, NB(diag)={nb_diag} ===")
    print(f"Blocs hors-diagonale analyses : {len(offdiag)} ({n_lowrank} bas-rang, {n_dense} repli dense)")
    print()
    print("--- Statistiques globales de rang (dense compte comme rang=NB) ---")
    print(f"  min={min(ranks)}  max={max(ranks)}  moyenne={statistics.mean(ranks):.2f}  "
          f"mediane={statistics.median(ranks)}  ecart-type={statistics.pstdev(ranks):.2f}")

    # Histogramme simple (10 classes entre min et max)
    lo, hi = min(ranks), max(ranks)
    if hi > lo:
        nbins = 10
        width = (hi - lo) / nbins
        bins = [0] * nbins
        for rk in ranks:
            idx = min(int((rk - lo) / width), nbins - 1)
            bins[idx] += 1
        print()
        print("--- Histogramme des rangs ---")
        for i, cnt in enumerate(bins):
            b_lo = lo + i * width
            b_hi = lo + (i + 1) * width
            bar = '#' * int(60 * cnt / max(bins))
            print(f"  [{b_lo:7.1f}, {b_hi:7.1f}) {cnt:5d} {bar}")

    # Rang moyen par distance |I-J| au bloc diagonal -- revele si le rang
    # depend de la distance (structure exploitable par une permutation
    # basee sur la distance/le clustering) ou s'il est homogene (rien a
    # gagner d'une reorganisation).
    by_dist = {}
    for (r, c, kind, rank) in offdiag:
        d = abs(r - c)
        by_dist.setdefault(d, []).append(effective_rank(kind, rank, nb_diag))
    print()
    print("--- Rang moyen par distance |I-J| au bloc diagonal (echantillon) ---")
    dists = sorted(by_dist.keys())
    step = max(1, len(dists) // 20)
    for d in dists[::step]:
        vals = by_dist[d]
        print(f"  |I-J|={d:3d} : moyenne={statistics.mean(vals):6.2f}  n={len(vals):4d}")

    # Rang moyen PAR bloc-ligne / bloc-colonne -- si certaines lignes/
    # colonnes ont un rang systematiquement plus eleve que d'autres,
    # c'est le signal exploitable par une permutation "trier par rang
    # croissant" (regrouper les blocs faciles/durs).
    row_ranks = {}
    for (r, c, kind, rank) in offdiag:
        row_ranks.setdefault(r, []).append(effective_rank(kind, rank, nb_diag))
    row_means = {r: statistics.mean(v) for r, v in row_ranks.items()}
    sorted_rows = sorted(row_means.items(), key=lambda kv: kv[1])
    print()
    print("--- Rang moyen par bloc-ligne (min -> max, echantillon) ---")
    n_show = min(10, len(sorted_rows))
    for r, m in sorted_rows[:n_show]:
        print(f"  bloc-ligne {r:4d} : rang moyen = {m:6.2f}  (le plus 'facile')")
    print("  ...")
    for r, m in sorted_rows[-n_show:]:
        print(f"  bloc-ligne {r:4d} : rang moyen = {m:6.2f}  (le plus 'dur')")

    row_vals = list(row_means.values())
    print()
    print("--- Dispersion du rang moyen ENTRE blocs-lignes (le signal qu'une permutation exploiterait) ---")
    print(f"  min={min(row_vals):.2f}  max={max(row_vals):.2f}  "
          f"ecart-type/moyenne (coef. variation) = {statistics.pstdev(row_vals)/statistics.mean(row_vals):.3f}")
    print("  (coefficient de variation proche de 0 => rangs quasi homogenes entre blocs-lignes,")
    print("   peu de marge pour une permutation ; proche de/au-dessus de 1 => forte heterogeneite,")
    print("   une permutation qui regroupe les blocs-lignes similaires a des chances d'aider.)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
