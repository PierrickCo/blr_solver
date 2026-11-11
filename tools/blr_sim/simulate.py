#!/usr/bin/env python3
"""Simulateur de repartition de la factorisation LDLt BLR+MPI.

Rejoue, sans rien calculer, la charge de chaque processus pour une grille
PR x PC et une regle de repartition des blocs quelconques, a partir du
travail estime de chaque mise a jour de Schur releve pendant une vraie
factorisation (variable d'environnement BLR_DUMP_WORK=1 : fichiers
blr_work_rank<id>.bin, un enregistrement (K, I, J, travail) par mise a jour
du bloc (I, J) a l'etape K ; I = J : bloc diagonal). Ce releve ne depend
pas de la grille du calcul d'origine.

Indicateurs (memes definitions que la trace BLR_TRACE_FACT=1) :
  desequilibre global      max des totaux par processus / moyenne - 1
  synchronisation          somme des max par etape / max des totaux - 1
  borne d'acceleration     P x moyenne / somme des max par etape
  diffusions par processus appels MPI_BCAST des panneaux (2 par bloc
                           aujourd'hui ; "groupees" : 2 par panneau et par
                           communicateur)

Repartitions :
  actuel    lignes cycliques, colonnes en cyclique decale (DIAG_SPREAD)
  cyclique  lignes et colonnes cycliques
  pondere   blocs-lignes et blocs-colonnes affectes aux lignes / colonnes
            de grille par poids decroissant (au moins charge), puis
            ameliorations locales sur la charge maximale

Usage :
  python3 simulate.py DOSSIER [--grilles 2x2,4x4,8x8] [--repartitions actuel,pondere]
"""
import argparse
import glob
import math
import os

import numpy as np

DT = np.dtype([("k", "<i4"), ("i", "<i4"), ("j", "<i4"), ("w", "<f8")])


def load(d):
    files = sorted(glob.glob(os.path.join(d, "blr_work_rank*.bin")))
    if not files:
        raise SystemExit(f"aucun fichier blr_work_rank*.bin dans {d}")
    a = np.concatenate([np.fromfile(f, dtype=DT) for f in files])
    return a, int(max(a["i"].max(), a["j"].max()))


def col_owner_current(nblk, pr, pc):
    g = math.gcd(pr, pc)
    per = (pr // g) * pc if g > 1 else 0
    j0 = np.arange(nblk)
    return (j0 + (j0 // per if per else 0)) % pc


def loads(wtot, rowmap, colmap, pr, pc):
    m = np.zeros((pr, pc))
    np.add.at(m, (rowmap[:, None].repeat(wtot.shape[1], 1), colmap[None, :].repeat(wtot.shape[0], 0)), wtot)
    return m


def weighted(wtot, pr, pc, sweeps=3):
    """Affectation gloutonne par poids decroissant, puis deplacements qui
    reduisent la charge maximale d'un processus."""
    nblk = wtot.shape[0]

    def greedy(weights, nbin):
        m = np.zeros(nblk, dtype=int)
        load = np.zeros(nbin)
        for b in np.argsort(-weights):
            p = int(np.argmin(load))
            m[b] = p
            load[p] += weights[b]
        return m

    rowmap = greedy(wtot.sum(1), pr)
    colmap = greedy(wtot.sum(0), pc)

    def improve(w, mp, other, n, nother):
        """Deplace chaque bloc (ligne de w) vers la ligne de grille qui
        minimise la charge maximale ; w[b, :] agrege par 'other'."""
        contrib = np.zeros((nblk, nother))
        np.add.at(contrib, (np.arange(nblk)[:, None].repeat(nblk, 1), other[None, :].repeat(nblk, 0)), w)
        base = np.zeros((n, nother))
        np.add.at(base, mp, contrib)
        for b in np.argsort(-w.sum(1)):
            base[mp[b]] -= contrib[b]
            rmax = base.max(1)
            order = np.argsort(-rmax)
            cand = (base + contrib[b]).max(1)                 # charge max de la ligne p si b y va
            others = np.where(np.arange(n) == order[0], rmax[order[1]] if n > 1 else 0.0, rmax[order[0]])
            p = int(np.argmin(np.maximum(cand, others)))
            mp[b] = p
            base[p] += contrib[b]

    for _ in range(sweeps):
        improve(wtot, rowmap, colmap, pr, pc)
        improve(wtot.T, colmap, rowmap, pc, pr)
    return rowmap, colmap


def metrics(a, nblk, rowmap, colmap, pr, pc):
    p = pr * pc
    proc = rowmap[a["i"] - 1] * pc + colmap[a["j"] - 1]
    w = np.zeros((p, nblk))
    np.add.at(w, (proc, a["k"] - 1), a["w"])
    tot = w.sum(1)
    wsync = w.max(0).sum()
    # diffusions : etape k, blocs I > k (lignes) et J > k (colonnes) ; 2 appels par bloc
    nrow = np.zeros((pr, nblk))
    ncol = np.zeros((pc, nblk))
    np.add.at(nrow, (rowmap, np.arange(nblk)), 1)
    np.add.at(ncol, (colmap, np.arange(nblk)), 1)
    # nombre de blocs d'indice > k dans chaque ligne / colonne de grille, pour k = 1..nblk
    rem_r = nrow[:, ::-1].cumsum(1)[:, ::-1] - nrow
    rem_c = ncol[:, ::-1].cumsum(1)[:, ::-1] - ncol
    calls = 2 * (rem_r.sum(1)[:, None] + rem_c.sum(1)[None, :])           # (pr, pc)
    packed = 2 * ((rem_r > 0).sum(1)[:, None] + (rem_c > 0).sum(1)[None, :])
    return dict(imb=tot.max() / tot.mean() - 1, sync=wsync / tot.max() - 1, bound=p * tot.mean() / wsync,
                calls=calls.mean(), packed=packed.mean(), wpercall=tot.mean() / max(calls.mean(), 1),
                mean=tot.mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dossier")
    ap.add_argument("--grilles", default="2x2,4x4,6x6,8x8")
    ap.add_argument("--repartitions", default="actuel,pondere")
    args = ap.parse_args()
    a, nblk = load(args.dossier)
    wtot = np.zeros((nblk, nblk))
    np.add.at(wtot, (a["i"] - 1, a["j"] - 1), a["w"])
    print(f"{len(a)} mises a jour, {nblk} blocs, travail total {a['w'].sum():.3e}")
    print(f"{'grille':>8} {'repartition':>11} {'deseq.':>7} {'synchro':>8} {'borne':>14} {'diffusions':>11} {'groupees':>9} "
          f"{'travail/diff.':>13}")
    for g in args.grilles.split(","):
        pr, pc = (int(x) for x in g.lower().split("x"))
        for rep in args.repartitions.split(","):
            if rep == "actuel":
                rowmap, colmap = np.arange(nblk) % pr, col_owner_current(nblk, pr, pc)
            elif rep == "cyclique":
                rowmap, colmap = np.arange(nblk) % pr, np.arange(nblk) % pc
            elif rep == "pondere":
                rowmap, colmap = weighted(wtot, pr, pc)
            else:
                raise SystemExit(f"repartition inconnue : {rep}")
            m = metrics(a, nblk, rowmap, colmap, pr, pc)
            print(f"{g:>8} {rep:>11} {100 * m['imb']:6.1f}% {100 * m['sync']:7.1f}% "
                  f"{m['bound']:7.1f}/{pr * pc:<6d} {m['calls']:11.0f} {m['packed']:9.0f} {m['wpercall']:13.3e}")


if __name__ == "__main__":
    main()
