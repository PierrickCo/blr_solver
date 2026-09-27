#!/usr/bin/env python3
"""
plot_blr_structure.py -- visualize how a BLR matrix's blocks are
compressed, from the CSVs written by example_blr_mpi_pivot_user.f90 (via
profiling_mod.f90) when run with --dump-structure=1.

Reads every blr_structure_<stage>_rank*.csv in the given directory
(columns: block_row,block_col,kind,rank -- one row per LOCAL block owned
by that process, kind is diagonal/dense/lowrank) for the two stages the
example writes ("assembly" and "factorization"), and draws a one-pixel-
per-block heatmap for each: diagonal blocks and dense (uncompressed)
off-diagonal blocks get fixed colors, low-rank blocks are colored by
their rank on a shared colormap -- so a side-by-side comparison shows
directly where the factorization's fill-in turned a compressed block
dense, or changed its rank. This is exactly the structure produced by
YOUR matrix (USER_GET_BLOCK in user_matrix_mod.f90), not a demo matrix.

Usage:
    python3 tools/plot_blr_structure.py [directory] [output.png]
    (defaults: directory='.', output='blr_structure.png')
"""
import csv
import glob
import os
import sys

STAGES = [('assembly', 'Fin assemblage'), ('factorization', 'Fin factorisation')]
COLOR_DIAGONAL = (0.15, 0.15, 0.15, 1.0)
COLOR_DENSE = (0.85, 0.15, 0.15, 1.0)


def load_stage(directory, stage):
    """Returns (nblk, grid) where grid[(row-1, col-1)] = (kind, rank),
    0-indexed for direct use as an image array -- or (0, {}) if no file
    for this stage was found."""
    grid = {}
    nblk = 0
    for path in glob.glob(os.path.join(directory, f'blr_structure_{stage}_rank*.csv')):
        with open(path, 'r', newline='') as f:
            for row in csv.DictReader(f):
                r, c = int(row['block_row']), int(row['block_col'])
                grid[(r - 1, c - 1)] = (row['kind'], int(row['rank']))
                nblk = max(nblk, r, c)
    return nblk, grid


def build_image(nblk, grid, rank_norm, cmap):
    import numpy as np
    img = np.ones((nblk, nblk, 4))
    img[:, :, 3] = 0.0  # blocks never assembled locally (shouldn't happen) stay transparent
    for (i, j), (kind, rank) in grid.items():
        if kind == 'diagonal':
            img[i, j] = COLOR_DIAGONAL
        elif kind == 'dense':
            img[i, j] = COLOR_DENSE
        else:
            img[i, j] = cmap(rank_norm(rank))
    return img


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    output = sys.argv[2] if len(sys.argv) > 2 else 'blr_structure.png'

    stage_data = [(label, load_stage(directory, stage)) for stage, label in STAGES]
    stage_data = [(label, nblk, grid) for label, (nblk, grid) in stage_data if nblk > 0]
    if not stage_data:
        print(f"Aucun blr_structure_<etape>_rank*.csv trouve dans {directory!r} -- "
              f"relancer l'exemple avec --dump-structure=1.")
        sys.exit(1)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib n'est pas installe. Installer avec :\n  pip install matplotlib")
        sys.exit(1)

    all_lowrank_ranks = [rank for _, _, grid in stage_data for kind, rank in grid.values()
                         if kind == 'lowrank']
    max_rank = max(all_lowrank_ranks, default=1)
    cmap = plt.get_cmap('viridis_r')
    rank_norm = mcolors.Normalize(vmin=0, vmax=max_rank)

    fig, axes = plt.subplots(1, len(stage_data), figsize=(6 * len(stage_data) + 1, 5.5), squeeze=False)
    axes = axes[0]

    for ax, (label, nblk, grid) in zip(axes, stage_data):
        img = build_image(nblk, grid, rank_norm, cmap)
        ax.imshow(img, interpolation='nearest')
        ax.set_title(f'{label} ({nblk}x{nblk} blocs)')
        ax.set_xlabel('Bloc colonne')
        ax.set_ylabel('Bloc ligne')

    handles = [
        mpatches.Patch(color=COLOR_DIAGONAL, label='diagonal (toujours dense)'),
        mpatches.Patch(color=COLOR_DENSE, label='dense (repli, non compresse)'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))

    sm = plt.cm.ScalarMappable(norm=rank_norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.85, pad=0.02)
    cbar.set_label('Rang (blocs bas-rang)')

    fig.suptitle('Structure de compression de la matrice BLR', fontsize=13)
    fig.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Ecrit {output} ({len(stage_data)} etape(s): " +
          ', '.join(f'{label} {nblk}x{nblk}' for label, nblk, _ in stage_data) + ").")


if __name__ == '__main__':
    main()
