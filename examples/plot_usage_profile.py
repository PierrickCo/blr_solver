#!/usr/bin/env python3
"""
plot_usage_profile.py -- plot the per-process RSS/CPU profile written by
example_blr_mpi_pivot.f90 when run with PROFILE_USAGE=1.

Reads every usage_profile_rank*.csv in the given directory (columns:
rank,stage,t_wall_s,rss_mb,cpu_pct -- one row per solver stage per MPI
process) and draws two stacked plots against elapsed wall time: memory
(RSS, Mo) and CPU usage (%, can exceed 100% with multiple MKL/OpenMP
threads), one line per process, with the solver's stage boundaries
(from rank 0) marked as vertical dashed lines.

Usage:
    python3 examples/plot_usage_profile.py [directory] [output.png]
    (defaults: directory='.', output='usage_profile.png')
"""
import csv
import glob
import os
import sys


def load_profiles(directory):
    profiles = {}  # rank -> list of (stage, t, rss_mb, cpu_pct)
    for path in sorted(glob.glob(os.path.join(directory, 'usage_profile_rank*.csv'))):
        with open(path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            rows = [(r['stage'], float(r['t_wall_s']), float(r['rss_mb']), float(r['cpu_pct']))
                    for r in reader]
        if not rows:
            continue
        base = os.path.basename(path)
        rank = int(base[len('usage_profile_rank'):-len('.csv')])
        profiles[rank] = rows
    return profiles


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    output = sys.argv[2] if len(sys.argv) > 2 else 'usage_profile.png'

    profiles = load_profiles(directory)
    if not profiles:
        print(f"Aucun usage_profile_rank*.csv trouve dans {directory!r} -- "
              f"relancer l'exemple avec PROFILE_USAGE=1.")
        sys.exit(1)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib n'est pas installe. Installer avec :\n  pip install matplotlib")
        sys.exit(1)

    fig, (ax_mem, ax_cpu) = plt.subplots(2, 1, sharex=True, figsize=(10, 7))

    for rank in sorted(profiles):
        rows = profiles[rank]
        t = [r[1] for r in rows]
        rss = [r[2] for r in rows]
        cpu = [r[3] for r in rows]
        ax_mem.plot(t, rss, marker='o', label=f'rang {rank}')
        # CPU% is an AVERAGE over the interval since the previous sample,
        # not an instantaneous value -- draw it as a step held over that
        # interval rather than interpolating linearly between points.
        ax_cpu.step(t, cpu, where='pre', marker='o', label=f'rang {rank}')

    ax_mem.set_ylabel('RSS (Mo)')
    ax_mem.set_title('Utilisation memoire par processus')
    ax_mem.grid(True, alpha=0.3)
    ax_mem.legend(fontsize=8, ncol=min(len(profiles), 4))

    ax_cpu.set_ylabel('CPU (%)')
    ax_cpu.set_xlabel('Temps ecoule (s)')
    ax_cpu.set_title('Utilisation CPU par processus (100% = 1 coeur pleinement utilise)')
    ax_cpu.grid(True, alpha=0.3)

    # Stage BOUNDARIES only (e.g. 'assemblage', not the fine-grained
    # 'assemblage@45%' samples PROFILE_USAGE also records inside that
    # stage) -- those carry the "@" separator and exist purely to give
    # the curves themselves more resolution; marking each with its own
    # vertical line would be ~20 lines per stage, illegible clutter.
    # Boundaries are taken from rank 0 (or the lowest rank present) as a
    # reference timeline -- other ranks may be missing one (e.g.
    # CHECK_ERROR is rank-0-only) but timing is roughly synchronized by
    # the solver's own MPI collectives. Fast stages (e.g. verification
    # passes on a small problem) can still land within a fraction of a
    # percent of each other -- individually labeling each would overlap
    # into an illegible smear, so nearby boundaries are merged into one
    # combined label.
    ref_rank = min(profiles)
    ref_rows = [r for r in profiles[ref_rank] if '@' not in r[0]]
    xlim = ax_mem.get_xlim()
    gap_threshold = 0.03 * (xlim[1] - xlim[0] or 1.0)
    clusters = []  # list of [t, [stage, ...]]
    for stage, t, _, _ in ref_rows:
        if clusters and t - clusters[-1][0] < gap_threshold:
            clusters[-1][1].append(stage)
        else:
            clusters.append([t, [stage]])
    for t, stages in clusters:
        for ax in (ax_mem, ax_cpu):
            ax.axvline(t, color='grey', linestyle='--', linewidth=0.7, alpha=0.6)
        ax_mem.text(t, ax_mem.get_ylim()[1], ' / '.join(stages), rotation=90,
                    fontsize=7, ha='right', va='top', color='grey')

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    print(f"Ecrit {output} ({len(profiles)} processus, "
          f"{len(profiles[ref_rank])} echantillons, {len(clusters)} etapes reperees).")


if __name__ == '__main__':
    main()
