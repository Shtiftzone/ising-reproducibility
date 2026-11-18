#!/usr/bin/env python3
"""Cross-temperature correlations of paired Ising time series.

Input: ROOT/size_L/conf-T-RUN.dat, with ten numeric columns:
  1 M, 2 E, 3 Chi_plus, 4 Chi_minus, then six cell counts.
For each row, EC=(Chi_plus-Chi_minus)/2; D=M+c*EC,
where c=2 (square) or c=1 (triangular).

Rows are compared by their recorded position in the original files.
The diagnostic tests for synchronized fluctuations between simulation
streams; failure to detect correlation does not prove independence.
ASCII source only.
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

PATTERN = re.compile(r'^conf-([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)-(\d+)\.dat$')
DEFAULT_SIZES = [64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072]
ALL_OBS = ('energy', 'm', 'abs_m', 'ec', 'abs_ec', 'd', 'abs_d')


def load_all(size_dir, lattice, expected_n):
    found = []
    for path in size_dir.glob('conf-*.dat'):
        match = PATTERN.fullmatch(path.name)
        if match:
            found.append((float(match.group(1)), int(match.group(2)), path))
    found.sort(key=lambda x: (x[0], x[1]))
    if not found:
        raise ValueError('No conf-T-RUN.dat files found')
    seen = set()
    temps, runs, datasets = [], [], []
    n = None
    for temp, run, path in found:
        if temp in seen:
            raise ValueError('Repeated temperature %.8f; multiple runs need explicit handling' % temp)
        seen.add(temp)
        raw = np.loadtxt(path, ndmin=2)
        if raw.ndim != 2 or raw.shape[1] != 10 or raw.shape[0] < 4:
            raise ValueError('Expected 10 columns and >=4 rows in %s; got %s' % (path, raw.shape))
        if not np.isfinite(raw).all():
            raise ValueError('Nonfinite values in ' + str(path))
        if n is None:
            n = raw.shape[0]
        if raw.shape[0] != n:
            raise ValueError('Unequal time series lengths in ' + str(size_dir))
        if expected_n and n != expected_n:
            raise ValueError('Expected %d measurements; got %d in %s' % (expected_n, n, path))
        m = raw[:, 0]
        energy = raw[:, 1]
        ec = 0.5 * (raw[:, 2] - raw[:, 3])
        d = m + (2.0 if lattice == 'square' else 1.0) * ec
        datasets.append({'energy':energy, 'm':m, 'abs_m':abs(m),
                         'ec':ec, 'abs_ec':abs(ec), 'd':d, 'abs_d':abs(d)})
        temps.append(temp)
        runs.append(run)
    return np.asarray(temps), runs, datasets


def get_standardized(datasets, observable):
    data = np.column_stack([ds[observable] for ds in datasets])
    centered = data - data.mean(axis=0)
    scale = np.linalg.norm(centered, axis=0)
    if np.any(scale <= 0):
        raise ValueError('Constant series for observable ' + observable)
    return centered / scale


def null_maxima(z, repeats, rng):
    n, k = z.shape
    upper = np.triu_indices(k, k=1)
    adj = (np.arange(k-1), np.arange(1, k))
    null_all = np.empty(repeats)
    null_adj = np.empty(repeats)
    shifted = np.empty_like(z)
    base = np.arange(n)
    for b in range(repeats):
        offsets = rng.integers(0, n, size=k)
        for j in range(k):
            shifted[:, j] = z[(base + offsets[j]) % n, j]
        cor = shifted.T @ shifted
        null_all[b] = np.max(np.abs(cor[upper]))
        null_adj[b] = np.max(np.abs(cor[adj]))
    return null_all, null_adj


def write_matrix(path, temps, cor):
    with path.open('w', newline='', encoding='ascii') as f:
        w = csv.writer(f)
        w.writerow(['temperature'] + ['%.6f' % t for t in temps])
        for t, row in zip(temps, cor):
            w.writerow(['%.6f' % t] + ['%.8f' % x for x in row])


def write_pairs(path, temps, runs, cor, cutoff_all, cutoff_adj):
    with path.open('w', newline='', encoding='ascii') as f:
        w = csv.writer(f)
        w.writerow(['T_i','T_j','run_i','run_j','pearson_r','abs_r','adjacent',
                    'exceeds_all_pairs_null_95','exceeds_adjacent_null_95'])
        for i in range(len(temps)):
            for j in range(i+1, len(temps)):
                r = float(cor[i,j]); near = (j == i+1)
                w.writerow(['%.6f' % temps[i], '%.6f' % temps[j], runs[i], runs[j],
                            '%.8f' % r, '%.8f' % abs(r), int(near),
                            int(abs(r) > cutoff_all), int(near and abs(r) > cutoff_adj)])


def plot_matrix(path, temps, cor, size, obs):
    fig, ax = plt.subplots(figsize=(8,7))
    im = ax.imshow(cor, cmap='coolwarm', vmin=-1, vmax=1, origin='lower',
                   interpolation='nearest')
    ids = np.unique(np.linspace(0, len(temps)-1, min(len(temps), 9), dtype=int))
    ax.set_xticks(ids, ['%.4f' % temps[i] for i in ids], rotation=45, ha='right')
    ax.set_yticks(ids, ['%.4f' % temps[i] for i in ids])
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Temperature')
    ax.set_title('L=%d: cross-temperature %s correlation' % (size, obs))
    fig.colorbar(im, ax=ax, label='Pearson r')
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--lattice', choices=('square','triangular'), required=True)
    p.add_argument('--sizes', nargs='+', type=int, default=DEFAULT_SIZES)
    p.add_argument('--observables', nargs='+', choices=ALL_OBS,
                   default=list(ALL_OBS))
    p.add_argument('--expected-n', type=int, default=200,
                   help='Set 0 to allow any consistent number of measurements')
    p.add_argument('--null-repeats', type=int, default=500)
    p.add_argument('--seed', type=int, default=20261104)
    args = p.parse_args()
    if args.null_repeats < 20:
        p.error('--null-repeats must be >= 20')
    rng = np.random.default_rng(args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for size in args.sizes:
        size_dir = args.root / ('size_%d' % size)
        if not size_dir.is_dir():
            print('L=%d: directory absent; skipping' % size)
            continue
        try:
            temps, runs, datasets = load_all(size_dir, args.lattice, args.expected_n)
        except (OSError, ValueError) as exc:
            print('L=%d: ERROR: %s' % (size, exc))
            continue
        if len(temps) < 2:
            print('L=%d: fewer than 2 temperatures; skipping' % size)
            continue
        out = args.output / ('size_%d' % size)
        out.mkdir(parents=True, exist_ok=True)
        for obs in args.observables:
            try:
                z = get_standardized(datasets, obs)
            except ValueError as exc:
                print('L=%d, %s: skipped: %s' % (size, obs, exc))
                continue
            cor = np.clip(z.T @ z, -1.0, 1.0)
            np.fill_diagonal(cor, 1.0)
            upper = np.abs(cor[np.triu_indices(len(temps), k=1)])
            adjacent = np.abs(np.diag(cor, k=1))
            null_all, null_adj = null_maxima(z, args.null_repeats, rng)
            max_all = float(upper.max())
            max_adj = float(adjacent.max())
            cut_all = float(np.quantile(null_all, 0.95))
            cut_adj = float(np.quantile(null_adj, 0.95))
            p_all = (1 + np.count_nonzero(null_all >= max_all)) / (args.null_repeats+1)
            p_adj = (1 + np.count_nonzero(null_adj >= max_adj)) / (args.null_repeats+1)
            write_matrix(out / ('correlations_%s.csv' % obs), temps, cor)
            write_pairs(out / ('pairs_%s.csv' % obs), temps, runs, cor, cut_all, cut_adj)
            plot_matrix(out / ('heatmap_%s.png' % obs), temps, cor, size, obs)
            rows.append([args.lattice,size,obs,len(temps),z.shape[0],max_all,cut_all,p_all,
                         max_adj,cut_adj,p_adj])
            print('L=%d %s: max |r|=%.3f (p=%.4f); adjacent max |r|=%.3f (p=%.4f)' %
                  (size,obs,max_all,p_all,max_adj,p_adj))
    if not rows:
        raise SystemExit('No results; verify root, filenames, and input format')
    with (args.output / 'summary.csv').open('w',newline='',encoding='ascii') as f:
        w = csv.writer(f)
        w.writerow(['lattice','L','observable','n_temperatures','n_measurements',
                    'max_abs_r','all_pairs_null_95','all_pairs_global_p',
                    'adjacent_max_abs_r','adjacent_null_95','adjacent_global_p'])
        w.writerows(rows)
    print('Output: ' + str(args.output))
    print('CAUTION: p-values compare synchronized ordering to independently shifted series;')
    print('they are exploratory diagnostics, not proof of independent simulation streams.')


if __name__ == '__main__':
    main()
