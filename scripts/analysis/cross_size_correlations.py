#!/usr/bin/env python3
"""Cross-size correlations of paired Ising Monte Carlo time series.

Input layout:
  ROOT/size_L/conf-T-RUN.dat

Each input file has ten numeric columns:
  1 M, 2 E, 3 Chi_plus, 4 Chi_minus, followed by six cell counts.
For each recorded configuration:
  EC = (Chi_plus - Chi_minus)/2
  D  = M + c*EC, c=2 for square and c=1 for triangular.

At each fixed temperature, this script compares the recorded time series
across different lattice sizes using the original row order. Pearson
correlations probe synchronized fluctuations, not correlations between
thermodynamic means or the full distribution of the observables.

Independent circular shifts of each size series provide exploratory null
calibration that preserves each series' temporal autocorrelation. Maxima
are calibrated over size pairs at fixed temperature and, separately,
across all temperatures (and optionally all observables). A lack of
significance does not establish independence, and shared-seed effects
may not be detectable as zero-lag linear correlation.

ASCII source only; no third-party dependencies beyond numpy/matplotlib.
"""

import argparse
import csv
import re
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

PATTERN = re.compile(
    r'^conf-([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)-(\d+)\.dat$'
)
DEFAULT_SIZES = [64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072]
ALL_OBS = ('energy', 'm', 'abs_m', 'ec', 'abs_ec', 'd', 'abs_d')


def index_files(root, sizes, run_filter, temp_digits):
    """Build index of temperature -> size -> (run, file)."""
    quant = Decimal(1).scaleb(-temp_digits)
    index = defaultdict(dict)
    for size in sizes:
        size_dir = root / ('size_%d' % size)
        if not size_dir.is_dir():
            print('L=%d: directory absent; skipping' % size)
            continue
        count = 0
        for path in size_dir.glob('conf-*.dat'):
            match = PATTERN.fullmatch(path.name)
            if not match:
                continue
            run = int(match.group(2))
            if run_filter is not None and run != run_filter:
                continue
            try:
                temp = Decimal(match.group(1)).quantize(quant)
            except InvalidOperation as exc:
                raise ValueError('Cannot parse temperature in %s: %s' % (path, exc))
            if size in index[temp]:
                prior = index[temp][size][1]
                raise ValueError(
                    'Multiple files for L=%d, T=%s: %s and %s. '
                    'Use --run to choose a single run.' % (size, temp, prior, path)
                )
            index[temp][size] = (run, path)
            count += 1
        print('L=%d: indexed %d files' % (size, count))
    return index


def read_observables(path, lattice, expected_n):
    raw = np.loadtxt(path, ndmin=2)
    if raw.ndim != 2 or raw.shape[1] != 10 or raw.shape[0] < 4:
        raise ValueError('Expected ten columns and at least four rows; got %s' % (raw.shape,))
    if not np.isfinite(raw).all():
        raise ValueError('Nonfinite numbers in file')
    if expected_n and raw.shape[0] != expected_n:
        raise ValueError('Expected %d measurements, got %d' % (expected_n, raw.shape[0]))
    m = raw[:, 0]
    energy = raw[:, 1]
    ec = 0.5 * (raw[:, 2] - raw[:, 3])
    d = m + (2.0 if lattice == 'square' else 1.0) * ec
    return {
        'energy': energy, 'm': m, 'abs_m': np.abs(m),
        'ec': ec, 'abs_ec': np.abs(ec), 'd': d, 'abs_d': np.abs(d)
    }


def load_temperature(entry, lattice, expected_n):
    """Read available L at this T; exclude invalid or unequal-length files."""
    accepted = []
    n = None
    for size in sorted(entry):
        run, path = entry[size]
        try:
            data = read_observables(path, lattice, expected_n)
        except (OSError, ValueError) as exc:
            print('WARNING: L=%d, %s: %s; skipping' % (size, path.name, exc))
            continue
        rows = len(data['m'])
        if n is None:
            n = rows
        if rows != n:
            print('WARNING: L=%d, %s: unequal series length %d vs %d; skipping' %
                  (size, path.name, rows, n))
            continue
        accepted.append((size, run, data))
    return accepted, n


def standardized_columns(accepted, obs):
    """Return standardized data and indices of nonconstant size series."""
    if not accepted:
        return None, []
    data = np.column_stack([entry[2][obs] for entry in accepted])
    centered = data - data.mean(axis=0)
    scale = np.linalg.norm(centered, axis=0)
    valid = np.flatnonzero(np.isfinite(scale) & (scale > 0))
    if len(valid) < 2:
        return None, valid
    return centered[:, valid] / scale[valid], valid


def null_maxima(z, offsets):
    """Max |r| under independent circular shifts, evaluated by FFT.

    offsets has shape (n_repeats, n_series). For columns i, j, the
    shifted cross-correlation uses lag offsets[:, j] - offsets[:, i].
    This agrees with direct circular shifting and dot products.
    """
    n, k = z.shape
    ii, jj = np.triu_indices(k, k=1)
    spectra = np.fft.rfft(z, axis=0)
    lag_cor = np.fft.irfft(
        np.conj(spectra[:, ii]) * spectra[:, jj], n=n, axis=0
    )
    lags = (offsets[:, jj] - offsets[:, ii]) % n
    vals = np.abs(lag_cor[lags, np.arange(len(ii))[None, :]])
    null_all = vals.max(axis=1)
    null_adj = vals[:, jj == ii+1].max(axis=1)
    return null_all, null_adj


def empirical_p(null, value):
    return (1 + np.count_nonzero(null >= value)) / (len(null) + 1)


def write_matrix(path, sizes, cor):
    with path.open('w', newline='', encoding='ascii') as f:
        w = csv.writer(f)
        w.writerow(['L'] + list(sizes))
        for size, row in zip(sizes, cor):
            w.writerow([size] + ['%.8f' % r for r in row])


def write_pairs(path, sizes, runs, cor, cutoff_all, cutoff_adj):
    with path.open('w', newline='', encoding='ascii') as f:
        w = csv.writer(f)
        w.writerow([
            'L_i', 'L_j', 'run_i', 'run_j', 'pearson_r', 'abs_r',
            'adjacent_available_sizes', 'exceeds_all_pairs_null_95',
            'exceeds_adjacent_null_95'
        ])
        for i in range(len(sizes)):
            for j in range(i+1, len(sizes)):
                r = float(cor[i, j])
                adjacent = (j == i+1)
                w.writerow([
                    sizes[i], sizes[j], runs[i], runs[j],
                    '%.8f' % r, '%.8f' % abs(r), int(adjacent),
                    int(abs(r) > cutoff_all),
                    int(adjacent and abs(r) > cutoff_adj)
                ])


def plot_matrix(path, sizes, cor, temp, obs):
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(
        cor, cmap='coolwarm', vmin=-1.0, vmax=1.0,
        origin='lower', interpolation='nearest'
    )
    ids = np.arange(len(sizes))
    ax.set_xticks(ids, [str(sizes[i]) for i in ids], rotation=45, ha='right')
    ax.set_yticks(ids, [str(sizes[i]) for i in ids])
    ax.set_xlabel('Lattice size L')
    ax.set_ylabel('Lattice size L')
    ax.set_title('T=%s: cross-size %s correlation' % (temp, obs))
    fig.colorbar(im, ax=ax, label='Pearson r')
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_csv(path, header, rows):
    with path.open('w', newline='', encoding='ascii') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def make_aggregate(repeats):
    return {
        'max_all': -1.0, 'max_adj': -1.0,
        'null_all': np.zeros(repeats), 'null_adj': np.zeros(repeats),
        'n_tests': 0, 'n_size_pairs': 0
    }


def update_aggregate(aggregate, obs_all, obs_adj, null_all, null_adj, count_pairs):
    aggregate['max_all'] = max(aggregate['max_all'], obs_all)
    aggregate['max_adj'] = max(aggregate['max_adj'], obs_adj)
    np.maximum(aggregate['null_all'], null_all, out=aggregate['null_all'])
    np.maximum(aggregate['null_adj'], null_adj, out=aggregate['null_adj'])
    aggregate['n_tests'] += 1
    aggregate['n_size_pairs'] += count_pairs


def summary_row(label, aggregate):
    return [
        label, aggregate['n_tests'], aggregate['n_size_pairs'],
        '%.8f' % aggregate['max_all'],
        '%.8f' % np.quantile(aggregate['null_all'], 0.95),
        '%.8f' % empirical_p(aggregate['null_all'], aggregate['max_all']),
        '%.8f' % aggregate['max_adj'],
        '%.8f' % np.quantile(aggregate['null_adj'], 0.95),
        '%.8f' % empirical_p(aggregate['null_adj'], aggregate['max_adj'])
    ]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--lattice', choices=('square', 'triangular'), required=True)
    p.add_argument('--sizes', nargs='+', type=int, default=DEFAULT_SIZES)
    p.add_argument('--observables', nargs='+', choices=ALL_OBS,
                   default=list(ALL_OBS))
    p.add_argument('--expected-n', type=int, default=200,
                   help='Set 0 to permit any common length at each T')
    p.add_argument('--run', type=int, default=None,
                   help='Choose one run if files from multiple runs exist')
    p.add_argument('--temp-digits', type=int, default=8,
                   help='Decimal places used to match temperature filenames')
    p.add_argument('--null-repeats', type=int, default=500)
    p.add_argument('--seed', type=int, default=20261112)
    p.add_argument('--plot-step', type=int, default=10,
                   help='Heatmaps every Nth T; 1=all, 0=disable plots')
    args = p.parse_args()
    if args.null_repeats < 20:
        p.error('--null-repeats must be >= 20')
    if args.plot_step < 0:
        p.error('--plot-step must be >= 0')
    if args.expected_n < 0:
        p.error('--expected-n must be >= 0')
    if not (0 <= args.temp_digits <= 12):
        p.error('--temp-digits must be in [0,12]')
    if len(set(args.sizes)) != len(args.sizes) or any(x <= 0 for x in args.sizes):
        p.error('--sizes must be positive and unique')
    if len(set(args.observables)) != len(args.observables):
        p.error('--observables must not contain duplicates')
    if not args.root.is_dir():
        p.error('Input root directory does not exist: %s' % args.root)
    rng = np.random.default_rng(args.seed)
    index = index_files(args.root, sorted(args.sizes), args.run, args.temp_digits)
    temps = sorted(index)
    if not temps:
        raise SystemExit('No matching input files found')
    args.output.mkdir(parents=True, exist_ok=True)

    by_obs = {obs: make_aggregate(args.null_repeats) for obs in args.observables}
    overall = make_aggregate(args.null_repeats)
    rows = []
    valid_temps = 0
    for t_index, temp in enumerate(temps):
        accepted, n = load_temperature(index[temp], args.lattice, args.expected_n)
        if len(accepted) < 2:
            print('T=%s: fewer than two usable sizes; skipping' % temp)
            continue
        valid_temps += 1
        sizes = [item[0] for item in accepted]
        runs = [item[1] for item in accepted]
        offsets_full = rng.integers(0, n, size=(args.null_repeats, len(sizes)))
        out = args.output / ('T_%s' % temp)
        out.mkdir(parents=True, exist_ok=True)
        for obs in args.observables:
            z, retained = standardized_columns(accepted, obs)
            if z is None:
                print('T=%s %s: fewer than two nonconstant series; skipping' % (temp, obs))
                continue
            if len(retained) < len(accepted):
                missing = sorted(set(range(len(accepted))) - set(retained))
                print('T=%s %s: constant series for sizes %s; excluding them' %
                      (temp, obs, [sizes[i] for i in missing]))
            sz = [sizes[i] for i in retained]
            run = [runs[i] for i in retained]
            cor = np.clip(z.T @ z, -1.0, 1.0)
            np.fill_diagonal(cor, 1.0)
            all_abs = np.abs(cor[np.triu_indices(len(sz), k=1)])
            adj_abs = np.abs(np.diag(cor, k=1))
            null_all, null_adj = null_maxima(z, offsets_full[:, retained])
            observed_all = float(all_abs.max())
            observed_adj = float(adj_abs.max())
            cutoff_all = float(np.quantile(null_all, 0.95))
            cutoff_adj = float(np.quantile(null_adj, 0.95))
            p_all = empirical_p(null_all, observed_all)
            p_adj = empirical_p(null_adj, observed_adj)
            write_matrix(out / ('correlations_%s.csv' % obs), sz, cor)
            write_pairs(out / ('pairs_%s.csv' % obs), sz, run, cor,
                        cutoff_all, cutoff_adj)
            if args.plot_step and (t_index % args.plot_step == 0 or t_index == len(temps)-1):
                plot_matrix(out / ('heatmap_%s.png' % obs), sz, cor, temp, obs)
            pair_count = len(all_abs)
            rows.append([
                args.lattice, str(temp), obs, len(sz), n, pair_count,
                '%.8f' % observed_all, '%.8f' % cutoff_all, '%.8f' % p_all,
                '%.8f' % observed_adj, '%.8f' % cutoff_adj, '%.8f' % p_adj
            ])
            update_aggregate(by_obs[obs], observed_all, observed_adj,
                             null_all, null_adj, pair_count)
            update_aggregate(overall, observed_all, observed_adj,
                             null_all, null_adj, pair_count)
            print('T=%s %s: %d sizes; max |r|=%.3f (p=%.4f); '
                  'neighbor max |r|=%.3f (p=%.4f)' %
                  (temp, obs, len(sz), observed_all, p_all, observed_adj, p_adj))

    if not rows:
        raise SystemExit('No results; verify root, files, and overlap of temperatures')
    write_csv(args.output / 'summary.csv', [
        'lattice', 'T', 'observable', 'n_sizes', 'n_measurements',
        'n_size_pairs', 'max_abs_r', 'all_pairs_null_95', 'all_pairs_p',
        'adjacent_max_abs_r', 'adjacent_null_95', 'adjacent_p'
    ], rows)
    global_rows = [summary_row(obs, by_obs[obs]) for obs in args.observables
                   if by_obs[obs]['n_tests'] > 0]
    write_csv(args.output / 'global_summary.csv', [
        'observable', 'n_temperature_tests', 'n_size_pair_tests',
        'max_abs_r', 'global_all_pairs_null_95', 'global_all_pairs_p',
        'max_adjacent_abs_r', 'global_adjacent_null_95', 'global_adjacent_p'
    ], global_rows)
    write_csv(args.output / 'overall_summary.csv', [
        'observable', 'n_temperature_observable_tests', 'n_size_pair_tests',
        'max_abs_r', 'global_all_pairs_null_95', 'global_all_pairs_p',
        'max_adjacent_abs_r', 'global_adjacent_null_95', 'global_adjacent_p'
    ], [summary_row('ALL', overall)])
    print('Analyzed %d temperatures; results: %s' % (valid_temps, args.output))
    print('global_summary.csv: maxima calibrated across all temperatures for each observable.')
    print('overall_summary.csv: maxima calibrated across all T, L pairs, and observables.')
    print('CAUTION: circular-shift tests detect zero-lag linear dependence only;')
    print('non-significant results are not evidence that streams are independent.')


if __name__ == '__main__':
    main()
