#!/usr/bin/env python3
"""Compute Ising temperature plots with nonoverlapping block-jackknife errors.

Expected input structure:
  RESULTS_DIR/T_<temperature>/size_<L>/mefile.txt
  RESULTS_DIR/T_<temperature>/size_<L>/epfile.txt
  RESULTS_DIR/T_<temperature>/size_<L>/enfile.txt

The files contain comma-separated measurements in recorded time order.
The first two columns of mefile are M, H; the first column of epfile
and enfile is the Euler characteristic of the positive/negative complex.

Observables: |M|, H, |(EC_plus - EC_minus)/2|, (EC_plus + EC_minus)/2.

For each (L,T), tau_int is estimated separately for each of the four
series, and a common nonoverlapping block size is selected from the
largest tau_int. The error is obtained by deleting one complete block
at a time. Autocorrelation estimates are exploratory at small sample
counts; inspect the diagnostics and repeat with fixed block sizes.

The script keeps the old CSV columns, including 'jackknife_se', so the
existing plotting script works without modification. Additional columns
and a separate diagnostics CSV report block sizes and autocorrelations.

Source code is ASCII-only.
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

TEMP_PATTERN = re.compile(r'^T_([-+]?\d+(?:\.\d+)?)$')

OBSERVABLES = (
    ('mag_abs', 'abs_magnetization', 'mag_vs_T.csv'),
    ('energy', 'energy', 'energy_vs_T.csv'),
    ('euler_sym_abs', 'abs_ec_sym', 'euler_sym_vs_T.csv'),
    ('ec_avg', 'ec_avg', 'ec_avg_vs_T.csv'),
)


def parse_temperature(dirname):
    match = TEMP_PATTERN.fullmatch(dirname)
    if not match:
        raise ValueError('Unrecognized temperature directory: ' + dirname)
    return float(match.group(1))


def read_float_matrix(path):
    """Read comma-separated numeric rows, rejecting invalid data."""
    rows = []
    with path.open('r', encoding='utf-8') as handle:
        for lineno, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                values = [float(item.strip()) for item in line.split(',')]
            except ValueError as exc:
                raise ValueError('%s:%d: expected comma-separated numbers' %
                                 (path, lineno)) from exc
            if rows and len(values) != len(rows[0]):
                raise ValueError('%s:%d: inconsistent number of columns' %
                                 (path, lineno))
            rows.append(values)
    if not rows:
        raise ValueError('No numeric measurements in ' + str(path))
    result = np.asarray(rows, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite measurements in ' + str(path))
    return result


def autocorrelation_time(values, window_factor=5.0, max_lag_fraction=0.25):
    """Estimate tau_int = 1/2 + sum_{lag=1}^W rho(lag).

    Uses a biased autocovariance estimate (denominator n) to avoid large
    noisy correlations at high lags. W is the first lag that meets the
    self-consistent window condition W >= window_factor * tau_int(W).
    Window search stops at floor(n * max_lag_fraction). The tau used to
    select blocks is never less than 0.5. This is a practical heuristic,
    not a guaranteed conservative error estimator.
    """
    x = np.asarray(values, dtype=float)
    n = len(x)
    if n < 4:
        raise ValueError('At least 4 samples required for tau_int')
    centered = x - x.mean()
    variance = float(np.dot(centered, centered) / n)
    if variance <= 0.0:
        return 0.5, 0, True, False

    max_lag = max(1, min(n - 1, int(n * max_lag_fraction)))
    # FFT autocovariance; zero padding avoids circular aliasing.
    fft_size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=fft_size)
    autocovariance = np.fft.irfft(spectrum * np.conj(spectrum), n=fft_size)
    rho = autocovariance[:max_lag + 1] / autocovariance[0]
    partial_sum = 0.0
    tau = 0.5
    for lag in range(1, max_lag + 1):
        partial_sum += float(rho[lag])
        tau = max(0.5, 0.5 + partial_sum)
        if lag >= window_factor * tau:
            return tau, lag, False, False
    return tau, max_lag, False, True


def select_block_size(n, tau_max, block_size, tau_multiplier, min_blocks):
    """Select equal-size blocks and flag auto targets exceeding the limit."""
    max_block_size = n // min_blocks
    if max_block_size < 1:
        raise ValueError('Only %d samples; need at least %d blocks' %
                         (n, min_blocks))

    capped = False
    if block_size:
        size = block_size
    else:
        target = max(1, int(math.ceil(tau_multiplier * tau_max)))
        capped = target > max_block_size
        # With short time series, the desired block length can conflict
        # with the minimum number of jackknife blocks. Cap and flag it.
        # This is a practical fallback, not proof of negligible bias.
        effective_target = min(target, max_block_size)
        divisors = [b for b in range(effective_target, max_block_size + 1)
                    if n % b == 0]
        # If no exact divisor is large enough, drop the trailing samples.
        size = divisors[0] if divisors else effective_target

    if size < 1 or size > max_block_size:
        raise ValueError('Block size %d gives fewer than %d blocks for n=%d' %
                         (size, min_blocks, n))
    n_blocks = n // size
    if n_blocks < min_blocks:
        raise ValueError('Not enough blocks for the requested minimum')
    return size, n_blocks, n_blocks * size, capped


def block_jackknife_mean_and_se(values, block_size, n_used):
    """Mean and delete-one-block jackknife SE from equally sized blocks."""
    x = np.asarray(values[:n_used], dtype=float)
    n_blocks = n_used // block_size
    blocks = x.reshape(n_blocks, block_size)
    block_sums = blocks.sum(axis=1)
    total = float(block_sums.sum())
    replicas = (total - block_sums) / (n_used - block_size)
    replicas_mean = float(replicas.mean())
    jackknife_var = (n_blocks - 1) / n_blocks * np.sum(
        (replicas - replicas_mean) ** 2
    )
    return float(total / n_used), float(np.sqrt(max(0.0, jackknife_var)))


def ordinary_se_of_mean(values):
    """Old delete-one-sample jackknife SE, equivalent to sample SD / sqrt(n)."""
    x = np.asarray(values, dtype=float)
    if len(x) < 2:
        return float('nan')
    return float(np.std(x, ddof=1) / np.sqrt(len(x)))


def load_temperature_observables(size_dir):
    me = read_float_matrix(size_dir / 'mefile.txt')
    ep = read_float_matrix(size_dir / 'epfile.txt')
    en = read_float_matrix(size_dir / 'enfile.txt')
    if me.shape[1] < 2 or ep.shape[1] < 1 or en.shape[1] < 1:
        raise ValueError('Expected >=2 columns in mefile.txt and '
                         '>=1 column in epfile.txt / enfile.txt in ' + str(size_dir))
    if not (me.shape[0] == ep.shape[0] == en.shape[0]):
        raise ValueError('Unequal numbers of recorded measurements in %s: '
                         'me=%d, ep=%d, en=%d (cannot align all rows)' %
                         (size_dir, len(me), len(ep), len(en)))
    m = me[:, 0]
    h = me[:, 1]
    ec_plus = ep[:, 0]
    ec_minus = en[:, 0]
    return {
        'mag_abs': np.abs(m),
        'energy': h,
        'euler_sym_abs': np.abs(0.5 * (ec_plus - ec_minus)),
        'ec_avg': 0.5 * (ec_plus + ec_minus),
    }


def run_analysis(args):
    results_dir = args.results_dir
    if not results_dir.is_dir():
        raise FileNotFoundError('Simulation directory does not exist: ' +
                                str(results_dir))
    out = (args.output_dir / args.lattice_type / args.representation
           / 'observables_vs_temperature' / ('size_%d' % args.size))
    out.mkdir(parents=True, exist_ok=True)

    rows = {key: [] for key, _, _ in OBSERVABLES}
    diagnostics = []
    temperature_dirs = sorted(
        (p for p in results_dir.iterdir()
         if p.is_dir() and TEMP_PATTERN.fullmatch(p.name)),
        key=lambda p: parse_temperature(p.name)
    )
    if not temperature_dirs:
        raise RuntimeError('No T_... directories in ' + str(results_dir))

    for temp_dir in temperature_dirs:
        temp = parse_temperature(temp_dir.name)
        size_dir = temp_dir / ('size_%d' % args.size)
        if not size_dir.is_dir():
            print('Skipping T=%.6f: missing %s' % (temp, size_dir))
            continue
        missing = [name for name in ('mefile.txt', 'epfile.txt', 'enfile.txt')
                   if not (size_dir / name).is_file()]
        if missing:
            print('Skipping T=%.6f: missing %s' % (temp, ', '.join(missing)))
            continue

        data = load_temperature_observables(size_dir)
        n = len(data['mag_abs'])
        if args.expected_n and n != args.expected_n:
            raise ValueError('T=%.6f, L=%d: expected %d samples, got %d' %
                             (temp, args.size, args.expected_n, n))
        if n < args.min_blocks:
            raise ValueError('T=%.6f, L=%d: only %d samples' %
                             (temp, args.size, n))

        estimates = {
            key: autocorrelation_time(
                data[key], window_factor=args.window_factor,
                max_lag_fraction=args.max_lag_fraction)
            for key, _, _ in OBSERVABLES
        }
        tau_max = max(estimates[key][0] for key, _, _ in OBSERVABLES)
        block_size, n_blocks, n_used, auto_block_capped = select_block_size(
            n, tau_max, args.block_size, args.tau_multiplier,
            args.min_blocks
        )
        if auto_block_capped:
            target = int(math.ceil(args.tau_multiplier * tau_max))
            print('WARNING: T=%.6f L=%d: auto block target %d exceeds '
                  'maximum %d with n=%d and min_blocks=%d; using block %d. '
                  'Inspect sensitivity to block size.' %
                  (temp, args.size, target, n // args.min_blocks, n,
                   args.min_blocks, block_size))
        dropped = n - n_used
        if dropped:
            print('WARNING: T=%.6f L=%d: dropping %d trailing samples '
                  'for equal-size blocks' % (temp, args.size, dropped))

        for key, observable, _ in OBSERVABLES:
            x = data[key]
            tau, window, constant, window_limited = estimates[key]
            mean, jk_se = block_jackknife_mean_and_se(x, block_size, n_used)
            iid_se = ordinary_se_of_mean(x[:n_used])
            ratio = jk_se / iid_se if iid_se > 0 else float('nan')
            row = {
                'lattice_type': args.lattice_type,
                'representation': args.representation,
                'L': args.size,
                'T': temp,
                'observable': observable,
                'mean': mean,
                'jackknife_se': jk_se,
                'n_samples': n,
                'n_used': n_used,
                'block_size': block_size,
                'n_blocks': n_blocks,
                'tau_int': tau,
                'autocorr_window': window,
            }
            rows[key].append(row)
            diagnostics.append({
                **row,
                'old_unblocked_se': iid_se,
                'block_to_unblocked_se_ratio': ratio,
                'dropped_samples': dropped,
                'constant_series': int(constant),
                'window_hit_max_lag': int(window_limited),
                'tau_max_across_observables': tau_max,
                'auto_block_capped': int(auto_block_capped),
                'target_block_size': (
                    args.block_size if args.block_size else
                    int(math.ceil(args.tau_multiplier * tau_max))
                ),
            })
        print('T=%.6f L=%d: n=%d, tau_max=%.3f, block=%d, blocks=%d' %
              (temp, args.size, n, tau_max, block_size, n_blocks))

    if not diagnostics:
        raise RuntimeError('No results; check input directory and lattice size')

    for key, _, filename in OBSERVABLES:
        path = out / filename
        df = pd.DataFrame(rows[key]).sort_values(['L', 'T'])
        df.to_csv(path, index=False, float_format='%.10g')
        print('Wrote ' + str(path))

    diag_path = out / 'autocorrelation_diagnostics.csv'
    diag_df = pd.DataFrame(diagnostics).sort_values(['L', 'T', 'observable'])
    diag_df.to_csv(diag_path, index=False, float_format='%.10g')
    print('Wrote ' + str(diag_path))
    print('Done. Review tau_int, block sizes, number of blocks, and SE ratios.')
    print('Note: block jackknife handles within-run time correlations; it does '
          'not test cross-temperature or cross-size dependence.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, required=True)
    parser.add_argument('--lattice-type', choices=('square', 'triangular'),
                        required=True)
    parser.add_argument('--representation', choices=('cell', 'vertex'),
                        required=True)
    parser.add_argument('--size', type=int, default=3072)
    parser.add_argument('--output-dir', type=Path,
                        default=Path('data/analysis_csv'))
    parser.add_argument('--expected-n', type=int, default=0,
                        help='Expected number of recorded samples; 0 disables check')
    parser.add_argument('--block-size', type=int, default=0,
                        help='Fixed block size; 0 selects from tau_int (default)')
    parser.add_argument('--tau-multiplier', type=float, default=4.0,
                        help='Minimum auto block length divided by max tau_int')
    parser.add_argument('--window-factor', type=float, default=5.0,
                        help='Self-consistent autocorrelation window factor')
    parser.add_argument('--max-lag-fraction', type=float, default=0.25,
                        help='Max autocorrelation lag as a fraction of n')
    parser.add_argument('--min-blocks', type=int, default=10,
                        help='Minimum number of nonoverlapping blocks')
    args = parser.parse_args()
    if args.size < 1 or args.expected_n < 0 or args.block_size < 0:
        parser.error('--size must be positive; --expected-n and --block-size >= 0')
    if args.tau_multiplier <= 0 or args.window_factor <= 0:
        parser.error('--tau-multiplier and --window-factor must be positive')
    if not (0 < args.max_lag_fraction < 1):
        parser.error('--max-lag-fraction must be between 0 and 1')
    if args.min_blocks < 3:
        parser.error('--min-blocks must be at least 3')
    run_analysis(args)


if __name__ == '__main__':
    main()
