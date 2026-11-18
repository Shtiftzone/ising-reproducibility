#!/usr/bin/env python3
"""Check measured-energy histogram overlap at neighboring temperatures.

Expected layout:
    ROOT/T_2.27000/size_3072/mefile.txt

mefile.txt must contain comma-separated columns M,E,..., as in the
organized Ising simulation outputs. Uses the SECOND column (total energy).

This is a diagnostic of recorded energy distributions, not an independent
validation of multihistogram reweighting. Histogram overlap depends on binning
and on how many samples were recorded. ASCII source only.
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TEMPERATURE_RE = re.compile(r"^T_([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$")
DEFAULT_SIZES = [64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072]


def read_energies(path, size):
    try:
        values = np.loadtxt(path, delimiter=",", usecols=(1,), ndmin=1)
    except ValueError as exc:
        raise ValueError("Cannot read second comma-separated column from %s: %s" % (path, exc))
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size < 4 or not np.isfinite(values).all():
        raise ValueError("Need >=4 finite energy values in %s" % path)
    return values / float(size * size)


def gather(root, size):
    series = []
    seen = set()
    for directory in root.iterdir():
        if not directory.is_dir():
            continue
        match = TEMPERATURE_RE.fullmatch(directory.name)
        if match is None:
            continue
        temperature = float(match.group(1))
        path = directory / ("size_%d" % size) / "mefile.txt"
        if not path.is_file():
            continue
        if temperature in seen:
            raise ValueError("Duplicate temperature %g at L=%d" % (temperature, size))
        seen.add(temperature)
        series.append((temperature, read_energies(path, size)))
    series.sort(key=lambda item: item[0])
    return series


def overlap_metrics(a, b, n_bins):
    combined = np.concatenate((a, b))
    low = float(combined.min())
    high = float(combined.max())
    if not high > low:
        # Both distributions are constant at exactly the same energy.
        return 1.0, 0.0, np.asarray([low - 0.5, high + 0.5])
    margin = max(1e-12, (high - low) * 1e-6)
    edges = np.linspace(low - margin, high + margin, n_bins + 1)
    counts_a, _ = np.histogram(a, bins=edges)
    counts_b, _ = np.histogram(b, bins=edges)
    probs_a = counts_a / counts_a.sum()
    probs_b = counts_b / counts_b.sum()
    overlap = float(np.minimum(probs_a, probs_b).sum())
    sd_a = float(np.std(a, ddof=1))
    sd_b = float(np.std(b, ddof=1))
    pooled_sd = np.sqrt((sd_a * sd_a + sd_b * sd_b) / 2.0)
    separation = float(abs(np.mean(a) - np.mean(b)) / pooled_sd) if pooled_sd > 0 else float("inf")
    return overlap, separation, edges


def plot_pair(path, size, ta, tb, a, b, ovl, separation, edges):
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.hist(a, bins=edges, density=True, alpha=0.48, label="T = %.6f (n=%d)" % (ta, a.size))
    ax.hist(b, bins=edges, density=True, alpha=0.48, label="T = %.6f (n=%d)" % (tb, b.size))
    ax.set_xlabel("Energy per spin, E / L^2")
    ax.set_ylabel("Probability density")
    ax.set_title("L=%d: adjacent-temperature energy histograms" % size)
    ax.text(0.02, 0.97, "Histogram overlap = %.3f\nMean separation = %.2f pooled SD" % (ovl, separation),
            ha="left", va="top", transform=ax.transAxes,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"})
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_overview(path, size, rows):
    tmid = [0.5 * (r["T_low"] + r["T_high"]) for r in rows]
    ovl = [r["hist_overlap"] for r in rows]
    fig, ax = plt.subplots(figsize=(7.4, 4))
    ax.plot(tmid, ovl, "o-", markersize=3, linewidth=1)
    ax.set_xlabel("Midpoint of adjacent measured temperatures")
    ax.set_ylabel("Histogram overlap coefficient")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("L=%d: energy histogram overlap" % size)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_csv(path, names, rows):
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Root with T_*/size_L/mefile.txt")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", nargs="+", type=int, default=DEFAULT_SIZES)
    parser.add_argument("--bins", type=int, default=25, help="Shared histogram bins for each pair")
    parser.add_argument("--plots-per-size", type=int, default=6,
                        help="Number of worst-overlap adjacent pairs to plot for each size")
    args = parser.parse_args()
    if not args.root.is_dir():
        parser.error("--root does not exist: %s" % args.root)
    if args.bins < 5 or args.plots_per_size < 0:
        parser.error("--bins must be >= 5 and --plots-per-size must be >= 0")
    args.output.mkdir(parents=True, exist_ok=True)

    all_rows = []
    size_rows = []
    for size in args.sizes:
        try:
            series = gather(args.root, size)
        except (ValueError, OSError) as exc:
            print("L=%d: ERROR: %s" % (size, exc))
            continue
        if len(series) < 2:
            print("L=%d: fewer than two temperatures with energy data, skipping" % size)
            continue
        rows = []
        for (ta, a), (tb, b) in zip(series[:-1], series[1:]):
            overlap, separation, edges = overlap_metrics(a, b, args.bins)
            rows.append({
                "L": size, "T_low": ta, "T_high": tb, "delta_T": tb - ta,
                "n_low": a.size, "n_high": b.size,
                "mean_E_per_spin_low": float(np.mean(a)),
                "mean_E_per_spin_high": float(np.mean(b)),
                "sd_E_per_spin_low": float(np.std(a, ddof=1)),
                "sd_E_per_spin_high": float(np.std(b, ddof=1)),
                "hist_overlap": overlap,
                "mean_separation_pooled_sd": separation,
            })
        all_rows.extend(rows)
        overlaps = np.array([r["hist_overlap"] for r in rows])
        size_rows.append({"L": size, "n_temperatures": len(series),
                          "n_pairs": len(rows), "min_overlap": float(np.min(overlaps)),
                          "median_overlap": float(np.median(overlaps)),
                          "fraction_below_0_1": float(np.mean(overlaps < 0.1)),
                          "fraction_below_0_2": float(np.mean(overlaps < 0.2))})
        plot_overview(args.output / ("overlap_vs_T_L%d.png" % size), size, rows)

        selected = sorted(range(len(rows)), key=lambda index: rows[index]["hist_overlap"])
        folder = args.output / ("size_%d" % size)
        folder.mkdir(exist_ok=True)
        for index in selected[:args.plots_per_size]:
            ta, a = series[index]
            tb, b = series[index + 1]
            ovl, sep, edges = overlap_metrics(a, b, args.bins)
            path = folder / ("worst_pair_%03d.png" % (index + 1))
            plot_pair(path, size, ta, tb, a, b, ovl, sep, edges)
        print("L=%d: %d pairs, min overlap=%.3f, median=%.3f; plots=%d" %
              (size, len(rows), float(overlaps.min()), float(np.median(overlaps)),
               min(len(rows), args.plots_per_size)))

    if not all_rows:
        raise SystemExit("No pairs found. Check input layout and --sizes.")
    fields = ["L", "T_low", "T_high", "delta_T", "n_low", "n_high",
              "mean_E_per_spin_low", "mean_E_per_spin_high",
              "sd_E_per_spin_low", "sd_E_per_spin_high",
              "hist_overlap", "mean_separation_pooled_sd"]
    write_csv(args.output / "adjacent_pairs.csv", fields, all_rows)
    write_csv(args.output / "summary_by_size.csv",
              ["L", "n_temperatures", "n_pairs", "min_overlap", "median_overlap",
               "fraction_below_0_1", "fraction_below_0_2"], size_rows)
    print("Output:", args.output)
    print("NOTE: Bin-dependent overlap is exploratory; low overlap is a warning,")
    print("not a proof that a multihistogram estimate is invalid.")


if __name__ == "__main__":
    main()
