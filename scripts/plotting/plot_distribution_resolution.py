#!/usr/bin/env python3
"""Plot finite-size resolution ratios for M, EC_sym, and D.

Input:
  data/analysis_csv/<lattice>/<representation>/distribution_resolution/
      entropy_effective_resolution.csv

Output:
  results/figures/<lattice>/<representation>/
      <lattice>_<representation>_resolution_ratios_vs_L.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RATIO_COLUMNS = (
    "unique_ratio_EC_over_M",
    "unique_ratio_D_over_M",
    "Keff_ratio_EC_over_M",
    "Keff_ratio_D_over_M",
)


def load_entropy_csv(path, lattice_type, representation):
    if not path.is_file():
        raise FileNotFoundError(f"Missing input file: {path}")

    df = pd.read_csv(path)
    required = {"lattice_type", "representation", "temperature", "L", *RATIO_COLUMNS}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if df.empty:
        raise ValueError(f"{path}: CSV is empty")

    for column, expected in (("lattice_type", lattice_type),
                             ("representation", representation)):
        if not df[column].eq(expected).all():
            raise ValueError(f"{path}: {column} does not match {expected!r}")

    numeric_columns = ["temperature", "L", *RATIO_COLUMNS]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="raise")
    if not np.isfinite(df[numeric_columns].to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: non-finite temperatures, sizes, or ratios")
    if (df["L"] <= 0).any() or df["L"].duplicated().any():
        raise ValueError(f"{path}: L must be positive and unique")

    temperatures = df["temperature"].to_numpy(dtype=float)
    if not np.allclose(temperatures, temperatures[0], rtol=0, atol=1e-9):
        raise ValueError(f"{path}: multiple temperatures present")

    return df.sort_values("L"), float(temperatures[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice-type", choices=["square", "triangular"], required=True)
    parser.add_argument("--representation", choices=["cell", "vertex"], required=True)
    parser.add_argument("--input-dir", type=Path, default=Path("data/analysis_csv"),
                        help="Base input directory (before lattice/representation).")
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures"),
                        help="Base output directory (before lattice/representation).")
    args = parser.parse_args()

    input_path = (args.input_dir / args.lattice_type / args.representation
                  / "distribution_resolution" / "entropy_effective_resolution.csv")
    output_dir = args.output_dir / args.lattice_type / args.representation
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = (output_dir /
                   f"{args.lattice_type}_{args.representation}_resolution_ratios_vs_L.png")

    df, temperature = load_entropy_csv(input_path, args.lattice_type, args.representation)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(df["L"], df["unique_ratio_EC_over_M"], marker="o", linestyle="--",
            label=r"$K_{EC}/K_M$")
    ax.plot(df["L"], df["unique_ratio_D_over_M"], marker="o", linestyle="--",
            label=r"$K_D/K_M$")
    ax.plot(df["L"], df["Keff_ratio_EC_over_M"], marker="o",
            label=r"$K_{\mathrm{eff},EC}/K_{\mathrm{eff},M}$")
    ax.plot(df["L"], df["Keff_ratio_D_over_M"], marker="o",
            label=r"$K_{\mathrm{eff},D}/K_{\mathrm{eff},M}$")

    ax.axhline(1.0, linestyle=":", linewidth=1)
    ax.set_xlabel(r"$L$")
    ax.set_ylabel("Resolution ratio relative to M")
    ax.set_title(rf"Relative finite-size resolution, $T={temperature:.5f}$")
    ax.set_xscale("log", base=2)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
