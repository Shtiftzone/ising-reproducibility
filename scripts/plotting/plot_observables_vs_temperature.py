#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PLOTS = {
    "square_mag_vs_T.csv": {
        "output": "2d_mag_vs_T_abs.png",
        "title": "Magnetization vs Temperature (2D)",
        "ylabel": r"$\langle |M| \rangle$",
    },
    "square_energy_vs_T.csv": {
        "output": "2d_energy_vs_T.png",
        "title": "Energy vs Temperature (2D)",
        "ylabel": r"$\langle E \rangle$",
    },
    "square_euler_sym_vs_T.csv": {
        "output": "2d_euler_sym_vs_T_abs.png",
        "title": r"Symmetrized Euler characteristic $EC_{\mathrm{sym}}$ vs Temperature (2D)",
        "ylabel": r"$\langle |EC_{\mathrm{sym}}| \rangle$",
    },
    "square_ec_avg_vs_T.csv": {
        "output": "2d_ec_avg_vs_T.png",
        "title": r"Average Euler characteristic $EC_{\mathrm{avg}}$ vs Temperature (2D)",
        "ylabel": r"$\langle EC_{\mathrm{avg}} \rangle$",
    },
    "triangular_mag_vs_T.csv": {
        "output": "tri_mag_vs_T_abs.png",
        "title": "Magnetization vs Temperature (Tri.)",
        "ylabel": r"$\langle |M| \rangle$",
    },
    "triangular_energy_vs_T.csv": {
        "output": "tri_energy_vs_T.png",
        "title": "Energy vs Temperature (Tri.)",
        "ylabel": r"$\langle E \rangle$",
    },
    "triangular_euler_sym_vs_T.csv": {
        "output": "tri_euler_sym_vs_T_abs.png",
        "title": r"Symmetrized Euler characteristic $EC_{\mathrm{sym}}$ vs Temperature (Tri.)",
        "ylabel": r"$\langle |EC_{\mathrm{sym}}| \rangle$",
    },
    "triangular_ec_avg_vs_T.csv": {
        "output": "tri_ec_avg_vs_T.png",
        "title": r"Average Euler characteristic $EC_{\mathrm{avg}}$ vs Temperature (Tri.)",
        "ylabel": r"$\langle EC_{\mathrm{avg}} \rangle$",
    },
}


def make_plot(csv_path: Path, output_path: Path, title: str, ylabel: str) -> None:
    df = pd.read_csv(csv_path)

    required = {"T", "mean", "jackknife_se"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{csv_path} is missing required columns: {sorted(missing)}"
        )

    df = df.sort_values("T")

    plt.figure(figsize=(6, 4))

    plt.errorbar(
        df["T"],
        df["mean"],
        yerr=df["jackknife_se"],
        fmt="o",
        ms=3,
        capsize=2,
        elinewidth=1,
        alpha=0.9,
    )

    plt.title(title)
    plt.xlabel(r"Temperature $T$")
    plt.ylabel(ylabel)

    plt.grid(
        True,
        linewidth=0.3,
        alpha=0.5,
    )

    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close()

    print(f"Wrote {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate temperature-dependence plots from "
            "figure-level CSV files."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/figure_csv"),
        help="Directory containing figure-level CSV files.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory for generated figures.",
    )

    args = parser.parse_args()

    if not args.input_dir.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist: {args.input_dir}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for csv_name, config in PLOTS.items():
        csv_path = args.input_dir / csv_name

        if not csv_path.is_file():
            print(f"Skipping missing file: {csv_path}")
            continue

        output_path = (
            args.output_dir
            / config["output"]
        )

        make_plot(
            csv_path,
            output_path,
            config["title"],
            config["ylabel"],
        )

    print("Done.")


if __name__ == "__main__":
    main()
