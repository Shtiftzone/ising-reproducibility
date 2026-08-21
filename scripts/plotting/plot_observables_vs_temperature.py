#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PLOTS = {
    "mag_vs_T.csv": {
        "output": "mag_vs_T_abs.png",
        "title": "Magnetization vs Temperature",
        "ylabel": r"$\langle |M| \rangle$",
    },
    "energy_vs_T.csv": {
        "output": "energy_vs_T.png",
        "title": "Energy vs Temperature",
        "ylabel": r"$\langle E \rangle$",
    },
    "euler_sym_vs_T.csv": {
        "output": "euler_sym_vs_T_abs.png",
        "title": (
            r"Symmetrized Euler characteristic "
            r"$EC_{\mathrm{sym}}$ vs Temperature"
        ),
        "ylabel": r"$\langle |EC_{\mathrm{sym}}| \rangle$",
    },
    "ec_avg_vs_T.csv": {
        "output": "ec_avg_vs_T.png",
        "title": (
            r"Average Euler characteristic "
            r"$EC_{\mathrm{avg}}$ vs Temperature"
        ),
        "ylabel": r"$\langle EC_{\mathrm{avg}} \rangle$",
    },
}


def make_plot(
    csv_path: Path,
    output_path: Path,
    title: str,
    ylabel: str,
) -> None:
    df = pd.read_csv(csv_path)

    required = {
        "T",
        "mean",
        "jackknife_se",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{csv_path} is missing required columns: "
            f"{sorted(missing)}"
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
        "--lattice-type",
        choices=["square", "triangular"],
        required=True,
    )

    parser.add_argument(
        "--representation",
        choices=["cell", "vertex"],
        required=True,
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/figure_csv"),
        help="Base directory containing figure-level CSV files.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory for generated figures.",
    )

    args = parser.parse_args()

    variant_input_dir = (
        args.input_dir
        / args.lattice_type
        / args.representation
    )

    if not variant_input_dir.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist: "
            f"{variant_input_dir}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    lattice_label = (
        "Square"
        if args.lattice_type == "square"
        else "Triangular"
    )

    representation_label = (
        "spin as cell"
        if args.representation == "cell"
        else "spin as vertex"
    )

    prefix = (
        f"{args.lattice_type}_"
        f"{args.representation}"
    )

    for csv_name, config in PLOTS.items():
        csv_path = (
            variant_input_dir
            / csv_name
        )

        if not csv_path.is_file():
            print(
                f"Skipping missing file: "
                f"{csv_path}"
            )
            continue

        output_path = (
            args.output_dir
            / f"{prefix}_{config['output']}"
        )

        title = (
            f"{lattice_label} Ising, "
            f"{representation_label}: "
            f"{config['title']}"
        )

        make_plot(
            csv_path,
            output_path,
            title,
            config["ylabel"],
        )

    print("Done.")


if __name__ == "__main__":
    main()