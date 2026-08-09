#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot KS distance and effective-resolution ratios "
            "from processed distribution-analysis CSV files."
        )
    )

    parser.add_argument(
        "--lattice-type",
        choices=["square", "triangular"],
        required=True,
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/analysis_csv"),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figures"),
    )

    args = parser.parse_args()

    input_dir = (
        args.input_dir
        / args.lattice_type
        / "distribution_resolution"
    )

    shape_path = (
        input_dir
        / "distribution_shape_tests.csv"
    )

    entropy_path = (
        input_dir
        / "entropy_effective_resolution.csv"
    )

    if not shape_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {shape_path}"
        )

    if not entropy_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {entropy_path}"
        )

    df_shape = pd.read_csv(shape_path)
    df_entropy = pd.read_csv(entropy_path)

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    temperature = float(
        df_shape["temperature"].iloc[0]
    )

    plt.figure(figsize=(7, 5))

    plt.plot(
        df_shape["L"],
        df_shape["KS_statistic_standardized"],
        marker="o",
    )

    plt.xlabel("L")
    plt.ylabel("KS distance")
    plt.title(
        rf"Distance between standardized distributions, "
        rf"$T={temperature:.5f}$"
    )

    plt.grid(
        True,
        alpha=0.3,
    )

    plt.tight_layout()

    ks_output = (
        args.output_dir
        / "KS_distance_standardized_vs_L.png"
    )

    plt.savefig(
        ks_output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(f"Wrote {ks_output}")

    plt.figure(figsize=(7, 5))

    plt.plot(
        df_entropy["L"],
        df_entropy["unique_ratio_EC_over_M"],
        marker="o",
        label=r"$K_{EC}/K_M$",
    )

    plt.plot(
        df_entropy["L"],
        df_entropy["Keff_ratio_EC_over_M"],
        marker="o",
        label=(
            r"$K_{\mathrm{eff},EC}/"
            r"K_{\mathrm{eff},M}$"
        ),
    )

    plt.axhline(
        1.0,
        linestyle="--",
        linewidth=1,
    )

    plt.xlabel("L")
    plt.ylabel("EC / M resolution ratio")
    plt.title(
        rf"Relative finite-size resolution, "
        rf"$T={temperature:.5f}$"
    )

    plt.legend()

    plt.xscale(
        "log",
        base=2,
    )

    plt.grid(
        True,
        alpha=0.3,
    )

    plt.tight_layout()

    resolution_output = (
        args.output_dir
        / "resolution_ratios_EC_over_M_vs_L.png"
    )

    plt.savefig(
        resolution_output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(f"Wrote {resolution_output}")
    print("Done.")


if __name__ == "__main__":
    main()
