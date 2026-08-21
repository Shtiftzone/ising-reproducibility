#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis, ks_2samp


DEFAULT_SIZES = [
    64, 96, 128, 192, 256, 384,
    512, 768, 1024, 1536, 2048, 3072,
]


def read_first_column(path: Path) -> np.ndarray:
    data = np.loadtxt(path, delimiter=",")

    if data.ndim == 1:
        return data.astype(float)

    return data[:, 0].astype(float)


def load_observables(size_dir: Path):
    me_path = size_dir / "mefile.txt"
    ep_path = size_dir / "epfile.txt"
    en_path = size_dir / "enfile.txt"

    for path in (me_path, ep_path, en_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    M = read_first_column(me_path)
    EP = read_first_column(ep_path)
    EN = read_first_column(en_path)

    n = min(len(M), len(EP), len(EN))

    if n == 0:
        raise ValueError(f"No samples found in {size_dir}")

    M = M[:n]
    EP = EP[:n]
    EN = EN[:n]

    EC_sym = (EN - EP) / 2.0

    return M, EC_sym


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)

    std = np.std(values, ddof=1)

    if not np.isfinite(std) or std == 0.0:
        return np.full(values.shape, np.nan, dtype=float)

    return (
        values - np.mean(values)
    ) / std


def entropy_effective_values(values: np.ndarray):
    values = np.asarray(values)

    unique_values, counts = np.unique(
        values,
        return_counts=True,
    )

    probabilities = counts / counts.sum()

    entropy = -np.sum(
        probabilities * np.log(probabilities)
    )

    k_eff = np.exp(entropy)

    return (
        int(len(unique_values)),
        float(entropy),
        float(k_eff),
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare magnetization and symmetrized Euler-characteristic "
            "distributions at a selected temperature."
        )
    )

    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Raw simulation results directory.",
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
        "--temperature",
        type=float,
        required=True,
        help=(
            "Temperature whose raw simulation distributions "
            "are analyzed."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/analysis_csv"),
    )

    parser.add_argument(
        "--sizes",
        nargs="+",
        type=int,
        default=DEFAULT_SIZES,
    )

    args = parser.parse_args()

    temperature_dir = (
        args.results_dir
        / f"T_{args.temperature:.5f}"
    )

    if not temperature_dir.is_dir():
        raise FileNotFoundError(
            f"Temperature directory does not exist: "
            f"{temperature_dir}"
        )

    output_dir = (
        args.output_dir
        / args.lattice_type
        / args.representation
        / "distribution_resolution"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    shape_records = []
    entropy_records = []

    for L in args.sizes:
        size_dir = (
            temperature_dir
            / f"size_{L}"
        )

        if not size_dir.is_dir():
            print(
                f"Skipping L={L}: "
                f"missing {size_dir}"
            )
            continue

        try:
            M, EC_sym = load_observables(size_dir)
        except FileNotFoundError as exc:
            print(
                f"Skipping L={L}: "
                f"missing {exc.filename}"
            )
            continue

        Z_M = standardize(M)
        Z_EC = standardize(EC_sym)

        Z_M = Z_M[np.isfinite(Z_M)]
        Z_EC = Z_EC[np.isfinite(Z_EC)]

        if len(Z_M) == 0 or len(Z_EC) == 0:
            print(
                f"Skipping L={L}: "
                "standardized distribution is empty"
            )
            continue

        ks_result = ks_2samp(
            Z_M,
            Z_EC,
        )

        shape_records.append({
            "lattice_type": args.lattice_type,
            "representation": args.representation,
            "temperature": args.temperature,
            "L": L,
            "n_samples": len(M),

            "M_mean": float(np.mean(M)),
            "M_std": float(np.std(M, ddof=1)),
            "M_skewness": float(
                skew(M, bias=False)
            ),
            "M_excess_kurtosis": float(
                kurtosis(
                    M,
                    fisher=True,
                    bias=False,
                )
            ),

            "EC_mean": float(np.mean(EC_sym)),
            "EC_std": float(
                np.std(EC_sym, ddof=1)
            ),
            "EC_skewness": float(
                skew(EC_sym, bias=False)
            ),
            "EC_excess_kurtosis": float(
                kurtosis(
                    EC_sym,
                    fisher=True,
                    bias=False,
                )
            ),

            "KS_statistic_standardized":
                float(ks_result.statistic),

            "KS_pvalue_standardized":
                float(ks_result.pvalue),
        })

        K_M, H_M, Keff_M = (
            entropy_effective_values(M)
        )

        K_EC, H_EC, Keff_EC = (
            entropy_effective_values(EC_sym)
        )

        entropy_records.append({
            "lattice_type": args.lattice_type,
            "representation": args.representation,
            "temperature": args.temperature,
            "L": L,
            "n_samples": len(M),

            "M_unique": K_M,
            "M_entropy": H_M,
            "M_Keff": Keff_M,

            "EC_unique": K_EC,
            "EC_entropy": H_EC,
            "EC_Keff": Keff_EC,

            "unique_ratio_EC_over_M": (
                K_EC / K_M
                if K_M > 0
                else np.nan
            ),

            "Keff_ratio_EC_over_M": (
                Keff_EC / Keff_M
                if Keff_M > 0
                else np.nan
            ),

            "M_unique_fraction": (
                K_M / len(M)
            ),

            "EC_unique_fraction": (
                K_EC / len(EC_sym)
            ),

            "M_Keff_fraction": (
                Keff_M / len(M)
            ),

            "EC_Keff_fraction": (
                Keff_EC / len(EC_sym)
            ),
        })

    if not shape_records:
        raise RuntimeError(
            "No distribution-shape results were produced."
        )

    if not entropy_records:
        raise RuntimeError(
            "No entropy-resolution results were produced."
        )

    df_shape = (
        pd.DataFrame(shape_records)
        .sort_values("L")
        .reset_index(drop=True)
    )

    df_entropy = (
        pd.DataFrame(entropy_records)
        .sort_values("L")
        .reset_index(drop=True)
    )

    shape_path = (
        output_dir
        / "distribution_shape_tests.csv"
    )

    entropy_path = (
        output_dir
        / "entropy_effective_resolution.csv"
    )

    df_shape.to_csv(
        shape_path,
        index=False,
        float_format="%.10g",
    )

    df_entropy.to_csv(
        entropy_path,
        index=False,
        float_format="%.10g",
    )

    print(f"Lattice type: {args.lattice_type}")
    print(f"Representation: {args.representation}")
    print(f"Wrote {shape_path}")
    print(f"Wrote {entropy_path}")
    print("Done.")


if __name__ == "__main__":
    main()