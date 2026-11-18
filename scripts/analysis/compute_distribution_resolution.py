#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis


DEFAULT_SIZES = [
    64, 96, 128, 192, 256, 384,
    512, 768, 1024, 1536, 2048, 3072,
]


def read_first_column(path: Path) -> np.ndarray:
    data = np.loadtxt(path, delimiter=",", ndmin=2)

    return data[:, 0].astype(float)


def load_observables(
    size_dir: Path,
    lattice_type: str,
):
    me_path = size_dir / "mefile.txt"
    ep_path = size_dir / "epfile.txt"
    en_path = size_dir / "enfile.txt"

    for path in (me_path, ep_path, en_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    M = read_first_column(me_path)
    EP = read_first_column(ep_path)
    EN = read_first_column(en_path)

    if not (len(M) == len(EP) == len(EN)):
        raise ValueError(
            f"Mismatched sample counts in {size_dir}: "
            f"M={len(M)}, EP={len(EP)}, EN={len(EN)}"
        )

    if len(M) == 0:
        raise ValueError(
            f"No samples found in {size_dir}"
        )

    EC_sym = (EN - EP) / 2.0

    if lattice_type == "square":
        # M = -2 EC_sym + D
        # D = F_0^{+4} - F_0^{-4}
        D = M + 2.0 * EC_sym

    elif lattice_type == "triangular":
        # M = -EC_sym + D
        # D = 1/2 (F_2^+ - F_2^-)
        D = M + EC_sym

    else:
        raise ValueError(
            f"Unknown lattice type: {lattice_type}"
        )

    return M, EC_sym, D


def entropy_effective_values(
    values: np.ndarray,
):
    values = np.asarray(values)

    unique_values, counts = np.unique(
        values,
        return_counts=True,
    )

    probabilities = (
        counts / counts.sum()
    )

    entropy = -np.sum(
        probabilities
        * np.log(probabilities)
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
            "Compare magnetization, symmetrized "
            "Euler-characteristic, and local "
            "domain-imbalance distributions "
            "at a selected temperature."
        )
    )

    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help=(
            "Raw simulation results directory."
        ),
    )

    parser.add_argument(
        "--lattice-type",
        choices=[
            "square",
            "triangular",
        ],
        required=True,
    )

    parser.add_argument(
        "--representation",
        choices=[
            "cell",
            "vertex",
        ],
        required=True,
    )

    parser.add_argument(
        "--temperature",
        type=float,
        required=True,
        help=(
            "Temperature whose raw simulation "
            "distributions are analyzed."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "data/analysis_csv"
        ),
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
            "Temperature directory does "
            f"not exist: {temperature_dir}"
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
            M, EC_sym, D = (
                load_observables(
                    size_dir,
                    args.lattice_type,
                )
            )

        except FileNotFoundError as exc:
            print(
                f"Skipping L={L}: "
                f"missing {exc}"
            )
            continue

        # --------------------------------------------------
        # Distribution-shape statistics
        # --------------------------------------------------

        shape_records.append({
            "lattice_type":
                args.lattice_type,

            "representation":
                args.representation,

            "temperature":
                args.temperature,

            "L":
                L,

            "n_samples":
                len(M),

            # Magnetization
            "M_mean":
                float(np.mean(M)),

            "M_std":
                float(
                    np.std(
                        M,
                        ddof=1,
                    )
                ),

            "M_skewness":
                float(
                    skew(
                        M,
                        bias=False,
                    )
                ),

            "M_excess_kurtosis":
                float(
                    kurtosis(
                        M,
                        fisher=True,
                        bias=False,
                    )
                ),

            # Symmetrized EC
            "EC_mean":
                float(
                    np.mean(EC_sym)
                ),

            "EC_std":
                float(
                    np.std(
                        EC_sym,
                        ddof=1,
                    )
                ),

            "EC_skewness":
                float(
                    skew(
                        EC_sym,
                        bias=False,
                    )
                ),

            "EC_excess_kurtosis":
                float(
                    kurtosis(
                        EC_sym,
                        fisher=True,
                        bias=False,
                    )
                ),

            # Domain imbalance D
            "D_mean":
                float(
                    np.mean(D)
                ),

            "D_std":
                float(
                    np.std(
                        D,
                        ddof=1,
                    )
                ),

            "D_skewness":
                float(
                    skew(
                        D,
                        bias=False,
                    )
                ),

            "D_excess_kurtosis":
                float(
                    kurtosis(
                        D,
                        fisher=True,
                        bias=False,
                    )
                ),

        })

        # --------------------------------------------------
        # Entropy / effective resolution
        # --------------------------------------------------

        K_M, H_M, Keff_M = (
            entropy_effective_values(M)
        )

        K_EC, H_EC, Keff_EC = (
            entropy_effective_values(
                EC_sym
            )
        )

        K_D, H_D, Keff_D = (
            entropy_effective_values(D)
        )

        entropy_records.append({
            "lattice_type":
                args.lattice_type,

            "representation":
                args.representation,

            "temperature":
                args.temperature,

            "L":
                L,

            "n_samples":
                len(M),

            # Magnetization
            "M_unique":
                K_M,

            "M_entropy":
                H_M,

            "M_Keff":
                Keff_M,

            # Symmetrized EC
            "EC_unique":
                K_EC,

            "EC_entropy":
                H_EC,

            "EC_Keff":
                Keff_EC,

            # Domain imbalance D
            "D_unique":
                K_D,

            "D_entropy":
                H_D,

            "D_Keff":
                Keff_D,

            # Ratios relative to M
            "unique_ratio_EC_over_M": (
                K_EC / K_M
                if K_M > 0
                else np.nan
            ),

            "unique_ratio_D_over_M": (
                K_D / K_M
                if K_M > 0
                else np.nan
            ),

            "Keff_ratio_EC_over_M": (
                Keff_EC / Keff_M
                if Keff_M > 0
                else np.nan
            ),

            "Keff_ratio_D_over_M": (
                Keff_D / Keff_M
                if Keff_M > 0
                else np.nan
            ),

            # Unique-value fractions
            "M_unique_fraction": (
                K_M / len(M)
            ),

            "EC_unique_fraction": (
                K_EC / len(EC_sym)
            ),

            "D_unique_fraction": (
                K_D / len(D)
            ),

            # Effective-resolution fractions
            "M_Keff_fraction": (
                Keff_M / len(M)
            ),

            "EC_Keff_fraction": (
                Keff_EC
                / len(EC_sym)
            ),

            "D_Keff_fraction": (
                Keff_D / len(D)
            ),
        })

    if not shape_records:
        raise RuntimeError(
            "No distribution-shape "
            "results were produced."
        )

    if not entropy_records:
        raise RuntimeError(
            "No entropy-resolution "
            "results were produced."
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
        / "distribution_shape.csv"
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

    print(
        f"Lattice type: "
        f"{args.lattice_type}"
    )

    print(
        f"Representation: "
        f"{args.representation}"
    )

    if args.lattice_type == "square":
        print(
            "D definition: "
            "D = M + 2 EC_sym"
        )

    elif args.lattice_type == "triangular":
        print(
            "D definition: "
            "D = M + EC_sym"
        )

    print(
        f"Wrote {shape_path}"
    )

    print(
        f"Wrote {entropy_path}"
    )

    print("Done.")


if __name__ == "__main__":
    main()