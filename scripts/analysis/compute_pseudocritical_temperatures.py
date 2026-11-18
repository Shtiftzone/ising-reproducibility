
#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SIZES = [
    64, 96, 128, 192, 256, 384,
    512, 768, 1024, 1536, 2048, 3072,
]

DEFAULT_OUTPUT_DIR = Path(
    "/data/home/mpx641/ising-reproducibility/data/analysis_csv"
)

COLUMN_INDEX = {
    "ec": 4,
    "m": 6,
    "d": 8,
}


def extract_pseudocritical_observables(
    sizes,
    input_dir,
    lattice_type,
    representation,
    variant,
):
    """
    Extract bootstrap pseudocritical temperatures and peak
    susceptibilities.

    Tc_b(L)    = argmax_T chi_b(T, L)
    chi_c_b(L) = max_T chi_b(T, L)

    Bootstrap IDs are preserved for paired analysis.
    """

    col_idx = COLUMN_INDEX[variant]

    bootstrap_records = []
    summary_records = []

    for L in sizes:
        file_path = input_dir / f"results_{L}.dat"

        if not file_path.is_file():
            print(f"Skipping L={L}: missing {file_path}")
            continue

        data = np.loadtxt(file_path)

        if data.ndim == 1:
            data = data.reshape(1, -1)

        if data.shape[1] <= col_idx:
            raise ValueError(
                f"{file_path} has {data.shape[1]} columns, "
                f"but {variant} requires column {col_idx}"
            )

        df = pd.DataFrame({
            "bootstrap": data[:, 0].astype(int),
            "T": data[:, 1],
            "chi": data[:, col_idx] / data[:, 1],
        })

        tc_values = []
        chi_c_values = []

        for bootstrap_id, df_b in df.groupby("bootstrap"):
            chi = df_b["chi"].to_numpy()
            temperatures = df_b["T"].to_numpy()

            if chi.size == 0:
                continue

            idx_max = np.argmax(chi)

            tc = float(temperatures[idx_max])
            chi_c = float(chi[idx_max])

            bootstrap_records.append({
                "lattice_type": lattice_type,
                "representation": representation,
                "variant": variant,
                "L": L,
                "bootstrap": int(bootstrap_id),
                "Tc": tc,
                "chi_c": chi_c,
            })

            tc_values.append(tc)
            chi_c_values.append(chi_c)

        tc_values = np.asarray(tc_values, dtype=float)
        chi_c_values = np.asarray(chi_c_values, dtype=float)

        if tc_values.size == 0:
            print(f"Skipping L={L}: no bootstrap maxima")
            continue

        summary_records.append({
            "lattice_type": lattice_type,
            "representation": representation,
            "variant": variant,
            "L": L,
            "Tc": float(tc_values.mean()),
            "Tc_std": (
                float(tc_values.std(ddof=1))
                if tc_values.size > 1
                else np.nan
            ),
            "chi_c": float(chi_c_values.mean()),
            "chi_c_std": (
                float(chi_c_values.std(ddof=1))
                if chi_c_values.size > 1
                else np.nan
            ),
            "N_boot": int(tc_values.size),
        })

    if not summary_records:
        raise RuntimeError(
            f"No pseudocritical observables found in {input_dir}"
        )

    bootstrap_df = (
        pd.DataFrame(bootstrap_records)
        .sort_values(["L", "bootstrap"])
        .reset_index(drop=True)
    )

    summary_df = (
        pd.DataFrame(summary_records)
        .sort_values("L")
        .reset_index(drop=True)
    )

    return bootstrap_df, summary_df


def run_variant(
    sizes,
    input_dir,
    output_dir,
    lattice_type,
    representation,
    variant,
):
    print()
    print(
        f"Extracting: {lattice_type}, "
        f"{representation}, {variant}"
    )

    bootstrap_df, summary_df = extract_pseudocritical_observables(
        sizes=sizes,
        input_dir=input_dir,
        lattice_type=lattice_type,
        representation=representation,
        variant=variant,
    )

    bootstrap_path = (
        output_dir / f"pseudocritical_bootstrap_{variant}.csv"
    )

    summary_path = (
        output_dir / f"pseudocritical_{variant}.csv"
    )

    bootstrap_df.to_csv(
        bootstrap_path,
        index=False,
        float_format="%.10g",
    )

    summary_df.to_csv(
        summary_path,
        index=False,
        float_format="%.10g",
    )

    print(f"Wrote {bootstrap_path}")
    print(f"Wrote {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract pseudocritical temperatures and peak "
            "susceptibilities from multihistogram bootstrap results."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing results_L.dat files.",
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
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
    )

    parser.add_argument(
        "--sizes",
        nargs="+",
        type=int,
        default=DEFAULT_SIZES,
    )

    args = parser.parse_args()

    if not args.input_dir.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist: {args.input_dir}"
        )

    output_dir = (
        args.output_dir
        / args.lattice_type
        / args.representation
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Lattice type: {args.lattice_type}")
    print(f"Representation: {args.representation}")
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")

    for variant in ("ec", "m", "d"):
        run_variant(
            args.sizes,
            args.input_dir,
            output_dir,
            args.lattice_type,
            args.representation,
            variant,
        )

    print()
    print("Done.")


if __name__ == "__main__":
    main()
