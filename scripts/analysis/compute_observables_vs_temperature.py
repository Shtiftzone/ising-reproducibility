#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


TEMP_PATTERN = re.compile(r"^T_([-+]?\d+(?:\.\d+)?)$")


def parse_temperature(dirname: str) -> float:
    match = TEMP_PATTERN.match(dirname)

    if not match:
        raise ValueError(
            f"Unrecognized temperature directory: {dirname}"
        )

    return float(match.group(1))


def read_float_matrix(path: Path) -> np.ndarray:
    """
    Read a comma-separated simulation output file into a 2D float array.
    Empty lines are ignored.
    """
    rows = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            values = [
                float(value.strip())
                for value in line.split(",")
            ]

            rows.append(values)

    if not rows:
        return np.empty((0, 0), dtype=float)

    n_columns = max(len(row) for row in rows)

    rows = [
        row + [np.nan] * (n_columns - len(row))
        for row in rows
    ]

    return np.asarray(rows, dtype=float)


def jackknife_se_of_mean(values: np.ndarray) -> float:
    """
    Jackknife standard error of the sample mean.
    """
    values = np.asarray(values, dtype=float)
    n = values.size

    if n < 2:
        return np.nan

    leave_one_out = (
        values.sum() - values
    ) / (n - 1)

    leave_one_out_mean = leave_one_out.mean()

    variance = (
        (n - 1)
        / n
        * np.sum(
            (leave_one_out - leave_one_out_mean) ** 2
        )
    )

    return float(np.sqrt(variance))


def aggregate_observable(
    values: np.ndarray,
    lattice_type: str,
    size: int,
    temperature: float,
    observable: str,
) -> dict:
    """
    Create one figure-level row containing the mean and jackknife error.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        raise ValueError(
            f"No finite values for {lattice_type}, "
            f"L={size}, T={temperature}, "
            f"observable={observable}"
        )

    return {
        "lattice_type": lattice_type,
        "L": size,
        "T": temperature,
        "observable": observable,
        "mean": float(np.mean(values)),
        "jackknife_se": jackknife_se_of_mean(values),
        "n_samples": int(values.size),
    }


def collect_lattice_observables(
    results_dir: Path,
    lattice_type: str,
    size: int,
) -> dict[str, pd.DataFrame]:
    """
    Compute the observables used in the temperature-dependence plots:

      <|M|>
      <E>
      <|EC_sym|>, where EC_sym = (chi_neg - chi_pos) / 2
      <EC_avg>,   where EC_avg = (chi_neg + chi_pos) / 2
    """
    if not results_dir.is_dir():
        raise FileNotFoundError(
            f"Simulation results directory does not exist: "
            f"{results_dir}"
        )

    rows = {
        "mag_abs": [],
        "energy": [],
        "euler_sym_abs": [],
        "ec_avg": [],
    }

    temperature_dirs = [
        path
        for path in results_dir.iterdir()
        if path.is_dir()
        and TEMP_PATTERN.match(path.name)
    ]

    temperature_dirs.sort(
        key=lambda path: parse_temperature(path.name)
    )

    if not temperature_dirs:
        raise RuntimeError(
            f"No temperature directories found in "
            f"{results_dir}"
        )

    for temperature_dir in temperature_dirs:
        temperature = parse_temperature(
            temperature_dir.name
        )

        size_dir = (
            temperature_dir
            / f"size_{size}"
        )

        if not size_dir.is_dir():
            print(
                f"Skipping {lattice_type}, "
                f"L={size}, "
                f"T={temperature:.5f}: "
                f"missing {size_dir}"
            )
            continue

        me_path = size_dir / "mefile.txt"
        ep_path = size_dir / "epfile.txt"
        en_path = size_dir / "enfile.txt"

        missing = [
            path.name
            for path in (
                me_path,
                ep_path,
                en_path,
            )
            if not path.is_file()
        ]

        if missing:
            print(
                f"Skipping {lattice_type}, "
                f"L={size}, "
                f"T={temperature:.5f}: "
                f"missing {', '.join(missing)}"
            )
            continue

        me = read_float_matrix(me_path)
        chi_pos = read_float_matrix(ep_path)
        chi_neg = read_float_matrix(en_path)

        if me.shape[1] < 2:
            raise ValueError(
                f"Expected at least two columns in "
                f"{me_path}"
            )

        if chi_pos.shape[1] < 1:
            raise ValueError(
                f"Expected at least one column in "
                f"{ep_path}"
            )

        if chi_neg.shape[1] < 1:
            raise ValueError(
                f"Expected at least one column in "
                f"{en_path}"
            )

        magnetization = me[:, 0]
        energy = me[:, 1]

        n_euler = min(
            len(chi_pos),
            len(chi_neg),
        )

        if n_euler == 0:
            raise ValueError(
                f"No Euler samples for "
                f"{lattice_type}, "
                f"L={size}, "
                f"T={temperature}"
            )

        chi_pos_values = chi_pos[:n_euler, 0]
        chi_neg_values = chi_neg[:n_euler, 0]

        ec_sym = (
            chi_neg_values
            - chi_pos_values
        ) / 2.0

        ec_avg = (
            chi_neg_values
            + chi_pos_values
        ) / 2.0

        rows["mag_abs"].append(
            aggregate_observable(
                np.abs(magnetization),
                lattice_type,
                size,
                temperature,
                "abs_magnetization",
            )
        )

        rows["energy"].append(
            aggregate_observable(
                energy,
                lattice_type,
                size,
                temperature,
                "energy",
            )
        )

        rows["euler_sym_abs"].append(
            aggregate_observable(
                np.abs(ec_sym),
                lattice_type,
                size,
                temperature,
                "abs_ec_sym",
            )
        )

        rows["ec_avg"].append(
            aggregate_observable(
                ec_avg,
                lattice_type,
                size,
                temperature,
                "ec_avg",
            )
        )

    result = {}

    for key, values in rows.items():
        dataframe = pd.DataFrame(values)

        if dataframe.empty:
            raise RuntimeError(
                f"No data produced for "
                f"{lattice_type}: {key}"
            )

        result[key] = (
            dataframe
            .sort_values(["L", "T"])
            .reset_index(drop=True)
        )

    return result


def write_lattice_csvs(
    observables: dict[str, pd.DataFrame],
    output_dir: Path,
    prefix: str,
) -> None:
    filenames = {
        "mag_abs":
            f"{prefix}_mag_vs_T.csv",

        "energy":
            f"{prefix}_energy_vs_T.csv",

        "euler_sym_abs":
            f"{prefix}_euler_sym_vs_T.csv",

        "ec_avg":
            f"{prefix}_ec_avg_vs_T.csv",
    }

    for key, dataframe in observables.items():
        path = output_dir / filenames[key]

        dataframe.to_csv(
            path,
            index=False,
            float_format="%.10g",
        )

        print(f"Wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute figure-level magnetization, energy, "
            "and Euler observables from square- and "
            "triangular-lattice Ising simulation outputs."
        )
    )

    parser.add_argument(
        "--square-results",
        type=Path,
        default=Path(
            "results/square_simulations"
        ),
        help=(
            "Square-lattice simulation "
            "results directory."
        ),
    )

    parser.add_argument(
        "--triangular-results",
        type=Path,
        default=Path(
            "results/triangular_simulations"
        ),
        help=(
            "Triangular-lattice simulation "
            "results directory."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "data/figure_csv"
        ),
        help=(
            "Directory for figure-level CSV files."
        ),
    )

    parser.add_argument(
        "--size",
        type=int,
        default=3072,
        help=(
            "Lattice size used for the plots "
            "(default: 3072)."
        ),
    )

    args = parser.parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "Computing square-lattice observables..."
    )

    square = collect_lattice_observables(
        args.square_results,
        lattice_type="square",
        size=args.size,
    )

    print()
    print(
        "Computing triangular-lattice observables..."
    )

    triangular = collect_lattice_observables(
        args.triangular_results,
        lattice_type="triangular",
        size=args.size,
    )

    print()

    write_lattice_csvs(
        square,
        args.output_dir,
        prefix="square",
    )

    write_lattice_csvs(
        triangular,
        args.output_dir,
        prefix="triangular",
    )

    print()
    print("Done.")


if __name__ == "__main__":
    main()
