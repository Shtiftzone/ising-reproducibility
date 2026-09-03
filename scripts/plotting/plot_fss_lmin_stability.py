#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import ScalarFormatter


SQUARE_TC_THEORY = 2.269185314213022
TRIANGULAR_TC_THEORY = 3.64095690651


def jittered_lmin(df_m, df_ec, factor):
    L_m = df_m["L_min"].to_numpy(dtype=float)
    L_ec = df_ec["L_min"].to_numpy(dtype=float)

    return (
        L_m * (1.0 - factor),
        L_ec * (1.0 + factor),
    )


def get_theory(lattice_type):
    if lattice_type == "square":
        return (
            SQUARE_TC_THEORY,
            r"Theory $T_c=\frac{2}{\ln(1+\sqrt{2})}$",
        )

    if lattice_type == "triangular":
        return (
            TRIANGULAR_TC_THEORY,
            r"Theory $T_c=\frac{4}{\ln 3}$",
        )

    raise ValueError(
        f"Unsupported lattice type: {lattice_type}"
    )


def validate_metadata(
    df,
    lattice_type,
    representation,
    variant,
):
    required = {
        "lattice_type",
        "representation",
        "variant",
        "L_min",
        "L_max",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Input CSV is missing columns: {sorted(missing)}"
        )

    lattice_values = (
        df["lattice_type"]
        .dropna()
        .unique()
    )

    representation_values = (
        df["representation"]
        .dropna()
        .unique()
    )

    variant_values = (
        df["variant"]
        .dropna()
        .unique()
    )

    if (
        len(lattice_values) != 1
        or lattice_values[0] != lattice_type
    ):
        raise ValueError(
            f"Input CSV does not correspond to "
            f"lattice_type='{lattice_type}'."
        )

    if (
        len(representation_values) != 1
        or representation_values[0] != representation
    ):
        raise ValueError(
            f"Input CSV does not correspond to "
            f"representation='{representation}'."
        )

    if (
        len(variant_values) != 1
        or variant_values[0] != variant
    ):
        raise ValueError(
            f"Input CSV does not correspond to "
            f"variant='{variant}'."
        )


def plot_nu(
    df_ec,
    df_m,
    output_path,
    title,
    jitter_factor,
):
    L_m, L_ec = jittered_lmin(
        df_m,
        df_ec,
        jitter_factor,
    )

    plt.figure(
        figsize=(9, 5.5)
    )

    plt.errorbar(
        L_m,
        df_m["nu"],
        yerr=df_m["dnu"],
        fmt="o-",
        capsize=4,
        elinewidth=1.2,
        markersize=5,
        linewidth=1.5,
        label="Magnetization",
    )

    plt.errorbar(
        L_ec,
        df_ec["nu"],
        yerr=df_ec["dnu"],
        fmt="s-",
        capsize=4,
        elinewidth=1.2,
        markersize=5,
        linewidth=1.5,
        label=r"$EC_{\mathrm{sym}}$",
    )

    plt.axhline(
        1.0,
        linestyle=":",
        linewidth=1.6,
        label=r"Theory $\nu=1$",
    )

    plt.xlabel(
        r"$L_{\min}$"
    )

    plt.ylabel(
        r"Fitted $\nu$"
    )

    plt.title(
        title
    )

    plt.xscale(
        "log"
    )

    plt.grid(
        True,
        alpha=0.25,
    )

    plt.legend(
        frameon=True,
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(
        f"Wrote {output_path}"
    )


def plot_tc(
    df_ec,
    df_m,
    output_path,
    title,
    tc_theory,
    jitter_factor,
    theory_label,
):
    L_m, L_ec = jittered_lmin(
        df_m,
        df_ec,
        jitter_factor,
    )

    plt.figure(
        figsize=(8, 5)
    )

    plt.errorbar(
        L_m,
        df_m["Tc_inf"],
        yerr=df_m["dTc_inf"],
        fmt="o-",
        capsize=3,
        label="Magnetization",
    )

    plt.errorbar(
        L_ec,
        df_ec["Tc_inf"],
        yerr=df_ec["dTc_inf"],
        fmt="s-",
        capsize=3,
        label=r"$EC_{\mathrm{sym}}$",
    )

    plt.axhline(
        tc_theory,
        linestyle=":",
        linewidth=1.5,
        label=theory_label,
    )

    plt.xlabel(
        r"$L_{\min}$"
    )

    plt.ylabel(
        r"Fitted $T_c$"
    )

    plt.title(
        title
    )

    plt.xscale(
        "log"
    )

    ax = plt.gca()

    formatter = ScalarFormatter(
        useOffset=False
    )

    formatter.set_scientific(
        False
    )

    ax.yaxis.set_major_formatter(
        formatter
    )

    plt.grid(
        True,
        alpha=0.3,
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(
        f"Wrote {output_path}"
    )


def load_pair(
    input_dir,
    suffix,
    lattice_type,
    representation,
):
    ec_path = (
        input_dir
        / f"lmin_ec_lmax_{suffix}.csv"
    )

    m_path = (
        input_dir
        / f"lmin_m_lmax_{suffix}.csv"
    )

    missing = []

    if not ec_path.is_file():
        missing.append(
            ec_path
        )

    if not m_path.is_file():
        missing.append(
            m_path
        )

    if missing:
        raise FileNotFoundError(
            "Missing input files: "
            + ", ".join(
                str(path)
                for path in missing
            )
        )

    df_ec = pd.read_csv(
        ec_path
    )

    df_m = pd.read_csv(
        m_path
    )

    validate_metadata(
        df_ec,
        lattice_type,
        representation,
        "ec",
    )

    validate_metadata(
        df_m,
        lattice_type,
        representation,
        "m",
    )

    return (
        df_ec,
        df_m,
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot stability of fitted Tc and nu "
            "against the minimum lattice size "
            "for one lattice/representation combination."
        )
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
        "--base-input-dir",
        type=Path,
        default=Path(
            "data/table_csv"
        ),
        help=(
            "Base directory containing "
            "<lattice>/<representation>/fixed_lmax."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "results/figures"
        ),
    )

    parser.add_argument(
        "--jitter-factor",
        type=float,
        default=0.01,
    )

    args = parser.parse_args()

    lattice_type = (
        args.lattice_type
    )

    representation = (
        args.representation
    )

    input_dir = (
        args.base_input_dir
        / lattice_type
        / representation
        / "fixed_lmax"
    )

    if not input_dir.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist: "
            f"{input_dir}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    tc_theory, theory_label = get_theory(
        lattice_type
    )

    print(
        f"Lattice type: "
        f"{lattice_type}"
    )

    print(
        f"Representation: "
        f"{representation}"
    )

    print(
        f"Input directory: "
        f"{input_dir}"
    )

    df_ec_max, df_m_max = load_pair(
        input_dir,
        "max",
        lattice_type,
        representation,
    )

    df_ec_second, df_m_second = load_pair(
        input_dir,
        "second",
        lattice_type,
        representation,
    )

    L_max_ec = int(
        df_ec_max["L_max"].iloc[0]
    )

    L_max_m = int(
        df_m_max["L_max"].iloc[0]
    )

    if L_max_ec != L_max_m:
        raise ValueError(
            f"EC and M have different L_max values "
            f"for the max scan: "
            f"{L_max_ec} vs {L_max_m}"
        )

    L_second_ec = int(
        df_ec_second["L_max"].iloc[0]
    )

    L_second_m = int(
        df_m_second["L_max"].iloc[0]
    )

    if L_second_ec != L_second_m:
        raise ValueError(
            f"EC and M have different L_max values "
            f"for the second scan: "
            f"{L_second_ec} vs {L_second_m}"
        )

    L_max = L_max_ec
    L_second = L_second_ec

    lattice_label = (
        "Square"
        if lattice_type == "square"
        else "Triangular"
    )

    representation_label = (
        "spin as cell"
        if representation == "cell"
        else "spin as vertex"
    )

    base_title = (
        f"{lattice_label} lattice, "
        f"{representation_label}"
    )

    prefix = (
        f"{lattice_type}_"
        f"{representation}"
    )

    plot_tc(
        df_ec_max,
        df_m_max,
        args.output_dir
        / f"Tc_vs_Lmin_{prefix}_Lmax_max.png",
        rf"{base_title}: fitted $T_c$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_max}$",
        tc_theory,
        args.jitter_factor,
        theory_label,
    )

    plot_tc(
        df_ec_second,
        df_m_second,
        args.output_dir
        / f"Tc_vs_Lmin_{prefix}_Lmax_second.png",
        rf"{base_title}: fitted $T_c$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_second}$",
        tc_theory,
        args.jitter_factor,
        theory_label,
    )

    plot_nu(
        df_ec_max,
        df_m_max,
        args.output_dir
        / f"nu_vs_Lmin_{prefix}_Lmax_max.png",
        rf"{base_title}: fitted $\nu$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_max}$",
        args.jitter_factor,
    )

    plot_nu(
        df_ec_second,
        df_m_second,
        args.output_dir
        / f"nu_vs_Lmin_{prefix}_Lmax_second.png",
        rf"{base_title}: fitted $\nu$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_second}$",
        args.jitter_factor,
    )

    print()
    print("Done.")


if __name__ == "__main__":
    main()
