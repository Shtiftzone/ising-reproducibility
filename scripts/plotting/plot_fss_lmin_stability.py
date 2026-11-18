#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import ScalarFormatter


SQUARE_TC_THEORY = 2.269185314213022
TRIANGULAR_TC_THEORY = 3.64095690651

plt.rcParams.update({
    "font.size": 20,
    "axes.titlesize": 24,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20,
    "legend.title_fontsize": 20,
    "figure.titlesize": 24,
    "figure.labelsize": 24,
})

COLOR_EC_SYM = "#0072B2"
COLOR_M = "#D55E00"
COLOR_D = "#009E73"

# Identical canvas size and axes margins for all four panels.
FIGURE_SIZE = (9, 5.5)
PLOT_MARGINS = dict(left=0.17, right=0.96, bottom=0.17, top=0.85)


def jittered_lmin(df_m, df_d, df_ec_sym, factor):
    L_m = df_m["L_min"].to_numpy(dtype=float)
    L_d = df_d["L_min"].to_numpy(dtype=float)
    L_ec_sym = df_ec_sym["L_min"].to_numpy(dtype=float)

    return (
        L_m * (1.0 - factor),
        L_d,
        L_ec_sym * (1.0 + factor),
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

    required_fit = {"Tc_inf", "dTc_inf", "nu", "dnu"}
    missing_fit = required_fit - set(df.columns)
    if missing_fit:
        raise ValueError(
            f"Input CSV is missing fit columns: {sorted(missing_fit)}"
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
    df_ec_sym,
    df_m,
    df_d,
    output_path,
    title,
    jitter_factor,
):
    L_m, L_d, L_ec_sym = jittered_lmin(
        df_m,
        df_d,
        df_ec_sym,
        jitter_factor,
    )

    plt.figure(figsize=FIGURE_SIZE)

    plt.errorbar(
        L_m,
        df_m["nu"],
        yerr=df_m["dnu"],
        fmt="o-",
        capsize=4,
        elinewidth=1.2,
        markersize=5,
        linewidth=1.5,
        color=COLOR_M,
        label=r"$m$",
    )

    plt.errorbar(
        L_d,
        df_d["nu"],
        yerr=df_d["dnu"],
        fmt="^-",
        capsize=4,
        elinewidth=1.2,
        markersize=5,
        linewidth=1.5,
        color=COLOR_D,
        label=r"$d$",
    )

    plt.errorbar(
        L_ec_sym,
        df_ec_sym["nu"],
        yerr=df_ec_sym["dnu"],
        fmt="s-",
        capsize=4,
        elinewidth=1.2,
        markersize=5,
        linewidth=1.5,
        color=COLOR_EC_SYM,
        label=r"$\mathrm{ec}_{\mathrm{sym}}$",
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

    plt.subplots_adjust(**PLOT_MARGINS)

    plt.savefig(
        output_path,
        dpi=300,
    )

    plt.close()

    print(
        f"Wrote {output_path}"
    )


def plot_tc(
    df_ec_sym,
    df_m,
    df_d,
    output_path,
    title,
    tc_theory,
    jitter_factor,
    theory_label,
):
    L_m, L_d, L_ec_sym = jittered_lmin(
        df_m,
        df_d,
        df_ec_sym,
        jitter_factor,
    )

    plt.figure(figsize=FIGURE_SIZE)

    plt.errorbar(
        L_m,
        df_m["Tc_inf"],
        yerr=df_m["dTc_inf"],
        fmt="o-",
        capsize=3,
        color=COLOR_M,
        label=r"$m$",
    )

    plt.errorbar(
        L_d,
        df_d["Tc_inf"],
        yerr=df_d["dTc_inf"],
        fmt="^-",
        capsize=3,
        color=COLOR_D,
        label=r"$d$",
    )

    plt.errorbar(
        L_ec_sym,
        df_ec_sym["Tc_inf"],
        yerr=df_ec_sym["dTc_inf"],
        fmt="s-",
        capsize=3,
        color=COLOR_EC_SYM,
        label=r"$\mathrm{ec}_{\mathrm{sym}}$",
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

    plt.subplots_adjust(**PLOT_MARGINS)

    plt.savefig(
        output_path,
        dpi=300,
    )

    plt.close()

    print(
        f"Wrote {output_path}"
    )


def load_triplet(
    input_dir,
    suffix,
    lattice_type,
    representation,
):
    # Keep existing filenames used by the analysis pipeline.
    ec_sym_path = (
        input_dir
        / f"lmin_ec_lmax_{suffix}.csv"
    )

    m_path = (
        input_dir
        / f"lmin_m_lmax_{suffix}.csv"
    )

    d_path = (
        input_dir
        / f"lmin_d_lmax_{suffix}.csv"
    )

    missing = []

    for path in (ec_sym_path, m_path, d_path):
        if not path.is_file():
            missing.append(path)

    if missing:
        raise FileNotFoundError(
            "Missing input files: "
            + ", ".join(
                str(path)
                for path in missing
            )
        )

    df_ec_sym = pd.read_csv(ec_sym_path)
    df_m = pd.read_csv(m_path)
    df_d = pd.read_csv(d_path)

    # Keep existing metadata identifiers used by the pipeline.
    validate_metadata(
        df_ec_sym,
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

    validate_metadata(
        df_d,
        lattice_type,
        representation,
        "d",
    )

    return (
        df_ec_sym,
        df_m,
        df_d,
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot stability of fitted Tc and nu "
            "against the minimum lattice size "
            "for ec_sym, m, and d observables."
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

    variant_output_dir = args.output_dir / lattice_type / representation
    variant_output_dir.mkdir(parents=True, exist_ok=True)

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

    df_ec_sym_max, df_m_max, df_d_max = load_triplet(
        input_dir,
        "max",
        lattice_type,
        representation,
    )

    df_ec_sym_second, df_m_second, df_d_second = load_triplet(
        input_dir,
        "second",
        lattice_type,
        representation,
    )

    L_max_values = {
        int(df_ec_sym_max["L_max"].iloc[0]),
        int(df_m_max["L_max"].iloc[0]),
        int(df_d_max["L_max"].iloc[0]),
    }

    if len(L_max_values) != 1:
        raise ValueError(
            "ec_sym, m, and d have different L_max values "
            f"for the max scan: {sorted(L_max_values)}"
        )

    L_second_values = {
        int(df_ec_sym_second["L_max"].iloc[0]),
        int(df_m_second["L_max"].iloc[0]),
        int(df_d_second["L_max"].iloc[0]),
    }

    if len(L_second_values) != 1:
        raise ValueError(
            "ec_sym, m, and d have different L_max values "
            f"for the second scan: {sorted(L_second_values)}"
        )

    L_max = next(iter(L_max_values))
    L_second = next(iter(L_second_values))

    lattice_label = (
        "Square"
        if lattice_type == "square"
        else "Triangular"
    )

    base_title = (
        f"{lattice_label} / {representation}"
    )

    prefix = (
        f"{lattice_type}_"
        f"{representation}"
    )

    plot_tc(
        df_ec_sym_max,
        df_m_max,
        df_d_max,
        variant_output_dir
        / f"Tc_vs_Lmin_{prefix}_Lmax_max.png",
        rf"{base_title}: $T_c$ vs $L_{{\min}}$ ($L_{{\max}}={L_max}$)",
        tc_theory,
        args.jitter_factor,
        theory_label,
    )

    plot_tc(
        df_ec_sym_second,
        df_m_second,
        df_d_second,
        variant_output_dir
        / f"Tc_vs_Lmin_{prefix}_Lmax_second.png",
        rf"{base_title}: $T_c$ vs $L_{{\min}}$ ($L_{{\max}}={L_second}$)",
        tc_theory,
        args.jitter_factor,
        theory_label,
    )

    plot_nu(
        df_ec_sym_max,
        df_m_max,
        df_d_max,
        variant_output_dir
        / f"nu_vs_Lmin_{prefix}_Lmax_max.png",
        rf"{base_title}: $\nu$ vs $L_{{\min}}$ ($L_{{\max}}={L_max}$)",
        args.jitter_factor,
    )

    plot_nu(
        df_ec_sym_second,
        df_m_second,
        df_d_second,
        variant_output_dir
        / f"nu_vs_Lmin_{prefix}_Lmax_second.png",
        rf"{base_title}: $\nu$ vs $L_{{\min}}$ ($L_{{\max}}={L_second}$)",
        args.jitter_factor,
    )

    print()
    print("Done.")


if __name__ == "__main__":
    main()