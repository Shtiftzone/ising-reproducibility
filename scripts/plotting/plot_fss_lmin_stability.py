#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SQUARE_TC_THEORY = 2.269185314213022
TRIANGULAR_TC_THEORY = 3.64095690651

def jittered_lmin(df_m, df_ec, factor):
    L_m = df_m["L_min"].to_numpy(dtype=float)
    L_ec = df_ec["L_min"].to_numpy(dtype=float)

    return (
        L_m * (1.0 - factor),
        L_ec * (1.0 + factor),
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

    plt.figure(figsize=(9, 5.5))

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

    plt.xlabel(r"$L_{\min}$")
    plt.ylabel(r"Fitted $\nu$")
    plt.title(title)
    plt.xscale("log")

    plt.grid(True, alpha=0.25)
    plt.legend(frameon=True)
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(f"Wrote {output_path}")


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

    plt.figure(figsize=(8, 5))

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

    if tc_theory is not None:
        plt.axhline(
            tc_theory,
            linestyle=":",
            linewidth=1.5,
            label=theory_label,
        )

    plt.xlabel(r"$L_{\min}$")
    plt.ylabel(r"Fitted $T_c$")
    plt.title(title)
    plt.xscale("log")

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(f"Wrote {output_path}")


def load_pair(input_dir, suffix):
    ec_path = input_dir / f"lmin_ec_lmax_{suffix}.csv"
    m_path = input_dir / f"lmin_m_lmax_{suffix}.csv"

    missing = []

    if not ec_path.is_file():
        missing.append(ec_path)

    if not m_path.is_file():
        missing.append(m_path)

    if missing:
        raise FileNotFoundError(
            "Missing input files: "
            + ", ".join(str(path) for path in missing)
        )

    return (
        pd.read_csv(ec_path),
        pd.read_csv(m_path),
    )


def plot_lattice(
    lattice_type,
    input_dir,
    output_dir,
    tc_theory,
    theory_label,
    jitter_factor,
):
    print()
    print(f"Processing lattice: {lattice_type}")

    try:
        df_ec_max, df_m_max = load_pair(
            input_dir,
            "max",
        )

        df_ec_second, df_m_second = load_pair(
            input_dir,
            "second",
        )

    except FileNotFoundError as exc:
        print(f"Skipping {lattice_type}: {exc}")
        return False

    L_max = int(df_ec_max["L_max"].iloc[0])
    L_second = int(
        df_ec_second["L_max"].iloc[0]
    )

    prefix = (
        "PBC"
        if lattice_type == "square"
        else "triangular"
    )

    plot_tc(
        df_ec_max,
        df_m_max,
        output_dir / f"Tc_vs_Lmin_{prefix}_Lmax_max.png",
        rf"Fitted $T_c$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_max}$",
        tc_theory,
        jitter_factor,
        theory_label,
    )

    plot_tc(
        df_ec_second,
        df_m_second,
        output_dir / f"Tc_vs_Lmin_{prefix}_Lmax_second.png",
        rf"Fitted $T_c$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_second}$",
        tc_theory,
        jitter_factor,
        theory_label,
    )

    plot_nu(
        df_ec_max,
        df_m_max,
        output_dir / f"nu_vs_Lmin_{prefix}_Lmax_max.png",
        rf"Fitted $\nu$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_max}$",
        jitter_factor,
    )

    plot_nu(
        df_ec_second,
        df_m_second,
        output_dir / f"nu_vs_Lmin_{prefix}_Lmax_second.png",
        rf"Fitted $\nu$ vs $L_{{\min}}$, "
        rf"$L_{{\max}}={L_second}$",
        jitter_factor,
    )

    return True


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot stability of fitted Tc and nu "
            "against the minimum lattice size."
        )
    )

    parser.add_argument(
        "--base-input-dir",
        type=Path,
        default=Path("data/table_csv"),
        help=(
            "Base directory containing square/fixed_lmax "
            "and triangular/fixed_lmax."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figures"),
    )

    parser.add_argument(
        "--jitter-factor",
        type=float,
        default=0.01,
    )

    args = parser.parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    generated_any = False

    square_input = (
        args.base_input_dir
        / "square"
        / "fixed_lmax"
    )

    generated_any |= plot_lattice(
        lattice_type="square",
        input_dir=square_input,
        output_dir=args.output_dir,
        tc_theory=SQUARE_TC_THEORY,
        theory_label=r"Theory $T_c=\frac{2}{\ln(1+\sqrt{2})}$",
        jitter_factor=args.jitter_factor,
    )

    triangular_input = (
        args.base_input_dir
        / "triangular"
        / "fixed_lmax"
    )

    generated_any |= plot_lattice(
        lattice_type="triangular",
        input_dir=triangular_input,
        output_dir=args.output_dir,
        tc_theory=TRIANGULAR_TC_THEORY,
        theory_label=None,
        jitter_factor=args.jitter_factor,
    )

    print()

    if generated_any:
        print("Done.")
    else:
        print("No plots generated: no complete input datasets were found.")


if __name__ == "__main__":
    main()
