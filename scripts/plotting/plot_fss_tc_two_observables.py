#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


DEFAULT_LMIN = {
    "square": {
        "cell": {
            "ec": 64,
            "m": 64,
        },
        "vertex": {
            "ec": 64,
            "m": 64,
        },
    },
    "triangular": {
        "cell": {
            "ec": 128,
            "m": 96,
        },
        "vertex": {
            "ec": 128,
            "m": 64,
        },
    },
}


def tc_scaling(L, Tc_inf, A, nu):
    return Tc_inf + A * L ** (-1.0 / nu)


def fit_tc(df, fix_nu=None):
    required = {"L", "Tc", "Tc_std"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Input dataframe is missing columns: {sorted(missing)}"
        )

    L = df["L"].to_numpy(dtype=float)
    Tc = df["Tc"].to_numpy(dtype=float)
    Tc_std = df["Tc_std"].to_numpy(dtype=float)

    valid = (
        np.isfinite(L)
        & np.isfinite(Tc)
        & np.isfinite(Tc_std)
        & (Tc_std > 0.0)
    )

    L = L[valid]
    Tc = Tc[valid]
    Tc_std = Tc_std[valid]

    if len(L) < 3:
        raise ValueError(
            "Not enough valid points for finite-size scaling fit."
        )

    order = np.argsort(L)
    L = L[order]
    Tc = Tc[order]
    Tc_std = Tc_std[order]

    if fix_nu is not None:

        def model(L, Tc_inf, A):
            return Tc_inf + A * L ** (-1.0 / fix_nu)

        popt, pcov = curve_fit(
            model,
            L,
            Tc,
            sigma=Tc_std,
            absolute_sigma=True,
            p0=[np.mean(Tc), 1.0],
            maxfev=20000,
        )

        Tc_inf, A = popt
        Tc_inf_err, A_err = np.sqrt(np.diag(pcov))

        nu = float(fix_nu)
        nu_err = 0.0

        def fitfun(L_values):
            return model(
                L_values,
                Tc_inf,
                A,
            )

    else:
        popt, pcov = curve_fit(
            tc_scaling,
            L,
            Tc,
            sigma=Tc_std,
            absolute_sigma=True,
            p0=[np.mean(Tc), 1.0, 1.0],
            maxfev=20000,
        )

        Tc_inf, A, nu = popt
        Tc_inf_err, A_err, nu_err = np.sqrt(
            np.diag(pcov)
        )

        def fitfun(L_values):
            return tc_scaling(
                L_values,
                Tc_inf,
                A,
                nu,
            )

    return {
        "L": L,
        "Tc": Tc,
        "Tc_std": Tc_std,
        "Tc_inf": float(Tc_inf),
        "A": float(A),
        "nu": float(nu),
        "Tc_inf_err": float(Tc_inf_err),
        "A_err": float(A_err),
        "nu_err": float(nu_err),
        "fitfun": fitfun,
    }


def get_default_lmin(
    lattice_type,
    representation,
    variant,
):
    try:
        value = DEFAULT_LMIN[
            lattice_type
        ][
            representation
        ][
            variant
        ]
    except KeyError:
        return None

    return value


def validate_metadata(
    df,
    lattice_type,
    representation,
    variant,
):
    if "lattice_type" in df.columns:
        values = df["lattice_type"].dropna().unique()

        if (
            len(values) != 1
            or values[0] != lattice_type
        ):
            raise ValueError(
                f"Input CSV does not correspond to "
                f"lattice_type='{lattice_type}'."
            )

    if "representation" in df.columns:
        values = df["representation"].dropna().unique()

        if (
            len(values) != 1
            or values[0] != representation
        ):
            raise ValueError(
                f"Input CSV does not correspond to "
                f"representation='{representation}'."
            )

    if "variant" in df.columns:
        values = df["variant"].dropna().unique()

        if (
            len(values) != 1
            or values[0] != variant
        ):
            raise ValueError(
                f"Input CSV does not correspond to "
                f"variant='{variant}'."
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot finite-size scaling of pseudocritical temperatures "
            "for Euler-characteristic and magnetization observables."
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
        default=Path("data/analysis_csv"),
        help=(
            "Base directory containing "
            "<lattice>/<representation>/ "
            "pseudocritical-temperature CSV files."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figures"),
    )

    parser.add_argument(
        "--lmin-ec",
        type=int,
        default=None,
        help=(
            "Minimum lattice size used for the EC fit. "
            "If omitted, the stored default for the selected "
            "lattice/representation is used."
        ),
    )

    parser.add_argument(
        "--lmin-m",
        type=int,
        default=None,
        help=(
            "Minimum lattice size used for the magnetization fit. "
            "If omitted, the stored default for the selected "
            "lattice/representation is used."
        ),
    )

    parser.add_argument(
        "--fix-nu-ec",
        type=float,
        default=None,
        help="Fix nu for the EC fit.",
    )

    parser.add_argument(
        "--fix-nu-m",
        type=float,
        default=None,
        help="Fix nu for the magnetization fit.",
    )

    parser.add_argument(
        "--nu-ref-x",
        type=float,
        default=1.0,
        help=(
            "Reference nu used only for the x-axis transformation "
            "x = L^(-1/nu_ref). Default: 1."
        ),
    )

    parser.add_argument(
        "--no-invert-x",
        action="store_true",
        help="Do not invert the x-axis.",
    )

    parser.add_argument(
        "--filename",
        type=str,
        default=None,
        help="Optional output filename.",
    )

    args = parser.parse_args()

    lattice_type = args.lattice_type
    representation = args.representation

    default_lmin_ec = get_default_lmin(
        lattice_type,
        representation,
        "ec",
    )

    default_lmin_m = get_default_lmin(
        lattice_type,
        representation,
        "m",
    )

    lmin_ec = (
        args.lmin_ec
        if args.lmin_ec is not None
        else default_lmin_ec
    )

    lmin_m = (
        args.lmin_m
        if args.lmin_m is not None
        else default_lmin_m
    )

    if lmin_ec is None:
        raise ValueError(
            "No default EC L_min is defined for "
            f"{lattice_type}/{representation}. "
            "Specify --lmin-ec explicitly."
        )

    if lmin_m is None:
        raise ValueError(
            "No default magnetization L_min is defined for "
            f"{lattice_type}/{representation}. "
            "Specify --lmin-m explicitly."
        )

    variant_input_dir = (
        args.input_dir
        / lattice_type
        / representation
    )

    ec_path = (
        variant_input_dir
        / "pseudocritical_ec.csv"
    )

    m_path = (
        variant_input_dir
        / "pseudocritical_m.csv"
    )

    if not ec_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {ec_path}"
        )

    if not m_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {m_path}"
        )

    df_ec = pd.read_csv(ec_path)
    df_m = pd.read_csv(m_path)

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

    df_ec_fit = (
        df_ec[df_ec["L"] >= lmin_ec]
        .copy()
        .sort_values("L")
    )

    df_m_fit = (
        df_m[df_m["L"] >= lmin_m]
        .copy()
        .sort_values("L")
    )

    print(
        f"Lattice type: {lattice_type}"
    )

    print(
        f"Representation: {representation}"
    )

    print(
        f"EC fit range: L >= {lmin_ec}"
    )

    print(
        f"M fit range:  L >= {lmin_m}"
    )

    fit_ec = fit_tc(
        df_ec_fit,
        fix_nu=args.fix_nu_ec,
    )

    fit_m = fit_tc(
        df_m_fit,
        fix_nu=args.fix_nu_m,
    )

    def x_of_L(L):
        return (
            np.asarray(L, dtype=float)
            ** (-1.0 / args.nu_ref_x)
        )

    L_min_plot = min(
        fit_ec["L"].min(),
        fit_m["L"].min(),
    )

    L_max_plot = max(
        fit_ec["L"].max(),
        fit_m["L"].max(),
    )

    L_grid = np.linspace(
        L_min_plot,
        L_max_plot,
        800,
    )

    x_grid = x_of_L(L_grid)

    x_ec = x_of_L(
        fit_ec["L"]
    )

    x_m = x_of_L(
        fit_m["L"]
    )

    Tc_fit_ec = fit_ec["fitfun"](
        L_grid
    )

    Tc_fit_m = fit_m["fitfun"](
        L_grid
    )

    plt.figure(
        figsize=(8.2, 6.0)
    )

    color_ec = "#0072B2"
    color_m = "#D55E00"

    plt.errorbar(
        x_ec,
        fit_ec["Tc"],
        yerr=fit_ec["Tc_std"],
        fmt="o",
        capsize=4,
        color=color_ec,
        label=r"$\mathrm{EC}$",
    )

    plt.errorbar(
        x_m,
        fit_m["Tc"],
        yerr=fit_m["Tc_std"],
        fmt="s",
        capsize=4,
        color=color_m,
        label=r"$M$",
    )

    plt.plot(
        x_grid,
        Tc_fit_ec,
        "-",
        color=color_ec,
        label=r"$\mathrm{EC}$ fit",
    )

    plt.plot(
        x_grid,
        Tc_fit_m,
        "-",
        color=color_m,
        label=r"$M$ fit",
    )

    if args.nu_ref_x == 1.0:
        plt.xlabel(
            r"$1/L$"
        )
    else:
        plt.xlabel(
            rf"$L^{{-1/{args.nu_ref_x:g}}}$"
        )

    plt.ylabel(
        r"$T_c(L)$"
    )

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

    plt.title(
        rf"{lattice_label} Ising, {representation_label}: "
        rf"$T_c(L)$ finite-size scaling"
    )

    plt.grid(True)

    if not args.no_invert_x:
        plt.gca().invert_xaxis()

    plt.gca().spines[
        "top"
    ].set_visible(False)

    plt.gca().spines[
        "right"
    ].set_visible(False)

    plt.legend(
        fontsize=9
    )

    plt.tight_layout()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if args.filename is None:
        filename = (
            f"{lattice_type}_"
            f"{representation}_"
            f"fss_Tc_EC_vs_M.png"
        )
    else:
        filename = args.filename

    output_path = (
        args.output_dir
        / filename
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print()
    print("==== FIT RESULTS (EC) ====")

    print(
        f"T_c(inf) = "
        f"{fit_ec['Tc_inf']:.6f} "
        f"+/- {fit_ec['Tc_inf_err']:.6f}"
    )

    print(
        f"A        = "
        f"{fit_ec['A']:.6f} "
        f"+/- {fit_ec['A_err']:.6f}"
    )

    print(
        f"nu       = "
        f"{fit_ec['nu']:.6f} "
        f"+/- {fit_ec['nu_err']:.6f}"
    )

    print()
    print("==== FIT RESULTS (M) ====")

    print(
        f"T_c(inf) = "
        f"{fit_m['Tc_inf']:.6f} "
        f"+/- {fit_m['Tc_inf_err']:.6f}"
    )

    print(
        f"A        = "
        f"{fit_m['A']:.6f} "
        f"+/- {fit_m['A_err']:.6f}"
    )

    print(
        f"nu       = "
        f"{fit_m['nu']:.6f} "
        f"+/- {fit_m['nu_err']:.6f}"
    )

    print()
    print(
        f"Saved plot as: {output_path}"
    )


if __name__ == "__main__":
    main()