#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


COLOR_EC_SYM = "#0072B2"
COLOR_M = "#D55E00"
COLOR_D = "#009E73"


DEFAULT_LMIN = {
    "square": {
        "cell": {
            "ec_sym": 64,
            "m": 64,
            "d": 64,
        },
        "vertex": {
            "ec_sym": 64,
            "m": 64,
            "d": 64,
        },
    },
    "triangular": {
        "cell": {
            "ec_sym": 96,
            "m": 64,
            "d": 64,
        },
        "vertex": {
            "ec_sym": 96,
            "m": 64,
            "d": 64,
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
        & (L > 0.0)
    )

    L = L[valid]
    Tc = Tc[valid]
    Tc_std = Tc_std[valid]

    n_parameters = 2 if fix_nu is not None else 3
    if len(L) <= n_parameters:
        raise ValueError(
            f"Need at least {n_parameters + 1} valid sizes for "
            f"a {n_parameters}-parameter finite-size scaling fit; "
            f"got {len(L)}."
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
            "for ec_sym, m, and d observables."
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
        "--lmin-ec-sym",
        "--lmin-ec",
        dest="lmin_ec_sym",
        type=int,
        default=None,
        help=(
            "Minimum lattice size used for the ec_sym fit. "
            "If omitted, the stored default for the selected "
            "lattice/representation is used."
        ),
    )

    parser.add_argument(
        "--lmin-m",
        type=int,
        default=None,
        help=(
            "Minimum lattice size used for the m fit. "
            "If omitted, the stored default for the selected "
            "lattice/representation is used."
        ),
    )

    parser.add_argument(
        "--lmin-d",
        type=int,
        default=None,
        help=(
            "Minimum lattice size used for the d fit. "
            "If omitted, the stored default for the selected "
            "lattice/representation is used."
        ),
    )

    parser.add_argument(
        "--fix-nu-ec-sym",
        "--fix-nu-ec",
        dest="fix_nu_ec_sym",
        type=float,
        default=None,
        help="Fix nu for the ec_sym fit.",
    )

    parser.add_argument(
        "--fix-nu-m",
        type=float,
        default=None,
        help="Fix nu for the m fit.",
    )

    parser.add_argument(
        "--fix-nu-d",
        type=float,
        default=None,
        help="Fix nu for the d fit.",
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

    if not np.isfinite(args.nu_ref_x) or args.nu_ref_x <= 0:
        parser.error("--nu-ref-x must be a positive finite number")
    for name in ("fix_nu_ec_sym", "fix_nu_m", "fix_nu_d"):
        value = getattr(args, name)
        if value is not None and (not np.isfinite(value) or value <= 0):
            parser.error(f"--{name.replace('_', '-')} must be positive and finite")

    lattice_type = args.lattice_type
    representation = args.representation

    default_lmin_ec_sym = get_default_lmin(
        lattice_type,
        representation,
        "ec_sym",
    )

    default_lmin_m = get_default_lmin(
        lattice_type,
        representation,
        "m",
    )

    default_lmin_d = get_default_lmin(
        lattice_type,
        representation,
        "d",
    )

    lmin_ec_sym = (
        args.lmin_ec_sym
        if args.lmin_ec_sym is not None
        else default_lmin_ec_sym
    )

    lmin_m = (
        args.lmin_m
        if args.lmin_m is not None
        else default_lmin_m
    )

    lmin_d = (
        args.lmin_d
        if args.lmin_d is not None
        else default_lmin_d
    )

    if lmin_ec_sym is None:
        raise ValueError(
            "No default ec_sym L_min is defined for "
            f"{lattice_type}/{representation}. "
            "Specify --lmin-ec-sym explicitly."
        )

    if lmin_m is None:
        raise ValueError(
            "No default m L_min is defined for "
            f"{lattice_type}/{representation}. "
            "Specify --lmin-m explicitly."
        )

    if lmin_d is None:
        raise ValueError(
            "No default d L_min is defined for "
            f"{lattice_type}/{representation}. "
            "Specify --lmin-d explicitly."
        )

    variant_input_dir = (
        args.input_dir
        / lattice_type
        / representation
    )

    # Read the current pseudocritical CSVs (ec, m, d).
    ec_sym_path = (
        variant_input_dir
        / "pseudocritical_ec.csv"
    )

    m_path = (
        variant_input_dir
        / "pseudocritical_m.csv"
    )

    d_path = (
        variant_input_dir
        / "pseudocritical_d.csv"
    )

    if not ec_sym_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {ec_sym_path}"
        )

    if not m_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {m_path}"
        )

    if not d_path.is_file():
        raise FileNotFoundError(
            f"Missing input file: {d_path}"
        )

    df_ec_sym = pd.read_csv(ec_sym_path)
    df_m = pd.read_csv(m_path)
    df_d = pd.read_csv(d_path)

    # Validate the current observable metadata.
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

    df_ec_sym_fit = (
        df_ec_sym[df_ec_sym["L"] >= lmin_ec_sym]
        .copy()
        .sort_values("L")
    )

    df_m_fit = (
        df_m[df_m["L"] >= lmin_m]
        .copy()
        .sort_values("L")
    )

    df_d_fit = (
        df_d[df_d["L"] >= lmin_d]
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
        f"ec_sym fit range: L >= {lmin_ec_sym}"
    )

    print(
        f"m fit range:      L >= {lmin_m}"
    )

    print(
        f"d fit range:      L >= {lmin_d}"
    )

    fit_ec_sym = fit_tc(
        df_ec_sym_fit,
        fix_nu=args.fix_nu_ec_sym,
    )

    fit_m = fit_tc(
        df_m_fit,
        fix_nu=args.fix_nu_m,
    )

    fit_d = fit_tc(
        df_d_fit,
        fix_nu=args.fix_nu_d,
    )

    def x_of_L(L):
        return (
            np.asarray(L, dtype=float)
            ** (-1.0 / args.nu_ref_x)
        )

    L_min_plot = min(
        fit_ec_sym["L"].min(),
        fit_m["L"].min(),
        fit_d["L"].min(),
    )

    L_max_plot = max(
        fit_ec_sym["L"].max(),
        fit_m["L"].max(),
        fit_d["L"].max(),
    )

    L_grid = np.linspace(
        L_min_plot,
        L_max_plot,
        800,
    )

    x_grid = x_of_L(L_grid)

    x_ec_sym = x_of_L(
        fit_ec_sym["L"]
    )

    x_m = x_of_L(
        fit_m["L"]
    )

    x_d = x_of_L(
        fit_d["L"]
    )

    Tc_fit_ec_sym = fit_ec_sym["fitfun"](
        L_grid
    )

    Tc_fit_m = fit_m["fitfun"](
        L_grid
    )

    Tc_fit_d = fit_d["fitfun"](
        L_grid
    )

    plt.figure(
        figsize=(8.2, 6.0)
    )

    plt.errorbar(
        x_ec_sym,
        fit_ec_sym["Tc"],
        yerr=fit_ec_sym["Tc_std"],
        fmt="o",
        capsize=4,
        color=COLOR_EC_SYM,
        label=r"$\mathrm{ec}_{\mathrm{sym}}$",
    )

    plt.errorbar(
        x_m,
        fit_m["Tc"],
        yerr=fit_m["Tc_std"],
        fmt="s",
        capsize=4,
        color=COLOR_M,
        label=r"$m$",
    )

    plt.errorbar(
        x_d,
        fit_d["Tc"],
        yerr=fit_d["Tc_std"],
        fmt="^",
        capsize=4,
        color=COLOR_D,
        label=r"$d$",
    )

    plt.plot(
        x_grid,
        Tc_fit_ec_sym,
        "-",
        color=COLOR_EC_SYM,
        label=r"$\mathrm{ec}_{\mathrm{sym}}$ fit",
    )

    plt.plot(
        x_grid,
        Tc_fit_m,
        "-",
        color=COLOR_M,
        label=r"$m$ fit",
    )

    plt.plot(
        x_grid,
        Tc_fit_d,
        "-",
        color=COLOR_D,
        label=r"$d$ fit",
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

    variant_output_dir = (
        args.output_dir
        / lattice_type
        / representation
    )
    variant_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if args.filename is None:
        filename = (
            f"{lattice_type}_"
            f"{representation}_"
            f"fss_Tc_ec_sym_vs_m_vs_d.png"
        )
    else:
        filename = args.filename

    output_path = (
        variant_output_dir
        / filename
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print()
    print("==== FIT RESULTS (ec_sym) ====")

    print(
        f"T_c(inf) = "
        f"{fit_ec_sym['Tc_inf']:.6f} "
        f"+/- {fit_ec_sym['Tc_inf_err']:.6f}"
    )

    print(
        f"A        = "
        f"{fit_ec_sym['A']:.6f} "
        f"+/- {fit_ec_sym['A_err']:.6f}"
    )

    print(
        f"nu       = "
        f"{fit_ec_sym['nu']:.6f} "
        f"+/- {fit_ec_sym['nu_err']:.6f}"
    )

    print()
    print("==== FIT RESULTS (m) ====")

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
    print("==== FIT RESULTS (d) ====")

    print(
        f"T_c(inf) = "
        f"{fit_d['Tc_inf']:.6f} "
        f"+/- {fit_d['Tc_inf_err']:.6f}"
    )

    print(
        f"A        = "
        f"{fit_d['A']:.6f} "
        f"+/- {fit_d['A_err']:.6f}"
    )

    print(
        f"nu       = "
        f"{fit_d['nu']:.6f} "
        f"+/- {fit_d['nu_err']:.6f}"
    )

    print()
    print(
        f"Saved plot as: {output_path}"
    )


if __name__ == "__main__":
    main()