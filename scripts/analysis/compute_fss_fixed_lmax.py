#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist


def tc_scaling(L, Tc_inf, A, nu):
    return Tc_inf + A * L ** (-1.0 / nu)


def fixed_lmax_scan(
    df,
    L_max,
    min_points=6,
    use_aicc=True,
):
    required = {
        "lattice_type",
        "representation",
        "variant",
        "L",
        "Tc",
        "Tc_std",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Input CSV is missing columns: {sorted(missing)}"
        )

    df = (
        df[df["L"] <= L_max]
        .copy()
        .sort_values("L")
        .reset_index(drop=True)
    )

    if df.empty:
        raise ValueError(
            f"No data available for L <= {L_max}"
        )

    lattice_types = df["lattice_type"].unique()
    representations = df["representation"].unique()
    variants = df["variant"].unique()

    if len(lattice_types) != 1:
        raise ValueError(
            "Input CSV must contain exactly one lattice type."
        )

    if len(representations) != 1:
        raise ValueError(
            "Input CSV must contain exactly one representation."
        )

    if len(variants) != 1:
        raise ValueError(
            "Input CSV must contain exactly one observable variant."
        )

    lattice_type = lattice_types[0]
    representation = representations[0]
    variant = variants[0]

    records = []
    prev = None

    for L_min in df["L"].to_numpy():
        df_fit = df[
            (df["L"] >= L_min)
            & (df["L"] <= L_max)
        ].copy()

        L = df_fit["L"].to_numpy(dtype=float)
        Tc = df_fit["Tc"].to_numpy(dtype=float)
        sigma = df_fit["Tc_std"].to_numpy(dtype=float)

        valid = (
            np.isfinite(L)
            & np.isfinite(Tc)
            & np.isfinite(sigma)
            & (sigma > 0.0)
        )

        L = L[valid]
        Tc = Tc[valid]
        sigma = sigma[valid]

        n = len(L)

        if n < min_points:
            continue

        try:
            popt, pcov = curve_fit(
                tc_scaling,
                L,
                Tc,
                sigma=sigma,
                absolute_sigma=True,
                p0=[Tc.mean(), 1.0, 1.0],
                maxfev=20000,
            )
        except Exception as exc:
            print(
                f"Fit failed for "
                f"{lattice_type}/{representation}, "
                f"variant={variant}, "
                f"L_min={L_min}, L_max={L_max}: {exc}"
            )
            continue

        Tc_inf, A, nu = popt
        dTc_inf, dA, dnu = np.sqrt(np.diag(pcov))

        Tc_fit = tc_scaling(L, *popt)

        chi2 = np.sum(
            ((Tc - Tc_fit) / sigma) ** 2
        )

        k = 3
        dof = n - k

        chi2_red = (
            chi2 / dof
            if dof > 0
            else np.nan
        )

        p_value = (
            1.0 - chi2_dist.cdf(chi2, dof)
            if dof > 0
            else np.nan
        )

        aic = chi2 + 2 * k

        if use_aicc and n > k + 1:
            aic += (
                2 * k * (k + 1)
                / (n - k - 1)
            )

        rec = {
            "lattice_type": lattice_type,
            "representation": representation,
            "variant": variant,
            "L_min": int(L_min),
            "L_max": int(L_max),
            "n_points": int(n),
            "Tc_inf": float(Tc_inf),
            "dTc_inf": float(dTc_inf),
            "A": float(A),
            "dA": float(dA),
            "nu": float(nu),
            "dnu": float(dnu),
            "chi2": float(chi2),
            "dof": int(dof),
            "chi2_red": float(chi2_red),
            "p_value": float(p_value),
            "AIC": float(aic),
        }

        if prev is None:
            rec["dnu_sigma"] = np.nan
            rec["dTc_inf_sigma"] = np.nan
        else:
            rec["dnu_sigma"] = (
                (nu - prev["nu"])
                / np.sqrt(
                    dnu**2
                    + prev["dnu"]**2
                )
            )

            rec["dTc_inf_sigma"] = (
                (Tc_inf - prev["Tc_inf"])
                / np.sqrt(
                    dTc_inf**2
                    + prev["dTc_inf"]**2
                )
            )

        records.append(rec)
        prev = rec

    if not records:
        raise RuntimeError(
            f"No successful fits for L_max={L_max}"
        )

    out = pd.DataFrame(records)

    out["delta_AIC"] = (
        out["AIC"] - out["AIC"].min()
    )

    return (
        out
        .sort_values("L_min")
        .reset_index(drop=True)
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Perform L_min finite-size-scaling scans "
            "for fixed values of L_max."
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
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/table_csv"),
    )

    parser.add_argument(
        "--min-points",
        type=int,
        default=6,
    )

    parser.add_argument(
        "--aic",
        action="store_true",
        help="Use AIC instead of AICc.",
    )

    args = parser.parse_args()

    input_dir = (
        args.input_dir
        / args.lattice_type
        / args.representation
    )

    output_dir = (
        args.output_dir
        / args.lattice_type
        / args.representation
        / "fixed_lmax"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    ec_path = input_dir / "pseudocritical_ec.csv"
    m_path = input_dir / "pseudocritical_m.csv"

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

    common_sizes = sorted(
        set(df_ec["L"].astype(int))
        & set(df_m["L"].astype(int))
    )

    if len(common_sizes) < 2:
        raise RuntimeError(
            "At least two common lattice sizes are required."
        )

    L_max = common_sizes[-1]
    L_second = common_sizes[-2]

    print(f"Lattice type: {args.lattice_type}")
    print(f"Representation: {args.representation}")
    print(f"L_max: {L_max}")
    print(f"Second L_max: {L_second}")
    print(
        "Information criterion: "
        + ("AIC" if args.aic else "AICc")
    )

    for label, current_Lmax in [
        ("max", L_max),
        ("second", L_second),
    ]:
        for variant, df in [
            ("ec", df_ec),
            ("m", df_m),
        ]:
            print()
            print(
                f"Running {variant} scan "
                f"with L_max={current_Lmax}"
            )

            result = fixed_lmax_scan(
                df=df,
                L_max=current_Lmax,
                min_points=args.min_points,
                use_aicc=not args.aic,
            )

            output_path = (
                output_dir
                / f"lmin_{variant}_lmax_{label}.csv"
            )

            result.to_csv(
                output_path,
                index=False,
                float_format="%.10g",
            )

            print(f"Wrote {output_path}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()