#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist


def tc_scaling(L, Tc_inf, A, nu):
    return Tc_inf + A * L ** (-1.0 / nu)


def lmin_aic_table(
    tc_df: pd.DataFrame,
    use_aicc: bool = True,
) -> pd.DataFrame:
    required = {
        "lattice_type",
        "representation",
        "variant",
        "L",
        "Tc",
        "Tc_std",
        "N_boot",
    }

    missing = required - set(tc_df.columns)

    if missing:
        raise ValueError(
            f"Input CSV is missing columns: {sorted(missing)}"
        )

    tc_df = (
        tc_df
        .sort_values("L")
        .reset_index(drop=True)
    )

    lattice_types = tc_df["lattice_type"].unique()
    representations = tc_df["representation"].unique()
    variants = tc_df["variant"].unique()

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

    sizes = tc_df["L"].to_numpy()

    if len(sizes) < 5:
        raise ValueError(
            "At least 5 lattice sizes are required."
        )

    records = []
    prev = None

    for L_min in sizes[:-4]:
        df_fit = tc_df[
            tc_df["L"] >= L_min
        ]

        L = df_fit["L"].to_numpy(dtype=float)
        Tc = df_fit["Tc"].to_numpy(dtype=float)
        sigma = df_fit["Tc_std"].to_numpy(dtype=float)

        valid = (
            np.isfinite(L)
            & np.isfinite(Tc)
            & np.isfinite(sigma)
            & (sigma > 0)
        )

        L = L[valid]
        Tc = Tc[valid]
        sigma = sigma[valid]

        n = len(L)

        if n < 4:
            print(
                f"Skipping L_min={L_min}: "
                f"only {n} valid points"
            )
            continue

        popt, pcov = curve_fit(
            tc_scaling,
            L,
            Tc,
            sigma=sigma,
            absolute_sigma=True,
            p0=[Tc.mean(), 1.0, 1.0],
            maxfev=20000,
        )

        Tc_inf, A, nu = popt
        dTc_inf, dA, dnu = np.sqrt(
            np.diag(pcov)
        )

        Tc_fit = tc_scaling(
            L,
            *popt,
        )

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
            1.0 - chi2_dist.cdf(
                chi2,
                dof,
            )
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
            "No successful L_min fits were produced."
        )

    out = pd.DataFrame(records)

    out["delta_AIC"] = (
        out["AIC"]
        - out["AIC"].min()
    )

    return (
        out
        .sort_values("L_min")
        .reset_index(drop=True)
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Perform finite-size-scaling L_min scans "
            "using pseudocritical-temperature CSV files."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help=(
            "CSV produced by "
            "compute_pseudocritical_temperatures.py."
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output CSV for the L_min scan.",
    )

    parser.add_argument(
        "--aic",
        action="store_true",
        help=(
            "Use ordinary AIC instead of AICc. "
            "AICc is used by default."
        ),
    )

    args = parser.parse_args()

    if not args.input.is_file():
        raise FileNotFoundError(
            f"Input file does not exist: {args.input}"
        )

    df = pd.read_csv(args.input)

    result = lmin_aic_table(
        df,
        use_aicc=not args.aic,
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    result.to_csv(
        args.output,
        index=False,
        float_format="%.10g",
    )

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()