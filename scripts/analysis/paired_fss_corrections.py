#!/usr/bin/env python3
"""Paired-bootstrap comparison of subleading FSS coefficients.

Model: Tc(L) = Tc_inf + a/L + b/L**2, with nu fixed to 1.
Within each lattice size, EC, M and D are paired by bootstrap ID.
Different lattice sizes are sampled independently.

Outputs include the complete three-parameter weighted FSS fit table,
bootstrap distributions of parameters, and paired comparisons of b.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist


VARIANTS = ("ec", "m", "d")
PARAMETERS = ("Tc_inf", "a", "b")
REPO_ROOT = Path(__file__).resolve().parents[2]


def read_inputs(input_dir, lmin, lmax):
    boots, summaries = {}, {}
    for variant in VARIANTS:
        bpath = input_dir / f"pseudocritical_bootstrap_{variant}.csv"
        spath = input_dir / f"pseudocritical_{variant}.csv"
        if not bpath.is_file() or not spath.is_file():
            missing = bpath if not bpath.is_file() else spath
            raise FileNotFoundError(f"Missing {missing}")

        b = pd.read_csv(bpath)
        s = pd.read_csv(spath)
        if not {"L", "bootstrap", "Tc"}.issubset(b.columns):
            raise ValueError(f"Missing L/bootstrap/Tc columns in {bpath}")
        if not {"L", "Tc", "Tc_std"}.issubset(s.columns):
            raise ValueError(f"Missing L/Tc/Tc_std columns in {spath}")

        b = b[(b.L >= lmin) & (b.L <= lmax)].copy()
        s = s[(s.L >= lmin) & (s.L <= lmax)].copy()
        if b.duplicated(["L", "bootstrap"]).any() or s.duplicated("L").any():
            raise ValueError(f"Duplicate (L, bootstrap) or L for {variant}")
        boots[variant] = b
        summaries[variant] = s.set_index("L")

    common_sizes = sorted(
        set.intersection(
            *(set(s.index) for s in summaries.values()),
            *(set(b.L) for b in boots.values()),
        )
    )
    if len(common_sizes) < 4:
        raise ValueError(f"Need at least 4 common lattice sizes, got {common_sizes}")

    paired = {}
    for L in common_sizes:
        frames = []
        for variant in VARIANTS:
            sub = boots[variant].loc[
                boots[variant].L == L, ["bootstrap", "Tc"]
            ]
            frames.append(sub.rename(columns={"Tc": variant}).set_index("bootstrap"))

        joined = pd.concat(frames, axis=1, join="inner").sort_index()
        joined = joined.replace([np.inf, -np.inf], np.nan).dropna()
        if len(joined) < 2:
            raise ValueError(f"L={L}: less than 2 common paired bootstrap IDs")
        counts = [int(sum(boots[v].L == L)) for v in VARIANTS]
        if any(n != len(joined) for n in counts):
            raise ValueError(
                f"L={L}: incomplete pairing: {counts}, common={len(joined)}"
            )
        paired[int(L)] = joined.loc[:, list(VARIANTS)].to_numpy(dtype=float)
        print(f"L={L}: {len(joined)} paired bootstrap samples")

    L = np.array(common_sizes, dtype=float)
    means = np.column_stack([
        summaries[v].loc[common_sizes, "Tc"].to_numpy(float)
        for v in VARIANTS
    ])
    sigma = np.column_stack([
        summaries[v].loc[common_sizes, "Tc_std"].to_numpy(float)
        for v in VARIANTS
    ])
    if (not np.isfinite(means).all() or not np.isfinite(sigma).all()
            or (sigma <= 0).any()):
        raise ValueError("Nonfinite Tc or invalid Tc_std in summary CSV")

    for j, variant in enumerate(VARIANTS):
        reconstructed = np.array([
            paired[int(l)][:, j].mean() for l in L
        ])
        if not np.allclose(reconstructed, means[:, j], atol=1e-8, rtol=0):
            raise ValueError(
                f"Summary Tc means disagree with paired bootstrap values for {variant}"
            )
    return L, paired, means, sigma


def compute_fits(L, paired, means, sigma, repetitions, seed):
    # Weighted linear least squares for the model with nu fixed to 1.
    X = np.column_stack((np.ones_like(L), 1.0 / L, 1.0 / L**2))
    fit_operators = []
    wls_std = np.empty((3, len(VARIANTS)), dtype=float)

    for j in range(len(VARIANTS)):
        weighted_X = X / sigma[:, j, None]
        if np.linalg.matrix_rank(weighted_X) != 3:
            raise ValueError("Fit matrix is rank deficient")

        weighted_pinv = np.linalg.pinv(weighted_X)
        fit_operators.append(weighted_pinv / sigma[:, j][None, :])

        # Absolute-sigma WLS covariance. Kept separately from bootstrap SEs.
        covariance = weighted_pinv @ weighted_pinv.T
        wls_std[:, j] = np.sqrt(np.diag(covariance))

    # Rows correspond to (Tc_inf, a, b); columns to (ec, m, d).
    baseline = np.column_stack([
        fit_operators[j] @ means[:, j] for j in range(len(VARIANTS))
    ])

    rng = np.random.default_rng(seed)
    draws = np.empty((repetitions, len(L), len(VARIANTS)), dtype=float)
    for i, l in enumerate(L.astype(int)):
        available = paired[l]
        idx = rng.integers(0, len(available), size=repetitions)
        draws[:, i, :] = available[idx, :]

    # Columns grouped by observable: (Tc_inf_ec, a_ec, b_ec, ...).
    fitted = np.column_stack([
        draws[:, :, j] @ fit_operators[j].T
        for j in range(len(VARIANTS))
    ])

    dof = len(L) - 3
    fitted_means = X @ baseline
    residuals = (means - fitted_means) / sigma
    chi2_values = np.sum(residuals**2, axis=0)
    p_values = chi2_dist.sf(chi2_values, dof)

    fit_quality = {
        "chi2": chi2_values,
        "dof": dof,
        "chi2_red": chi2_values / dof,
        "p_value": p_values,
        "wls_std": wls_std,
    }
    return baseline, fitted, fit_quality


def save_outputs(output_dir, L, baseline, fitted, fit_quality,
                 lattice_type, representation, ci):
    output_dir.mkdir(parents=True, exist_ok=True)
    n = len(fitted)
    records = {"replicate": np.arange(n, dtype=int)}

    for j, v in enumerate(VARIANTS):
        for i, p in enumerate(PARAMETERS):
            records[f"{p}_{v}"] = fitted[:, j * len(PARAMETERS) + i]

    b_ec = records["b_ec"]
    b_m = records["b_m"]
    b_d = records["b_d"]
    records["delta_b_ec_m"] = b_ec - b_m
    records["delta_b_ec_d"] = b_ec - b_d
    records["delta_abs_b_ec_m"] = np.abs(b_ec) - np.abs(b_m)
    records["delta_abs_b_ec_d"] = np.abs(b_ec) - np.abs(b_d)

    pd.DataFrame(records).to_csv(
        output_dir / "paired_fss_bootstrap.csv", index=False,
        float_format="%.12g"
    )

    # Complete table of central FSS fits and uncertainties for the paper.
    fit_rows = []
    for j, v in enumerate(VARIANTS):
        row = {
            "lattice_type": lattice_type,
            "representation": representation,
            "variant": v,
            "L_min": int(L.min()),
            "L_max": int(L.max()),
            "n_sizes": len(L),
        }
        for i, p in enumerate(PARAMETERS):
            samples = records[f"{p}_{v}"]
            row[p] = float(baseline[i, j])
            row[f"{p}_bootstrap_std"] = float(samples.std(ddof=1))
            row[f"{p}_wls_std"] = float(fit_quality["wls_std"][i, j])
            row[f"{p}_ci_low"] = float(np.percentile(samples, (100 - ci) / 2))
            row[f"{p}_ci_high"] = float(np.percentile(samples, 100 - (100 - ci) / 2))
        row.update({
            "chi2": float(fit_quality["chi2"][j]),
            "dof": int(fit_quality["dof"]),
            "chi2_red": float(fit_quality["chi2_red"][j]),
            "p_value": float(fit_quality["p_value"][j]),
        })
        fit_rows.append(row)

    fit_table = pd.DataFrame(fit_rows)
    fit_table.to_csv(
        output_dir / "fss_subleading_fit_table.csv", index=False,
        float_format="%.12g"
    )

    # Preserve the old b-comparison summary format and output filename.
    summary_rows = []
    for v in VARIANTS:
        arr = records[f"b_{v}"]
        j = VARIANTS.index(v)
        summary_rows.append({
            "quantity": f"b_{v}",
            "fit_to_mean": baseline[2, j],
            "bootstrap_mean": arr.mean(),
            "bootstrap_std": arr.std(ddof=1),
            "ci_low": np.percentile(arr, (100 - ci) / 2),
            "ci_high": np.percentile(arr, 100 - (100 - ci) / 2),
            "fraction_above_zero": (arr > 0).mean(),
        })

    for label in (
        "delta_b_ec_m", "delta_b_ec_d",
        "delta_abs_b_ec_m", "delta_abs_b_ec_d",
    ):
        arr = records[label]
        if label == "delta_b_ec_m":
            fitted_center = baseline[2, 0] - baseline[2, 1]
        elif label == "delta_b_ec_d":
            fitted_center = baseline[2, 0] - baseline[2, 2]
        elif label == "delta_abs_b_ec_m":
            fitted_center = abs(baseline[2, 0]) - abs(baseline[2, 1])
        else:
            fitted_center = abs(baseline[2, 0]) - abs(baseline[2, 2])

        summary_rows.append({
            "quantity": label,
            "fit_to_mean": fitted_center,
            "bootstrap_mean": arr.mean(),
            "bootstrap_std": arr.std(ddof=1),
            "ci_low": np.percentile(arr, (100 - ci) / 2),
            "ci_high": np.percentile(arr, 100 - (100 - ci) / 2),
            "fraction_above_zero": (arr > 0).mean(),
        })

    summary = pd.DataFrame(summary_rows)
    summary.insert(0, "lattice_type", lattice_type)
    summary.insert(1, "representation", representation)
    summary.insert(2, "L_min", int(L.min()))
    summary.insert(3, "L_max", int(L.max()))
    summary.insert(4, "n_sizes", len(L))
    summary.to_csv(
        output_dir / "paired_fss_summary.csv", index=False,
        float_format="%.12g"
    )

    for label, xlabel in [
        ("delta_b_ec_m", r"$b_{EC_{sym}}-b_M$"),
        ("delta_b_ec_d", r"$b_{EC_{sym}}-b_D$"),
        ("delta_abs_b_ec_m", r"$|b_{EC_{sym}}|-|b_M|$"),
        ("delta_abs_b_ec_d", r"$|b_{EC_{sym}}|-|b_D|$"),
    ]:
        arr = records[label]
        low, high = np.percentile(arr, [(100 - ci) / 2, 100 - (100 - ci) / 2])
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(arr, bins=50, edgecolor="white", linewidth=0.4)
        ax.axvline(0, color="black", linestyle="--", label="Zero")
        ax.axvline(low, color="tab:red", linestyle=":", label=f"{ci:g}% percentile CI")
        ax.axvline(high, color="tab:red", linestyle=":")
        ax.axvline(arr.mean(), color="tab:green", label="Bootstrap mean")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Number of resampled fits")
        ax.set_title(
            f"{lattice_type.capitalize()} / {representation}   "
            f"L={int(L.min())}...{int(L.max())}"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / f"{label}.png", dpi=200)
        plt.close(fig)

    print("\nSubleading FSS fit table:")
    cols = [
        "variant", "Tc_inf", "Tc_inf_bootstrap_std", "a",
        "a_bootstrap_std", "b", "b_bootstrap_std", "chi2_red", "p_value"
    ]
    print(fit_table[cols].to_string(index=False, float_format=lambda x: f"{x:.8g}"))
    print("\nPaired b comparisons:")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print(f"\nSaved CSVs and plots to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice-type", required=True, choices=["square", "triangular"])
    parser.add_argument("--representation", default="vertex", choices=["vertex", "cell"])
    parser.add_argument("--input-base", type=Path, default=REPO_ROOT / "data/analysis_csv")
    parser.add_argument(
        "--output-base", type=Path,
        default=REPO_ROOT / "data/analysis_csv/fss_corrections"
    )
    parser.add_argument("--lmin", type=int, default=64)
    parser.add_argument("--lmax", type=int, default=3072)
    parser.add_argument("--n-resamples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=123456)
    parser.add_argument("--ci", type=float, default=95.0)
    args = parser.parse_args()
    if args.n_resamples < 100 or not (0 < args.ci < 100):
        parser.error("Require at least 100 resamples and CI in (0,100)")
    if args.lmin > args.lmax:
        parser.error("--lmin must not exceed --lmax")

    indir = args.input_base / args.lattice_type / args.representation
    outdir = args.output_base / args.lattice_type / args.representation
    L, paired, means, sigma = read_inputs(indir, args.lmin, args.lmax)
    baseline, fitted, fit_quality = compute_fits(
        L, paired, means, sigma, args.n_resamples, args.seed
    )
    save_outputs(
        outdir, L, baseline, fitted, fit_quality,
        args.lattice_type, args.representation, args.ci
    )


if __name__ == "__main__":
    main()
