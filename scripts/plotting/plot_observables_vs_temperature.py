#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MultipleLocator, ScalarFormatter


PLOTS = {
    "mag_vs_T.csv": {
        "output": "mag_vs_T_abs.png",
        "title": "Magnetization vs Temperature",
        "ylabel": r"$\langle |M| \rangle$",
        "tick_group": "order_parameter",
    },
    "energy_vs_T.csv": {
        "output": "energy_vs_T.png",
        "title": "Energy vs Temperature",
        "ylabel": r"$\langle E \rangle$",
        "tick_group": "energy_like",
    },
    "euler_sym_vs_T.csv": {
        "output": "euler_sym_vs_T_abs.png",
        "title": r"$EC_{\mathrm{sym}}$ vs Temperature",
        "ylabel": r"$\langle |EC_{\mathrm{sym}}| \rangle$",
        "tick_group": "order_parameter",
    },
    "ec_avg_vs_T.csv": {
        "output": "ec_avg_vs_T.png",
        "title": r"$EC_{\mathrm{avg}}$ vs Temperature",
        "ylabel": r"$\langle EC_{\mathrm{avg}} \rangle$",
        "tick_group": "energy_like",
    },
}


class FixedDecimalScalarFormatter(ScalarFormatter):
    """
    Scientific formatter with a fixed number of decimal places.
    """

    def __init__(self, decimals: int):
        self.decimals = decimals

        super().__init__(
            useOffset=False,
            useMathText=True,
        )

    def _set_format(self):
        self.format = f"%1.{self.decimals}f"

        if self._usetex or self._useMathText:
            self.format = (
                r"$\mathdefault{%s}$"
                % self.format
            )


def format_y_axis(
    ax,
    values: pd.Series,
    errors: pd.Series,
    tick_group: str,
) -> None:
    """
    Format related observables consistently.

    M and EC_sym:
        integer labels after scientific scaling,
        e.g. 0, 1, 2, 3 x 10^5.

    E and EC_avg:
        one decimal place after scientific scaling,
        e.g. -1.8, -1.6, -1.4 x 10^7.

    The tick spacing is automatically reduced until
    at least two major ticks fall inside the data range.
    """

    lower = (values - errors).min()
    upper = (values + errors).max()

    max_abs = max(
        abs(lower),
        abs(upper),
    )

    if max_abs == 0:
        return

    exponent = math.floor(
        math.log10(max_abs)
    )

    scale = 10.0 ** exponent

    if tick_group == "order_parameter":
        decimals = 0
        step = scale

    elif tick_group == "energy_like":
        decimals = 1
        step = 0.2 * scale

    else:
        raise ValueError(
            f"Unknown tick group: {tick_group}"
        )

    # Reduce the tick spacing until at least two major ticks
    # fall inside the actual data range.
    while True:
        first_tick = (
            math.ceil(lower / step)
            * step
        )

        last_tick = (
            math.floor(upper / step)
            * step
        )

        if last_tick >= first_tick:
            n_ticks = (
                int(
                    round(
                        (last_tick - first_tick)
                        / step
                    )
                )
                + 1
            )
        else:
            n_ticks = 0

        if n_ticks >= 2:
            break

        step /= 2.0

    ax.yaxis.set_major_locator(
        MultipleLocator(step)
    )

    formatter = FixedDecimalScalarFormatter(
        decimals=decimals
    )

    formatter.set_scientific(True)

    formatter.set_powerlimits(
        (exponent, exponent)
    )

    ax.yaxis.set_major_formatter(
        formatter
    )


def make_plot(
    csv_path: Path,
    output_path: Path,
    title: str,
    ylabel: str,
    tick_group: str,
) -> None:
    df = pd.read_csv(csv_path)

    required = {
        "T",
        "mean",
        "jackknife_se",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{csv_path} is missing required columns: "
            f"{sorted(missing)}"
        )

    df = df.sort_values("T")

    fig, ax = plt.subplots(
        figsize=(6, 4)
    )

    ax.errorbar(
        df["T"],
        df["mean"],
        yerr=df["jackknife_se"],
        fmt="o",
        ms=3,
        capsize=2,
        elinewidth=1,
        alpha=0.9,
    )

    ax.set_title(title)

    ax.set_xlabel(
        r"Temperature $T$"
    )

    ax.set_ylabel(
        ylabel
    )

    format_y_axis(
        ax,
        df["mean"],
        df["jackknife_se"],
        tick_group,
    )

    ax.grid(
        True,
        linewidth=0.3,
        alpha=0.5,
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"Wrote {output_path}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate temperature-dependence plots from "
            "figure-level CSV files."
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
        "--input-dir",
        type=Path,
        default=Path(
            "data/figure_csv"
        ),
        help=(
            "Base directory containing "
            "figure-level CSV files."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "results/figures"
        ),
        help=(
            "Base directory for generated figures."
        ),
    )

    args = parser.parse_args()

    variant_input_dir = (
        args.input_dir
        / args.lattice_type
        / args.representation
    )

    if not variant_input_dir.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist: "
            f"{variant_input_dir}"
        )

    variant_output_dir = (
        args.output_dir
        / args.lattice_type
        / args.representation
    )

    variant_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    lattice_label = (
        "Sq."
        if args.lattice_type == "square"
        else "Tri."
    )

    for csv_name, config in PLOTS.items():
        csv_path = (
            variant_input_dir
            / csv_name
        )

        if not csv_path.is_file():
            print(
                f"Skipping missing file: "
                f"{csv_path}"
            )
            continue

        output_path = (
            variant_output_dir
            / config["output"]
        )

        title = (
            f"{config['title']} "
            f"({lattice_label})"
        )

        make_plot(
            csv_path=csv_path,
            output_path=output_path,
            title=title,
            ylabel=config["ylabel"],
            tick_group=config["tick_group"],
        )

    print("Done.")


if __name__ == "__main__":
    main()