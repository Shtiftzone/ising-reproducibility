# Software requirements

This repository contains C, C++, Bash, and Python code for two-dimensional Ising simulations, multihistogram reweighting, statistical analysis, diagnostics, and figure generation.

## System requirements

The workflow is intended for a Unix-like environment, preferably Linux, with:

- Bash
- GNU Make
- A C compiler supporting C11 (GCC is recommended)
- G++ for the C++ multihistogram code and the autocorrelation diagnostics
- Python 3 and `pip`
- `awk`, `sed`, `find`, `sort`, `head`, `paste`, `wc`, and other standard command-line utilities
- GNU `time`, for the `-v` option used by the simulation runners

Some shell scripts use GNU-specific features, such as `sort -V`, and the simulation runners invoke `time -v`. Consequently, the scripts are not guaranteed to run unchanged on every Unix-like system (for example, macOS with its default BSD utilities).

The provided Bash runners do not depend on a particular job scheduler. Large production runs can be submitted to a cluster scheduler separately.

## C and C++ compilation

### Simulation code

The simulation code is located in `src/simulation/`. Its Makefile explicitly uses the following C compilation options:

```text
-std=c11 -O3 -Wall -pedantic
```

It compiles the square- and triangular-lattice programs for both the `cell` and `vertex` constructions, for the supported lattice sizes.

From the repository root:

```bash
make -C src/simulation CC=gcc
```

The explicit `CC=gcc` selects GCC. Without this override, GNU Make may use the compiler configured as `cc` in its environment or built-in defaults.

### Multihistogram code

The multihistogram source is located in `src/multihistogram/` and is compiled with `g++`:

```bash
make -C src/multihistogram
```

This produces the executable:

```text
src/multihistogram/Rw
```

The multihistogram Makefile uses `-O2` and defines `SSE2` and `AVX2` preprocessor macros. It does **not** explicitly select a C++ language standard using a `-std=...` flag: the language mode is therefore determined by the compiler's default settings.

The presence of architecture-related macros makes it particularly important to check compatibility when building on a different platform or CPU architecture. The repository does not currently provide a documented, platform-independent build configuration.

## Bundled C/C++ dependencies

The simulation includes RANLUX source files under:

```text
src/simulation/ranlux-3.4/
```

The multihistogram component also includes its required RANLUX implementation (for example, `ranlxs.cpp`). These sources are built locally; a separate RANLUX installation is not required.

The compilation rules do not specify additional third-party C/C++ libraries beyond the bundled source and the standard toolchain/runtime libraries.

## Python dependencies

Python is used for statistical analysis, diagnostics, finite-size-scaling fits, and figure generation.

From the repository root, install the dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

The repository currently pins these versions:

| Package | Version |
|---|---|
| NumPy | `2.4.6` |
| pandas | `3.0.5` |
| Matplotlib | `3.11.1` |
| SciPy | `1.17.1` |

These are the versions specified in `requirements.txt`; they are not a guarantee that every operating system or Python version is supported. Use a Python interpreter for which `pip` can install all pinned packages.

A virtual environment is recommended to avoid conflicts with system-wide packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Additional requirements for autocorrelation diagnostics

The scripts:

- `scripts/analysis/diagnose_autocorr.py`
- `scripts/analysis/analyze_block_bootstrap.py`

compile a small temporary C++ driver and call the repository's existing `AutoCorr` implementation. They therefore require `g++` at **execution time**, not just during the initial project build.

The source and header are located at:

```text
src/multihistogram/Autocorr.cpp
src/multihistogram/Autocorr.h
```

For these scripts, pass the source explicitly using:

```bash
--autocorr-source src/multihistogram/Autocorr.cpp
```

The scripts search for the matching header in the source directory by default. Neither script runs the multihistogram calculation; both are diagnostics of previously recorded Monte Carlo time series.

## Reproducing figures versus regenerating raw data

There are two distinct levels of reproduction:

1. **Figures from processed CSVs.** Install the Python dependencies, obtain the required processed CSVs from the external data archive (once available) or generate them locally, and place them in the documented `data/analysis_csv/` and `data/table_csv/` directories. The plotting scripts can then be run without rebuilding or rerunning the C/C++ simulation and multihistogram programs.
2. **Full simulation and analysis workflow.** Build both C/C++ components, run simulations, prepare multihistogram inputs, perform reweighting, and run the analysis scripts. This requires substantially more computing time and storage.

The Git repository contains source code, scripts, documentation, dependency specifications, and configuration files, including `data/analysis_parameters/reweighting_windows.csv`. Raw simulation data, prepared multihistogram input, reweighting outputs, processed CSVs, tables, and generated figures are not tracked in Git. These data products are intended to be distributed separately through an external data archive. Until the archive is available, they must be generated locally or supplied from existing datasets.

For commands and expected input layouts, see:

- [Simulation workflow](simulation_workflow.md)
- [Multihistogram workflow](multihistogram_workflow.md)
- [Analysis workflow](analysis_workflow.md)
