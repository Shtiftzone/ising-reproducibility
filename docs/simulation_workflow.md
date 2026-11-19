# Simulation workflow

## Overview

This document describes the simulation workflow for the two-dimensional Ising model on square and triangular lattices.

Two constructions are available for each lattice:

- `cell`: spins represented by cells
- `vertex`: spins represented by vertices

The simulation code supports all four lattice/construction combinations. The analysis presented in the associated publication focuses on `square/cell` and `triangular/vertex`.

The source code is located in `src/simulation/`, and the execution scripts are in `scripts/simulation/`. The scripts are scheduler-independent and can be run in a Unix-like environment.

## Building the simulation executables

From the repository root:

```bash
cd src/simulation
make
cd ../..
```

The Makefile builds the following executables for each supported lattice size:

```text
ising2d-cell-L      # square lattice, spin-as-cell
ising2d-vertex-L    # square lattice, spin-as-vertex
isingtr-cell-L      # triangular lattice, spin-as-cell
isingtr-vertex-L    # triangular lattice, spin-as-vertex
```

Here, `L` is the linear lattice size of the periodic `L x L` lattice. The supported sizes are:

```text
64 96 128 192 256 384 512 768 1024 1536 2048 3072
```

## Running the simulations

The simulation runners accept the arguments:

```text
<executables_dir> <results_dir> <construction> [T_start] [T_end] [T_step] [Nconf]
```

The first three arguments are required. The optional arguments specify the temperature interval, its spacing, and the number of recorded configurations per temperature and lattice size.

Both runners process all supported lattice sizes for each temperature in the requested range. If an executable is missing, the corresponding lattice size is skipped.

### Square lattice

Cell construction (used in the publication):

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
    src/simulation \
    results/square_cell_simulations \
    cell
```

Vertex construction:

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
    src/simulation \
    results/square_vertex_simulations \
    vertex
```

The default square-lattice parameters are:

```text
T_start = 2.22
T_end   = 2.32
T_step  = 0.0005
Nconf   = 200
```

To override the defaults, append the four optional arguments, for example:

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
    src/simulation results/square_cell_simulations cell \
    2.22 2.32 0.0005 200
```

### Triangular lattice

Vertex construction (used in the publication):

```bash
bash scripts/simulation/run_triangular_2d_simulations.sh \
    src/simulation \
    results/triangular_vertex_simulations \
    vertex
```

Cell construction:

```bash
bash scripts/simulation/run_triangular_2d_simulations.sh \
    src/simulation \
    results/triangular_cell_simulations \
    cell
```

The default triangular-lattice parameters are:

```text
T_start = 3.59
T_end   = 3.72
T_step  = 0.0005
Nconf   = 200
```

To override the defaults, append the four optional arguments, for example:

```bash
bash scripts/simulation/run_triangular_2d_simulations.sh \
    src/simulation results/triangular_vertex_simulations vertex \
    3.59 3.72 0.0005 200
```

## Monte Carlo updates and measurement spacing

The simulation code implements a hybrid Metropolis-Wolff update procedure. One update iteration consists of one sequential Metropolis sweep followed by 15 Wolff cluster updates.

The following values are defined in `src/simulation/ising.c`:

```text
Ntherm  = 2000
Nwolff  = 15
Ndecorr = 4   # C-program default if no final argument is supplied
```

The thermalization stage performs `Ntherm` update iterations before the first measurement. During measurement, `Ndecorr` update iterations are performed between consecutive recorded configurations.

**Both supplied Bash simulation runners explicitly pass `Ndecorr = 10` to the C executable**, overriding its internal default of 4. Thus, simulations run through these scripts use **10 update iterations between recorded configurations**, consistent with the corrected description in the associated paper.

The number of saved configurations is specified by `Nconf` (200 by default). The C-program default of 4 applies only when the executable is invoked without an explicit `Ndecorr` argument.

## Output structure

Results are organized by lattice, construction, temperature, and lattice size:

```text
results/<simulation_type>/
    seed.txt
    simulation_logs_*.txt
    T_<temperature>/
        size_<L>/
            mefile.txt
            bpfile.txt
            bnfile.txt
            epfile.txt
            enfile.txt
            fpfile.txt
            fnfile.txt
```

Temperature directory names use five decimal places. Examples include:

```text
results/square_cell_simulations/T_2.22000/size_128/
results/square_vertex_simulations/T_2.22000/size_128/
results/triangular_cell_simulations/T_3.59000/size_128/
results/triangular_vertex_simulations/T_3.59000/size_128/
```

The measurement files contain:

| File | Contents |
|---|---|
| `mefile.txt` | Magnetization and energy |
| `bpfile.txt` | Betti numbers of positive-spin complexes |
| `bnfile.txt` | Betti numbers of negative-spin complexes |
| `epfile.txt` | Euler characteristic of positive-spin complexes |
| `enfile.txt` | Euler characteristic of negative-spin complexes |
| `fpfile.txt` | Cell counts of positive-spin complexes |
| `fnfile.txt` | Cell counts of negative-spin complexes |

The multihistogram input preparation step uses `mefile.txt`, `epfile.txt`, `enfile.txt`, `fpfile.txt`, and `fnfile.txt`. The Betti-number files are generated but are not used by that step.

## Random seed handling

Each simulation runner creates `seed.txt` in its selected results directory if the file does not already exist. The initial value written by the runners is:

```text
123456
```

The C executable reads the seed file at the start of each invocation and initializes its RANLUX random-number generator with that value. The current implementation **does not update `seed.txt` with the final generator state** after a simulation.

Consequently, invocations using the same seed file are initialized with the same seed. The workflow does **not** assign independent seeds to different temperatures or lattice sizes, and statistical independence between their samples must not be inferred from the seed-handling procedure.

The description above documents the behavior of the code used to produce the simulation results; it does not imply that independent random-number streams were generated for different runs.

## Execution logs and computational cost

The Bash runners append timing information to:

```text
results/square_cell_simulations/simulation_logs_square_cell.txt
results/triangular_vertex_simulations/simulation_logs_triangular_vertex.txt
```

The same naming pattern applies to the other constructions. The runners execute temperatures and lattice sizes sequentially; cluster-level parallelization can be supplied externally.

Full production runs, especially at the largest lattice sizes, can be computationally expensive.

## Generated data and subsequent analysis

Full simulation outputs are generated locally and are not included in the source-code repository. Large datasets are distributed separately.

Raw measurements are used for multihistogram input preparation, temperature-dependent observables, distribution-resolution analysis, and statistical diagnostics. See:

- [Multihistogram workflow](multihistogram_workflow.md)
- [Analysis and figure reproduction](analysis_workflow.md)
