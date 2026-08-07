# Simulation workflow

## Overview

This repository contains portable scripts for running two-dimensional Ising simulations on:

1. the square lattice
2. the triangular lattice

The simulation executables are built from the source code in:

```text
src/simulation/
```

The generated simulation outputs are not tracked by Git. They should be treated as generated data and archived separately for release if needed.

## Building the simulation executables

From the repository root:

```bash
cd src/simulation
make
cd ../..
```

This builds the executables used in the analysis:

```text
ising2d-L    # square lattice
isingtr-L    # triangular lattice
```

where `L` is the linear lattice size of the `L x L` grid.

The default lattice sizes are:

```text
64 96 128 192 256 384 512 768 1024 1536 2048 3072
```

## Square-lattice simulations

The square-lattice simulations can be run with:

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
  src/simulation \
  results/square_simulations \
  2.22 \
  2.32 \
  0.0005 \
  200
```

Arguments:

```text
1. executable directory
2. output results directory
3. starting temperature
4. final temperature
5. temperature step
6. number of configurations
```

The default square-lattice temperature range is:

```text
T_start = 2.22
T_end   = 2.32
T_step  = 0.0005
Nconf   = 200
```

## Triangular-lattice simulations

The triangular-lattice simulations can be run with:

```bash
bash scripts/simulation/run_triangular_2d_simulations.sh \
  src/simulation \
  results/triangular_simulations \
  3.59 \
  3.72 \
  0.0005 \
  200
```

Arguments:

```text
1. executable directory
2. output results directory
3. starting temperature
4. final temperature
5. temperature step
6. number of configurations
```

The default triangular-lattice temperature range is:

```text
T_start = 3.59
T_end   = 3.72
T_step  = 0.0005
Nconf   = 200
```

## Output structure

Both simulation scripts write outputs in the same directory structure:

```text
results/<simulation_type>/T_TTTTT/size_L/
```

For example:

```text
results/square_simulations/T_2.22000/size_128/
results/triangular_simulations/T_3.59000/size_128/
```

For each temperature and lattice size, the simulation writes:

```text
mefile.txt
bpfile.txt
bnfile.txt
epfile.txt
enfile.txt
fpfile.txt
fnfile.txt
```

The multihistogram input preparation workflow currently uses:

```text
mefile.txt
epfile.txt
enfile.txt
fpfile.txt
fnfile.txt
```

The files `bpfile.txt` and `bnfile.txt` are produced by the simulation but are not currently used by the multihistogram input preparation script.

## Seed handling

Each simulation script creates a local seed file in the selected results directory if one does not already exist:

```text
results/<simulation_type>/seed.txt
```

The default initial seed is:

```text
123456
```

For a final reproducibility release, the exact seed files used for production runs should be included in the data release or documented explicitly.

## Generated data

The generated simulation outputs are not tracked by Git.

They should be archived separately, together with metadata describing:

```text
lattice type
lattice size L
temperature T
number of configurations
seed file
simulation executable version / Git commit
date of generation
machine or cluster environment
```
