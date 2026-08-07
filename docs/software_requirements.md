# Software requirements

This repository contains C, C++, Bash, and Python code used for two-dimensional Ising model simulations, multihistogram reweighting, and figure generation.

## System requirements

The workflow is intended for a Unix-like environment with:

* Bash
* GNU Make
* GCC
* G++
* awk
* standard Unix command-line tools

The scripts are written as portable Bash scripts and do not require a specific job scheduler.

## C/C++ requirements

The two-dimensional Ising simulation code is written in C.

The multihistogram reweighting code is written in C++.

The default builds use:

* `gcc` for the simulation code
* `g++` for the multihistogram code
* C11 for the simulation code
* C++14 for the multihistogram code

## External C/C++ dependencies

### RANLUX

A local copy of RANLUX used by the simulation code is included in:

```text id="qcrw8n"
src/simulation/ranlux-3.4/
```

The multihistogram source also contains the RANLUX files required by that component.

No additional external C/C++ libraries are required by the current two-dimensional workflow.

## Python requirements

Python is used for post-processing and figure generation.

Install the Python dependencies from the repository root with:

```bash id="17m7me"
python -m pip install -r requirements.txt
```

The current Python requirements are:

* numpy
* pandas
* matplotlib
* scipy

## Building the simulation code

From the repository root:

```bash id="d09ywn"
cd src/simulation
make
cd ../..
```

This builds the square- and triangular-lattice simulation executables.

## Building the multihistogram code

From the repository root:

```bash id="rjq2pk"
cd src/multihistogram
make
cd ../..
```

This builds the multihistogram executable `Rw`.

## Generated data

The simulation and multihistogram workflows generate output files locally.

Large simulation outputs and processed datasets are distributed separately from the source-code repository.
