# Ising reproducibility

This repository contains the code and workflow used for two-dimensional Ising model simulations, multihistogram reweighting, and figure generation.

The analysis includes simulations on:

* square lattices
* triangular lattices

The repository is organized as a reproducibility package containing source code, portable workflow scripts, documentation, and lightweight data products. Large generated datasets are distributed separately.

## Associated publication

This repository accompanies:

**[Paper title]**
[Authors]

Paper: **[DOI / arXiv link to be added]**

The version of the repository associated with the published analysis will be archived as a persistent release.

## Repository structure

```text
src/
  simulation/          C simulation code
  multihistogram/      C++ multihistogram reweighting code

scripts/
  simulation/          simulation workflow scripts
  multihistogram/      temperature selection, input preparation, and reweighting scripts

docs/
  simulation_workflow.md
  multihistogram_workflow.md
  software_requirements.md

data/                  directories populated by the analysis workflow
results/               generated analysis outputs

requirements.txt       Python dependencies
```

The `data/` and `results/` directories are populated by the workflow and are not part of the source-code repository, except for lightweight data products explicitly included for figure reproduction.

## Workflow overview

The analysis follows the pipeline:

```text
Build simulation executables
        ¡
Run square- or triangular-lattice simulations
        ¡
Generate lattice-specific reweighting temperature files
        ¡
Prepare multihistogram input data
        ¡
Build the multihistogram executable
        ¡
Run multihistogram reweighting
        ¡
Prepare figure data
        ¡
Generate figures
```

Square- and triangular-lattice datasets are kept separate throughout the workflow.

## Software requirements

The workflow uses C, C++, Bash, and Python in a Unix-like environment.

Detailed requirements and build information are provided in:

```text
docs/software_requirements.md
```

Python dependencies can be installed with:

```bash
python -m pip install -r requirements.txt
```

## Quick start

### 1. Build the simulation executables

```bash
cd src/simulation
make
cd ../..
```

### 2. Run the simulations

Square lattice:

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
  src/simulation \
  results/square_simulations \
  2.22 \
  2.32 \
  0.0005 \
  200
```

Triangular lattice:

```bash
bash scripts/simulation/run_triangular_2d_simulations.sh \
  src/simulation \
  results/triangular_simulations \
  3.59 \
  3.72 \
  0.0005 \
  200
```

These commands run all lattice sizes defined in the simulation workflow. Full production runs are computationally expensive and can be parallelized using the scheduler or execution environment available to the user.

### 3. Generate the reweighting temperature files

Square lattice:

```bash
bash scripts/multihistogram/generate_square_reweighting_temperatures.sh
```

Triangular lattice:

```bash
bash scripts/multihistogram/generate_triangular_reweighting_temperatures.sh
```

This creates:

```text
data/temperatures/square/
data/temperatures/triangular/
```

The reweighting windows are explicitly defined by lattice size in the corresponding scripts.

### 4. Prepare multihistogram input data

Square lattice:

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures/square \
  results/square_simulations \
  data/eul2d/square
```

Triangular lattice:

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures/triangular \
  results/triangular_simulations \
  data/eul2d/triangular
```

### 5. Build the multihistogram executable

```bash
cd src/multihistogram
make
cd ../..
```

### 6. Run multihistogram reweighting

Square lattice:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/square \
  data/eul2d/square \
  results/multihistogram/square
```

Triangular lattice:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/triangular \
  data/eul2d/triangular \
  results/multihistogram/triangular
```

The default reweighting refinement factor is `50`.

## Detailed documentation

More detailed descriptions of each stage are available in:

* `docs/software_requirements.md`
* `docs/simulation_workflow.md`
* `docs/multihistogram_workflow.md`

## Data availability

Large simulation outputs, prepared multihistogram inputs, and multihistogram results are distributed separately from the source-code repository.

Data release: **[DOI / persistent data repository link to be added]**

After obtaining the accompanying dataset, its directory structure corresponds to the paths used by the workflow:

```text
results/square_simulations/
results/triangular_simulations/

data/eul2d/square/
data/eul2d/triangular/

results/multihistogram/square/
results/multihistogram/triangular/
```

Lightweight data products required to reproduce the final figures will be included directly in the repository.

## Computational cost

Full production simulations span multiple lattice sizes and temperatures and are intended for parallel execution on computational infrastructure.

Reference runtime and machine information:

```text
Hardware / cluster: [to be added]
Compiler versions:  [to be added]
Simulation runtime: [to be added]
MH runtime:         [to be added]
```

These values describe the reference production environment and are not requirements for running the workflow.

## Reproducibility

The scripts in this repository define the parameters, directory structure, and processing steps used in the analysis.

Simulation runs use the seed handling implemented in the simulation code and workflow scripts. The repository preserves the random-number-generation behavior used for the production analysis.

The workflow is scheduler-independent; cluster-specific parallelization can be added without changing the analysis steps.

## Figure reproduction

Scripts and lightweight data products for reproducing the figures will be added after the final analysis outputs are available.

The intended workflow is:

```text
multihistogram outputs
        ¡
figure-level data
        ¡
paper figures
```

## Citation

Machine-readable citation information will be provided in:

```text
CITATION.cff
```

## License

Licensing information will be provided in:

```text
LICENSE
```
