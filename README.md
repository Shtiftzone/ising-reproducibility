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
  docs/analysis_workflow.md

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

From the repository root:

    cd src/simulation
    make
    cd ../..

This builds construction-specific executables of the form:

    ising2d-cell-L
    ising2d-vertex-L
    isingtr-cell-L
    isingtr-vertex-L

where L is the lattice size.

The two available constructions are:

    cell
    vertex

### 2. Run the simulations

The simulation runner scripts take the construction as an explicit argument: cell or vertex.

To reproduce the full workflow for both constructions, run both constructions for each lattice type.

Square lattice, cell construction:

    bash scripts/simulation/run_square_2d_simulations.sh \
      src/simulation \
      results/square_cell_simulations \
      cell

Square lattice, vertex construction:

    bash scripts/simulation/run_square_2d_simulations.sh \
      src/simulation \
      results/square_vertex_simulations \
      vertex

Triangular lattice, cell construction:

    bash scripts/simulation/run_triangular_2d_simulations.sh \
      src/simulation \
      results/triangular_cell_simulations \
      cell

Triangular lattice, vertex construction:

    bash scripts/simulation/run_triangular_2d_simulations.sh \
      src/simulation \
      results/triangular_vertex_simulations \
      vertex

The default simulation parameters are:

Square lattice:

    T_start = 2.22
    T_end   = 2.32
    T_step  = 0.0005
    Nconf   = 200

Triangular lattice:

    T_start = 3.62
    T_end   = 3.72
    T_step  = 0.0005
    Nconf   = 200

These defaults can be overridden by passing the optional arguments explicitly.

For example:

    bash scripts/simulation/run_square_2d_simulations.sh \
      src/simulation \
      results/square_cell_simulations \
      cell \
      2.22 \
      2.32 \
      0.0005 \
      200

The same argument pattern applies to the triangular-lattice runner.

Full production runs are computationally expensive and can be parallelized using the scheduler or execution environment available to the user.

### 3. Generate the reweighting temperature files

Square lattice:

    bash scripts/multihistogram/generate_square_reweighting_temperatures.sh

Triangular lattice:

    bash scripts/multihistogram/generate_triangular_reweighting_temperatures.sh

This creates:

    data/temperatures/square/
    data/temperatures/triangular/

The reweighting windows are explicitly defined by lattice size in the corresponding scripts.

### 4. Prepare multihistogram input data

Prepare the multihistogram input separately for each lattice type and construction.

Square lattice, cell construction:

    bash scripts/multihistogram/prepare_eul2d_input.sh \
      data/temperatures/square \
      results/square_cell_simulations \
      data/eul2d/square/cell

Square lattice, vertex construction:

    bash scripts/multihistogram/prepare_eul2d_input.sh \
      data/temperatures/square \
      results/square_vertex_simulations \
      data/eul2d/square/vertex

Triangular lattice, cell construction:

    bash scripts/multihistogram/prepare_eul2d_input.sh \
      data/temperatures/triangular \
      results/triangular_cell_simulations \
      data/eul2d/triangular/cell

Triangular lattice, vertex construction:

    bash scripts/multihistogram/prepare_eul2d_input.sh \
      data/temperatures/triangular \
      results/triangular_vertex_simulations \
      data/eul2d/triangular/vertex

### 5. Build the multihistogram executable

From the repository root:

    cd src/multihistogram
    make
    cd ../..

### 6. Run multihistogram reweighting

Run the multihistogram reweighting separately for each lattice type and construction.

Square lattice, cell construction:

    bash scripts/multihistogram/run_multihistogram_all.sh \
      src/multihistogram \
      data/temperatures/square \
      data/eul2d/square/cell \
      results/multihistogram/square/cell

Square lattice, vertex construction:

    bash scripts/multihistogram/run_multihistogram_all.sh \
      src/multihistogram \
      data/temperatures/square \
      data/eul2d/square/vertex \
      results/multihistogram/square/vertex

Triangular lattice, cell construction:

    bash scripts/multihistogram/run_multihistogram_all.sh \
      src/multihistogram \
      data/temperatures/triangular \
      data/eul2d/triangular/cell \
      results/multihistogram/triangular/cell

Triangular lattice, vertex construction:

    bash scripts/multihistogram/run_multihistogram_all.sh \
      src/multihistogram \
      data/temperatures/triangular \
      data/eul2d/triangular/vertex \
      results/multihistogram/triangular/vertex

The default reweighting refinement factor is 50.

A different refinement factor can be passed as the final argument. For example:

    bash scripts/multihistogram/run_multihistogram_all.sh \
      src/multihistogram \
      data/temperatures/square \
      data/eul2d/square/cell \
      results/multihistogram/square/cell \
      100
## Detailed documentation

More detailed descriptions of each stage are available in:

* `docs/software_requirements.md`
* `docs/simulation_workflow.md`
* `docs/multihistogram_workflow.md`
* `docs/analysis_workflow.md`

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
