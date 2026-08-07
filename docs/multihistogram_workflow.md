# Multihistogram workflow

## Overview

The multihistogram analysis requires configuration files organized by lattice size.

The expected input structure is:

data/eul2d/size_L/

where `L` is the lattice size of the `L x L` grid.

Each `size_L` directory contains files named:

conf-T-IDX.dat

where:

- `T` is the simulation temperature formatted with four decimal places
- `IDX` is the integer index assigned to that temperature

Example:

conf-2.2600-0.dat

## Temperature files

Temperature files are named:

temperatures_L.txt

Each line has the format:

T IDX

Example:

2.2600 0
2.2605 1
2.2610 2

The temperature windows used here were selected manually around the relevant peak region for each lattice size. They should not be interpreted as automatically computed FWHM intervals.

The temperature files can be generated with:

bash scripts/multihistogram/generate_reweighting_temperatures.sh data/temperatures

## Preparing multihistogram input data

The input files can be prepared with:

bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures \
  /path/to/simulation/results \
  data/eul2d

The script reads simulation outputs from directories of the form:

/path/to/simulation/results/T_TTTT/size_L/

For each temperature and lattice size, it combines:

- mefile.txt
- epfile.txt
- enfile.txt
- fpfile.txt
- fnfile.txt

The output files are written to:

data/eul2d/size_L/conf-T-IDX.dat

## Output columns

Each row in `conf-T-IDX.dat` contains ten columns:

1. first value from `mefile.txt`
2. second value from `mefile.txt`
3. value from `epfile.txt`
4. value from `enfile.txt`
5. first value from `fpfile.txt`
6. second value from `fpfile.txt`
7. third value from `fpfile.txt`
8. first value from `fnfile.txt`
9. second value from `fnfile.txt`
10. third value from `fnfile.txt`

## Notes

The generated `.dat` files are not tracked by Git. They should be treated as generated data and archived separately if needed for release.

## Reweighting grid refinement

The multihistogram runner uses a refinement factor to define the number of output temperature points.

By default:

POINTS = (N_temperatures - 1) * 50 + 1

where `N_temperatures` is the number of temperatures in `temperatures_L.txt`.

This means that each interval between consecutive input temperatures is divided into 50 smaller intervals. For example, if the input spacing is 0.0005, the output reweighting spacing is approximately 0.0005 / 50 = 0.00001.

The refinement factor can be changed by passing a fourth argument:

bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures \
  results/multihistogram \
  100

The output files are named deterministically as:

results_L.dat
output2_L.dat

where `L` is the lattice size.
