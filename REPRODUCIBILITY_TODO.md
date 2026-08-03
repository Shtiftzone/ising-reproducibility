# Reproducibility TODO

## Goal

Prepare the code, data, documentation, and release materials needed to independently reproduce the Ising multihistogram analysis, following the TELOS-style reproducibility workflow.

## Release components

- [ ] Paper / preprint
- [ ] Analysis workflow repository
- [ ] Data release
- [ ] Software Availability Statement
- [ ] Data Availability Statement

## Reproducibility levels

- Level 1: reproduce figures and tables from released CSV files.
- Level 2: reproduce released CSV files from processed simulation data.
- Level 3: reproduce processed simulation data from raw simulation output.
- Level 4: reproduce raw simulation output by rerunning simulations.

## Results to reproduce

| ID | Result/Figure/Table | What it shows | Input data | Script | Output CSV | Status |
|---|---|---|---|---|---|---|
| R1 | Pseudo-critical temperature estimates | Tc(L) estimates for different lattice sizes | TODO | TODO | data/figure_csv/tc_estimates.csv | TODO |
| R2 | Multihistogram observable curves | Reweighted observables as functions of temperature | TODO | TODO | data/figure_csv/multihistogram_curves.csv | TODO |
| R3 | Specific heat peaks | Peak position and height of specific heat | TODO | TODO | data/figure_csv/specific_heat_peaks.csv | TODO |
| R4 | Susceptibility peaks | Peak position and height of susceptibility | TODO | TODO | data/figure_csv/susceptibility_peaks.csv | TODO |
| R5 | Binder cumulant crossings | Binder cumulant crossing temperatures | TODO | TODO | data/figure_csv/binder_crossings.csv | TODO |
| R6 | Final critical temperature extrapolation | Extrapolation of Tc(L) to infinite volume | TODO | TODO | data/figure_csv/tc_extrapolation.csv | TODO |
| R7 | Final estimates table | Final numerical estimates quoted in the paper | TODO | TODO | data/figure_csv/final_estimates.csv | TODO |

## Data release

- [ ] Identify all raw simulation outputs.
- [ ] Include raw simulation outputs in the data release, unless impossible or explicitly justified.
- [ ] Prepare processed data files.
- [ ] Prepare CSV files for all plotted data.
- [ ] Prepare CSV files for all tabulated data.
- [ ] Prepare CSV files for all fit parameters quoted in the paper.
- [ ] Create simulation metadata file.
- [ ] Document all data formats in README or docs/data_description.md.
- [ ] Check that the archive contains no unwanted files, e.g. .DS_Store, temporary files, logs.
- [ ] Publish data release on Zenodo/OSF/equivalent.
- [ ] Obtain DOI for data release.
- [ ] Add data release DOI to README and paper.

## Analysis workflow release

- [ ] Add current simulation code or document where it lives.
- [ ] Add current multihistogram analysis code.
- [ ] Add current plotting code.
- [ ] Add script to check required input files.
- [ ] Add script to reproduce processed data.
- [ ] Add script to reproduce figure CSV files.
- [ ] Add script to reproduce figures.
- [ ] Add script to reproduce tables.
- [ ] Remove hard-coded absolute paths.
- [ ] Replace personal/local paths with relative paths.
- [ ] Add command-line arguments where needed.
- [ ] Add requirements.txt or environment.yml.
- [ ] Record final commit hash used for the paper.
- [ ] Tag final version, e.g. v1.0-paper.
- [ ] Create GitHub release.
- [ ] Archive analysis workflow on Zenodo or equivalent.
- [ ] Obtain DOI for analysis workflow.
- [ ] Add analysis workflow DOI to README and paper.

## Documentation

- [ ] Update README with setup instructions.
- [ ] Update README with exact reproduction commands.
- [ ] Document repository structure.
- [ ] Document raw data location.
- [ ] Document processed data files and columns.
- [ ] Document figure CSV files and columns.
- [ ] Document software dependencies.
- [ ] Document HPC environment if needed.
- [ ] Document random seeds and simulation parameters.
- [ ] Add Software Availability Statement draft.
- [ ] Add Data Availability Statement draft.
- [ ] Add citation information.
- [ ] Decide license.

## Cross-checks

- [ ] Fresh clone test: clone repository into a new directory.
- [ ] Install dependencies from scratch.
- [ ] Run data check script.
- [ ] Reproduce figures from figure CSV files.
- [ ] Reproduce figure CSV files from processed data.
- [ ] Ask another collaborator to cross-check the analysis workflow.
- [ ] Ask another collaborator to cross-check data release contents.
- [ ] Confirm that all paper figures/tables correspond to released CSV files.
- [ ] Confirm that paper contains Software Availability Statement.
- [ ] Confirm that paper contains Data Availability Statement.
- [ ] Confirm that data DOI is included in paper.
- [ ] Confirm that analysis workflow DOI is included in paper.