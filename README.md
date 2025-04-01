# eHarmonize
#### Contact: Alyssa Zhu <alyssahz@usc.edu>

[ENIGMA](https://enigma.ini.usc.edu/) Harmonization (`eharmonize`) is a python-based package that harmonizes provided data to included lifespan reference curves. The package is currently set up for outputs of the ENIGMA-DTI pipeline.

Please see wiki for additional details and information. 

Current version (Release v0.0.0): [![DOI](https://zenodo.org/badge/622720163.svg)](https://doi.org/10.5281/zenodo.15116824)

## Installation 

As this package has not been added to PyPI, please follow the instructions below carefully.

One-step process:
`pip install git+https://github.com/ahzhu/eharmonize.git`

Two-step process:
1. To clone the repository:
`git clone https://github.com/ahzhu/eharmonize.git`

2. After downloading or cloning (in the same directory that `git clone` was run in):
`pip install ./eharmonize`

N.B. The above installation commands will require administrator/root privileges of the system. If you don't have administrator/root privileges, please try adding the `--user` flag.

Example:
`pip install --user ./eharmonize`

## Usages

### Display available subcommands
`eharmonize --help`

### Display flags/options for a specific subcommand
`eharmonize [subcommand] --help`

### Use subcommand
`eharmonize [subcommand] [--flags inputs]`

## Subcommands

### harmonize-fa

The **harmonize-fa** subcommand takes in FA spreadsheets from the ENIGMA-DTI pipeline and exports FA values harmonized to the included reference curves.

The input CSV should have columns for `Age` and `Sex`. Optional columns include `SITE` and `Dx`. If `SITE` is not included, all data will be assumed to have come from a single site/scanner. If `Dx` is not included, all data will be assumed to have come from controls. In the `Dx` column, subjects should be denoted as either `case` or `control`. The harmonization parameters will be calculated based on the control group and then applied to the cases. 

N.B. We currently do not support harmonization for case-only studies. If such data is detected, the program will self-quit.

Required Inputs:
* `--incsv CSV`             File path of the CSV with the necessary covariates
                            and FA measures 
* `--outdir DIR`            Directory to write out outputs to

Optional Inputs:
* `--reference version`    Version number of desired reference 
                           [default: v0.1] 
* `--rerun`                If used, will overwrite previous outputs 

Example usage:
 
`eharmonize harmonize-fa --incsv enigma_FA_spreadsheet.csv --outdir /path/to/write/to`
`eharmonize harmonize-fa --incsv enigma_FA_spreadsheet.csv --outdir /path/to/write/to --reference v0.0 --rerun`

### apply-harmonization

The **apply-harmonization** subcommand will apply an existing harmonization model to new data. The input CSV should follow the same requirements as in **harmonize-fa**. The `--model` and `--log` inputs are outputs from **harmonize-fa**. 

Required Inputs:
* `--incsv CSV`            File path of the CSV file with the necessary covariates 
                           and FA measures 
* `--model MODEL`          Model file as output by initial harmonization 
* `--log JSON`             JSON file from initial harmonization pant IDs to include
* `--outdir DIR`           Directory to write out outputs to

Optional Inputs: 
* `--metric`               [default: FA]

Example usage:
 
`eharmonize apply-harmonization --incsv enigma_FA_new_spreadsheet.csv --model --log --outdir`
`eharmonize apply-harmonization --incsv enigma_FA_new_spreadsheet.csv --model --log --outdir --metric FA`

### harmonize-dti

In Development

## Acknowledgments 

We would like to thank the study coordinators and participants of the public datasets used to create the DTI reference curves.

* PING ([Jernigan et al. 2016](https://doi.org/10.1016/j.neuroimage.2015.04.057))
* HCP-Development ([Somerville et al. 2018](https://doi.org/10.1016/j.neuroimage.2018.08.050))
* ABCD ([Hagler et al. 2019](https://doi.org/10.1016/j.neuroimage.2019.116091))
* SLIM ([Qiu et al.](http://dx.doi.org/10.15387/fcp_indi.retro.slim))
* HCP ([Van Essen et al. 2013](https://doi.org/10.1016/j.neuroimage.2013.05.041))
* Cam-CAN ([Shafto et al. 2014](https://doi.org/10.1186/s12883-014-0204-1))
* PPMI ([Marek et al. 2011](https://doi.org/10.1016/j.pneurobio.2011.09.005))
* UK Biobank ([Miller et al. 2016](https://doi.org/10.1038/nn.4393); Application Number 11559)
* OASIS3 ([LaMontagne et al. 2019](https://doi.org/10.1101/2019.12.13.19014902))
* ADNI3 ([Weiner et al. 2016](https://doi.org/10.1016/j.jalz.2016.10.006))
