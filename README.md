GEX
================

Generate XICs and mass spectra for theoretical charge states of target
protein sequences with or without isotopic depletion or
carbamidomethylation.

## Installation

Install from GitHub:

``` r
remotes::install_github("davidsbutcher/GEX")
```

## Usage

The main functionality of the package is split into two functions,
`make_every_XIC` and `make_every_spectrum`. The output of
`make_every_XIC` must be used for `make_every_spectrum`.

### Arguments to `make_every_XIC`

#### Mandatory

  - `rawFileDir` Character vector, length 1. Directory containing the
    raw file to be searched. Can be any level above the raw file in the
    directory structure.

  - `rawFileName` Character vector, length 1. Full name of the raw file
    to be searched including extension.

  - `targetSeqData` Character vector, length 1. Full or relative path to
    the csv file containing the target sequences. The csv must have
    columns specifying name, proteoform sequence, PTM name, PTM formula
    to add and PTM formula to subtract. See
    `input/example_sequences.csv` for an example. Column names should be
    `UNIPROTKB`, `ProteoformSequence`, `PTMname`, `FormulaToAdd`, and
    `FormulaToSubtract` unless specified otherwise using optional
    arguments.

#### Optional

  - `outputDir` Character vector, length 1. Directory to place output
    results. A new directory will be created in `outputDir` providing
    the date and time of the analysis, raw file name, number of
    sequences searched, and whether depleted isotopes were used.

  - `mz_range` Numeric vector, length 2. Range of m/z values to use for
    filtering theoretical isotopic distributions. Should most likely be
    set to m/z range of MS1 scans in raw file. Defaults to 600 - 2000.

  - `use_depleted_isotopes` Boolean value (TRUE or FALSE). Determines
    whether depleted isotopic abundances are used to generate
    theoretical isotopic distributions. Defaults to FALSE.

  - `mass_range` Numeric vector, length 2. Mass range of proteoforms
    from input csv to include in search. This range is compared to
    monoisotopic masses calculated from the proteoform sequence.
    Defaults to 0 - 100,000.

  - `abund_cutoff` Numeric vector, length 1. Percentage cutoff for
    relative abundance of theoretical isotopologues - all theoretical
    isotopologues with relative abundance (compared to most abundant
    isotopologue at 100%) lower than this value are truncated prior to
    XIC generation. Defaults to 5.

  - `XIC_tol` Numeric vector, length 1. Tolerance (in ppm) to use for
    generation of XICs from theoretical isotopic distributions. Defaults
    to 2. Can be set higher for non-complex samples.

  - `use_IAA` Boolean value (TRUE or FALSE). Determines whether
    cysteines in target sequences should be considered
    carbamidomethylated when generating chemical formulas prior to
    calculating theoretical isotopologue m/z values. Note that this is
    an all-or-nothing parameter, with no option for partial
    carbamidomethylation. Defaults to FALSE.

  - `target_charges` Numeric vector, length 2. Specified what range of
    charge states to use for calculation of theoretical isotopic
    distributions. Defaults to 1 - 50.

  - `target_col_name` Character vector, length 1. Allows the name of the
    target sequences name column from the input csv to be specified.
    Defaults to `UNIPROTKB`.

  - `target_sequence_col_name` Character vector, length 1. Allows the
    name of the target sequences proteoform sequence column from the
    input csv to be specified. Defaults to `ProteoformSequence`.

  - `PTMname_col_name` Character vector, length 1. Allows the name of
    the target sequences PTM name column from the input csv to be
    specified. Defaults to `PTMname`.

  - `PTMformula_col_name1` Character vector, length 1. Allows the name
    of the target sequences PTM formula to add column from the input csv
    to be specified. Defaults to `FormulaToAdd`.

  - `PTMformula_col_name2` Character vector, length 1. Allows the name
    of the target sequences PTM formula to subtract column from the
    input csv to be specified. Defaults to `FormulaToSubtract`.

### Arguments to `make_every_spectrum`

The only mandatory argument to `make_every_spectrum` is the output from
`make_every_XIC`, which should be the first argument (see example
below).

#### Optional

  - `mz_window` Numeric vector, length 1. Size of m/z window to use in
    generation of zoomed mass spectra. Defaults to 3.

### Example run with minimal arguments

``` r
library(GEX)

XIC_output <-
   make_every_XIC(
      rawFileDir = "C:/Users/dbutcher/Documents/raw_files",
      rawFileName = "example_file.raw",
      targetSeqData = "input/example_sequences.csv",
      outputDir = "C:/Users/dbutcher/Documents/GEX_output"
   )

GEX_timer <-
   make_every_spectrum(
      XIC_output
   )
```
