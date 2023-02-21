# assignments

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/aberHRML/assignments/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aberHRML/assignments/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/aberHRML/assignments/branch/master/graph/badge.svg)](https://app.codecov.io/gh/aberHRML/assignments?branch=master)
[![license](https://img.shields.io/badge/license-GNU%20GPL%20v3.0-blue.svg)](https://github.com/aberHRML/assignments/blob/master/DESCRIPTION)
[![GitHub release](https://img.shields.io/github/release/aberHRML/assignments.svg)](https://GitHub.com/aberHRML/assignments/releases/)
<!-- badges: end -->
 
> An R package for automated molecular formula assignment of ultra-high resolution ESI-MS based metabolomics data

### Overview

This R package provides an automated molecular formula assignment approach for electrospray ionisation ultra-high resolution mass spectrometry (ESI-HRMS) metabolomics data. This includes data from direct and flow injection/infustion (FIE-HRMS) fingerprinting as well as liquid chromatography mass spectrometry (LC-HRMS) profiling. The approach includes correlation analysis, relationship calculation, molecular formula generation and selection and graphical component selection based on adducts, isotopes and transformations.

### Installation

The `assignments` package can be installed from GitHub using the
following:

``` r
remotes::install_github('aberHRML/assignments')
```

### Learn more

The package documentation can be browsed online at
<https://aberHRML.github.io/assignments/>.

If this is your first time using `assignments` see the
[vignette](https://aberhrml.github.io/assignments/articles/assignments.html) for information on how to get started.

If you believe you’ve found a bug in `assignments`, please file a bug (and, if possible, a [reproducible example](https://reprex.tidyverse.org)) at
<https://github.com/aberHRML/assignments/issues>.
