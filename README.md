MultiAssayExperiment
==============

[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/MultiAssayExperiment.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MultiAssayExperiment)
[![Platforms](http://www.bioconductor.org/shields/availability/release/MultiAssayExperiment.svg)](https://www.bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html#archives)
[![Travis Build Status](https://travis-ci.org/waldronlab/MultiAssayExperiment.svg?branch=RELEASE_3_8)](https://travis-ci.org/waldronlab/MultiAssayExperiment)
[![Build status](https://ci.appveyor.com/api/projects/status/rf25e9h995wnto7n?svg=true)](https://ci.appveyor.com/project/LiNk-NY/multiassayexperiment-94gjw/branch/RELEASE_3_8)
[![wercker status](https://app.wercker.com/status/2aa523f23142715771256b85187d7bcb/s/ "wercker status")](https://app.wercker.com/project/byKey/2aa523f23142715771256b85187d7bcb)
[![Coverage Status](https://codecov.io/gh/waldronlab/MultiAssayExperiment/branch/RELEASE_3_8/graphs/badge.svg?branch=RELEASE_3_8)](https://codecov.io/github/waldronlab/MultiAssayExperiment?branch=RELEASE_3_8)
[![Downloads](http://www.bioconductor.org/shields/downloads/MultiAssayExperiment.svg)](https://bioconductor.org/packages/stats/bioc/MultiAssayExperiment)

## Overview

`MultiAssayExperiment` provides an integrative representation of multi-'omics
data on a set of samples. You can choose to represent your data using this
infrastructure or use pre-processed datasets from `curatedTCGAData`.

## Introduction

We recommend installing the stable release version of MultiAssayExperiment in
Bioconductor. This can be done using `BiocManager`:

```
if (!require("BiocManager"))
    install.packages("BiocManager")

library(BiocManager)

install("MultiAssayExperiment")
```

## Documentation website

The `MultiAssayExperiment` quick reference, vignettes, and function help pages
are available by browsing to the
[package website](http://waldronlab.io/MultiAssayExperiment).

## Cheatsheet

We provide a reference cheat sheet for a summary of available functionality
in `MultiAssayExperiment`.

<a href="https://github.com/waldronlab/cheatsheets/blob/master/MultiAssayExperiment_QuickRef.pdf"><img src="https://raw.githubusercontent.com/waldronlab/cheatsheets/master/pngs/MultiAssayExperiment_QuickRef.png" width="989" height="1091"/></a>

## Ready-to-use datasets from The Cancer Genome Atlas (TCGA)

For easy-to-use and analysis-ready `MultiAssayExperiment` objects, use the
`curatedTCGAData` experiment data package:

```
install("curatedTCGAData")
```

## How to make working with TCGA data easier

`TCGAutils` is a handy package for working with TCGA data represented as
`MultiAssayExperiment` objects. It is highly recommended to use `TCGAutils` for
TCGA identifier manipulation, sample identification and separation, and
much more.

To install `TCGAutils`, run the following:

```
install("TCGAutils")
```

## Getting help

If you find a bug or error within `MultiAssayExperiment`,
`curatedTCGAData`, or `TCGAutils`, please provide a reproducible example in the
respective GitHub issues pages:
'http://github.com/waldronlab/<package_name>/issues'.

If you have a question about using any of the packages listed above, create
a post on https://support.bioconductor.org/.

## The `MultiAssayExperiment` Bioconductor Special Interest Group

This group meets remotely to discuss this project approximately every 3 weeks.
If you are interested, please join the
[MultiAssayExperiment Google Group](https://groups.google.com/forum/#!forum/biocmultiassay)
and see the
[calendar](https://www.google.com/calendar/embed?src=9ar0qc8mpkv6b9intgmdcdf0ss%40group.calendar.google.com&ctz=America/New_York)
of upcoming meetings.

## The Fine Print

### Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect
all people who contribute through reporting issues, posting feature
requests, updating documentation, submitting pull requests or patches,
and other activities.

We are committed to making participation in this project a
harassment-free experience for everyone, regardless of level of
experience, gender, gender identity and expression, sexual
orientation, disability, personal appearance, body size, race, age, or
religion.

Examples of unacceptable behavior by participants include the use of
sexual language or imagery, derogatory comments or personal attacks,
trolling, public or private harassment, insults, or other
disrespectful conduct.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct. Project
maintainers who do not follow the Code of Conduct may be removed from
the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior
may be reported by opening an issue or contacting one or more of the
project maintainers.

This Code of Conduct is adapted from the [Contributor
Covenant](http://contributor-covenant.org), version 1.0.0, available
at
[http://contributor-covenant.org/version/1/0/0/](http://contributor-covenant.org/version/1/0/0/)
