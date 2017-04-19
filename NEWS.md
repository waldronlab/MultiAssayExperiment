# MultiAssayExperiment 1.1.52

* `reduce` removed and broken up into `mergeReplicates` and `intersectColumns`
* Additional helper introduced: `intersectRows`

# MultiAssayExperiment 1.1.49

## New features

* `pData` deprecated in favor of `colData`
* Quick start vignette now available

## Bug fixes and minor improvements

* Fixed API function link
* Removed coercion to old `RangedRaggedAssay` class
* Improved `listToMap`

# MultiAssayExperiment 1.1.44

## Bug fixes and minor improvements

* Renamed `PrepMultiAssay` to `prepMultiAssay` (lower `p` following convention)

# MultiAssayExperiment 1.1.43

## New features

* The `MultiAssayExperiment` quickstart guide vignette added
* Deprecation of the `RangedRaggedAssay` class. Use `RaggedExperiment` package
instead.
* `reduce` function simplified
* `mapFrom` convenience argument added to the `c,MultiAssayExperiment-method`
* `assay` and `assays` methods have been revised to conform to 
`SummarizedExperiment` standards

## Bug fixes and minor improvements

* `API()` now points to the correct web document
* `ExperimentList` constructor no longer coerces `GRangesList` to `RangedRaggedAssay`
* Documentation changes: consolidate man pages for `experiments`, `experiments<-`,
`sampleMap`, and `sampleMap<-`
* removal of internal `getHits` method, simplified helper function in place
* `prepMultiAssay` helper now returns a `list` with names corresponding to the
`MultiAssayExperiment` constructor function

# MultiAssayExperiment 1.1.27

## New features

* `c` method implemented for experiments with 1:1 sample matches in `pData` rows
* `MultiAssayExperiment` show method improved
* Double bracket `[[` extracts single experiment (replacement also included)
* Internal `getHits` methods removed and refactored `subsetByRows`
* `subsetBypData` available
* `rearrange` method now supports "wide" format outputs

## Bug fixes and minor improvements

* Updates to HDF5 vignette
* More examples to documentation
* Numerous bug fixes
* `mapToList` uses the more efficient `splitAsList` function

# MultiAssayExperiment 1.1.17

## New features

* `upsetSamples` function implemented

# MultiAssayExperiment 1.1.16

## New features

* Implement `shape` argument for `rearrange` function: `wide` now available

## Bug fixes and minor improvements

* Updated vignettes: `DelayedMatrix` & PRAD `MultiAssayExperiment` object

# MultiAssayExperiment 1.1.15

## New features

* `disjoin` method for `RangedRaggedAssay` 

## Bug fixes and minor improvements

* `show` method for `RangedRaggedAssay` abbreviated. No longer summarizes data with `assay`
* Documentation changes for `reduce` and `disjoin` 

# MultiAssayExperiment 1.1.12

## New features 

* `gather`/`collect` function name changed to `rearrange`
* `clusterSex` now `clusterOn`, works with characteristic of choice

# MultiAssayExperiment 1.1.11

## New features

* Renamed `gather` function to `collect`

# MultiAssayExperiment 1.1.10

## New features 

* Double bracket method for MultiAssayExperiment available

# MultiAssayExperiment 1.1.9

## New features

* `clusterSex` function available for clustering gender from expression data

## Bug fixes and minor improvements

* Improvements to documentation

# MultiAssayExperiment 1.1.6

## New features

* Added an example `HNSC` dataset

## Bug fixes and minor improvements

* Improve documentation of `assay` method for the `RangedRaggedAssay`
* Bug fixes for `assay` method
* Removed method pollution for other Bioconductor classes
* `assay` method only shows numeric or character data

# MultiAssayExperiment 1.1.2

## New features

* `extract` method renamed to `gather`
* `gather` allows for inclusion of pData columns
* `gather` method supports common classes; creates a "tidy" DataFrame with
pData rownames, `ExperimentList` rownames, `ExperimentList` columns,
assay names, and optional pData columns

## Bug fixes and minor improvements

* Fix `assay` arguments for the `RangedRaggedAssay` method
* Subsetting by column now arranges `sampleMap` in proper order

# MultiAssayExperiment 1.1.1

## New features

* MultiAssayExperiment now in release!
* `extract` method not available for common classes - creates `tidy` data.frame
from data

## Bug fixes and minor improvements

* Documentation updated with new roxygen version

# MultiAssayExperiment 0.101.49

## New features

* Example section added to vignette for converting data frames to Bioconductor
objects

## Minor improvements

* A proper `dimnames` method added to `MultiAssayExperiment`

# MultiAssayExperiment 0.101.45

## New features

* `dimnames` method added to `RangedRaggedAssay`

## Bug fixes and minor improvements

* Improved `RangedRaggedAssay` rowname construction
* Improved `show` method for the `RangedRaggedAssay` class

# MultiAssayExperiment 0.101.44

## New features

* `$` (DollarSign) method available for `MultiAssayExperiment` to access `pData`
column

# MultiAssayExperiment 0.101.43

## Bug fixes and minor improvements

* Improved `MultiAssayExperiment` constructor now handles stray assays, colnames,
pData rownames, and `sampleMap` rows

# MultiAssayExperiment 0.101.42

## New features

* `metadata<-` set method now available for the `MultiAssayExperiment`

## Bug fixes and minor improvements

* `metadata` argument available in the `MultiAssayExperiment` constructor function
* Fix bug when subsetting for unmatched samples/colnames (drop = FALSE) in constructor

# MultiAssayExperiment 0.101.40

## New features

* Improved `MultiAssayExperiment` constructor with renamed argument "experiments"
for the `ExperimentList` or `list` input.

# MultiAssayExperiment 0.101.39

## New features

* `updateObject` method now available for old instances of the `MultiAssayExperiment`
* Users with invalid `MultiAssayExperiments` should update and re-serialize them

## Bug fixes and minor improvements

* `drop` argument now works as intended when using `List` inherited objects

# MultiAssayExperiment 0.101.38

## New features

* `complete.cases` method available for the `MultiAssayExperiment` class

# MultiAssayExperiment 0.101.37

## New features

* `sampleMap` column names renamed to __assay__ (prev. "assayname"), __primary__,
and __colname__ (prev. "assay")
* New vignete available for creating `MultiAssayExperiment` objects with TCGA data

## Bug fixes and minor improvements

* Added informative error to `MultiAssayExperiment` constructor
* Improved `show` method display
* Removed warning message when only `ExperimentList` argument provided

# MultiAssayExperiment 0.101.36

* `Elist` class renamed to `ExperimentList`
* `ExperimentList` constructor is homonymous
* `ExperimentList` accessor now called `experiments`
* `ExperimentList` replacement method is now `experiments<-`
* Updated vignettes to reflect change of names

# MultiAssayExperiment 0.101.34

* `assay` method for `RangedRaggedAssay` works on inner metadata columns now
* vignette examples available for `HDF5Array` package
* Improved outline for main vignette

# MultiAssayExperiment 0.99.202

## New Features

* `assay` method available for `RangedRaggedAssay` and other classes.
Created to obtain raw data from certain classes (see `?assay,(class),ANY-method`).

# MultiAssayExperiment 0.99.194

## New Features

* Subsetting by non-character i (#108)
* `PrepMultiAssay` helper function now available to aid in creating object (#122)
* Vignette now building (#125)
* Preliminary assay method for RangedRaggedAssay

## Bug fixes and minor improvements

* superflous `subset` function removed
* `sampleMap` uses character vectors instead of Rle
* Elist order consistent when subsetting
* mapToList preserves list order

# MultiAssayExperiment 0.99.14

* NEWS file is now live!
* Package now in `Bioc-devel`!
* More to come!

## New Features

* Replacement method for `colnames` now available for the `RangedRaggedAssay` class.

## Bug fixes and minor improvements

* Documentation additions
