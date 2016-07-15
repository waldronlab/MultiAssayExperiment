# MultiAssayExperiment 0.101.39

## New features

* `updateObject` method now available for old instances of the `MultiAssayExperiment`
* Users with invalid `MultiAssayExperiments` should update and re-serialize them

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
