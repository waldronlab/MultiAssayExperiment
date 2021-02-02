# Changes in version 1.18.0

## New features

* `renameColname` and `renamePrimary` provide renaming facilities for column
names in experiments and `rownames` in the `colData`, respectively
* Users can now rename some or all the column names in experiments using
`colnames(mae) <- value`
* When replacing `colData` or `experiments` (including `[[<-`), new `rownames`
and `colnames` (respectively) are checked against existing values and an error
is given when none match
* `splitAssay` allows users to separate / split columns across assays
* `makeHitList` is a facilitator function to create the logical lists that
are required as input to `splitAssay`

## Bug fixes and minor improvements

* Updated the constructor function to auto-populate `rownames` in `colData`
when it is missing (@LTLA, #287)
* The metadata now includes names of dropped experiments

# Changes in version 1.16.0

## New features

* Coercion methods from `list`/`List` to `MultiAssayExperiment` method now
available.

## Bug fixes and minor improvements

* Provide more details in documentation for `mergeReplicates`
* Improved documentation for accessor function return values, helper function
examples (@llrs, #281)
* Fixed bug when using `longFormat` with character assay matrices
(@jonocarroll, #282)

# Changes in version 1.14.0

## New features

* `exportClass` creates a number of `.csv` data files for exporting data
* Allow vector input `i` for selecting assays in `longFormat` (@lgatto, #266)
* Updates to 'Using `MultiAssayExperiment` with `DelayedMatrix`' vignette

## Bug fixes and minor improvements

* Warn when `colData` rownames and `ExperimentList` colnames are empty
(@LTLA #262)
* Add informative error message for `ExperimentList` (@lgatto, #265)
* Informative warning when dropping `ExperimentList` element columns (@lwaldron)
* Fixes to constructor functions, `MultiAssayExperiment` and
`MatchedAssayExperiment` (@lgatto, #267 #268, @lwaldron)
* Add warning when `j` in `mae[i, j, k]` is longer than `colData` rows
* Strict argument matching between generic and methods
* Updates due to `class(matrix())`
* `UpsetSamples` more robust to differences in names between split `sampleMap`
and `names(ExperimentList)` (@jonocarroll, #269)
* Refactored and improved `UpsetSamples`
* `ExperimentList` propagation of `mcols` and `metadata` (@vobencha, #270)
* Enforcement of `validObject` with replacement methods `colData` and
`sampleMap` (@vobencha, #271)

# Changes in version 1.12.0

## Bug fixes and minor improvements

* Improvements to the main vignette, `MultiAssayExperiment` class
schematic now included (@mtmorgan, #261)
* Updated documentation for the `upsetSamples` function
* Update code to use `splitAsList` from `S4Vectors` (@hpages)
* Fixed bug with metadata disappearing from `ExperimentList` when replacing it
inside a `MultiAssayExperiment` object (@lawremi, #259)
* Fixed the formatting of the NEWS file

# Changes in version 1.10.0

## New features

* `getWithColData` now allows easy extraction of experiments
(such as `SummarizedExperiment`) with associated `colData` data
* Single bracket replace method implemented for `MultiAssayExperiment` assays

## Bug fixes and minor improvements

* `isEmpty` method fixed for `ExperimentList`s containing matrices
* `MultiAssayExperiment` now inherits from the standard `Annotated` virtual
class
* `c` method better distinguishes between `list` and `ExperimentList` inputs
* Improvements on `.getHits` internal method for obtaining correct queries on
row metadata
* Subsetting mechanism tweaked to do nothing when subsetting by `NULL` rows
compared to empty rows (i.e., `character(0L)`)
* Improved README.md

# Changes in version 1.8.0

## New features

* The single bracket replacement method `[<-` added to support assignment
of assay datasets
* Users can now rename experiments in a MultiAssayExperiment with
`names(x) <- value`
* `replicated` and `mergeReplicates` functions have been refactored and
improved
* combining MultiAssayExperiments now possible with `c` function

## Bug fixes and minor improvements

* `wideFormat` function improvements and bug fixes with name indicator
subsetting
* `BiocGenerics:::replaceSlots` used instead of replace methods
* Added tests for `anyReplicated`, `c`, and `names<-` functions
* Unit tests added for replacement method testing

# Changes in version 1.7.14

## New features

* Subsetting `MultiAssayExperiment` by a `list` or `List` type class now
returns experiments in the input order for rows, columns, and assays

## Bug fixes and minor improvements

* Renamed objects in examples for brevity and descriptiveness
* Updated `importFrom` directives
* Internal sample names check now only works on non-empty `colnames`
* Various documentation improvements
* `listToMap` re-written for efficiency
* Various improvements to subsetting mechanism
* `subsetByAssay` bug fixed when using an integer index (@vjcitn, #)

# Changes in version 5.108

## New features

* `DataFrame` now exported for users (@DarioS, #242)
* `c` is smarter at matching `colnames` with `primary` names and creating a
`sampleMap`

## Bug fixes and minor improvements

* Added an `isEmpty` method for `ExperimentList` to account for an empty matrix
* Documentation improvements to `MultiAssayExperiment-class` and
`MultiAssayExperiment-helpers`
* `c` internals improved

# Changes in version 5.102

## New features

* The `MatchedAssayExperiment` constructor function now works either a
single `MultiAssayExperiment` or the essential components of one.

# Changes in version 5.101

## New features

* Renamed `duplicated` function to `replicated`
* Added coercion functions from `List` and `list` to `ExperimentList`
* Improve speed of reshape functions from previous change (`wideFormat`)

## Bug fixes and minor improvements

* Explicitly check for `DataFrame` in `ExperimentList`
* Fixed documentation warnings for inexact links
* Fix subsetting order in bracket method (`[`)
* Minor vignette changes
* Supply a collapse character for `wideFormat` column names

# Changes in version 1.5.65

## Bug fixes and minor improvements

* `upsetSamples` does not munge experiment names with special characters when
`check.names = FALSE` (by default keeps hyphens, underscores, etc.). A
`nameFilter` functional argument allows operations such as `substr` on the
experiment names. (@vjcitn, #231)

# Changes in version 1.5.64

## New features

* Remove `clusterOn` function and move to `Bioconductor/MultiOmicQC` package
on GitHub

# Changes in version 1.5.63

## New features

* duplicated has been deprecated, use `replicated` and `anyReplicated`
* removed dependencies on `tidyr` and `reshape2`

## Bug fixes and minor improvements

* Updates to `prepMultiAssay`
* Enhancements to the main vignette
* New format for NEWS section

# Changes in version 1.5.38

## New features

* Moved the API shiny function to waldronlab/MultiAssayShiny package
* Reduced imports (removed shinydashboard and shiny)
* Method requirement checks for classes are practical using `try()`
* Deprecated methods removed: `pData`
* Deprecated class removed: `RangedRaggedAssay`
* Assay-selective subsetting implemented via `list`/`List` class subsettors

## Bug fixes and minor improvements

* updated `duplicated` function now returns FALSE for non-duplicated samples
* Improved `ExperimentList` constructor now handles multiple `key = value`
entries
* Removed `updateObject` before giving warning
* Removed old RTCGAToolbox example vignette
* Official manuscript citation added
* Improved examples (removed `ExpressionSet` legacy objects)
* Improved test scripts

# Changes in version 1.1.37

## New features

* `MatchedAssayExperiment` subclass added for matched samples in all assays
* Supply mini ACC dataset `data(miniACC)`
* Provide reference table for methods in package, see vignettes
* Merge with GitHub development version

## Bug fixes and minor improvements

* ensure assay column in sampleMap is a factor
* rearrange long DataFrame correctly
* remove support for RangedRaggedAssay - deprecate
* `drop = FALSE` in single column subset of colData
* default sampleMap representation as empty DataFrame with colnames
* added combine `c` vignette section for adding experiments to an existing

# Changes in version 1.1.59

## New features

* `rearrange` is now broken up into `longFormat` and `wideFormat` functions.
* Helper functions now have a dedicated man page, see: `?'MultiAssayExperiment-helpers'`

## Bug fixes and minor improvements

* A subset can affect the order of rows (previously it didn't)
* `rownames` are exclusively used to create `longFormat` `DataFrame`s
* The `longFormat,ExperimentList-method` now returns a long DataFrame
* Minor improvements to tests

# Changes in version 1.1.52

## New features

* `reduce` removed and broken up into `mergeReplicates` and `intersectColumns`
* Additional helper introduced: `intersectRows`

# Changes in version 1.1.49

## New features

* `pData` deprecated in favor of `colData`
* Quick start vignette now available

## Bug fixes and minor improvements

* Fixed API function link
* Removed coercion to old `RangedRaggedAssay` class
* Improved `listToMap`

# Changes in version 1.1.44

## Bug fixes and minor improvements

* Renamed `PrepMultiAssay` to `prepMultiAssay` (lower `p` following convention)

# Changes in version 1.1.43

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

# Changes in version 1.1.27

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

# Changes in version 1.1.17

## New features

* `upsetSamples` function implemented

# Changes in version 1.1.16

## New features

* Implement `shape` argument for `rearrange` function: `wide` now available

## Bug fixes and minor improvements

* Updated vignettes: `DelayedMatrix` & PRAD `MultiAssayExperiment` object

# Changes in version 1.1.15

## New features

* `disjoin` method for `RangedRaggedAssay`

## Bug fixes and minor improvements

* `show` method for `RangedRaggedAssay` abbreviated. No longer summarizes data with `assay`
* Documentation changes for `reduce` and `disjoin`

# Changes in version 1.1.12

## New features

* `gather`/`collect` function name changed to `rearrange`
* `clusterSex` now `clusterOn`, works with characteristic of choice

# Changes in version 1.1.11

## New features

* Renamed `gather` function to `collect`

# Changes in version 1.1.10

## New features

* Double bracket method for MultiAssayExperiment available

# Changes in version 1.1.9

## New features

* `clusterSex` function available for clustering gender from expression data

## Bug fixes and minor improvements

* Improvements to documentation

# Changes in version 1.1.6

## New features

* Added an example `HNSC` dataset

## Bug fixes and minor improvements

* Improve documentation of `assay` method for the `RangedRaggedAssay`
* Bug fixes for `assay` method
* Removed method pollution for other Bioconductor classes
* `assay` method only shows numeric or character data

# Changes in version 1.1.2

## New features

* `extract` method renamed to `gather`
* `gather` allows for inclusion of pData columns
* `gather` method supports common classes; creates a "tidy" DataFrame with
pData rownames, `ExperimentList` rownames, `ExperimentList` columns,
assay names, and optional pData columns

## Bug fixes and minor improvements

* Fix `assay` arguments for the `RangedRaggedAssay` method
* Subsetting by column now arranges `sampleMap` in proper order

# Changes in version 1.1.1

## New features

* MultiAssayExperiment now in release!
* `extract` method not available for common classes - creates `tidy` data.frame
from data

## Bug fixes and minor improvements

* Documentation updated with new roxygen version

# Changes in version 101.49

## New features

* Example section added to vignette for converting data frames to Bioconductor
objects

### Minor improvements

* A proper `dimnames` method added to `MultiAssayExperiment`

# Changes in version 101.45

## New features

* `dimnames` method added to `RangedRaggedAssay`

## Bug fixes and minor improvements

* Improved `RangedRaggedAssay` rowname construction
* Improved `show` method for the `RangedRaggedAssay` class

# Changes in version 101.44

## New features

* `$` (DollarSign) method available for `MultiAssayExperiment` to access `pData`
column

# Changes in version 101.43

## Bug fixes and minor improvements

* Improved `MultiAssayExperiment` constructor now handles stray assays, colnames,
pData rownames, and `sampleMap` rows

# Changes in version 101.42

## New features

* `metadata<-` set method now available for the `MultiAssayExperiment`

## Bug fixes and minor improvements

* `metadata` argument available in the `MultiAssayExperiment` constructor function
* Fix bug when subsetting for unmatched samples/colnames (drop = FALSE) in constructor

# Changes in version 101.40

## New features

* Improved `MultiAssayExperiment` constructor with renamed argument "experiments"
for the `ExperimentList` or `list` input.

# Changes in version 101.39

## New features

* `updateObject` method now available for old instances of the `MultiAssayExperiment`
* Users with invalid `MultiAssayExperiments` should update and re-serialize them

## Bug fixes and minor improvements

* `drop` argument now works as intended when using `List` inherited objects

# Changes in version 101.38

## New features

* `complete.cases` method available for the `MultiAssayExperiment` class

# Changes in version 101.37

## New features

* `sampleMap` column names renamed to __assay__ (prev. "assayname"), __primary__,
and __colname__ (prev. "assay")
* New vignete available for creating `MultiAssayExperiment` objects with TCGA data

## Bug fixes and minor improvements

* Added informative error to `MultiAssayExperiment` constructor
* Improved `show` method display
* Removed warning message when only `ExperimentList` argument provided

# Changes in version 101.36

* `Elist` class renamed to `ExperimentList`
* `ExperimentList` constructor is homonymous
* `ExperimentList` accessor now called `experiments`
* `ExperimentList` replacement method is now `experiments<-`
* Updated vignettes to reflect change of names

# Changes in version 101.34

* `assay` method for `RangedRaggedAssay` works on inner metadata columns now
* vignette examples available for `HDF5Array` package
* Improved outline for main vignette

# Changes in version 9.202

## New features

* `assay` method available for `RangedRaggedAssay` and other classes.
Created to obtain raw data from certain classes (see `?assay,(class),ANY-method`).

# Changes in version 9.194

## New features

* Subsetting by non-character i (#108)
* `PrepMultiAssay` helper function now available to aid in creating object (#122)
* Vignette now building (#125)
* Preliminary assay method for RangedRaggedAssay

## Bug fixes and minor improvements

* superflous `subset` function removed
* `sampleMap` uses character vectors instead of Rle
* Elist order consistent when subsetting
* mapToList preserves list order

# Changes in version 0.99.14

* NEWS file is now live!
* Package now in `Bioc-devel`!
* More to come!

## New features

* Replacement method for `colnames` now available for the `RangedRaggedAssay` class.

## Bug fixes and minor improvements

* Documentation additions
