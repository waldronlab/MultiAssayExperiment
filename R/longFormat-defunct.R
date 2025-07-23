#' @name longFormat-defunct
#'
#' @title Defunct longFormat method
#'
#' @description The `longFormat` method is defunct and will be removed in the
#'   next release. Please use the `longForm` method instead.
#'
#' @details The `longFormat` "ANY" class method, works with classes such as
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} and
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class] as
#' well as `matrix` to provide a consistent long and skinny
#' [`DataFrame`][S4Vectors::DataFrame-class].
#'
#' @section longFormat:
#' The 'longFormat' method takes data from the [`ExperimentList`]
#' in a `MultiAssayExperiment` and returns a uniform
#' `DataFrame`. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' columns of the MultiAssayExperiment colData named by `colDataCols` character
#' vector argument. (`MultiAssayExperiment` method only). The `i` argument
#' allows the user to specify the assay value for the
#' `SummarizedExperiment` assay function's `i` argument.
#'
#' @param object Any supported class object
#'
#' @param colDataCols A `character`, `logical`, or `numeric`
#'     index for `colData` columns to be included
#'
#' @param i longFormat: The i-th assay in
#'     `SummarizedExperiment`-like objects. A vector input is
#'     supported in the case that the `SummarizedExperiment` object(s) has more
#'     than one assay (default 1L),
#'     renameColname: Either a `numeric` or `character` index
#'     indicating the assay whose colnames are to be renamed
#'
#' @param ... Additional arguments. See details for more information.
#'
#' @importFrom BiocBaseUtils lifeCycle
#'
#' @return A `DataFrame` with selected `colDataCols`, if any.
#'
#' @aliases longFormat
#' @export
setGeneric(
    "longFormat",
    function(object, colDataCols = NULL, i = 1L, ...)
        standardGeneric("longFormat")
)

#' @rdname longFormat-defunct
#' @exportMethod longFormat
setMethod(
    "longFormat", "MultiAssayExperiment",
    function(object, colDataCols = NULL, i = 1L, ...) {
        lifeCycle(
            "longForm", package = "MultiAssayExperiment", title = "longFormat",
            cycle = "defunct"
        )
    }
)

#' @rdname longFormat-defunct
#' @exportMethod longFormat
setMethod("longFormat", "ANY", function(object, colDataCols, i = 1L, ...) {
    lifeCycle(
        "longForm", package = "MultiAssayExperiment", title = "longFormat",
        cycle = "defunct"
    )
})

#' @rdname longFormat-defunct
#' @exportMethod longFormat
setMethod(
    "longFormat", "ExperimentList",
    function(object, colDataCols, i = 1L, ...) {
        lifeCycle(
            "longForm", package = "MultiAssayExperiment", title = "longFormat",
            cycle = "defunct"
        )
    }
)
