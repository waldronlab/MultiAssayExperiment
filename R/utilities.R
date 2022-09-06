.msg <-
    function(..., width=getOption("width"), wrap. = TRUE)
    ## Use this helper to format all error / warning / message text
{
    txt <- paste(...)
    if (wrap.) {
        txt <- strwrap(txt, width=width, exdent=2)
        paste(txt, collapse="\n")
    } else {
        txt
    }
}

.warning <-
    function(..., call.=FALSE, immediate.=FALSE)
{
    warning(.msg(...), call.=call., immediate.=immediate.)
    invisible(TRUE)
}

## stats::setNames
.setNames <- function(object = nm, nm) {
    names(object) <- nm
    object
}

.normarg_verbose <- function(verbose) {
    if (!BiocBaseUtils::isTRUEorFALSE(verbose, na.ok = TRUE))
        stop(S4Vectors::wmsg("'verbose' must be FALSE, TRUE, or NA"))
    if (is.na(verbose))
        verbose <-
            getOption("DelayedArray.verbose.block.processing", default = FALSE)
    verbose
}
