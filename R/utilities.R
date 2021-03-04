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

