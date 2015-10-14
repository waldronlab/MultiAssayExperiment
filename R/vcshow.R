
vcshow = function (object) 
{
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    cat("A", o_class, "object with", o_len, ifelse(o_len == 1L, 
        "experiment", "experiments named:\n"))
    classes = sapply(elist(object), class)
    names(o_names) = classes
    cat(selectSome(o_names), sep = ", ")
    cat("\n")
    cat("use elist() to obtain list of experiment instances,\n")
    cat("   masterPheno() for phenotype data frame, sampleMap() for sample\n")
    cat("   availability list.\n")
}

