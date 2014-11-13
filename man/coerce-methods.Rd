\name{coerce-methods}
\docType{methods}
\alias{coerce}
\alias{exs2se,ExpressionSet-method}
\alias{coerce,ExpressionSet,exs2se-method}
\title{ convert ExpressionSet to SummarizedExperiment }
\description{
convert ExpressionSet to SummarizedExperiment
}
\section{Methods}{
\describe{

\item{\code{SummarizedExperiment(assays = "ExpressionSet", ...}}{
There are many ways of establishing mappable feature names in an
ExpressionSet instance.  This constructor strives for reasonable generality.
Additional arguments that can be supplied are
\preformatted{
     assayname = "exprs",    # for naming SimpleList element in SummarizedExperiment
     fngetter = 
           function(z) rownames(exprs(z)),   # extract usable feature names
     annDbGetter = 
          function(z) { 
              clnanno = sub(".db", "", annotation(z)) 
              stopifnot(require(paste0(annotation(z), ".db"), character.only=TRUE) ) 
              get(paste0(annotation(z), ".db"))   # obtain resource for mapping feature names to coordinates
              }, 
     probekeytype = "PROBEID",   # chipDb field to use
     duphandler = function(z) {    # action to take to identify and process duplicated features
          if (any(isd <- duplicated(z[,"PROBEID"]))) 
              return(z[!isd,,drop=FALSE]) 
          z 
          }, 
     signIsStrand = TRUE,   # verify that signs of addresses define strand
     ucsdChrnames = TRUE    # prefix 'chr' to chromosome token
  }  

}

}}
\examples{
library(Biobase)
data(sample.ExpressionSet)
exs2se(sample.ExpressionSet)
}
\keyword{methods}
