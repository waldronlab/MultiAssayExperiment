\name{key2GRanges}
\alias{key2GRanges}
\title{key2GRanges}
\description{Convert an identifier (such as gene symbols) to a GRanges object, using
a Bioconductor .db annotation database.}
\usage{
key2GRanges(x, db = "org.Hs.eg.db", keytype = "SYMBOL", duphandler = function(z) {
    if (any(isd <- duplicated(z[, keytype]))) 
        return(z[!isd, , drop = FALSE])
    z
}
 signIsStrand = TRUE, ucsdChrnames = TRUE)
}
\arguments{
  \item{x}{character vector of the keys to convert}
  \item{db}{Name of the Bioconductor .db annotation package to use}
  \item{keytype}{Must be a valid key type for the .db package, e.g. see keytypes(org.Hs.eg.db)}
  \item{duphandler}{
}
  \item{signIsStrand}{If TRUE, negative values are annotated as being on the - strand. }
  \item{ucsdChrnames}{If TRUE, use chr1 etc.}
}

\value{Returns a GRanges object}

\author{multiAssay SIG}




\examples{
key2GRanges(c("TP53", "BRCA1"))
key2GRanges(c("1", "10"), keytype="ENTREZID")
}
