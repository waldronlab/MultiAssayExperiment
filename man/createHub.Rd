\name{createHub}
\alias{createHub}
\title{
Create a loadedHub from objects in memory.
}
\description{
As simple a function as possible for converting omic data objects in
memory (e.g. ExpressionSet and SummarizedExperiment) to a loadedHub object.
}
\usage{
createHub(masterpheno, objlist, drop = FALSE, samplemaps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{masterpheno}{
A data.frame that should contain phenotype data for *all* samples in any assay of
the experiment.  rownames should be sample IDs and columns are phenotypes.
}
  \item{objlist}{
A list of data objects (one for each assay type used) with informative
names, e.g. list(mrna=eset1, copynumber=summarizedexpt2)
}
  \item{drop}{
Logical, if TRUE, drop samples from individual assays that are not
mapped to masterpheno.  If FALSE, stop with an error if any sample is
not mapped to masterpheno.
}
  \item{samplemaps}{
In the future, will provide a way to have different sample identifiers
for each assay.
}
}
\details{
A helpful function for creating loadedHub objects.
}
\value{
An object of class loadedHub.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
multiassay SIG
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
if(interactive()){
ov = dir(system.file("extdata/tcga_ov",
   package="biocMultiAssay"), full=TRUE, pattern="\\\\.rda$")
drop = grep("pheno", ov)
if (length(drop)>0) {
  pdpath=ov[drop]
  ov=ov[-drop]
  }
#
# informal labels for constituents
#
tags = c("ov CNV gistic", "ov methy 450k", "ov affy", "ov agilent", "ov mirna", "ov RNA-seq")
ovlist <- lapply(ov, function(x) get(load(x)))
names(ovlist) <- tags
lovhub2 <- createHub(masterpheno=pData(ovlist[[2]]), objlist=ovlist, drop=TRUE)
lovhub2
object.size(lovhub2)
}
}
