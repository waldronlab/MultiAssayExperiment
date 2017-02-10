#' Head-Neck Squamous Cell Carcinoma (HNSC)
#' A \code{\link{MultiAssayExperiment}} object containing a random sample
#' of 20 patients with copy number data. For example purposes only.
#'
#' @format A \code{MultiAssayExperiment} with 2 experiments, each containing
#' 400 rows and 20 columns:
#' \describe{
#'     \item{CNASNP}{Copy Number Alteration}
#'     \item{CNVSNP}{Copy Number Variation}
#' }
#' @source The Cancer Genome Atlas
"hnscSample"

## Script for downloading and extracting HNSC data
system(paste("aws s3 cp s3://multiassayexperiments/hnscMAEO.rds",
       "/home/$USER/Downloads/"))
hnsc <- readRDS(file.path(Sys.getenv("HOME"), "Downloads", "hnscMAEO.rds"))
hnsc <- hnsc[1:20, 1:11, c("CNASNP", "CNVSNP")]
save(hnsc, file = "inst/extdata/hnscSample.rda", compress = "bzip2")
