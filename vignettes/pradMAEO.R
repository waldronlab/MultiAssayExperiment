#Title: pradMAEO.R
#Author: Lucas Schiffer
#Date: January 7, 2016

#install bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()

#install RTCGAToolbox
BiocInstaller::biocLite("LiNk-NY/RTCGAToolbox")
library(RTCGAToolbox)

#install biocMultiAssay
BiocInstaller::biocLite("LiNk-NY/biocMultiAssay")
library(biocMultiAssay)

#background reading
?RTCGAToolbox

#cohort list
getFirehoseDatasets()

#obtain some of the demographic data
#prostate adenocarcinoma = PRAD
prad <- getFirehoseData("PRAD", runDate = "20150402")
bcPRAD <- extract(prad, type = NULL, clinical = TRUE)
bcPRAD

#begin formating data into a MultiAssayExperiment object
suppressPackageStartupMessages(library(biocMultiAssay))
empty <- MultiAssayExperiment()
empty