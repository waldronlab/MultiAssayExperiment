library(RTCGAToolbox)
library(pipeR)

lapply(object@elist, FUN = function(assay) {colnames(assay)} ) %>>% rapply(f = intersect, x = .)

lapply(locomm@elist, FUN = function(assay) {colnames(assay)} ) %>>% Map(f = intersect, x = .)
lapply(lovhub@elist, FUN = function(assay) {colnames(assay)} ) %>>% Map(f = intersect, x = .) %>>% lapply(length) 

lapply(lovhub@elist, FUN = function(assay) {colnames(assay)} ) %>>% 
	lapply(seq_along(.), FUN = function(samples, i) { intersect(samples[i], samples[i+1]) }, samples = .)
 
samps <- lapply(lovhub@elist, FUN = function(assay) {colnames(assay)} ) 
a <- c()
for(i in seq_len(length(samps)-1)){
 if(i == 1){
 a <- cbind(a, intersect(samps[i], samps[i+1]))
} else { a <- intersect(a, samps[[i+1]])
} }
