if(all(grepl("[A-Z].[0-9]{2}.[0-9]{4}", sapply(object@elist, sampleNames))))

allidsfound <- sort(Reduce(union, lapply(lovhub@elist, sampleNames)))
idmap <- matrix(NA, ncol = length(lovhub@elist), nrow = length(allidsfound))
rownames(idmap) <- allidsfound

samps <- lapply(lovhub@elist, sampleNames)

for(i in seq_along(samps)){
idmap[, i] <- samps[[i]][match(rownames(idmap), samps[[i]])]
}

Reduce(intersect, lapply(lovhub@elist, sampleNames)) 
Reduce(union, lapply(lovhub@elist, sampleNames))

##------------------
## TODO: Create a mapping list for each assay (col1 = patientID, col2 = sampleName)
##------------------



library(RTCGAToolbox)
bcIDR(samps[[1]], position = 3)


