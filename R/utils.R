
allidsfound <- sort(Reduce(union, lapply(lovhub@elist, sampleNames)))
if(any(grepl("\\.", allidsfound))){
 	allidsfound <- gsub("\\.", "-", allidsfound)	
}
# Create a matrix map 
idmap <- matrix(NA, ncol = length(lovhub@elist), nrow = length(allidsfound))
rownames(idmap) <- allidsfound
colnames(idmap) <- names(lovhub@elist)

samps <- lapply(lovhub@elist, sampleNames)
samps <- lapply(samps, FUN = function(code) {gsub("\\.", "-", code)})	

for(i in seq_along(samps)){
idmap[, i] <- samps[[i]][match(rownames(idmap), samps[[i]])]
}


Reduce(intersect, lapply(lovhub@elist, sampleNames)) 
Reduce(union, lapply(lovhub@elist, sampleNames))

# Create a list map
idmap2 <- lapply(seq_along(lovhub@elist), 
	FUN = function(experiment) { data.frame(cbind(barcodeID = allidsfound, sampleID = lovhub@elist[[experiment]]$alt_sample_name[match(allidsfound, samps[[experiment]])] )) } )
names(idmap2) <- names(lovhub@elist)

idmap3 <- lapply(seq_along(ids), 
FUN = function(expt) { data.frame(cbind(barcodeID = allidsfound, sampleID = ids[[expt]][match(bcIDR(allidsfound),bcIDR(ids[[expt]]))] )) } ) 
names(idmap3) <- names(ids)

test1 <- samps[[1]]
test2 <- lovhub@elist[[1]]$alt_sample_name
masterpheno <- pData(lovhub@elist[[1]])

findIDs <- function(ehub){
	sampnames <- lapply(ehub@elist, FUN = function(expt) { 
						if(length(pData(expt)) == 0){
							return(sampleNames(expt))
						}else{
							return(expt$alt_sample_name)
	}}
	)
	if(any(rapply(sampnames, f = function(bcode) { grepl("\\.", bcode) } ) )){
		inx <- which(lapply(sampnames, function(expt) { which(any(grepl("\\.", expt))) } ) == 1)
		sampnames[inx] <- lapply(sampnames[inx], FUN = function(expt) { gsub("\\.", "-", expt) } )
	}
	return(sampnames)
}

bb <- findIDs(lovhub)
library(RTCGAToolbox)
ids <- lapply(bb, function(x) na.omit(x))

lapply(idmap3,FUN =  function(expt) dim(na.omit(expt)))

apply(idmap, 2, function(x) { x %in% bcIDR(rownames(masterpheno)) } )

.validMap <- function(mapobj) {


return(logi)
}


