### ==============================================
### multiAssayExperiment object
### ----------------------------------------------


#' An integrative MultiAssay class for experiment data
#' 
#' @slot elist A list of data across different types of assays 
#' @slot masterPheno A data.frame of all clinical data available across experiments
#' @slot sampleMap A list of translatable identifiers of samples and participants
#' @slot metadata Additional data describing the \code{\linkS4class{MultiAssayExperiment}} class 
setClass("MultiAssayExperiment", representation(
												elist="Assays", 
												masterPheno = "data.frame",
												sampleMap = "list", 
												metadata = "ANY"), 
		 prototype(
				   elist = SummarizedExperiment::Assays()
				   )
		 )

.validMultiAssayExperiment <- function(object){
c(.checkMasterPheno(object), 
.checkSampleMap(object),
.checkElist(object))
}

setValidity2("MultiAssayExperiment", .validMultiAssayExperiment)

#' Show method for MultiAssayExperiment class
#' 
#' @param object A \code{\linkS4class{MultiAssayExperiment}} 
#' @return Returns a list of contents for the MultiAssayExperiment
# setMethod("show", "MultiAssayExperiment", function(object) {
# 		  objdim <- lapply(seq_along(object@elist), FUN = function(j, expt) {	
# 						   dd <- matrix(NA, nrow = length(expt), ncol = 2)
# 						   if(any(is(expt[j], "data.frame"), is(expt[j], "matrix"))){
# 							   dimmat <- matrix(c(dim(expt[j])[1], dim(expt[j])[2]), ncol = 2)
# 							   colnames(dimmat) <- c("Features", "Samples")
# 							   dd <- rbind(dd, dimmat)
# 						   }
# 	} , expt = object@elist)
# print(objdim)
# })


### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))
setMethod("sampleMap", "MultiAssayExperiment", function(x)
getElement(x, "sampleMap"))

setGeneric("elist", function(x) standardGeneric("elist"))
setMethod("elist", "MultiAssayExperiment", function(x)
getElement(x, "elist"))

setGeneric("masterPheno", function(x) standardGeneric("masterPheno"))
setMethod("masterPheno", "MultiAssayExperiment", function(x)
getElement(x, "masterPheno"))


