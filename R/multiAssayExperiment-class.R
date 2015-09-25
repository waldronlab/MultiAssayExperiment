### ==============================================
### multiAssayExperiment object
### ----------------------------------------------


#' An integrative MultiAssay class for experiment data
#' 
#' @slot elist A list of data across different types of assays 
#' @slot masterPheno A data.frame of all clinical data available across experiments
#' @slot sampleMap A list of translatable identifiers of samples and participants
#' @slot metadata Additional data describing the \code{\linkS4class{MultiAssayExperiment}} class 
setClass("MultiAssayExperiment", representation(elist="list", masterPheno = "data.frame",
	sampleMap = "list", metadata = "ANY"))

### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))
setMethod("sampleMap", "MultiAssayExperiment", function(x)
getElement(x, "sampleMap"))

setMethod("elist", "MultiAssayExperiment", function(x)
getElement(x, "elist"))

setMethod("masterPheno", "MultiAssayExperiment", function(x)
getElement(x, "masterPheno"))


