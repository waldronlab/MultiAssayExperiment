library(rhdf5)
library(HDF5Array)

exMatrix <- matrix(rnorm(10e4), ncol = 20)
rownames(exMatrix) <- paste0("GENE", seq_len(nrow(exMatrix)))
colnames(exMatrix) <- paste0("SampleID", seq_len(ncol(exMatrix)))

format(object.size(exMatrix), unit = "Mb")
head(exMatrix)

h5save(exMatrix,
       file =
         "~/Documents/GitHub/MultiAssayExperiment/inst/extdata/exMatrix.h5")
H5close()

fpath <- "~/Documents/GitHub/MultiAssayExperiment/inst/extdata/"
DelayedArray(HDF5Dataset(file = file.path(fpath, "exMatrix.h5"),
                         name = "exMatrix"))

# file.remove(file.path(fpath, "exMatrix.h5"))

