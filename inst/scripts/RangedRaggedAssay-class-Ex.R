## Load sample MultiAssayExperiment Copy Number data
dataFile <- system.file("extdata", "hnscSample.rds",
    package = "MultiAssayExperiment")

hnscSample <- readRDS(dataFile)

## assay method for a RangedRaggedAssay
assay(hnscSample[[1L]], mcolname = "Segment_Mean")[1:5, 1:3]

hnscSample[[2]] <- disjoin(hnscSample[[2L]])

matrices <- assay(hnscSample, mcolname = "Segment_Mean")

length(matrices)
class(matrices)
