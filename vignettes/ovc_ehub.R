
## ----, echo=FALSE, message=FALSE-----------------------------------------
library(biocMultiAssay)


## ------------------------------------------------------------------------
cn <- read.table(system.file("extdata/tcga_ov/raw", "all_lesions.conf_99.txt", package="biocMultiAssay"), 
                 row.names=1, header=TRUE, sep="\t", as.is=TRUE)


## ------------------------------------------------------------------------
chr <- sub(":.+", "", cn$Region.Limits)
head(chr)


## ------------------------------------------------------------------------
head(cn$Region.Limits, 1)
ranges <- sub(".+?:", "", cn$Region.Limits)
head(ranges, 1)
ranges <- sub("\\(.+", "", ranges)
head(ranges, 1)
ranges <- data.frame(t(do.call(cbind, strsplit(ranges, "-"))), stringsAsFactors=FALSE)
ranges[, 1] <- as.integer(ranges[, 1])
ranges[, 2] <- as.integer(ranges[, 2])
colnames(ranges) <- c("start", "end")
head(ranges, 1)


## ----, message=FALSE-----------------------------------------------------
library(GenomicRanges)
gr <- GRanges(seqnames=Rle(chr),
        ranges=IRanges(start=ranges$start, end=ranges$end),
        q.values=cn$q.values)
genome(gr) <- "hg19"
gr


## ------------------------------------------------------------------------
cna.dat <- cn[, -1:-8]
colnames(cna.dat) <- sub("\\.[0-9]*[ABCD].+", "", colnames(cna.dat))
cna.dat$X <- NULL
table(sapply(cna.dat, class))
cna.dat <- as.matrix(cna.dat)
class(cna.dat)
dim(cna.dat)
cna.dat[1:3, 1:3]


## ------------------------------------------------------------------------
gistic.calls <- SummarizedExperiment(assays=SimpleList(cna.dat), 
                                  rowData=gr)
colnames(gistic.calls) <- colnames(cna.dat)
gistic.calls


## ------------------------------------------------------------------------
cn <- read.table(system.file("extdata/tcga_ov/raw", "focal_input.seg.txt", 
                             package="biocMultiAssay"),
                 header=TRUE, sep="\t", as.is=TRUE)
cn$id <- make.names(sub("\\-[0-9]*[ABCD].+", "", cn[, 1]))


## ------------------------------------------------------------------------
grl <- lapply(unique(cn$id), function(x){
  cn1 <- cn[cn$id==x, ]
  gr <- GRanges(seqnames=Rle(paste0("chr", cn1$Chromosome)),
          ranges=IRanges(start=cn1$Start.bp, end=cn1$End.bp),
          Num.Markers=cn1$Num.Markers,
          Seg.CN=cn1$Seg.CN)
  genome(gr) <- "hg19"
  gr
})
names(grl) <- unique(cn$id)
grl <- GRangesList(grl)


## ------------------------------------------------------------------------
if (!file.exists("all_data_by_genes.txt.gz")) {
download.file("https://dl.dropboxusercontent.com/u/15152544/TCGA/all_data_by_genes.txt.gz",
              destfile="all_data_by_genes.txt.gz", method="wget")
}
cn <- read.table(gzfile("all_data_by_genes.txt.gz"), as.is=TRUE, row.names=1, sep="\t", header=TRUE)
colnames(cn) <- sub("\\.[0-9]*[ABCD].+", "", colnames(cn))
cn <- cn[, -1:-2]


## ------------------------------------------------------------------------
rowd <- key2GRanges(rownames(cn))
hasrowd = match(names(rowd), rownames(cn), nomatch=0)
ex = as.matrix(cn[hasrowd,])
stopifnot(identical(names(rowd), rownames(ex)))
gistic.cn <- SummarizedExperiment(assays=SimpleList(ex), 
                                  rowData=rowd)
colnames(gistic.cn) <- colnames(ex)
gistic.cn


## ------------------------------------------------------------------------
load(system.file("extdata/tcga_ov","TCGA_eset.rda",package="biocMultiAssay"))
TCGA_eset


## ------------------------------------------------------------------------
tmp <- TCGA_eset
featureNames(tmp) <- as.character(featureData(tmp)$probeset)
tmp@annotation <- "hgu133a"
TCGA_se <- exs2se(tmp)

