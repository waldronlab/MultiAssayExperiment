
R version 3.1.2 (2014-10-31) -- "Pumpkin Helmet"
Copyright (C) 2014 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin13.4.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 
> ##Retrieve NCI-60 MS-based protein data summarized by gene in CSV file 
> ##from http://129.187.44.58:7070/NCI60/main/download
> ##Three data types were available and are returned separately in CSV files:
> ##proteomes, deep proteomes and kinomes
> 
> 
> require(affy)
Loading required package: affy
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following object is masked from ‘package:stats’:

    xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, as.vector, cbind, colnames,
    do.call, duplicated, eval, evalq, Filter, Find, get, intersect,
    is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax,
    pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce, rep.int,
    rownames, sapply, setdiff, sort, table, tapply, union, unique,
    unlist, unsplit

Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.

> require(biocMultiAssay)
Loading required package: biocMultiAssay
Loading required package: GenomicRanges
Loading required package: S4Vectors
Loading required package: stats4
Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: dplyr

Attaching package: ‘dplyr’

The following objects are masked from ‘package:GenomicRanges’:

    intersect, setdiff, union

The following object is masked from ‘package:GenomeInfoDb’:

    intersect

The following objects are masked from ‘package:IRanges’:

    collapse, desc, intersect, setdiff, slice, union

The following object is masked from ‘package:S4Vectors’:

    rename

The following objects are masked from ‘package:BiocGenerics’:

    intersect, setdiff, union

The following object is masked from ‘package:stats’:

    filter

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

> require(data.table)
Loading required package: data.table
Warning message:
In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
  there is no package called ‘data.table’
> 
> 
> ##input NCI-60 MS-based protein data
> load("nci60_prot_MS_log10.rda")
> 
> summarize <- function(x,y,z){
+  cnames <- y
+  dat <- x
+  outfile.name <- z
+ ##select protein data quantified by iBAQ approach
+ cnames <- gsub(" ",".",cnames)
+ keep.cols <- c("Fasta.headers",cnames[grepl("iBAQ",cnames)==T])
+ dat <- dat[, match(keep.cols,cnames)]
+ 
+ ids <- unlist(lapply(dat[,1],
+  function(x)(strsplit(strsplit(x,"Gene_Symbol=")[[1]][2]," ")[[1]][1])))
+ dat[,1] <- ids
+ 
+ ##remove entries without valid HGNC symbols 
+ dat <- dat[is.na(dat[,1])==F,]
+ dat <- dat[is.na(match(dat[,1],"-"))==T,]
+ ##null values correspond to undetected proteins 
+ dat[dat==0] <- NA
+ 
+ dat <- data.table(dat)
+ setkey(dat, "Fasta.headers")
+ ##average rows with identical HGNC symbols:
+ dat <- dat[, lapply(.SD, mean), by="Fasta.headers"]
+ 
+ setnames(dat, colnames(dat), gsub("Fasta.headers","id",keep.cols))
+ 
+ ##Create ExpressionSet
+ eset <- ExpressionSet(assayData=as.matrix(dat))
+ featureNames(eset) <- dat$id
+ 
+ ##Write matrix 
+ write.csv(dat, row.names=FALSE, file=outfile.name)
+ }
> 
> summarize(dat,cnames,"nci60_prot_MS_log10.csv")
Error in summarize(dat, cnames, "nci60_prot_MS_log10.csv") : 
  could not find function "data.table"
Execution halted
