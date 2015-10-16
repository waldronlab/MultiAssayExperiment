---
title: "biocMultiAssayToyExample"
author: "Marcel Ramos, Levi Waldron"
date: "September 4, 2015"
output: html_document
---

NOTE: This README.md is a compiled version of the biocMultiAssayToyExample.Rmd vignette.

# Generate toy data

In this example we have 4 patients, and a bit of metadata on them:

```r
masPheno <- data.frame(sex=c("M", "F", "M", "F"),
						  age=38:41,
						  row.names=c("Jack", "Jill", "Bob", "Barbara"))
masPheno
```

```
##         sex age
## Jack      M  38
## Jill      F  39
## Bob       M  40
## Barbara   F  41
```

We have three matrix-like datasets.  First let's say expression data:


```r
library(Biobase)
(arraydat <- matrix(seq(101, 108), ncol=4, dimnames=list(c("ENST00000294241", "ENST00000355076"), c("array1", "array2", "array3", "array4"))))
```

```
##                 array1 array2 array3 array4
## ENST00000294241    101    103    105    107
## ENST00000355076    102    104    106    108
```

```r
arraypdat <- as(data.frame(slope53=rnorm(4), row.names=c("array1", "array2", "array3", "array4")), "AnnotatedDataFrame")
exprdat <- ExpressionSet(assayData=arraydat, phenoData=arraypdat)
exprdat
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 2 features, 4 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2 array3 array4
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:
```

The following map matches masterPheno sample names to exprdata sample names.  Note that row orders aren't initially matched up.

```r
(exprmap <- data.frame(master=rownames(masPheno)[c(1, 2, 4, 3)], assay=c("array1", "array2", "array3", "array4")))
```

```
##    master  assay
## 1    Jack array1
## 2    Jill array2
## 3 Barbara array3
## 4     Bob array4
```

Now methylation data.  It uses gene identifiers also, but measures a partially overlapping set of genes.  For fun, let's store this as a simple matrix. Also, it contains a replicate for one of the patients.


```r
(methyldat <- matrix(1:10, ncol=5, 
                     dimnames=list(c("ENST00000355076", "ENST00000383706"),
                                   c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5"))))
```

```
##                 methyl1 methyl2 methyl3 methyl4 methyl5
## ENST00000355076       1       3       5       7       9
## ENST00000383706       2       4       6       8      10
```

The following map matches masterPheno sample names to methyldat sample names.


```r
(methylmap <- data.frame(master = c("Jack", "Jack", "Jill", "Barbara", "Bob"),
                        assay = c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5")))
```

```
##    master   assay
## 1    Jack methyl1
## 2    Jack methyl2
## 3    Jill methyl3
## 4 Barbara methyl4
## 5     Bob methyl5
```

Now we have a microRNA platform, which has no common identifiers.  It is also missing data for Jill.  Just for fun, let's use the same sample naming convention as we did for arrays.


```r
(microdat <- matrix(201:212, ncol=3, 
                    dimnames=list(c("hsa-miR-21", "hsa-miR-191", "hsa-miR-148a", "hsa-miR148b"), 
                                  c("micro1", "micro2", "micro3"))))
```

```
##              micro1 micro2 micro3
## hsa-miR-21      201    205    209
## hsa-miR-191     202    206    210
## hsa-miR-148a    203    207    211
## hsa-miR148b     204    208    212
```

And the following map matches masterPheno sample names to microdat sample names.

```r
(micromap <- data.frame(master = c("Jack", "Barbara", "Bob"),
                        assay = c("micro1", "micro2", "micro3")))
```

```
##    master  assay
## 1    Jack micro1
## 2 Barbara micro2
## 3     Bob micro3
```

Let's include a `GRangesList`:  


```r
library(GenomicRanges)
gr1 <-
  GRanges(seqnames = "chr3", ranges = IRanges(58000000, 59502360), #completely encompasses ENST00000355076
          strand = "+", score = 5L, GC = 0.45)
gr2 <-
  GRanges(seqnames = c("chr3", "chr3"),
          ranges = IRanges(c(58493000, 3), width=9000), #first is within ENST0000035076
          strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
gr3 <-
  GRanges(seqnames = c("chr1", "chr2"),
          ranges = IRanges(c(1, 4), c(3, 9)),
          strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
names(grl) <- c("snparray1", "snparray2", "snparray3")
grl
```

```
## GRangesList object of length 3:
## $snparray1 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames               ranges strand |     score        GC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]     chr3 [58000000, 59502360]      + |         5      0.45
## 
## $snparray2 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames               ranges strand | score  GC
##   [1]     chr3 [58493000, 58501999]      + |     3 0.3
##   [2]     chr3 [       3,     9002]      - |     4 0.5
## 
## $snparray3 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames ranges strand | score  GC
##   [1]     chr1 [1, 3]      - |     6 0.4
##   [2]     chr2 [4, 9]      - |     2 0.1
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

The following `data.frame` matches masterPheno sample to the `GRangesList`: 


```r
(rangemap <- data.frame(master = c("Jack", "Jill", "Jill"), 
						assay = c("snparray1", "snparray2", "snparray3")))
```

```
##   master     assay
## 1   Jack snparray1
## 2   Jill snparray2
## 3   Jill snparray3
```

Adding the new `RangedSummarizedExperiment` class:

Create a `GenomicRangesList` object for the `RangedSummarizedExperiment`:


```r
gg1 <-
  GRanges(seqnames = "chr2", ranges = IRanges(48000000, 49502360),
          strand = "-", score = 8L, GC = 0.95)
gg2 <-
  GRanges(seqnames = c("chr1", "chr2"),
          ranges = IRanges(c(48493000, 30000), width=9000),
          strand = c("-", "+"), score = 6:7, GC = c(0.8, 1.0))
gg3 <-
  GRanges(seqnames = c("chr1", "chr2"),
          ranges = IRanges(c(1, 4), c(3, 9)),
          strand = c("+", "+"), score = c(8L, 4L), GC = c(1.4, 1.1))
gg4 <- 
  GRanges(seqnames = c("chr1", "chr2"), 
          ranges = IRanges(c(1, 6), c(3, 9)),
          strand = c("+", "-"), score = 2.8, GC = c(0.3, 0.5))
grl2 <- GRangesList("grange1" = gg1, "grange2" = gg2, "grange3" = gg3, "grange4" = gg4)
names(grl2) <- c("mysnparray1", "mysnparray2", "mysnparray3", "mysnparray4")
```


```r
library(SummarizedExperiment)
nrows <- 4
ncols <- 4
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 2),
                     row.names=c("Jack", "Jill", "Bob", "Barbara"))
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=grl2, colData=colData)
```


```r
(rangemap2 <- data.frame(master = c("Jack", "Jill", "Bob", "Barbara"), 
                        assay = c("mysnparray1", "mysnparray2", "mysnparray3", "mysnparray4")))
```

```
##    master       assay
## 1    Jack mysnparray1
## 2    Jill mysnparray2
## 3     Bob mysnparray3
## 4 Barbara mysnparray4
```

# Setup for creating the `MultiAssayExperiment` object
Create an ID map for all available experiments.  Names required, and must be identical to names of `elist`.


```r
idmap <- list(exprmap, methylmap, micromap, rangemap, rangemap2)
names(idmap) <- c("Affy", "Methyl 450k", "Mirna", "CNV gistic", "CNV gistic2")
idmap
```

```
## $Affy
##    master  assay
## 1    Jack array1
## 2    Jill array2
## 3 Barbara array3
## 4     Bob array4
## 
## $`Methyl 450k`
##    master   assay
## 1    Jack methyl1
## 2    Jack methyl2
## 3    Jill methyl3
## 4 Barbara methyl4
## 5     Bob methyl5
## 
## $Mirna
##    master  assay
## 1    Jack micro1
## 2 Barbara micro2
## 3     Bob micro3
## 
## $`CNV gistic`
##   master     assay
## 1   Jack snparray1
## 2   Jill snparray2
## 3   Jill snparray3
## 
## $`CNV gistic2`
##    master       assay
## 1    Jack mysnparray1
## 2    Jill mysnparray2
## 3     Bob mysnparray3
## 4 Barbara mysnparray4
```

Allowing for the possibility of ID maps entered as dataframes, convert to conventional list: 


```r
library(reshape2) # reshape2 used for example only
dfmap <- melt(idmap, id.var = c("master", "assay"))
library(biocMultiAssay)
toListMap(dfmap)
```

```
## $assay
##     master       assay
## 1     Jack      array1
## 2     Jill      array2
## 3  Barbara      array3
## 4      Bob      array4
## 5     Jack     methyl1
## 6     Jack     methyl2
## 7     Jill     methyl3
## 8  Barbara     methyl4
## 9      Bob     methyl5
## 10    Jack      micro1
## 11 Barbara      micro2
## 12     Bob      micro3
## 13    Jack   snparray1
## 14    Jill   snparray2
## 15    Jill   snparray3
## 16    Jack mysnparray1
## 17    Jill mysnparray2
## 18     Bob mysnparray3
## 19 Barbara mysnparray4
## 
## $L1
##     master          L1
## 1     Jack        Affy
## 2     Jill        Affy
## 3  Barbara        Affy
## 4      Bob        Affy
## 5     Jack Methyl 450k
## 6     Jack Methyl 450k
## 7     Jill Methyl 450k
## 8  Barbara Methyl 450k
## 9      Bob Methyl 450k
## 10    Jack       Mirna
## 11 Barbara       Mirna
## 12     Bob       Mirna
## 13    Jack  CNV gistic
## 14    Jill  CNV gistic
## 15    Jill  CNV gistic
## 16    Jack CNV gistic2
## 17    Jill CNV gistic2
## 18     Bob CNV gistic2
## 19 Barbara CNV gistic2
```

Create an named list of experiments `objlist` for the MultiAssay function


```r
objlist <- list("Affy" = exprdat, "Methyl 450k" = methyldat, "Mirna" = microdat, "CNV gistic" = grl, "CNV gistic2" = rse)
```

# Create a `multiAssayExperiment` class object


```r
myMultiAssay <- MultiAssay(objlist, masPheno, idmap)
myMultiAssay
```

```
## A MultiAssayExperiment object with 5 
##  listed experiments and their user-defined names: 
##    Affy 
##    Methyl 450k 
##    Mirna 
##    CNV gistic 
##    CNV gistic2 
## A "masterPheno" slot:
##         sex age
## Jack      M  38
## Jill      F  39
## Bob       M  40
## Barbara   F  41
## A "sampleMap" slot:
## $Affy
##    master  assay
## 1    Jack array1
## 2    Jill array2
## 3 Barbara array3
## 4     Bob array4
## 
## $`Methyl 450k`
##    master   assay
## 1    Jack methyl1
## 2    Jack methyl2
## 3    Jill methyl3
## 4 Barbara methyl4
## 5     Bob methyl5
## 
## $Mirna
##    master  assay
## 1    Jack micro1
## 2 Barbara micro2
## 3     Bob micro3
## 
## $`CNV gistic`
##   master     assay
## 1   Jack snparray1
## 2   Jill snparray2
## 3   Jill snparray3
## 
## $`CNV gistic2`
##    master       assay
## 1    Jack mysnparray1
## 2    Jill mysnparray2
## 3     Bob mysnparray3
## 4 Barbara mysnparray4
## 
## An "elist" slot:
## $Affy
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 2 features, 4 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2 array3 array4
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $`Methyl 450k`
##                 methyl1 methyl2 methyl3 methyl4 methyl5
## ENST00000355076       1       3       5       7       9
## ENST00000383706       2       4       6       8      10
## 
## $Mirna
##              micro1 micro2 micro3
## hsa-miR-21      201    205    209
## hsa-miR-191     202    206    210
## hsa-miR-148a    203    207    211
## hsa-miR148b     204    208    212
## 
## $`CNV gistic`
## GRangesList object of length 3:
## $snparray1 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames               ranges strand |     score        GC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]     chr3 [58000000, 59502360]      + |         5      0.45
## 
## $snparray2 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames               ranges strand | score  GC
##   [1]     chr3 [58493000, 58501999]      + |     3 0.3
##   [2]     chr3 [       3,     9002]      - |     4 0.5
## 
## $snparray3 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames ranges strand | score  GC
##   [1]     chr1 [1, 3]      - |     6 0.4
##   [2]     chr2 [4, 9]      - |     2 0.1
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic2`
## class: RangedSummarizedExperiment 
## dim: 4 4 
## metadata(0):
## assays(1): counts
## rownames(4): mysnparray1 mysnparray2 mysnparray3 mysnparray4
## rowRanges metadata column names(0):
## colnames(4): Jack Jill Bob Barbara
## colData names(1): Treatment
## 
## A "metadata" slot:
## NULL
```

## Subsetting by Sample


```r
logicID <- identifyBySample(myMultiAssay, 1:2)
logicID
```

```
## $Affy
## [1]  TRUE  TRUE FALSE FALSE
## 
## $`Methyl 450k`
## [1]  TRUE  TRUE  TRUE FALSE FALSE
## 
## $Mirna
## [1]  TRUE FALSE FALSE
## 
## $`CNV gistic`
## [1] TRUE TRUE TRUE
## 
## $`CNV gistic2`
## [1]  TRUE  TRUE FALSE FALSE
```

```r
subMultiAssay <- subsetBySample(myMultiAssay, logicID)
```

```
## Warning in methods:::.selectDotsMethod(classes, .MTable, .AllMTable):
## multiple direct matches: "SimpleList", "list"; using the first of these
```

```r
as.list(elist(subMultiAssay))
```

```
## $Affy
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 2 features, 2 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $`Methyl 450k`
##                 methyl1 methyl2 methyl3
## ENST00000355076       1       3       5
## ENST00000383706       2       4       6
## 
## $Mirna
##              micro1
## hsa-miR-21      201
## hsa-miR-191     202
## hsa-miR-148a    203
## hsa-miR148b     204
## 
## $`CNV gistic`
## GRangesList object of length 3:
## $snparray1 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames               ranges strand |     score        GC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]     chr3 [58000000, 59502360]      + |         5      0.45
## 
## $snparray2 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames               ranges strand | score  GC
##   [1]     chr3 [58493000, 58501999]      + |     3 0.3
##   [2]     chr3 [       3,     9002]      - |     4 0.5
## 
## $snparray3 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames ranges strand | score  GC
##   [1]     chr1 [1, 3]      - |     6 0.4
##   [2]     chr2 [4, 9]      - |     2 0.1
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic2`
## class: RangedSummarizedExperiment 
## dim: 2 4 
## metadata(0):
## assays(1): counts
## rownames(2): mysnparray1 mysnparray2
## rowRanges metadata column names(0):
## colnames(4): Jack Jill Bob Barbara
## colData names(1): Treatment
```

```r
subsetByAssay(subMultiAssay, c(TRUE, FALSE, FALSE, FALSE), drop = TRUE)[[1]] %>% exprs 
```

```
##                 array1 array2
## ENST00000294241    101    103
## ENST00000355076    102    104
```

```r
subMultiAssay
```

```
## A MultiAssayExperiment object with 5 
##  listed experiments and their user-defined names: 
##    Affy 
##    Methyl 450k 
##    Mirna 
##    CNV gistic 
##    CNV gistic2 
## A "masterPheno" slot:
##      sex age
## Jack   M  38
## Jill   F  39
## A "sampleMap" slot:
## $Affy
## [1] "array1" "array2"
## 
## $`Methyl 450k`
## [1] "methyl1" "methyl2" "methyl3"
## 
## $Mirna
## [1] "micro1"
## 
## $`CNV gistic`
## [1] "snparray1" "snparray2" "snparray3"
## 
## $`CNV gistic2`
## [1] "mysnparray1" "mysnparray2"
## 
## An "elist" slot:
## $Affy
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 2 features, 2 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $`Methyl 450k`
##                 methyl1 methyl2 methyl3
## ENST00000355076       1       3       5
## ENST00000383706       2       4       6
## 
## $Mirna
##              micro1
## hsa-miR-21      201
## hsa-miR-191     202
## hsa-miR-148a    203
## hsa-miR148b     204
## 
## $`CNV gistic`
## GRangesList object of length 3:
## $snparray1 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames               ranges strand |     score        GC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]     chr3 [58000000, 59502360]      + |         5      0.45
## 
## $snparray2 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames               ranges strand | score  GC
##   [1]     chr3 [58493000, 58501999]      + |     3 0.3
##   [2]     chr3 [       3,     9002]      - |     4 0.5
## 
## $snparray3 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames ranges strand | score  GC
##   [1]     chr1 [1, 3]      - |     6 0.4
##   [2]     chr2 [4, 9]      - |     2 0.1
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic2`
## class: RangedSummarizedExperiment 
## dim: 2 4 
## metadata(0):
## assays(1): counts
## rownames(2): mysnparray1 mysnparray2
## rowRanges metadata column names(0):
## colnames(4): Jack Jill Bob Barbara
## colData names(1): Treatment
## 
## A "metadata" slot:
## NULL
```

Endogenous operation, returns a MultiAssay object containing elist of length 1, map of length 1, and masterPheno for only Jack, Barbara, and Bob.  The "Mirna" argument is used to index the `elist` object using `[`, so could also be `integer` or `logical`:


```r
subsetByAssay(myMultiAssay, "Mirna", drop=FALSE)
```

```
## A MultiAssayExperiment object with 1 
##  listed experiment and their user-defined name: 
##    Mirna 
## A "masterPheno" slot:
##         sex age
## Jack      M  38
## Jill      F  39
## Bob       M  40
## Barbara   F  41
## A "sampleMap" slot:
## $Mirna
##    master  assay
## 1    Jack micro1
## 2 Barbara micro2
## 3     Bob micro3
## 
## An "elist" slot:
## $Mirna
##              micro1 micro2 micro3
## hsa-miR-21      201    205    209
## hsa-miR-191     202    206    210
## hsa-miR-148a    203    207    211
## hsa-miR148b     204    208    212
## 
## A "metadata" slot:
## NULL
```

Not endogenous, returns just the matrix contained in `elist$Mirna`


```r
subsetByAssay(myMultiAssay, "Mirna", drop=TRUE)
```

```
## $Mirna
##              micro1 micro2 micro3
## hsa-miR-21      201    205    209
## hsa-miR-191     202    206    210
## hsa-miR-148a    203    207    211
## hsa-miR148b     204    208    212
```

## Subsetting by Feature 

This operation returns a `MultiAssayExperiment` class, with any `elist` element not containing the feature having zero rows.

Until we make subsetting by a gene ID work on ranges, the following will return c(TRUE, TRUE, FALSE, FALSE). When we implement subsetting by ranges, it should return c(TRUE, TRUE, FALSE, TRUE). 

Returns c(TRUE, TRUE, FALSE, FALSE):

```r
identifyByFeature(myMultiAssay, "ENST00000355076")  
```

```
## [1]  TRUE  TRUE FALSE FALSE FALSE
```

Returns MultiAssayExperiment where `Affy` and `Methyl 450k` contain only ENST0000035076 row, and "Mirna" and "CNV gistic" have zero rows:


```r
featSubsetted0 <- subsetByFeature(myMultiAssay, "ENST00000355076")
class(featSubsetted0)
```

```
## [1] "MultiAssayExperiment"
## attr(,"package")
## [1] "biocMultiAssay"
```

```r
class(elist(featSubsetted0))
```

```
## [1] "SimpleList"
## attr(,"package")
## [1] "S4Vectors"
```

```r
as.list(elist(featSubsetted0))
```

```
## $Affy
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 1 features, 4 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2 array3 array4
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $`Methyl 450k`
##                 methyl1 methyl2 methyl3 methyl4 methyl5
## ENST00000355076       1       3       5       7       9
## 
## $Mirna
## [1] "micro1" "micro2" "micro3"
## 
## $`CNV gistic`
## $`CNV gistic`$snparray1
## GRangesList object of length 0:
## <0 elements>
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic`$snparray2
## GRangesList object of length 0:
## <0 elements>
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic`$snparray3
## GRangesList object of length 0:
## <0 elements>
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## 
## $`CNV gistic2`
## class: RangedSummarizedExperiment 
## dim: 0 4 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowRanges metadata column names(0):
## colnames(4): Jack Jill Bob Barbara
## colData names(1): Treatment
```

In the following, `Affy` ExpressionSet keeps both rows but with their order reversed, and `Methyl 450k` keeps only its second row.


```r
featSubsetted <- subsetByFeature(myMultiAssay, c("ENST00000355076", "ENST00000294241"))
as.list(elist(featSubsetted))
```

```
## $Affy
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 2 features, 4 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2 array3 array4
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $`Methyl 450k`
##                 methyl1 methyl2 methyl3 methyl4 methyl5
## ENST00000355076       1       3       5       7       9
## 
## $Mirna
## [1] "micro1" "micro2" "micro3"
## 
## $`CNV gistic`
## $`CNV gistic`$snparray1
## GRangesList object of length 0:
## <0 elements>
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic`$snparray2
## GRangesList object of length 0:
## <0 elements>
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic`$snparray3
## GRangesList object of length 0:
## <0 elements>
## 
## -------
## seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## 
## $`CNV gistic2`
## class: RangedSummarizedExperiment 
## dim: 0 4 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowRanges metadata column names(0):
## colnames(4): Jack Jill Bob Barbara
## colData names(1): Treatment
```

## Identify assays that contain any / all of a vector of features

Note that the output of this function could be used as the input for `subsetByAssay`.

```r
identifyByFeature(myMultiAssay, c("ENST00000355076", "ENST00000294241"), requireall=FALSE) 
```

```
## [1]  TRUE  TRUE FALSE FALSE FALSE
```

```r
identifyByFeature(myMultiAssay, c("ENST00000355076", "ENST00000294241"), requireall=TRUE) 
```

```
## [1]  TRUE FALSE FALSE FALSE FALSE
```

### Feature extraction by Ranges

Use `GenomicRanges::subsetByOverlaps`, passing its arguments to allow flexible types of subsetting. 


```r
rangeSubset <- GRanges(seqnames = c("chr1"), strand = c("-", "+", "-"), ranges = IRanges(start = c(1, 4, 6), width = 3))
subsetted <- subsetByFeature(myMultiAssay, rangeSubset)
as.list(elist(subsetted))
```

```
## $Affy
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 0 features, 4 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: array1 array2 array3 array4
##   varLabels: slope53
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $`Methyl 450k`
##      methyl1 methyl2 methyl3 methyl4 methyl5
## 
## $Mirna
##      micro1 micro2 micro3
## 
## $`CNV gistic`
## $`CNV gistic`$snparray1
## GRanges object with 0 ranges and 2 metadata columns:
##    seqnames    ranges strand |     score        GC
##       <Rle> <IRanges>  <Rle> | <integer> <numeric>
##   -------
##   seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic`$snparray2
## GRanges object with 0 ranges and 2 metadata columns:
##    seqnames    ranges strand |     score        GC
##       <Rle> <IRanges>  <Rle> | <integer> <numeric>
##   -------
##   seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## $`CNV gistic`$snparray3
## GRanges object with 1 range and 2 metadata columns:
##       seqnames    ranges strand |     score        GC
##          <Rle> <IRanges>  <Rle> | <integer> <numeric>
##   [1]     chr1    [1, 3]      - |         6       0.4
##   -------
##   seqinfo: 3 sequences from an unspecified genome; no seqlengths
## 
## 
## $`CNV gistic2`
## class: RangedSummarizedExperiment 
## dim: 0 4 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowRanges metadata column names(0):
## colnames(4): Jack Jill Bob Barbara
## colData names(1): Treatment
```

TODO: Pass more arguments to `subsetByOverlaps`

