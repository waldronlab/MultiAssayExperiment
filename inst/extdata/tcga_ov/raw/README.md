# Copy number by GISTIC2

## Relevant files:

*  **all_lesions.conf_99.txt** - data file with ranges
*  **SNP6 Copy number analysis (GISTIC2).pdf** - full GISTIC2 documentation output

## Methods Install [firehose_get command-line tool from the
broad](https://confluence.broadinstitute.org/display/GDAC/Download),
then:

    firehose_get -tasks gistic analyses latest ov
    cp analyses__2014_07_15/OV/20140715/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2014071500.0.0/all_lesions.conf_99.txt .

## Notes from GISTIC output

Columns 1-9 of all_lesions.conf_99.txt present the data about the
significant regions as follows:

1. Unique Name: A name assigned to identify the region.
2. Descriptor: The genomic descriptor of that region.
3. Wide Peak Limits: The 'wide peak' boundaries most likely to contain the targeted genes. These are listed in genomic coordinates and marker (or probe) indices.
4. Peak Limits: The boundaries of the region of maximal amplification or deletion.
5. Region Limits: The boundaries of the entire significant region of amplification or deletion.
6. Q values: The Q value of the peak region.
7. Residual Q values: The Q value of the peak region after removing ('peeling off') amplifications or deletions that overlap other, more significant peak regions in the same chromosome.
8. Broad or Focal: Identifies whether the region reaches significance due primarily to broad events (called 'broad'), focal events (called 'focal'), or independently significant broad and focal events (called 'both').
9. Amplitude Threshold: Key giving the meaning of values in the subsequent columns associated with each sample.

# Default mRNA

## Relevant files

*  **OV.medianexp.txt** - One row per gene symbol
*  [Download link](https://dl.dropboxusercontent.com/u/15152544/TCGA/OV.medianexp.txt)

## Methods

    firehose_get -tasks mrna stddata latest ov
    tar xvfz stddata__2014_09_02/OV/20140902/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_3.2014090200.0.0.tar.gz
    cp stddata__2014_09_02/OV/20140902/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_3.2014090200.0.0/OV.medianexp.txt .

## Demo data

    Hybridization REF       TCGA-01-0628-11A-01R-0361-03    TCGA-01-0630-11A-01R-0361-03
    Composite Element REF   Signal  Signal
    C9orf152        5.68406070638714        5.74059738582655
    ELMO2   6.94072114738553        6.9554576454791
    RPS11   12.112295206898 11.6037710967151

# Methylation 27K data

## Relevant files:

*  OV.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt
*  [Download link](https://dl.dropboxusercontent.com/u/15152544/TCGA/OV.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt)

## Methods
    firehose_get -tasks methylation stddata latest ov
    tar xfz stddata__2014_09_02/OV/20140902/gdac.broadinstitute.org_OV.Methylation_Preprocess.Level_3.2014090200.0.0.tar.gz
    cp stddata__2014_09_02/OV/20140902/gdac.broadinstitute.org_OV.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014090200.0.0/OV.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt .

## Demo data

    Hybridization REF       TCGA-01-0628-11A-01D-0383-05    TCGA-01-0628-11A-01D-0383-05    TCGA-01-0628-11A-01D-0383-05    TCGA-01-0628-11A-01D-0383-05
    Composite Element REF   Beta_value      Gene_Symbol     Chromosome      Genomic_Coordinate
    cg00000292      0.799408581806102       ATP2A1  16      28890100
    cg00002426      0.339004437825915       SLMAP   3       57743543
    cg00003994      0.0281193010321092      MEOX2   7       15725862

# Custom Affy for OVC

## Relevant files: 

*  **TCGA_eset.rda** - RMA-processed ExpressionSet with phenoData
      from HT-HG_U133A CEL files.  Just going to put this here because
      it's necessary for the known sample mix-ups, but I don't think in
      general we'll be processing CEL files.
*  [Download link](https://dl.dropboxusercontent.com/u/15152544/TCGA/TCGA_eset.rda)

# Multiple assays by RTCGAToolbox

I ran RTCGAToolbox with all data options TRUE (except for
miRNASeq_Gene, which generated an error and is reported to the
developer).

## Relevant files:
*  **ovres.rda** - object of class "FirehoseData"
*  [Download link](https://dl.dropboxusercontent.com/u/15152544/TCGA/ovres.rda)

## Methods

    if( !require(RTCGAToolbox) ){
        library(devtools)
        install_github("mksamur/RTCGAToolbox")
        library(RTCGAToolbox)
    }
    all.dates <- getFirehoseRunningDates()
    an.dates <- getFirehoseAnalyzeDates()
    ovres <- getFirehoseData("OV", runDate=all.dates[1], gistic2_Date=an.dates[1],
                         RNAseq_Gene=TRUE, Clinic=TRUE, miRNASeq_Gene=FALSE,
                         RNAseq2_Gene_Norm=TRUE, CNA_SNP=TRUE, CNV_SNP=TRUE, 
                         CNA_Seq=TRUE, CNA_CGH=TRUE,  Methylation=TRUE, 
                         Mutation=TRUE, mRNA_Array=TRUE, miRNA_Array=TRUE, RPPA=TRUE)
    save(ovres, file="ovres.rda", compress="bzip2")

    > slotNames(ovres)
        [1] "Dataset"         "Clinical"        "RNASeqGene"      "RNASeq2GeneNorm"
        [5] "miRNASeqGene"    "CNASNP"          "CNVSNP"          "CNAseq"         
        [9] "CNACGH"          "Methylation"     "mRNAArray"       "miRNAArray"     
       [13] "RPPAArray"       "Mutations"       "GISTIC"         
