# Copy number by GISTIC2

* Relevant files:
    + **all_lesions.conf_99.txt** - data file with ranges
    + **SNP6 Copy number analysis (GISTIC2).pdf** - full GISTIC2 documentation output

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

* Relevant files:
    + **OV.medianexp.txt** - One row per gene symbol
    + Available from: https://dl.dropboxusercontent.com/u/15152544/TCGA/OV.medianexp.txt
    + R: download.file("https://dl.dropboxusercontent.com/u/15152544/TCGA/OV.medianexp.txt", destfile="OV.medianexp.txt", method="wget")

## Methods
    firehose_get -tasks mrna stddata latest ov
    tar xvfz stddata__2014_09_02/OV/20140902/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_3.2014090200.0.0.tar.gz
    cp stddata__2014_09_02/OV/20140902/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_3.2014090200.0.0/OV.medianexp.txt .

# Custom Affy for OVC

* Relevant files: 
    + **TCGA_eset.rda** - fRMA-processed ExpressionSet
      from HT-HG_U133A CEL files.  Just going to put this here because
      it's necessary for the known sample mix-ups, but I don't think in
      general we'll be processing CEL files.
    + Available from: https://dl.dropboxusercontent.com/u/15152544/TCGA/TCGA_eset.rda
    + R: download.file("https://dl.dropboxusercontent.com/u/15152544/TCGA/TCGA_eset.rda", destfile="TCGA_eset.rda", method="wget")
