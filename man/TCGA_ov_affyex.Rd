\name{TCGA_ov_affyex}
\alias{TCGA_ov_affyex}
\docType{data}
\title{
TCGA level 3 expression data
}
\description{
TCGA level 3 expression data
}
\usage{data(TCGA_ov_affyex)}
\format{
  The format is: \cr
Formal class 'SummarizedExperiment' [package "GenomicRanges"] with 4 slots \cr
  ..@ exptData:Formal class 'SimpleList' [package "IRanges"] with 4 slots \cr
  .. .. ..@ listData       :List of 1 \cr
  .. .. .. ..$ initExptData:Formal class 'MIAME' [package "Biobase"] with 13 slots \cr
  .. .. .. .. .. ..@ name             : chr "Integrated genomic analyses of ovarian carcinoma. Nature 2011, 474:609-615." \cr
  .. .. .. .. .. ..@ lab              : chr "Cancer Genome Atlas Research Network 2011" \cr
  .. .. .. .. .. ..@ contact          : chr "" \cr
  .. .. .. .. .. ..@ title            : chr "Integrated genomic analyses of ovarian carcinoma." \cr
  .. .. .. .. .. ..@ abstract         : chr "A catalogue of molecular aberrations that cause ovarian cancer is critical for developing and deploying therapies that will imp"| __truncated__ \cr
  .. .. .. .. .. ..@ url              : chr "" \cr
  .. .. .. .. .. ..@ pubMedIds        : chr "21720365" \cr
  .. .. .. .. .. ..@ samples          : list() \cr
  .. .. .. .. .. ..@ hybridizations   : list() \cr
  .. .. .. .. .. ..@ normControls     : list() \cr
  .. .. .. .. .. ..@ preprocessing    :List of 1 \cr
  .. .. .. .. .. .. ..$ : chr "rma" \cr
  .. .. .. .. .. ..@ other            :List of 7 \cr
  .. .. .. .. .. .. ..$ platform_title       : chr "[HT_HG-U133A] Affymetrix HT Human Genome U133A Array" \cr
  .. .. .. .. .. .. ..$ platform_shorttitle  : chr "Affymetrix HT_HG-U133A" \cr
  .. .. .. .. .. .. ..$ platform_summary     : chr "hthgu133a" \cr
  .. .. .. .. .. .. ..$ platform_manufacturer: chr "Affymetrix" \cr
  .. .. .. .. .. .. ..$ platform_distribution: chr "commercial" \cr
  .. .. .. .. .. .. ..$ platform_accession   : chr "GPL3921" \cr
  .. .. .. .. .. .. ..$ warnings             :List of 1 \cr
  .. .. .. .. .. .. .. ..$ warnings: chr "The following samples are likely from specimens also used in GSE26712: TCGA.13.0725, TCGA.13.0885, TCGA.13.0887, TCGA.13.0890, "| __truncated__ \cr
  .. .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots \cr
  .. .. .. .. .. .. .. ..@ .Data:List of 2 \cr
  .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 0 0 \cr
  .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0 \cr
  .. .. ..@ elementType    : chr "ANY" \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ rowData :Formal class 'GRanges' [package "GenomicRanges"] with 6 slots \cr
  .. .. ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots \cr
  .. .. .. .. ..@ values         : Factor w/ 24 levels "chr1","chr10",..: 2 4 15 16 4 16 7 12 9 13 ... \cr
  .. .. .. .. ..@ lengths        : int [1:10785] 1 1 1 1 2 1 1 2 1 1 ... \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. .. .. ..@ start          : int [1:12791] 52559169 9220304 43088127 137842560 53701240 125549925 151531861 67493367 69685127 219128852 ... \cr
  .. .. .. .. ..@ width          : int [1:12791] 86267 48255 28750 8670 14173 77947 14416 53708 185851 6042 ... \cr
  .. .. .. .. ..@ NAMES          : chr [1:12791] "220951_s_at" "217757_at" "219488_at" "221131_at" ... \cr
  .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots \cr
  .. .. .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 2 1 2 1 2 1 2 1 2 1 ... \cr
  .. .. .. .. ..@ lengths        : int [1:6235] 5 2 3 2 2 1 2 1 8 1 ... \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. .. .. ..@ rownames       : NULL \cr
  .. .. .. .. ..@ nrows          : int 12791 \cr
  .. .. .. .. ..@ listData       : Named list() \cr
  .. .. .. .. ..@ elementType    : chr "ANY" \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots \cr
  .. .. .. .. ..@ seqnames   : chr [1:24] "chr1" "chr10" "chr11" "chr12" ... \cr
  .. .. .. .. ..@ seqlengths : int [1:24] NA NA NA NA NA NA NA NA NA NA ... \cr
  .. .. .. .. ..@ is_circular: logi [1:24] NA NA NA NA NA NA ... \cr
  .. .. .. .. ..@ genome     : chr [1:24] NA NA NA NA ... \cr
  .. .. ..@ metadata       :List of 2 \cr
  .. .. .. ..$ dropped   : chr [1:198] "216418_at" "217211_at" "216814_at" "217237_at" ... \cr
  .. .. .. ..$ nmultiaddr: int 7312 \cr
  ..@ colData :Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. ..@ rownames       : NULL \cr
  .. .. ..@ nrows          : int 578 \cr
  .. .. ..@ listData       :List of 29 \cr
  .. .. .. ..$ alt_sample_name                    : chr [1:578] "TCGA-20-0987-01A-02R-0434-01" "TCGA-23-1031-01A-01R-0434-01" "TCGA-24-0979-01A-01R-0434-01" "TCGA-23-1117-01A-02R-0434-01" ... \cr
  .. .. .. ..$ unique_patient_ID                  : chr [1:578] "TCGA-20-0987" "TCGA-23-1031" "TCGA-24-0979" "TCGA-23-1117" ... \cr
  .. .. .. ..$ sample_type                        : chr [1:578] "tumor" "tumor" "tumor" "tumor" ... \cr
  .. .. .. ..$ histological_type                  : chr [1:578] "ser" "ser" "ser" "ser" ... \cr
  .. .. .. ..$ primarysite                        : chr [1:578] "ov" "ov" "ov" "ov" ... \cr
  .. .. .. ..$ arrayedsite                        : logi [1:578] NA NA NA NA NA NA ... \cr
  .. .. .. ..$ summarygrade                       : chr [1:578] "high" "high" "high" "high" ... \cr
  .. .. .. ..$ summarystage                       : chr [1:578] "late" "late" "late" "late" ... \cr
  .. .. .. ..$ tumorstage                         : int [1:578] 3 4 4 3 4 3 3 4 3 3 ... \cr
  .. .. .. ..$ substage                           : chr [1:578] "c" NA NA "c" ... \cr
  .. .. .. ..$ grade                              : int [1:578] 3 3 3 3 3 2 3 3 3 3 ... \cr
  .. .. .. ..$ age_at_initial_pathologic_diagnosis: int [1:578] 61 60 53 42 45 78 74 73 45 45 ... \cr
  .. .. .. ..$ pltx                               : chr [1:578] "y" "y" "y" "y" ... \cr
  .. .. .. ..$ tax                                : chr [1:578] "y" "y" "y" "n" ... \cr
  .. .. .. ..$ neo                                : chr [1:578] "n" "n" "n" "n" ... \cr
  .. .. .. ..$ days_to_tumor_recurrence           : int [1:578] 442 574 428 418 1446 61 870 68 2616 797 ... \cr
  .. .. .. ..$ recurrence_status                  : chr [1:578] "recurrence" "norecurrence" "recurrence" "recurrence" ... \cr
  .. .. .. ..$ days_to_death                      : int [1:578] 701 574 1264 1013 1446 61 789 84 2616 816 ... \cr
  .. .. .. ..$ vital_status                       : chr [1:578] "deceased" "deceased" "deceased" "deceased" ... \cr
  .. .. .. ..$ os_binary                          : logi [1:578] NA NA NA NA NA NA ... \cr
  .. .. .. ..$ relapse_binary                     : logi [1:578] NA NA NA NA NA NA ... \cr
  .. .. .. ..$ site_of_tumor_first_recurrence     : chr [1:578] "locoregional" NA "locoregional" "locoregional" ... \cr
  .. .. .. ..$ primary_therapy_outcome_success    : chr [1:578] NA "partialresponse" "completeresponse" "completeresponse" ... \cr
  .. .. .. ..$ debulking                          : chr [1:578] "optimal" "optimal" "suboptimal" "optimal" ... \cr
  .. .. .. ..$ percent_normal_cells               : int [1:578] 16 8 0 0 0 0 12 15 NA 0 ... \cr
  .. .. .. ..$ percent_stromal_cells              : int [1:578] 1 5 25 NA 0 36 2 5 NA NA ... \cr
  .. .. .. ..$ percent_tumor_cells                : int [1:578] 72 80 72 NA NA 64 84 72 NA NA ... \cr
  .. .. .. ..$ batch                              : int [1:578] 12 12 12 12 12 12 12 12 12 12 ... \cr
  .. .. .. ..$ uncurated_author_metadata          : chr [1:578] "age_at_initial_pathologic_diagnosis: 61///anatomic_organ_subdivision: NA///bcr_patient_uuid: 34f285ab-6c31-4865-a0a7-a57af3567d"| __truncated__ "age_at_initial_pathologic_diagnosis: 60///anatomic_organ_subdivision: NA///bcr_patient_uuid: c7c3fe4a-327a-4b7e-a259-989b729712"| __truncated__ "age_at_initial_pathologic_diagnosis: 53///anatomic_organ_subdivision: NA///bcr_patient_uuid: 9a30dee8-962a-49db-a85a-c39515c91b"| __truncated__ "age_at_initial_pathologic_diagnosis: 42///anatomic_organ_subdivision: NA///bcr_patient_uuid: a88b7e66-5f12-4023-a7e2-fcfbd1f259"| __truncated__ ... \cr
  .. .. ..@ elementType    : chr "ANY" \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ assays  :Reference class 'ShallowSimpleListAssays' [package "GenomicRanges"] with 1 fields \cr
  .. ..$ data:Formal class 'SimpleList' [package "IRanges"] with 4 slots \cr
  .. .. .. ..@ listData       :List of 1 \cr
  .. .. .. .. ..$ exprs: num [1:12791, 1:578] 2.92 10.35 4.21 3.32 4.61 ... \cr
  .. .. .. ..@ elementType    : chr "ANY" \cr
  .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. ..@ metadata       : list() \cr
  .. ..and 12 methods, \cr
}
\details{
Expression data from TCGA level 3
}
\source{
NCI TCGA
}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(TCGA_ov_affyex)
}
\keyword{datasets}
