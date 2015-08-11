
.db_index <- function(db, tablename) {
#
# need some bulletproofing on column name availability
#
  message("indexing")
  idx <- c(
     "CREATE INDEX sampidx on %%TAB%%(sampleID);",
     "CREATE INDEX featidx on %%TAB%%(feature);",
     "CREATE INDEX atypeidx on %%TAB%%(assaytype);")
  idx <- gsub("%%TAB%%", tablename, idx)
  for (i in idx) {
     message(i)
     dbSendQuery(db, i)
  }
  db
}

.createOvarianHub = function() {
 library(biocMultiAssay)
 ov = dir(system.file("extdata/tcga_ov",
   package="biocMultiAssay"), full=TRUE, pattern="\\.rda$")
 drop = grep("pheno", ov)
 if (length(drop)>0) {
   pdpath=ov[drop]
   ov=ov[-drop]
  }
 #
 # informal labels for constituents, CamelCase
 #
 tags = c("OvarCnvGistic", "OvarMethy450k", "OvarAffyExpr",
    "OvarAgilent", "OvarMiRna", "OvarRNASeq")
 ovlist <- lapply(ov, function(x) get(load(x)))
 names(ovlist) <- tags
 createMA(masterpheno=pData(ovlist[[2]]), objlist=ovlist, drop=TRUE)
}

.eset2longdf <- function(es, tumortype="ovarian", assaytype, ...) {
  fn = featureNames(es)
  nsamp = ncol(es)
  ids = sampleNames(es)
  ids = rep(ids, each=length(fn))
  fn = rep(fn, nsamp)
  data.frame(sampleID=ids, feature=fn, value=as.numeric(exprs(es)),
         tumortype=tumortype, assaytype=assaytype, stringsAsFactors=FALSE)
}

.db_create_long <-
  function(dbname, hub, tablename, masterInd=2) # dbname is filename of SQLite database
{
  stopifnot(is(hub, "MultiAssayExperiment"))  # could weaken
  library(RSQLite)
  db <- dbConnect(dbDriver("SQLite"), dbname)
  nass <- length(hub@elist)
  alldf = lapply(1:nass, function(x) .eset2longdf(hub@elist[[x]],
      assaytype=names(hub@elist)[x]))
  alld = do.call(rbind, alldf)
  mphen = pData(hub@elist[[masterInd]])
  dbWriteTable(db, tablename, alld)
  dbWriteTable(db, "masterPheno", mphen)
  db
}

#ohub = .createOvarianHub()
#OvarCon <- .db_create_long("OvarLong.sqlite", ohub, "OvarLong")
#OvarCon <- .db_index(OvarCon, "OvarLong")

OvarHub = function() {
   library(RSQLite)
   library(dplyr)
   con <- dbConnect(dbDriver("SQLite"), Sys.getenv("OVARLONG_PATH"))
   src_sql("sqlite", con, path=Sys.getenv("OVARLONG_PATH"),
         info = dbGetInfo(con))
}

OvarHub = function() {
   library(RSQLite)
   library(dplyr)
   # need to prevent reconnecting
   if (!exists(".ovarHubCon")) {
   .ovarHubCon <<-
        dbConnect(dbDriver("SQLite"), Sys.getenv("OVARLONG_PATH"))
   message("global instance of .ovarHubCon created")
   }
   src_sql("sqlite", .ovarHubCon, path=Sys.getenv("OVARLONG_PATH"),
         info = dbGetInfo(.ovarHubCon))
}

OvMasterPheno = function() tbl(OvarHub(), "masterPheno")
OvarLong = function() tbl(OvarHub(), "OvarLong")

tabulateSamples = function(tgen=OvarLong,
  featureToCheck="A2M") {  # how to get an exemplar?
  stopifnot(is(featureToCheck,"character"))
  stopifnot(length(featureToCheck)==1)
  tgen() %>% dplyr::select(sampleID, assaytype) %>%
         filter(feature==featureToCheck) %>%
         group_by(assaytype) %>% summarize(n=n())
}

tabulateFeatures = function(tgen=OvarLong,
   sampleToCheck="TCGA.04.1331") {
  stopifnot(is(sampleToCheck,"character"))
  stopifnot(length(sampleToCheck)==1)
  tgen() %>% dplyr::select(sampleID, assaytype) %>%
          filter(sampleID==sampleToCheck) %>%
          group_by(assaytype) %>% summarize(n=n())
}

