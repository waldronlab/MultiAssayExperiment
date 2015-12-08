# 
# .db_index <- function(db, tablename) {
# #
# # need some bulletproofing on column name availability
# #
#   message("indexing")
#   idx <- c(
#      "CREATE INDEX sampidx on %%TAB%%(sampleID);",
#      "CREATE INDEX featidx on %%TAB%%(feature);",
#      "CREATE INDEX atypeidx on %%TAB%%(assaytype);")
#   idx <- gsub("%%TAB%%", tablename, idx)
#   for (i in idx) {
#      message(i)
#      DBI::dbSendQuery(db, i)
#   }
#   db
# }
# 
# .createOvarianHub = function() {
#  ov = dir(system.file("extdata/tcga_ov",
#    package="biocMultiAssay"), full.names=TRUE, pattern="\\.rda$")
#  drop = grep("pheno", ov)
#  if (length(drop)>0) {
#    pdpath=ov[drop]
#    ov=ov[-drop]
#   }
#  #
#  # informal labels for constituents, CamelCase
#  #
#  tags = c("OvarCnvGistic", "OvarMethy450k", "OvarAffyExpr",
#     "OvarAgilent", "OvarMiRna", "OvarRNASeq")
#  ovlist <- lapply(ov, function(x) get(load(x)))
#  names(ovlist) <- tags
#  MultiAssay(masterPheno=Biobase::pData(ovlist[[2]]), explist=ovlist) # drop=TRUE unused
# }
# 
# .eset2longdf <- function(es, tumortype="ovarian", assaytype, ...) {
#   fn = Biobase::featureNames(es)
#   nsamp = ncol(es)
#   ids = Biobase::sampleNames(es)
#   ids = rep(ids, each=length(fn))
#   fn = rep(fn, nsamp)
#   data.frame(sampleID=ids, feature=fn, value=as.numeric(Biobase::exprs(es)),
#          tumortype=tumortype, assaytype=assaytype, stringsAsFactors=FALSE)
# }
# 
# .db_create_long <-
#   function(dbname, hub, tablename, masterInd=2) # dbname is filename of SQLite database
# {
#   stopifnot(is(hub, "MultiAssayExperiment"))  # could weaken
#   db <- RSQLite::dbConnect(DBI::dbDriver("SQLite"), dbname)
#   nass <- length(hub@Elist)
#   alldf = lapply(1:nass, function(x) .eset2longdf(hub@Elist[[x]],
#       assaytype=names(hub@Elist)[x]))
#   alld = do.call(rbind, alldf)
#   mphen = Biobase::pData(hub@Elist[[masterInd]])
#   RSQLite::dbWriteTable(db, tablename, alld)
#   RSQLite::dbWriteTable(db, "masterPheno", mphen)
#   db
# }
# 
# #ohub = .createOvarianHub()
# #OvarCon <- .db_create_long("OvarLong.sqlite", ohub, "OvarLong")
# #OvarCon <- .db_index(OvarCon, "OvarLong")
# 
# OvarHub = function() {
#    con <- RSQLite::dbConnect(DBI::dbDriver("SQLite"), Sys.getenv("OVARLONG_PATH"))
#    dplyr::src_sql("sqlite", con, path=Sys.getenv("OVARLONG_PATH"),
#          info = RSQLite::dbGetInfo(con))
# }
# 
# OvarHub = function() {
#    # need to prevent reconnecting
#    if (!exists(".ovarHubCon")) {
#    .ovarHubCon <<-
#      RSQLite::dbConnect(DBI::dbDriver("SQLite"), Sys.getenv("OVARLONG_PATH"))
#    message("global instance of .ovarHubCon created")
#    }
#    dplyr::src_sql("sqlite", .ovarHubCon, path=Sys.getenv("OVARLONG_PATH"),
#          info = RSQLite::dbGetInfo(.ovarHubCon))
# }
# 
# OvMasterPheno = function() dplyr::tbl(OvarHub(), "masterPheno")
# OvarLong = function() dplyr::tbl(OvarHub(), "OvarLong")
# 
# tabulateSamples = function(tgen=OvarLong,
#   featureToCheck="A2M") {  # how to get an exemplar?
#   requireNamespace(package = "dplyr", quietly = TRUE)
#   stopifnot(is(featureToCheck,"character"))
#   stopifnot(length(featureToCheck)==1)
#   tgen() %>% select(sampleID, assaytype) %>%
#          filter(feature==featureToCheck) %>%
#          group_by(assaytype) %>% summarize(n=n())
# }
# 
# tabulateFeatures = function(tgen=OvarLong,
#   sampleToCheck="TCGA.04.1331") {
#   requireNamespace(package = "dplyr", quietly = TRUE)
#   stopifnot(is(sampleToCheck,"character"))
#   stopifnot(length(sampleToCheck)==1)
#   tgen() %>% select(sampleID, assaytype) %>%
#           filter(sampleID==sampleToCheck) %>%
#           group_by(assaytype) %>% summarize(n=n())
# }
# 
