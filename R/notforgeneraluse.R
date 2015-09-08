# I guess eventually this won't be exported, but need it just for now.
downloadData <- function(datadir){
  s3url <- "http://s3.amazonaws.com/cancerhub/"
  con=url(paste0(s3url, "filelist.txt"))
  filedf <- read.table(con, as.is=TRUE)
  fileurls <- filedf[, 4]
  filedests <- file.path(datadir, fileurls)
  fileurls <- paste0(s3url, fileurls)
  for (i in 1:length(fileurls)){
    if(!file.exists(filedests[i])){
      dir.create(dirname(filedests[i]), showWarnings=FALSE)
      download.file(url=fileurls[i], destfile=filedests[i])
    }}
}