\name{downloadData}
\alias{downloadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download all the data that used to be in extdata/.
}
\description{
Download all the data that used to be in extdata/, from the Amazon S3 bucket.
}
\usage{
downloadData(datadir="downloaddir")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datadir}{
datadir is the directory to where the S3 bucket will be dumped.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
none
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Levi Waldron
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (datadir) 
{
    s3url <- "http://s3.amazonaws.com/cancerhub/"
    con = url(paste0(s3url, "filelist.txt"))
    filedf <- read.table(con, as.is = TRUE)
    fileurls <- filedf[, 4]
    filedests <- file.path(datadir, fileurls)
    fileurls <- paste0(s3url, fileurls)
    for (i in 1:length(fileurls)) {
        if (!file.exists(filedests[i])) {
            dir.create(dirname(filedests[i]), showWarnings = FALSE)
            download.file(url = fileurls[i], destfile = filedests[i])
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
