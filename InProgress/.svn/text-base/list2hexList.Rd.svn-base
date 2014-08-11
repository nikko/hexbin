\name{list2hexList}
\alias{list2hexList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
list2hexList(binlst)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binlst}{ ~~Describe \code{binlst} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(binlst){
  if(length(binlst) < 2) stop(" need at least 2 hex bin objects")
  if(!all(unlist(lapply(binlst, is, "hexbin"))))
    stop("All Elements of list must be hexbin objects")
  if(!bnds.check(binlst))
    stop("All bin objects in list need the same xbnds and ybnds")
  if(!xbins.check(binlst))
    stop("All bin objects in list need the same number of bins")
  if(!shape.check(binlst))
    stop("All bin objects in list need the same shape parameter")
  mx <- max(unlist(lapply(binlst,function(h)max(h@count))))
  mn <- min(unlist(lapply(binlst,function(h)min(h@count))))           
  xbins <- binlst[[1]]@xbins
  xbnds <- binlst[[1]]@xbnds
  ybnds <- binlst[[1]]@ybnds
  shape <- binlst[[1]]@shape
  hl <- new("hexbinList",n=length(binlst),hbins=binlst, Xbnds=xbnds,
            Ybnds=ybnds, Xbins=xbins, Shape=shape,
            Bnames=names(binlst), CntBnds=c(mn,mx))
  hl
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
