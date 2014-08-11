\name{hexList}
\alias{hexList}
\alias{hexbinList-class}
\alias{coerce,list,hexbinList-method}
\title{Conditional Bivariate Binning into Hexagon Cells }
\description{
  Creates a list of \code{\link{hexbin}} objects.  Basic components are
  a cell id and a count of points falling in each occupied cell.
  Basic methods are \code{\link[methods]{show}()}, \code{plot()} %(\link{plot.hexbin})
  and \code{\link{summary}()}, but also \code{\link{erode}}.
  % .. \code{\link{smooth.hexbin}}
}
\usage{
hexList(x, y = NULL, given = NULL, xbins = 30, shape = 1,
        xbnds = NULL, ybnds = NULL, xlab = NULL, ylab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{given}{ ~~Describe \code{given} here~~ }
  \item{xbins}{ ~~Describe \code{xbins} here~~ }
  \item{shape}{ ~~Describe \code{shape} here~~ }
  \item{xbnds}{ ~~Describe \code{xbnds} here~~ }
  \item{ybnds}{ ~~Describe \code{ybnds} here~~ }
  \item{xlab}{ ~~Describe \code{xlab} here~~ }
  \item{ylab}{ ~~Describe \code{ylab} here~~ }
}
\details{
  There is also a \code{\link[methods:as]{coerce}} method to produce
  \code{hexbinList} objects from \code{\link{list}}s.
  %% i.e., \code{as(list(....), "hexbinList")} should work
}
\value{

  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\author{Nicholas Lewin-Koh}
% \note{ ~~further notes~~ }

\seealso{\code{\link{hexbin}}, \code{\link{hdiffplot}} }
\examples{

}
\keyword{dplot}
\keyword{misc}
