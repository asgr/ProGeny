\name{progenyUpdateIMF}
\alias{progenyUpdateIMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrate IMF Between Mass Limits
}
\description{
This function only needs to be run when generating a ProGeny output for a new combination of Isochrones and/or IMF.
}
\usage{
progenyUpdateIMF(Iso, IMFfunc, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see ParsecIso in \code{\link{progenyIsoDownload}}).
}
  \item{IMFfunc}{
Function; the IMF to be analysed. \code{ProSpect} comes with IMF_Chabrier and IMF_Kroupa, which are good starting points. Must accept a \option{masslow}/\option{massmax} arguments as input.
}
  \item{\dots}{
Other arguments that can be accepted by the input \option{IMFfunc}.
}
}
\details{
TBD
}
\value{

\item{lo}{The low mass limit for the target bin.}
\item{hi}{The high mass limit for the target bin.}
\item{IMFint}{The integrated number counts inside the target IMF bin.}
}
\references{
Nothing yet
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyMakeSSP}}, \code{\link{IMF_evo}}
}
\examples{
## Nothing yet
}
