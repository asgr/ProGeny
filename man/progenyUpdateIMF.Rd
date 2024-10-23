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
This is a low level function that generally only need to be run if manually generating a single spectrum with \code{\link{progenyIso2Spec}}. In comparison (for simpler use) \code{\link{progenyMakeSSP}} runs this internally each time it is executed. This adds a small amount of overhead when making a lot of SSPs, but reduces user errors.
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
\dontrun{
#load iso
Iso  = read.fst('path/to/Iso.fst', as.data.table=TRUE)

#Integrate the IMF
IMFint = progenyUpdateIMF(Iso, IMF_Chabrier)
}
}
