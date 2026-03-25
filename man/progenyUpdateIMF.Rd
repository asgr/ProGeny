\name{progenyUpdateIMF}
\alias{progenyUpdateIMF}
\title{
Integrate IMF Between Mass Limits
}
\description{
This function only needs to be run when generating a ProGeny output for a new combination of Isochrones and/or IMF.
}
\usage{
progenyUpdateIMF(Iso, IMFfunc, max_width = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see ParsecIso in \code{\link{progenyIsoDownload}}).
}
  \item{IMFfunc}{
Function; the IMF to be analysed. \code{ProGeny} comes with \code{\link{IMF_Chabrier}}, \code{\link{IMF_Kroupa}} and \code{\link{IMF_Salpeter}}, which are good starting points. Must accept \option{masslow}/\option{massmax} arguments as input.
}
  \item{max_width}{
Numeric scalar; if NULL this does nothing, if set to numeric scalar then this is the upper limit allowed in any given bin width during IMF integration. This prevents pathological behavious when you get big gaps in Mini masses during phase transitions. Since bins are processed in log space, a value of 0.1/0.2 would be a factor 1.25/1.58 (these are sensible limits to use based on good quality ProGeny isochrones.
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
Robotham & Bellstedt, 2025, RASTI, 4, 19 \cr
Bellstedt & Robotham, 2025, MNRAS, 540, 2703
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyMakeSSP}}, \code{\link{IMF_evo}}
}
\examples{
\dontrun{
#Load isochrone:
Iso  = read.fst('path/to/Iso.fst', as.data.table=TRUE)

#Integrate the IMF:
IMFint = progenyUpdateIMF(Iso, IMF_Chabrier)
}
}
