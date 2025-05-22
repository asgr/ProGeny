\name{runShinyProGeny}
\alias{runShinyProGeny}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive Tool to Build SSPs
}
\description{
A simple interactive tool to interactively build the various inputs needed to create a new SSP. Users get graphical feedback, making it clear what isochrone, spectral atmospheres, interpolation and IMF are being used.
}
\usage{
runShinyProGeny()
}
%- maybe also 'usage' for other objects documented here.
\details{
Run the the side effect of launching the interactive Shiny UI.
}
\value{
On each tab it is possible to return the specified outputs product alone (isochrones, atmospheres or interpolation grids), or in the last tab you can return the full SSP.
}
\references{
Robotham & Bellstedt, 2025, RASTI \cr
Bellstedt & Robotham, 2025, MNRAS
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyAtmosLoad}}, \code{\link{progenyMakeSSP}}, \code{\link{IMF}}
}
\examples{
# Just run it!
\dontrun{
output = runShinyProGeny()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
