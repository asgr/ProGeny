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
runShinyProGeny(app_mode = 'user',
  iso_path = "~/Google Drive/My Drive/ProGeny_isochrone/",
  atmos_path = "~/Google Drive/My Drive/ProGeny_atmos/", cores_server = 4)
}

\arguments{
  \item{app_mode}{
Character scalar; either 'user' or 'server'. The 'former 'user' option lets you return outputs at various stages, and lets you dynamically specify data locations for your own computer and choose the number of cores to run on. The 'server' option is more restricted to reflect how the Shiny App should run on a server.
}
  \item{iso_path}{
Character scalar; path to the isochrones to be used only in 'server' \option{mode}. Defaults just happen to work on ASGR's machine, so almost certainly need to be changed.
}
  \item{atmos_path}{
Character scalar; path to the atmospheres to be used only in 'server' \option{mode}. Defaults just happen to work on ASGR's machine, so almost certainly need to be changed.
}
  \item{cores_server}{
Integer scalar; the number of server cores to run on.
}
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
