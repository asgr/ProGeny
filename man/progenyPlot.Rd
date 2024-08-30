\name{progenyPlot}
\alias{progenyPlot}
\alias{progenyAtmosPlot}
\alias{progenyIsoPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Progeny Plotting
}
\description{
Create useful diagnostic plots of ProGeny related objects.
}
\usage{
progenyAtmosPlot(Spec_combine, add = FALSE, do_base = TRUE, do_extend = TRUE,
  do_hot = TRUE, do_AGB = TRUE, do_white = TRUE,
  col = c("black", "grey", "blue", "red", "darkgreen"), pch = 16, cex = 1,
  xlim = c(2000, 4e+05), ylim = c(-1, 9), log = "x", ...)

progenyIsoPlot(Iso, add = FALSE, sampleN = 1e4, seed = 666, zsel = 'Lum', zunit = 'Lsol',
  zstretch = 'log', col = hcl.colors(101, "geyser"), pch = 16, cex = 0.5,
  xlim = c(2e3,4e5), ylim = c(-1,9), log = 'x', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Spec_combine}{
List; the required spectral atmspheres. This must be loaded via \code{\link{progenyAtmosLoad}}, where the inputs for that come via \code{\link{progenyAtmosDownload}}.
}
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see ParsecIso in \code{\link{progenyIsoDownload}}).
}
  \item{add}{
Logical; should information be added to the current plot device (TRUE) or should a fresh plot be made (FALSE).
}
  \item{sampleN}{
Numeric/character scalar; the sparse sampling of the \option{Iso} being plotted. Since these are large you generally do not want to plot all points. Setting to 'all' will plot all points.
}
  \item{seed}{
Integer scalar; random seed to use when sparse sampling with \option{sampleN}.
}
  \item{zsel}{
Character scalar; the column to use for z axis colouring.
}
  \item{zunit}{
Character scalar; the units for \option{zsel} (used for plot labelling).
}
  \item{zstretch}{
Character scalar; the stretch to apply to \option{zsel}, either 'lin' (for linear) or 'log' (for log).
}
  \item{do_base}{
Logical; should the 'base' spectral atmospheres be plotted?
}
  \item{do_extend}{
Logical; should the 'extend' spectral atmospheres be plotted?
}
  \item{do_hot}{
Logical; should the 'hot' spectral atmospheres be plotted?
}
  \item{do_AGB}{
Logical; should the 'AGB' spectral atmospheres be plotted?
}
  \item{do_white}{
Logical; should the 'white' spectral atmospheres be plotted?
}
  \item{col}{
Various; either colours vectors (for atmospheres) or colour palette (for isochrones) to use.
}
  \item{pch}{
The point 'pch'.
}
  \item{cex}{
The point 'cex'.
}
  \item{xlim}{
Vector; range of x data to display.
}
  \item{ylim}{
Vector; range of y data to display.
}
  \item{log}{
Character scalar; log axis arguments to be passed to plot. E.g. use 'x', 'y', 'xy' or 'yx' as appropriate. '' assumes no logging of any axes.
}
  \item{\dots}{
Additional arguments to pass to the \code{\link{points}} plotting function.
}
}
\details{
TBD
}
\value{
Nothing, run for the side effect of generating plot (or plot points).
}
\references{
Nothing yet
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyMakeSSP}}
}
\examples{
## Nothing yet
}
