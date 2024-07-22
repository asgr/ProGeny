\name{progenyDownload}
\alias{progenyIsoDownload}
\alias{progenyAtmosDownload}
\alias{progenyAtmosLoad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download and Internally Load Isochrones and Stellar Atmospheres.
}
\description{
Functions to discover, download and load the various isochrone and stellar atmosphere libraries needed to generate SSPs.
}
\usage{
progenyIsoDownload(Iso = "avail", destpath = "",
  URL = "https://tinyurl.com/progeny_isochrone/")

progenyAtmosDownload(atmos = "avail", destpath = "",
  URL = "https://tinyurl.com/progeny_atmos/")

progenyAtmosLoad(destpath= '', base = 'combine_PHOENIX_husser',
  extend = 'combine_PHOENIX_allard', hot = 'combine_TMAP_werner',
  AGB = 'combine_AGB_lancon', white = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
Character scalar; either 'avail' which will then display all the isochrones currently available at \option{URL} and other key information, or the name of the target isochrone (that must be available!). Note all isochrones are saved as FST files.
}
  \item{atmos}{
Character scalar; either 'avail' which will then display all the stellar atmospheres currently available at \option{URL} and other key information, or the name of the target stellar atmosphere (that must be available!). Note all stellar atmospheres are saved as multi-extension FITS files.
}
  \item{destpath}{
Character scalar; where the target ProGeny isochrones / atmospheres should be downloaded to.
}
  \item{URL}{
Character scalar; location of the ProGeny isochrones / atmospheres libraries. This should not need to be altered by the user.
}
  \item{base}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the 'base' stellar atmosphere library to use. The 'base' library will be the primary (i.e. preferred) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. For \code{progenyAtmosPlot} should the 'base' grid of Teff versus logG be plotted?
}
  \item{extend}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'extend' stellar atmosphere library to use. The 'extend' library will be the secondary (i.e. used when \option{base} does not have a good match) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. For \code{progenyAtmosPlot} should the 'extend' grid of Teff versus logG be plotted?
}
  \item{hot}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'hot' stellar atmosphere library to use. The 'hot' library will be the preferred option for the hot star spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. For \code{progenyAtmosPlot} should the 'hot' grid of Teff versus logG be plotted?
}
  \item{AGB}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'AGB' stellar atmosphere library to use. The 'AGB' library will be the preferred option for the AGB spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}. For \code{progenyAtmosPlot} should the 'AGB' grid of Teff versus logG be plotts?
}
  \item{white}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'white' dwarf stellar atmosphere library to use. The 'white' dwarf library will be the preferred option for the white dwarf spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}. For \code{progenyAtmosPlot} should the 'white' grid of Teff versus logG be plotted?
}
}
\details{
TBD
}
\value{
\code{progenyIsoDownload} and \code{progenyAtmosDownload}. If \option{stellpop} = 'avail' then the contents of avail.txt will be displayed (using the default method). Otherwise, the target SSP rda/FITS file is downloaded and the destination path (+ filename) is returned.

\code{progenyAtmosLoad} creates the list sturcture used by \code{ProGeny} for combinations of stellar atomospheres. This will be the \option{Spec_combine} input for other functions.

\code{progenyAtmosPlot} creates a plot of all the Teff versus logG atmosphere grid components. This is really to check for isochrone versus spectral atmosphere coverage, where a particular target isochrone can easily be overplotted. Ideally you want all parts of an isochrone to fall inside the grid. Anything entirely outside will need to be approximated with a pure black body (which is rarely ideal).
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
