\name{progenyDownload}
\alias{progenyIsoDownload}
\alias{progenyAtmosDownload}
\alias{progenyAtmosLoad}
\alias{progenyIsoFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download and Internally Load Isochrones and Stellar Atmospheres.
}
\description{
Functions to discover, download and load the various isochrone and stellar atmosphere libraries needed to generate SSPs.
}
\usage{
progenyIsoDownload(URL = "https://tinyurl.com/progeny-isochrone/", ...)

progenyAtmosDownload(URL = "https://tinyurl.com/progeny-atmos/", ...)

progenyAtmosLoad(destpath= '', base = 'combine_C3K_conroy',
  extend = 'combine_PHOENIX_allard', hot = 'combine_TMAP_werner',
  AGB = 'combine_AGB_lancon', white = NULL, wavegrid = NULL, ...)

progenyIsoFormat(Iso)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{URL}{
Character scalar; location of the ProGeny isochrones / atmospheres libraries. This should not need to be altered by the user.
}
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones.
}
  \item{destpath}{
Character scalar; folder path where the target atomospheres will be laoded from. They all need to be in the same folder.
}
  \item{base}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the 'base' stellar atmosphere library to use (without the '.fits' extension). The 'base' library will be the primary (i.e. preferred) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{extend}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'extend' stellar atmosphere library to use (without the '.fits' extension). The 'extend' library will be the secondary (i.e. used when \option{base} does not have a good match) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{hot}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'hot' stellar atmosphere library to use (without the '.fits' extension). The 'hot' library will be the preferred option for the hot star spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{AGB}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'AGB' stellar atmosphere library to use (without the '.fits' extension). The 'AGB' library will be the preferred option for the AGB spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
  \item{white}{
Character scalar/logical; for \code{progenyAtmosLoad}  name of the 'white' dwarf stellar atmosphere library to use (without the '.fits' extension). The 'white' dwarf library will be the preferred option for the white dwarf spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
  \item{wavegrid}{
Numeric vector; the target spectral wavelength grid in Angstroms. If NULL then all inputs must have the same native wavelength grid or bad things will happen. Basic checks our made though (that they are at least the same length). If not NULL then all stellar atmosphere spectra are rebinned using the \code{ProSpect} \code{specReBin} function, so \code{ProSpect} must also be isntalled in this case.
}
  \item{\dots}{
For download functions, additional arguments to pass to \code{\link{browseURL}}. For \code{progenyAtmosLoad}, additional arguments to pass to \code{ProSpect} \code{specReBin}. Only relevant if \option{wavegrid} is not NULL.
}
}
\details{
The download functions are currently just re-directs to the DropBox page hosting the data. The direct download functions that I use in \code{ProSpect} (\code{speclib_download}) do not really work properly anymore because DropBox changed the way it create shared folders and links. Oh well, this will have to do for now!
}
\value{
\code{progenyIsoDownload} and \code{progenyAtmosDownload} just take users to the relevant DropBox link currently.

\code{progenyAtmosLoad} creates the list sturcture used by \code{ProGeny} for combinations of stellar atomospheres. This will be the \option{Spec_combine} input for other functions.

\code{progenyAtmosPlot} creates a plot of all the Teff versus logG atmosphere grid components. This is really to check for isochrone versus spectral atmosphere coverage, where a particular target isochrone can easily be overplotted. Ideally you want all parts of an isochrone to fall inside the grid. Anything entirely outside will need to be approximated with a pure black body (which is rarely ideal).

\code{progenyIsoFormat} will correctly sort the columns of the isochrone: should be logZ, logAge, Mini ordered [slowest to fastest changing]). Returns the new Iso as the output.
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
