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

progenyAtmosLoad(destpath= '',
  base = 'combine_PHOENIX_alpha+0.0', extend = 'combine_basel_cut', hot = 'combine_hot',
  AGB = 'combine_AGB', white = 'combine_white')
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
Character scalar; name of the 'base' stellar atmosphere library to use. The 'base' library will be the primary (i.e. preferred) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{extend}{
Character scalar; name of the 'extend' stellar atmosphere library to use. The 'extend' library will be the secondary (i.e. used when \option{base} does not have a good match) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{hot}{
Character scalar; name of the 'hot' stellar atmosphere library to use. The 'hot' library will be the preferred option for the hot star spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{AGB}{
Character scalar; name of the 'AGB' stellar atmosphere library to use. The 'AGB' library will be the preferred option for the AGB spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
  \item{white}{
Character scalar; name of the 'white' dwarf stellar atmosphere library to use. The 'white' dwarf library will be the preferred option for the white dwarf spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmopsheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
}
\details{
TBD
}
\value{
\code{progenyIsoDownload} and \code{progenyAtmosDownload}. If \option{stellpop} = 'avail' then the contents of avail.txt will be displayed (using the default method). Otherwise, the target SSP rda/FITS file is downloaded and the destination path (+ filename) is returned.

\code{progenyAtmosLoad} creates the list sturcture used by \code{ProGeny} for combinations of stellar atomospheres. This will be the \option{Spec_combine} input for other functions.
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
