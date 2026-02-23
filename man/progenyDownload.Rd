\name{progenyDownload}
\alias{progenyDownload}
\alias{progenyIsoDownload}
\alias{progenyAtmosDownload}
\alias{progenyAtmosLoad}
\alias{progenyIsoFormat}
\title{
Download and Internally Load Isochrones and Stellar Atmospheres.
}
\description{
Functions to discover, download and load the various isochrone and stellar atmosphere libraries needed to generate SSPs.
}
\usage{
progenyIsoDownload(URL = "https://drive.google.com/drive/folders/1SnUUDgXGOiRZgcUZUKnWpG0WN0Lzg1Zs?usp=sharing", ...)

progenyAtmosDownload(URL = "https://drive.google.com/drive/folders/1oYS1JBQzP56cvTXuvDaw8b1ckwIZSnzw?usp=sharing", ...)

progenyAtmosLoad(destpath = '', base = 'combine_C3K_conroy',
  extend = 'combine_PHOENIX_allard', hot = 'combine_OB_PoWR',
  AGB = 'combine_AGB_lancon', white = 'combine_TMAP_werner',
  WR = "combine_WNE_PoWR", wavegrid = NULL, cores = 8, ...)

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
Character scalar; folder path where the target atmospheres will be loaded from. They all need to be in the same folder.
}
  \item{base}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the base stellar atmosphere library to use (without the '.fits' extension). The 'base' library will be the primary (i.e. preferred) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmospheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{extend}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the extend stellar atmosphere library to use (without the '.fits' extension). The 'extend' library will be the secondary (i.e. used when \option{base} does not have a good match) \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmospheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{hot}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the hot stellar atmosphere library to use (without the '.fits' extension). The 'hot' library will be the preferred option for the hot star spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmospheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}.
}
  \item{AGB}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the AGB stellar atmosphere library to use (without the '.fits' extension). The 'AGB' library will be the preferred option for the AGB spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmospheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
  \item{white}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the white dwarf stellar atmosphere library to use (without the '.fits' extension). The 'white' dwarf library will be the preferred option for the white dwarf spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmospheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
  \item{WR}{
Character scalar/logical; for \code{progenyAtmosLoad} name of the Wolf-Rayet stellar atmosphere library to use (without the '.fits' extension). The 'WR' library will be the preferred option for the Wolf-Rayet spectra \option{Teff}, \option{logG} and \option{logZ} solution when interpolating the atmospheres onto the isochrone grids. See \code{\link{progenyInterpGrid}}. This can be flagged in the isochrones with the \option{label} column, see \code{\link{progenyInterpBest}}.
}
  \item{wavegrid}{
Numeric vector; the target spectral wavelength grid in Angstroms. If NULL then all inputs must have the same native wavelength grid or bad things will happen. Basic checks are made though (that they are at least the same length). If not NULL then all stellar atmosphere spectra are rebinned using the \code{ProSpect} \code{specReBin} function, so \code{ProSpect} must also be installed in this case.
}
  \item{cores}{
Integer scalar; the number of cores to run on.
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

\code{progenyAtmosLoad} creates the list structure used by \code{ProGeny} for combinations of stellar atmospheres. This will be the \option{Spec_combine} input for other functions.

\code{progenyIsoFormat} will correctly sort the columns of the isochrone: should be logZ, logAge, Mini ordered [slowest to fastest changing]). Returns the new Iso as the output.
}
\references{
Robotham & Bellstedt, 2025, RASTI, 4, 19 \cr
Bellstedt & Robotham, 2025, MNRAS, 540, 2703
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyMakeSSP}}
}
\examples{
\dontrun{
progenyIsoDownload()
Iso = read.fst('path/to/Iso.fst')
Iso = progenyIsoFormat(Iso)

progenyAtmosDownload()
SpecCombine = progenyAtmosLoad('path/to/atmos/')

#We can plot some atmospheres coloured by temperature:

Temp_vec = seq(2000,50000,by=1000)
col_vec = hcl.colors(length(Temp_vec), 'Zissou 1', rev=TRUE)

magplot(NA, NA, type='n', log='xy', xlim=c(3e2,1e5), ylim=c(1e-10,1),
  xlab='Wavelength / Ang', ylab = expression(L[lambda]))
for(i in seq_along(Temp_vec)){
  loc = SpecCombine$base$info[Teff==Temp_vec[i] & logG == 4 & logZ == 0,which=T]
  if(length(loc) > 0){
    lines(SpecCombine$base$wave, SpecCombine$base$spec[loc,], col=col_vec[i])
  }
}

magbar('topright', range=c(2000,50000), col=col_vec, title='Temp / K', titleshift=0.5)
}
}
