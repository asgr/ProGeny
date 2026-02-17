\name{progenyIso2Spec}
\alias{progenyIso2Spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a Specific Age and Z Isochrone and IMF to a Stacked Spectrum
}
\description{
Given a full Isochrone table and IMF, this will subset by a target Age and Z and stacked all the best matching spectral atmospheres to make a single galaxy spectrum. Running this with all Age and Z combinations is basically how we construct a full SSP (with \code{\link{progenyMakeSSP}}).
}
\usage{
progenyIso2Spec(logAge = 8.4, logZ = 0, logA = NULL, Iso, IMFint, Spec_combine,
  Interp_combine, interp = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logAge}{
Numeric scalar; target log10 Age in Gyrs.
}
  \item{logZ}{
Numeric scalar; target log10 Z (relative to solar, where \option{logZ} = 0 means solar).
}
  \item{logA}{
Numeric scalar; target log10 Alpha abundance (relative to solar, where \option{logA} = 0 means solar abundance). If NULL it is ignored. This is only relevant when \option{Iso} contains a logA column (otherwise it is ignored). Note isochrone to spectral matching will have already been carried out when making \option{Interp_combine} (which can also use and interpolate \option{logA}), so the interpolation we are talking about here is just between isochrone \option{logA}. In practice all current isochrones have fixed abundance, but you could combine together the varying abundance BaSTI ischrones (-0.2, 0 and +0.4 dex options) to allow interpolation between values (+0.2 dex, say).
}
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see ParsecIso in \code{\link{progenyIsoDownload}}).
}
  \item{IMFint}{
A pre-integrated target IMF (made using \option{progenyUpdateIMF}, where you need the IMFint list element).
}
  \item{Spec_combine}{
List; the required spectral atmspheres. This must be loaded via \code{\link{progenyAtmosLoad}}, where the inputs for that come via \code{\link{progenyAtmosDownload}}.
}
  \item{Interp_combine}{
List; information that allows us to map Teff, logG, logZ information from the target \option{Iso} with the provided \option{Spec_combine}. This is made via \code{\link{progenyInterpGrid_All}}.
}
  \item{interp}{
Logical; should the spectrum be interpolated? This is relevant if you wish to create metallicity outputs that are not on the same grid as the input isochrone. TRUE means a weighted log-linear interpolation is carried out using the surrounding 4 values of logAge and logZ, FALSE means the nearest exact grid values are used regardless of what is requested.
}
}
\details{
This is a lower level function that can create a single realisation spectra for a target age and metallicity. The SSP generating function \code{\link{progenyMakeSSP}} basically loops over this function for the grid of all possible ages and metallicities available in the given isochrone.
}
\value{
Numeric vector; the stacked spectrum. Units will be Lsol/Ang.
}
\references{
Robotham & Bellstedt, 2025, RASTI, Vol, Page \cr
Bellstedt & Robotham, 2025, MNRAS, Vol, Page
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyMakeSSP}}
}
\examples{
\dontrun{
#Load isochrone:
Iso  = read.fst('path/to/Iso.fst', as.data.table=TRUE)

#Load stellar atmosheres:
Spec_combine = progenyAtmosLoad('path/to/atmos/')

#Interpolate stellar atmosheres onto isochrone grid:
Interp_combine = progenyInterpGrid_All(Iso=Iso, Spec_combine=Spec_combine)

#Select best stellar atmoshere library:
Iso = progenyInterpBest(Iso=Iso, Interp_combine=Interp_combine)

#Integrate the IMF:
IMFint = progenyUpdateIMF(Iso, IMF_Chabrier)

#Generate a single spectrum at Age = 1 Gyr and Z = Zsol:
Spec = progenyIso2Spec(logAge=9, logZ=0, Iso=Iso, IMFint=IMFint$IMFint,
  Spec_combine=Spec_combine, Interp_combine=Interp_combine)

#Plot it:
magplot(Spec_combine$base$wave, Spec, log='xy', type='l', xlim=c(1e3,1e5),
  xlab="Wavelength / Ang", ylab="Lsun / Ang (for 1 Msun SF)")
}
}
