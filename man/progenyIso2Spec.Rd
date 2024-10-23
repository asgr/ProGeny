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
progenyIso2Spec(logAge = 8.4, logZ = 0, Iso, IMFint, Spec_combine, Interp_combine)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logAge}{
Numeric scalar; target log10 Age in Gyrs.
}
  \item{logZ}{
Numeric scalar; target log10 Z (relative to solar, where \option{logZ} = 0 means solar).
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
}
\details{
This is a lower level function that can create a single realisation spectra for a target age and metallicity. The SSP generating function \code{\link{progenyMakeSSP}} basically loops over this function for the grid of all possible ages and metallicities available in the given isochrone.
}
\value{
Numeric vector; the stacked spectrum. Units will be Lsol/Ang.
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
