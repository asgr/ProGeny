\name{progenyMakeSSP}
\alias{progenyMakeSSP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High Level SSP Creator
}
\description{
The main interface via which users will likely use ProGeny to create SSPs. This can either create SSPs on the exact logAge and logZ grid that exists within the provided isochrone, or be interplated onto an arbitrary logAge and logZ grid.
}
\usage{
progenyMakeSSP(Iso, IMFfunc = IMF_Chabrier, ..., rem_frac = 'get', Spec_combine,
  Interp_combine, logAge_steps = NULL, logZ_steps = NULL, Mini_range = NULL, Zsol = 0.02,
  cores = 8, Labels = list(Zlab = "Metallicity", Agelab = "Time since ZAM / Yrs",
  Wavelab = "Wavelength / Ang", Lumlab = "Lsun / Ang (for 1 Msun SF)",
  LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see e.g. MistIso in \code{\link{progenyIsoDownload}}).
}
  \item{IMFfunc}{
Function; the IMF to be analysed. \code{ProSpect} comes with IMF_Chabrier and IMF_Kroupa, which are good starting points. Must accept a \option{masslow}/\option{massmax} arguments as input.
}
  \item{\dots}{
Other arguments that can be accepted by the input \option{IMFfunc}.
}
  \item{rem_frac}{
Character/numeric scalar; the fraction of stellar mass that remans locked up in remnants after a star leaves the main sequence. If 'get' then this is estimated direclty from the isochrones, but this only works if the isochrones at least partially follow the mass loss phases of the star (MIST and PARSEC). Other the value needs to be specified directly (usually somewhere between 0.5 - 0.6 in practice). This is necessary for BaSTI_New, which terminates at the point of mass loss (so this phase is not well tracked).
}
  \item{Spec_combine}{
List; the required spectral atmspheres. This should probably be loaded via \code{\link{progenyAtmosLoad}}, where the inputs for that come via \code{\link{progenyAtmosDownload}}.
}
  \item{Interp_combine}{
List; information that allows us to map Teff, logG, logZ information from the target \option{Iso} with the provided \option{Spec_combine}. This is made via \code{\link{progenyInterpGrid_All}}.
}
  \item{logAge_steps}{
Numeric vector; if provided it is the log10(Age/yr) values that you want the output SSP to be interpolated at. If this argument is NULL we use the sorted unique logAge values from the \option{Iso} provided (the default). This means you can generate SSPs with fairly arbitrary age grids, with the caveat the \option{logAge_steps} must be fully inside the logAge extremes in the provided \option{Iso} (we can only interpolate, not extrapolate).
}
  \item{logZ_steps}{
Numeric vector; if provided it is the log10(Z/Zsol) values that you want the output SSP to be interpolated at. If this argument is NULL we use the sorted unique logZ values from the \option{Iso} provided (the default). This means you can generate SSPs with fairly arbitrary metallicity grids, with the caveat the \option{logZ_steps} must be fully inside the logZ extremes in the provided \option{Iso} (we can only interpolate, not extrapolate).
}
  \item{Mini_range}{
Numeric vector; optional. Lower and upper ranges of \option{Iso} Mini to use. This can be useful if wanting to combine mass ranges from different isochrone to make a better combined SSP (say one is better for low mass stars, and another better for high mass stars). If provided it must be length 2.
}
  \item{Zsol}{
Numeric scalar; value of the Sun's solar metallicity (to convert logZ to Z).
}
  \item{cores}{
Integer scalar; how many cores to run on.
}
  \item{Labels}{
Character vector; labels for the \code{ProSpect} format speclib (SSP). Probably leave as they are, unless you really know what you are doing.
}
}
\details{
For best accuracy the \option{logAge_steps} and \option{logZ_steps} arguments should not be used, since spectral interpolation is never perfect. Internally the interpolation is done in log-linear age and metallicity bilinear weighting, since this does the best job of approximating a target value of logAge and logZ during testing.
}
\value{
A \code{ProSpect} compatible speclib (SSP) format looks like a list containing:

\item{Z}{Numeric vector of available template metallicities (ascending order).}
\item{Age}{Numeric vector of available template ages in Yrs (ascending order).}
\item{AgeBins}{Numeric vector of available template age bin limits in Yrs (ascending order, must be 1 longer than "Age" vector to cover lower and upper limits).}
\item{AgeWeights}{Numeric vector of available template age weights (bin sizes) in Yrs (same length as "Age").}
\item{Wave}{Numeric vector of spectral wavelength in Angstroms (ascending order).}
\item{Labels}{List of the labels required for plotting (see \code{BC03}).}
\item{Zspec}{List containing the spectra in units Lsun / Ang (for 1 Msun SF). The list must correspond to the 'Z' vector metallicities. Each list item is a numeric matrix of spectra with Nrow equal to the number of age bins (corresponding with 'Age') and Ncol the number of wavelength element (corresponding with 'Wave').}
\item{Zevo}{List containing the mean evolutionary tracks of the SSPs (as determined by the route isochrones). The list must correspond to the 'Z' vector metallicities. Each list item is a data.frame with columns 'SMstar' (fraction of mass in luminous stars) 'SMgas' (fraction of mass in gas) 'SMtot' (fraction of total mass, always 1 in practice) 'SFR' (star formation rate, always 0 in practice) 'SMrem' (fraction of mass in stellar remnants, e.g. black holes etc). Nrow should be equal to the number of age bins (corresponding with 'Age').}
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

#Note above you might want to set label_AGB and label_white to restrict how AGB
#and white dwarf atmospheres are added.

#Finally we can generate a full ProSpect ready SSP:

SSP = progenyMakeSSP(Iso=Iso, IMFfunc=IMF_Chabrier, Spec_combine=Spec_combine,
  Interp_combine=Interp_combine)

#We can check this with the ProSpect function:

library(ProSpect)
speclib_check(SSP)
}
}
