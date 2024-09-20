\name{progenyMakeSSP}
\alias{progenyMakeSSP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High Level SSP Creator
}
\description{
The main interface via which users will likely use ProGeny to create SSPs.
}
\usage{
progenyMakeSSP(Iso, IMFfunc, ..., rem_frac = 'get', Spec_combine, Interp_combine,
  Zsol = 0.02, cores = 8, Labels = list(Zlab = "Metallicity",
  Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
  Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see ParsecIso in \code{\link{progenyIsoDownload}}).
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
TBD
}
\value{
A \code{ProSpect} compatible speclib (SSP) format looks like a list containing:

Z: Numeric vector of available template metallicities (ascending order).
Age: Numeric vector of available template ages in Yrs (ascending order).
AgeBins: Numeric vector of available template age bin limits in Yrs (ascending order, must be 1 longer than "Age" vector to cover lower and upper limits).
AgeWeights: Numeric vector of available template age weights (bin sizes) in Yrs (same length as "Age").
Wave: Numeric vector of spectral wavelength in Angstroms (ascending order).
Labels: List of the labels required for plotting (see \code{BC03}).
Zspec: List containing the spectra in units Lsun / Ang (for 1 Msun SF). The list must correspond to the 'Z' vector metallicities. Each list item is a numeric matrix of spectra with Nrow equal to the number of age bins (corresponding with 'Age') and Ncol the number of wavelength element (corresponding with 'Wave').
Zevo: List containing the mean evolutionary tracks of the SSPs (as determined by the route isochrones). The list must correspond to the 'Z' vector metallicities. Each list item is a data.frame with columns 'SMstar' (fraction of mass in luminous stars) 'SMgas' (fraction of mass in gas) 'SMtot' (fraction of total mass, always 1 in practice) 'SFR' (star formation rate, always 0 in practice) 'SMrem' (fraction of mass in stellar remnants, e.g. black holes etc). Nrow should be equal to the number of age bins (corresponding with 'Age').
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
