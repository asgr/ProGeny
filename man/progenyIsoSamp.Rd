\name{progenyIsoSamp}
\alias{progenyIsoSamp}
\title{
Sample a Subset of Stellar Spectra
}
\description{
A helper function to create small libraries of representative stellar spectra given various observing constraints. This was originally made to generate templates for stellar fitting in WAVES.
}
\usage{
progenyIsoSamp(Iso, IMFfunc, ..., Nsamp = 10000, mag_lim = NULL, rad_lim = NULL,
  Lum_weight_func = NULL, logAge_weight_func = NULL, logZ_weight_func = NULL,
  label_use = seq(-1, 9), mag_Lsol_1kpc = 14.65, seed = 666)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see e.g. MistIso in \code{\link{progenyIsoDownload}}).
}
  \item{IMFfunc}{
Function; the IMF to be analysed. \code{ProGeny} comes with \code{\link{IMF_Chabrier}}, \code{\link{IMF_Kroupa}} and \code{\link{IMF_Salpeter}}, which are good starting points. Must accept \option{masslow}/\option{massmax} arguments as input.
}
  \item{\dots}{
Other arguments that can be accepted by the input \option{IMFfunc}.
}
  \item{Nsamp}{
Integer scalar; the number of samples to make.
}
  \item{mag_lim}{
Numeric scalar; approximate faint AB mag limit of selection (see \option{mag_Lsol_1kpc}).
}
  \item{rad_lim}{
Numeric scalar; approximate radial limit of selection (kpc). This should be sensible for the Milky Way if provided probably, so a few kpc at most.
}
  \item{Lum_weight_func}{
Function; optional additional weighting to apply based on the isochrone luminosity (in linear Lsol units).
}
  \item{logAge_weight_func}{
Function; optional additional weighting to apply based on the isochrone age (in log10(Age/Yr) units).
}
  \item{logZ_weight_func}{
Function; optional additional weighting to apply based on the isochrone metallicity (in log10(Z/Zsol) units).
}
  \item{label_use}{
Integer vector; subset of isochrone phase labels to use. E.g. if you only want main sequence stars specify as needed here (this would be 0 for MIST isochrones).
}
  \item{mag_Lsol_1kpc}{
Numeric scalar; rought scaling for the nominal filter implied by \option{mag_lim}. The default is the magnitude of our own Sun at a distance of 1 kpc in the r-band, so if using an optical filter this can be left as is. For other conversions apply the difference between the absolute AB magnitude of the Sun in the r-band (~4.65) to the new band as found here \href{https://mips.as.arizona.edu/~cnaw/sun.html}{https://mips.as.arizona.edu/~cnaw/sun.html}.
}
  \item{seed}{
Integer scalar; random seed to use for sampling.
}
}
\details{
Mostly useful for making small libraries of stellar spectra to represent different populations of stars.
}
\value{
The correctly sampled isochrone, which probably then needs passing through \code{\link{progenyInterpGrid_All}} and \code{\link{progenyInterpBest}} to actually match the subset isochrone with the stellar spectral atmospheres.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyInterpGrid_All}}, \code{\link{progenyInterpBest}}
}
\examples{
\dontrun{
#Load isochrone:
Iso  = read.fst('path/to/Iso.fst', as.data.table=TRUE)

#Sample isochrone (roughly mimic the Milk Way):
MistIso_sub = progenyIsoSamp(Iso, IMF_Chabrier, mag_lim=20, rad_lim=4,
  logAge_weight_func = function(x){dnorm(x,mean=9.6, sd=0.3)},
  logZ_weight_func = function(x){dnorm(x,mean=0,sd=0.3)},
  Nsamp=1e3)
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
