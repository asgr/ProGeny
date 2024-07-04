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
TBD
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
## Nothing yet
}
