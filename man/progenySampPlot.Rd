\name{progenySampPlot}
\alias{progenySampPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of Wavelength Sampling
}
\description{
Generate useful plots to check wavelength sampling for different stellar spectra libraries and SPLs/SSPs etc.
}
\usage{
progenySampPlot(wave, z = 0, v = NULL, h = NULL, add = FALSE, xlim = "auto",
  ylim = "auto", xlab = "Wavelength / Ang", ylab = "Sampling Resolution", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; wavelength in units of \option{xlab} (i.e. we assume this is Angstrom, but change the \option{xlab} as appropriate).
}
  \item{z}{
Numeric scalar; redshift to put wavelength grid at. This is useful when trying to see at what redshift you become SSP spectral resolution limited for a given instrument etc.
}
  \item{v}{
Numeric vector; vertical (wavelength) reference lines to add. These will be drawn with lty = 3 (dotted) and be red.
}
  \item{h}{
Numeric vector; horizontal (sampling resolution) reference lines to add. These will be drawn with lty = 3 (dotted) and be blue.
}
  \item{add}{
Logical; if TRUE then line is added, if FALSE then a new plot is created.
}
  \item{xlim}{
Numeroc vector; range of x data to display with format c(low_xlim, high_xlim) . Default of 'auto' will span the full range.
}
  \item{ylim}{
Numeroc vector; range of y data to display with format c(low_ylim, high_ylim) . Default of 'auto' will span 0 to the max sampling resolution.
}
  \item{xlab}{
Character scalar; x axis name.
}
  \item{ylab}{
Character scalar; y axis name.
}
  \item{\dots}{
Other arguments to pass to \code{\link{magplot}} (\option{add} = FALSE) or \code{\link{lines}} (\option{add} = TRUE).
}
}
\details{
Produces useful sampling resolution plots. This can be estimated directly from a wavelength grid, and is notably not the same thing as computing the spectral resolution (which requires deeper instrumental knowledge or high quality spectra to estimate).

As a warning, many people conflate the idea of spectral resolution and sampling resolution. Much like images, spectra have both a PSF component with a FWHM (usually called a spectral resolution) and a pixel sampling resolution (usually called the sampling resolution). Normally the pixel sampling will oversample the spectral sampling by a factor of a few, but this is not always the case. Minimally you need a SSP to at least target the spectral resolutions (since there is not much information below this scale), but it is often a good idea to target the sampling resolution too. Typical ratios are 1 (under-sampled) 3 (well sampled) and 5+ (over-sampled). E.g. 4MOST LRS targets at least an over-sampling factor of 2.5 at all wavelengths (but is often well above this).

Spectral resolution is closely related to the line spread function (LSF), and varies strongly with wavelength usually. \code{ProSpect} has functionality to use knowledge of the LSF when doing spectral fitting (read more there).
}
\value{
As well as the plot, outputs a data.frame with columns:

\item{wave}{Wavelength; same units as input via \option{wave}. Optionally shifted if \option{z} > 0.}
\item{res}{Sampling resolution.}

Note the number of rows will be 1 fewer than the length of the input \option{wave} because we can only compute the sampling between elements. This also means the values are different even when \option{z} = 0 because we take the mean the \option{wave} values where either the first or last element of the vector is ignored (soeffectively everything is offset by 1 and the extremes are ignored).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{magplot}}
}
\examples{
wave_grid_lin = seq(3e3, 1e4, by=1)
wave_grid_log = 10^seq(log10(3e3), log10(1e4), by=0.0001)

progenySampPlot(wave_grid_lin, col='red')
progenySampPlot(wave_grid_log, col='blue', add=TRUE)
}
\concept{ resolution }

