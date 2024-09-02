\name{IMF_evo}
\alias{IMF_evo}
\alias{IMF_Kroupa_evo}
\alias{IMF_Lacey_evo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evoloving Kroupa Initial Mass Functions (IMF)
}
\description{
An evolving form of the differential Kroupa IMF. This can be used to encode a evolving IMF in a target simple stellar population (SSP) library with ProGeny.
}
\usage{
IMF_Kroupa_evo(mass, Age = 0, Age_lim = c(0, 13.8), alpha1_lim = c(0.3, 0.3),
  alpha2_lim = c(1.3, 1.3), alpha3_lim = c(2.3, 2), masslow_lim = c(0.01, 0.01),
  mass1_lim = c(0.08, 0.08), mass2_lim = c(0.5, 0.5), massmax_lim = c(150, 150),
  Lookback_Age = 0, massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25,
  method = 'linear', ...)

IMF_Lacey_evo(mass, Age = 0, Age_lim = c(0, 10), alpha1_lim = c(0.4, 2),
  alpha2_lim = c(0.4, 2), alpha3_lim = c(2.3, 2), masslow_lim = c(0.1, 0.1),
  mass1_lim = c(0.1, 0.1), mass2_lim = c(1, 1), massmax_lim = c(100, 100),
  Lookback_Age = 0, massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25,
  method = 'constant', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mass}{
Numeric vector or data.frame/list; the target stellar mass to be evaluated (Msol). If a numeric vector then it will be the value of the IMF at the mass specfied. If a data.frame/list it must have components mass$lo and mass$hi (same length), and the output is the IMF integral between these mass limits (either numeric integral if \option{massmult} = FALSE, or mass integral if \option{massmult} = TRUE).
}
  \item{Age}{
Numeric scalar; the linear age of the population in Gyrs. This basically informs the interpolation relative the range defined by \option{Age_lim}. The \code{\link{approxfun}} function is used to calculate the interpolation between \option{Age} given the vector of \option{Age_lim}  (which is treat as \option{x}). and the parameter of interest (which is treat as \option{y}).
}
  \item{Age_lim}{
Numeric vector; the lower and upper ages to interpolate between. \option{Age} values outside of this range are bounded to the \option{Age_lim} limits specified.
}
  \item{alpha1_lim}{
Numeric vector; the low mass slope (masslow < mass < mass1) of the Kroupa IMF at \option{Age_lim}. We interpolate over this domain using the specified \option{Age}.
}
  \item{alpha2_lim}{
Numeric vector; the mid mass slope (mass1 < mass < mass2) of the Kroupa IMF at \option{Age_lim}. We interpolate over this domain using the specified \option{Age}.
}
  \item{alpha3_lim}{
Numeric vector; the high mass slope (mass2 < mass < massmax) of the Kroupa IMF at \option{Age_lim}. We interpolate over this domain using the specified \option{Age}. This is the only evolving parameter by default, since we generally expect the IMF to become top-heavy with age (if anything). The lower value for older stars mean a flatter IMF, i.e. more mass is formed in massive stars.
}
  \item{masslow_lim}{
Numeric vector; the extreme low mass limit (Msol), where mass < masslow = 0, of the Kroupa IMF at \option{Age_lim}. \option{massform} is then the mass formed between \option{masslow} and \option{massmax}. We interpolate over this domain using the specified \option{Age}.
}
  \item{mass1_lim}{
Numeric vector; the low mass limit (Msol), of the Kroupa IMF at \option{Age_lim}. We interpolate over this domain using the specified \option{Age}. See \option{alpha1_lim}.
}
  \item{mass2_lim}{
Numeric vector; the mid mass limit (Msol), of the Kroupa IMF at \option{Age_lim}. We interpolate over this domain using the specified \option{Age}. See \option{alpha2_lim}.
}
  \item{massmax_lim}{
Numeric vector; the max mass limit (Msol), where mass > massmax = 0, of the Kroupa IMF at \option{Age_lim}. \option{massform} is then the mass formed between \option{masslow} and \option{massmax}. We interpolate over this domain using the specified \option{Age}.
}
  \item{Lookback_Age}{
Numeric scalar' the age of the Universe at the observation of interest in Gyrs. The default of 0 means we are observing the galaxy 'today', or at least at very low redshift.
}
  \item{massform}{
Numeric scalar; the mass formed between \option{masslow} and \option{massmax}, i.e. what the IMF should integrate to (Msol). Since the default is 1, this means the IMF returned is a true PDF and can be used as such.
}
  \item{massmult}{
Should the returns value be dN/dM (FALSE) or M.dN/dM (TRUE)? dN/dM is useful when evaluating the number of stars formed within limits, whilst M.dN/dM is used for evaluating the mass of stars formed within limits.
}
  \item{rel.tol}{
Numeric scaler; relative accuracy requested. See \code{\link{integrate}}.
}
  \item{method}{
Specifies the \code{\link{approxfun}} interpolation method to be used. Choices are "linear" or "constant".
}
  \item{\dots}{
Additional arguments to pass into \code{\link{approxfun}}.
}
}
\details{
These are potentially useful evolving forms of the IMF. We include a simple linearly evolving form fo the Kroupa IMF, where by default we interpolate linearly and smoothly along the age in Gyrs. We also have a step wise version of the Lacey evolving IMF. In more quiescent star formation this is determined to be Kennicutt 1983 like, and when star formation rate becomes bursty (so above a lookback of about 10 Gyrs, which is redshift 1.7) we transition sharply to a top-heavy single power law form (since \option{method} = 'constant'). For \code{IMF_Kroupa_evo} at 0 Gyrs (the default) we get the canonical Kroupa IMF shape, and at the earliest times in the Universe (13.8 Gyrs) we have a top-heavy (more massive stars) IMF.

Note for galaxies observed at different epochs the correct looback age needs to be accounted for. The simplest way to deal with this is to leave the \option{Age_lim} alone, and add the lookback age of the Universe as the \option{lookback_Age} argument. This is because if IMF evolves, it is probably relative to the age of the actual Universe, not relative to the age of stellar populations themselves (i.e a 10 Gyr stellar population today probably had a similar IMF to a 6 Gyr stellar population that formed at a lookback time of 4 Gyrs ago). In detail this does mean that you need to create a distinct evolving IMF SSP for every target observation, but in practice creating the library of SSPs that are in Gyr steps would suffice for most use cases (i.e. 0 < redshift < 0.1 can probably all use a common SSP).
}
\value{
If \option{massmult} = FALSE then the return is dN/dM for the specified \option{mass}, if \option{massmult} = TRUE then the return is M.dN/dM.
}
\references{
Kroupa, 2001, MNRAS, 322, 231
Lacey et al, 2016, MNRAS, 462, 3854
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyUpdateIMF}}
}
\examples{
magcurve(IMF_Kroupa_evo(x, Age=0), 0.01, 150, log='xy', lty=1, lwd=2, xlab = 'Star Mass / Msol', ylab='Relative Number Density')
curve(IMF_Kroupa_evo(x, Age=6.9), 0.01, 150, log='xy', lty=2, add=TRUE, lwd=2)
curve(IMF_Kroupa_evo(x, Age=13.8), 0.01, 150, log='xy', lty=3, add=TRUE, lwd=2)

curve(IMF_Lacey_evo(x, Age=0), 0.01, 150, log='xy', lty=1, add=TRUE, col='red', lwd=2)
curve(IMF_Lacey_evo(x, Age=6.9), 0.01, 150, log='xy', lty=2, add=TRUE, col='red', lwd=2)
curve(IMF_Lacey_evo(x, Age=13.8), 0.01, 150, log='xy', lty=3, add=TRUE, col='red', lwd=2)

legend('topright', legend=c('Kroupa Evo', 'Lacey Evo'), pch=16, col=c('black', 'red'))
legend('bottomleft', legend=c('0 Gyrs', '6.9 Gyrs', '13.8 Gyrs'), lty=1:3, lwd=2)
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
