\name{progenyParamWarp}
\alias{progenyParamWarp}
\alias{progenyWarpInterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Warp [x,y] Source Data to Target
}
\description{
A common issue when comparing isochrones is various parameters might behave similarly but not identically. This function allows as to warp the x-axis of source data to create optimal overlap in source and target y-axis.
}
\usage{
progenyParamWarp(x_src, y_src, x_tar, y_tar, smooth = FALSE, check_order = FALSE,
  open.end = TRUE, open.begin = FALSE, ...)

progenyWarpInterp(x_src, y_src, x_tar, y_tar, ParamWarp_out, wt = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_src}{
Numeric vector / matrix / data.frame; x source data we want to warp.
}
  \item{y_src}{
Numeric vector / matrix / data.frame; y source data we want to maximially overlap with \option{t_tar} after warping. Must be same length as \option{x_src} if vector. If matrix then rows must match to \option{x_src} and columns are the multi-variate parameters we want to warp against.
}
  \item{x_tar}{
Numeric vector / matrix / data.frame; x target data.
}
  \item{y_tar}{
Numeric vector / matrix / data.frame; y target data we want to maximially overlap with \option{y_src} after warping. Must be same length as \option{x_tar}. If matrix then rows must match to \option{x_tar} and columns are the multi-variate parameters we want to warp against.
}
  \item{smooth}{
Logical; should smoothing be applied to the inputs? Usually some smoothing is a good thing, but features can also be lost which might make warping harder (probably worth trying both options).
}
  \item{check_order}{
Logical; should \option{x_src} and \option{x_tar} be checked for sorting. If already sorted (either ascending or descending is fine) leave this as FALSE since it will run faster.
}
  \item{open.end}{
Logical; can the end points be open (TRUE) or should they be forced to align (FALSE). For isochrones this usually need to be set to TRUE.
}
  \item{open.begin}{
Logical; can the beginning points be open (TRUE) or should they be forced to align (FALSE). For isochrones this generally this does not need to be changed from FALSE.
}
  \item{ParamWarp_out}{
List, the output of \code{progenyParamWarp}.
}
  \item{wt}{
Numeric scalar; how biased towards the target parameters the output is. \option{wt} = 0 would return the src curve, and \option{wt} = 1 would return the tar curve.
}
  \item{\dots}{
Other arguments to pass to \code{dtw}.
}
}
\details{
This function uses the \code{dtw} package and function to uncover the optimal warping (or re-mapping) between source and target data. The end result are two warping functions that will convert source x data to target. The warping is only done in the x-dimension, and the y-dimension is kept intact.

The main use case for this function is to allow isochrone properties (like Mini and Lum) to be distorted between isochrones so properties can be compared/used. A good option for \option{y_src} and \option{y_tar} could be a combination of input like log10(Lum), log10(Teff), logG. This makes us sensitive to most of the evoluationary features when interpolating between isochrone ages.

Note source inputs (\option{x_src} and \option{y_src}) do no need to be the same length as the target inputs (\option{x_tar} and \option{y_tar}).

Witht the output \code{warp_src2tar} and \code{warp_tar2src} function, the \option{wt} parameter describes how biased to the RHS the output is, so in \code{warp_src2tar} \option{wt} = 1 would give purely tar mapping, whereas \option{wt} = 0 would simply return the input x.

For \code{progenyWarpInterp} the \option{x_src} and \option{x_tar} should be the same parameter as passed to \code{progenyParamWarp} (e.g. \option{Mini} is typical), but \option{y_src} and \option{y_tar} can be totally different. In this sense we can use the computed x-dimension warp to create a new Teff, Lum or logG curve. \option{y_src} and \option{y_tar} can be matrices / data.frames, meaning we can conveniently warp multiple parameters (Teff, Lum or logG), which is useful for isochrones.
}
\value{
  \item{src}{Data.frame; combined (cbind) \option{x_src} and \option{y_src}.}
  \item{tar}{Data.frame; combined (cbind) \option{x_tar} and \option{y_tar}.}
  \item{warp_src2tar}{Function to convert \option{x_src} to \option{x_tar} in a manner that produces the maximum overlap of \option{y_src} and \option{y_tar}.}
  \item{warp_tar2src}{Function to convert \option{x_tar} to \option{x_src} in a manner that produces the maximum overlap of \option{y_src} and \option{y_tar}.}
  \item{dtw_alignment}{Full output of \code{dtw}.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyTrackInterp}}
}
\examples{
library(magicaxis)

# Source curve
x_src = seq(0, 10, length.out = 400)
y_src = sin(x_src) + 0.2 * cos(3 * x_src)

# Nonlinear warp + noise for target
x_tar = seq(0, 10, len = 800)
x_tar_warped = 10 * ( (x_tar / 10)^1.3 )  # monotone distortion
ysp_src = smooth.spline(x_src, y_src, spar = 0.6)
y_tar   = predict(ysp_src, x_tar_warped)$y + rnorm(length(x_tar), sd = 0.05)

fit = progenyParamWarp(x_src, y_src, x_tar, y_tar)

# Original curves
magplot(x_src, y_src, type = "l", col = "blue", lwd = 2,
     main = "Source and target (original)",
     xlab = "x", ylab = "y")
lines(x_tar, y_tar, col = "red", lwd = 2)
lines(fit$warp_src2tar(x_src), y_src, col='blue', lty=2, lwd=2)
legend("bottomleft", legend = c("Source", "Target", "Source Warp"),
       col = c("blue", "red", "blue"), lty = c(1,1,2), bty = "n")

magplot(fit$warp_src2tar, xlim=c(0,10), xlab='x_src', ylab='x_tar')
abline(0,1,lty=2)

\dontrun{
#this is an example with the pre-made ProGeny MIST isochrone, were we want to interpolate between ages:
library(fst)
library(foreach)
library(data.table)
library(magicaxis)

Iso_MIST = read.fst('/path/to/ProGeny_isochrone/MistIso.fst', as.data.table = TRUE)

# Imagine we have logAge = 9 and 9.2 available:
src = Iso_MIST[abs(logAge - 9) < 1e-4 & logZ==0, list(Mini, logL=log10(Lum), logT=log10(Teff), logG)]
tar = Iso_MIST[abs(logAge - 9.2) < 1e-4 & logZ==0, list(Mini, logL=log10(Lum), logT=log10(Teff), logG)]
# And we wish to create logAge = 9.1:
cmp = Iso_MIST[abs(logAge - 9.1) < 1e-4 & logZ==0, list(Mini, logL=log10(Lum), logT=log10(Teff), logG)]

#Fit it:
fit = progenyParamWarp(x_src=src, x_tar=tar, smooth=FALSE)

#Plot some results:
par(mfrow = c(3, 1))
par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, src$logL, type = "l", col = "blue", xlab = "Mini / Msol", ylab = "Lum / Lsol", unlog='y')
lines(tar$Mini, tar$logL, col = "red")
lines(fit$warp_src2tar(src$Mini), src$logL, col='blue', lty=2)
legend("topleft", legend = c("Source", "Target", "Source Warp"), col = c("blue", "red", "blue"), lty = c(1,1,2), bty = "n")

par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, src$logT, type = "l", col = "blue", xlab = "Mini / Msol", ylab = "Teff / K", unlog='y')
lines(tar$Mini, tar$logT, col = "red")
lines(fit$warp_src2tar(src$Mini), src$logT, col='blue', lty=2)
legend("topleft", legend = c("Source", "Target", "Source Warp"), col = c("blue", "red", "blue"), lty = c(1,1,2), bty = "n")

par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, src$logG, type = "l", col = "blue", xlab = "Mini / Msol", ylab = "logG")
lines(tar$Mini, tar$logG, col = "red")
lines(fit$warp_src2tar(src$Mini), src$logG, col='blue', lty=2)
legend("topleft", legend = c("Source", "Target", "Source Warp"), col = c("blue", "red", "blue"), lty = c(1,1,2), bty = "n")

#Now interpolate to our target age:
wt = 0.5 # 0.5 since we are comparing to an isochrone halfway in logAge space
temp_warp = progenyWarpInterp(x_src = src, x_tar = tar, ParamWarp_out = fit, wt=wt)

par(mfrow = c(4, 1))
par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, src$logL, type = "l", col = "blue", xlab = "Mini / Msol", ylab = "Lum / Lsol", unlog='y')
lines(tar$Mini, tar$logL, col = "red")
lines(cmp$Mini, cmp$logL, col = "darkgreen")
lines(temp_warp$Mini, temp_warp$logL, col="darkgreen", lty=2)
legend("topleft", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")

par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, src$logT, type = "l", col = "blue", xlab = "Mini / Msol", ylab = "Teff / K", unlog='y')
lines(tar$Mini, tar$logT, col = "red")
lines(cmp$Mini, cmp$logT, col = "darkgreen")
lines(temp_warp$Mini, temp_warp$logT, col="darkgreen", lty=2)
legend("topleft", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")

par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, src$logG, type = "l", col = "blue", xlab = "Mini / Msol", ylab = "logG")
lines(tar$Mini, tar$logG, col = "red")
lines(cmp$Mini, cmp$logG, col = "darkgreen")
lines(temp_warp$Mini, temp_warp$logG, col="darkgreen", lty=2)
legend("topleft", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")

# Mass mapping is a special case, because in general we just need to scale it:
mass_scale = max(temp_warp$Mini)/max(src$Mini)
par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$Mini, Iso_MIST[abs(logAge - 9) < 1e-4 & logZ==0,Mass], type = "l", col = "blue", xlab = "Mini / Msol", ylab = "Mass / Msol")
lines(tar$Mini, Iso_MIST[abs(logAge - 9.2) < 1e-4 & logZ==0,Mass], col = "red")
lines(cmp$Mini, Iso_MIST[abs(logAge - 9.1) < 1e-4 & logZ==0,Mass], col = "darkgreen")
lines(src$Mini*mass_scale, Iso_MIST[abs(logAge - 9) < 1e-4 & logZ==0,Mass]*mass_scale, col="darkgreen", lty=2)
legend("topleft", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")

# Look at the classic Teff versus logG space:
par(mfrow = c(3, 1))
par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$logT, src$logG, type = "l", col = "blue", xlab = "Teff / K", ylab = "logG", unlog='x')
lines(tar$logT, tar$logG, col = "red")
lines(cmp$logT, cmp$logG, col = "darkgreen")
lines(temp_warp$logT, temp_warp$logG, col="darkgreen", lty=2)
legend("bottomright", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")

# Definitely not perfect, but not too bad:
par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$logT, src$logG, type = "l", col = "blue", xlab = "Teff / K", ylab = "logG", xlim=c(3.7,3.9), ylim=c(3,4))
lines(tar$logT, tar$logG, col = "red")
lines(cmp$logT, cmp$logG, col = "darkgreen")
lines(temp_warp$logT, temp_warp$logG, col="darkgreen", lty=2)
legend("bottomright", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")

# TP AGB region is complex to interpolate too:
par(mar=c(3.1,3.1,1.1,1.1))
magplot(src$logT, src$logG, type = "l", col = "blue", xlab = "Teff / K", ylab = "logG", xlim=c(3.4,3.6), ylim=c(-1,0))
lines(tar$logT, tar$logG, col = "red")
lines(cmp$logT, cmp$logG, col = "darkgreen")
lines(temp_warp$logT, temp_warp$logG, col="darkgreen", lty=2)
legend("bottomright", legend = c("Source", "Target", "Compare", "Param Warp"), col = c("blue", "red", 'darkgreen', "darkgreen"), lty = c(1,1,1,2), bty = "n")
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
