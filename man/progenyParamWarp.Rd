\name{progenyParamWarp}
\alias{progenyParamWarp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Warp [x,y] Source Data to Target
}
\description{
A common issue when comparing isochrones is various parameters might behave similarly but not identically. This function allows as to warp the x-axis of source data to create optimal overlap in source and target y-axis.
}
\usage{
progenyParamWarp(x_src, y_src, x_tar, y_tar, check_order = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_src}{
Numeric vector; x source data we want to warp.
}
  \item{y_src}{
Numeric vector; y source data we want to maximially overlap with \option{t_tar} after warping. Must be same length as \option{x_src}.
}
  \item{x_tar}{
Numeric vector; x target data.
}
  \item{y_tar}{
Numeric vector; y target data we want to maximially overlap with \option{y_src} after warping. Must be same length as \option{x_tar}.
}
  \item{check_order}{
Logical; should \option{x_src} and \option{x_tar} be checked for sorting
}
  \item{\dots}{
Other arguments to pass to \code{dtw}.
}
}
\details{
This function uses the \code{dtw} package and function to uncover the optimal warping (or re-mapping) between source and target data. The end result are two warping functions that will convert source x data to target. The warping is only done in the x-dimension, and the y-dimension is kept intact.

The main use case for this function is to allow isochrone properties (like Mini and Lum) to be distorted between isochrones so properties can be compared/used.

Note source inputs (\option{x_src} and \option{y_src}) do no need to be the same length as the target inputs (\option{x_tar} and \option{y_tar}).
}
\value{
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
x_tar = seq(0, 10, length.out = 100)
x_tar_warped = 10 * ( (x_tar / 10)^1.3 )  # monotone distortion
ysp_src = smooth.spline(x_src, y_src, spar = 0.6)
y_tar   = predict(ysp_src, x_tar_warped)$y + rnorm(length(x_tar), sd = 0.05)

fit <- warp_by_curvature_dtw(x_src, y_src, x_tar, y_tar)

# Original curves
magplot(x_src, y_src, type = "l", col = "blue", lwd = 2,
     main = "Source and target (original)",
     xlab = "x", ylab = "y")
lines(x_tar, y_tar, col = "red", lwd = 2)
lines(fit$warp_src2tar(x_src), y_src, col='blue', lty=2, lwd=2)
legend("bottomleft", legend = c("Source", "Target", "Source Warp"),
       col = c("blue", "red", "blue"), lty = c(1,1,2), bty = "n")

magcurve(fit$warp_src2tar(x), xlim=c(0,10), xlab='x_src', ylab='x_tar')
abline(0,1,lty=2)
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
