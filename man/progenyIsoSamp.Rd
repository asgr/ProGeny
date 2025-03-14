\name{progenyIsoSamp}
\alias{progenyIsoSamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
progenyIsoSamp(Iso, IMFfunc, ..., Nsamp = 10000, Lum_weight_func = NULL, logAge_weight_func = NULL, logZ_weight_func = NULL, label_use = seq(-1, 9), seed = 666)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Iso}{
%%     ~~Describe \code{Iso} here~~
}
  \item{IMFfunc}{
%%     ~~Describe \code{IMFfunc} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{Nsamp}{
%%     ~~Describe \code{Nsamp} here~~
}
  \item{Lum_weight_func}{
%%     ~~Describe \code{Lum_weight_func} here~~
}
  \item{logAge_weight_func}{
%%     ~~Describe \code{logAge_weight_func} here~~
}
  \item{logZ_weight_func}{
%%     ~~Describe \code{logZ_weight_func} here~~
}
  \item{label_use}{
%%     ~~Describe \code{label_use} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (Iso, IMFfunc, ..., Nsamp = 10000, Lum_weight_func = NULL, 
    logAge_weight_func = NULL, logZ_weight_func = NULL, label_use = seq(-1, 
        9), seed = 666) 
{
    set.seed(seed)
    setDT(Iso)
    Iso_temp = copy(Iso)
    setkeyv(Iso_temp, c("logZ", "logAge", "Mini"))
    logZ = logAge = Mini = Mass = NULL
    Iso_temp[, `:=`(IMFint, 0)]
    message("Integrating the IMF")
    Iso_temp[Lum > 1e-06, `:=`(IMFint, progenyUpdateIMF(Iso_temp[Lum > 
        1e-06, ], IMFfunc, ...)$IMFint)]
    message("Applying weights")
    if (!is.null(Lum_weight_func)) {
        Iso_temp[, `:=`(IMFint, IMFint * Lum_weight_func(Lum))]
    }
    if (!is.null(logAge_weight_func)) {
        Iso_temp[, `:=`(IMFint, IMFint * logAge_weight_func(logAge))]
    }
    if (!is.null(logZ_weight_func)) {
        Iso_temp[, `:=`(IMFint, IMFint * logZ_weight_func(logZ))]
    }
    message("Creating inversion function")
    Iso_temp[IMFint == 0, `:=`(IMFint, 1e-300)]
    Iso_temp = Iso_temp[label \%in\% label_use, ]
    temp_cumsum = cumsum(Iso_temp$IMFint)
    temp_cumsum = temp_cumsum/max(temp_cumsum)
    sel_dup = which(duplicated(temp_cumsum)) - 1L
    if (length(sel_dup) > 0) {
        temp_cumsum = temp_cumsum[-sel_dup]
        Iso_temp = Iso_temp[-sel_dup, ]
    }
    temp_map = approxfun(temp_cumsum, 1:length(temp_cumsum), 
        yleft = 1, yright = length(temp_cumsum))
    message("Sampling")
    sel = ceiling(temp_map(runif(Nsamp)))
    return(Iso_temp[sel, ])
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
