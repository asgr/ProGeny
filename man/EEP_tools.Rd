\name{EEP_tools}
\alias{identify_primary_eeps}
\alias{metric_distance}
\alias{build_eep_track}
\alias{plot_eep}
\title{Identify and Interpolate Equivalent Evolutionary Points (EEPs)}
\description{
Utilities to (i) identify primary Equivalent Evolutionary Points (EEPs) on a stellar evolution \code{track},
(ii) compute a cumulative metric distance along regions of a track, (iii) build an interpolated EEP
representation of a track, and (iv) plot the original and EEP-based Hertzsprung--Russell diagram (HRD) and other diagnostic plots.
These functions are designed to work with MESA-like tracks whose columns are mapped via \code{cols}.
}

\section{Shared Inputs}{
All functions expect a data frame \code{track} with at least the following columns (or mapped equivalents):
\itemize{
  \item \code{star_age}: stellar age in years (numeric)
  \item \code{log_L}: log$_{10}$ of luminosity in solar units (numeric)
  \item \code{log_Teff}: log$_{10}$ of effective temperature in Kelvin (numeric)
  \item \code{log_center_T}: log$_{10}$ of central temperature (numeric)
  \item \code{center_h1}: central H mass fraction (numeric)
  \item \code{center_he4}: central He mass fraction (numeric)
}
Optional columns, if present (or if their names are provided in \code{cols}), enable additional phase detections:
\itemize{
  \item \code{log_LH}: log$_{10}$ H-burning luminosity (for precise ZAMS detection)
  \item \code{center_gamma}: central Coulomb coupling parameter Gamma_c (for WDCS detection)
  \item \code{star_mass}: total stellar mass (for envelope-based PostAGB detection)
  \item \code{m_shell_H}, \code{m_shell_He}: shell mass coordinates (for TP-AGB detection)
  \item \code{center_c12}: central carbon mass fraction (for core C-burning detection)
  \item \code{he_core_mass}: He-core mass (for envelope mass inference)
}
The \code{cols} argument maps the semantic names used by these functions to the actual column names in your \code{track}.
}

\usage{
identify_primary_eeps(
  track,
  cols = list(
    star_age = "star_age", logL = "log_L", logTeff = "log_Teff",
    logTc = "log_center_T", Xc = "center_h1", Yc = "center_he4",
    logLH = "log_LH", Gamma_c = "center_gamma", M_star = "star_mass",
    M_Hshell = "m_shell_H", M_Heshell = "m_shell_He", Xc_C = "center_c12",
    M_core = "he_core_mass"
  ),
  params = list(
    prems_logTc = 5.0, zams_lburn_frac = 0.999, zams_dXc_max = 0.0015,
    iams_Xc = 0.3, tams_Xc = 1e-12, rgbtip_Yc = 0.01, zacheb_Yc = 0.03,
    tacheb_Yc = 1e-4, tpagb_Yc = 1e-6, tpagb_delta_M = 0.1, cburn_Xc_X = 1e-4,
    postagb_env_frac = 0.2, wdcs_Gamma_max = 100
  )
)

metric_distance(
  track, ind,
  vars = c("logL", "logTeff"), weights = NULL,
  cols = list(
    star_age = "star_age", logL = "log_L", logTeff = "log_Teff",
    logTc = "log_center_T", Xc = "center_h1", Yc = "center_he4",
    logLH = "log_LH", Gamma_c = "center_gamma", M_star = "star_mass",
    M_Hshell = "m_shell_H", M_Heshell = "m_shell_He", Xc_C = "center_c12",
    M_core = "he_core_mass"
  )
)

build_eep_track(
  track, primary_idx, n_secondary_between = 50,
  metric_vars = c("logL", "logTeff"), metric_weights = NULL,
  cols = list(
    star_age = "star_age", logL = "log_L", logTeff = "log_Teff",
    logTc = "log_center_T", Xc = "center_h1", Yc = "center_he4",
    logLH = "log_LH", Gamma_c = "center_gamma", M_star = "star_mass",
    M_Hshell = "m_shell_H", M_Heshell = "m_shell_He", Xc_C = "center_c12",
    M_core = "he_core_mass"
  )
)

plot_eep(
  track, eep_track,
  cols = list(x = "log_Teff", y = "log_L"),
  xlab = cols$x, ylab = cols$y,
  show_secondary = TRUE, ...
)
}

\arguments{
  \item{track}{A \code{data.frame} with stellar evolution track data. Columns are mapped via \code{cols}.}
  \item{cols}{A named \code{list} mapping semantic names used by these functions (e.g., \code{"logL"}, \code{"logTeff"})
              to column names present in \code{track}. Entries may be set to \code{NULL} to indicate absence.}
  \item{params}{For \code{identify_primary_eeps}, a named \code{list} of thresholds controlling phase picks:
  \describe{
    \item{\code{prems_logTc}}{log$_{10}$ central temperature threshold for \strong{PreMS} start.}
    \item{\code{zams_lburn_frac}}{Required fraction of luminosity from H-burning at ZAMS (if \code{log_LH} available).}
    \item{\code{zams_dXc_max}}{Maximum allowed central-H depletion prior to ZAMS (fallback).}
    \item{\code{iams_Xc}}{Central H abundance for \strong{IAMS}.}
    \item{\code{tams_Xc}}{Central H abundance for \strong{TAMS}.}
    \item{\code{rgbtip_Yc}}{Allowed drop in central He prior to choosing \strong{RGBTip}.}
    \item{\code{zacheb_Yc}}{Allowed drop in central He prior to \strong{ZACHeB}.}
    \item{\code{tacheb_Yc}}{Central He threshold for \strong{TACHeB}.}
    \item{\code{tpagb_Yc}}{Central He threshold used in TP-AGB detection.}
    \item{\code{tpagb_delta_M}}{Shell mass separation criterion (H vs He) for TP-AGB.}
    \item{\code{cburn_Xc_X}}{Central carbon threshold for end of core C burning.}
    \item{\code{postagb_env_frac}}{Envelope mass fraction threshold for PostAGB.}
    \item{\code{wdcs_Gamma_max}}{Maximum \eqn{\Gamma_c} for entry onto WD cooling sequence.}
  }}
  \item{ind}{For \code{metric_distance}, integer indices defining the segment along \code{track}.}
  \item{vars}{For \code{metric_distance}, character vector of variables (mapped via \code{cols}) used to build the metric.}
  \item{weights}{Optional numeric weights (same length as \code{vars}); defaults to equal weights.}
  \item{primary_idx}{For \code{build_eep_track}, a named integer vector of primary EEP indices, as returned by \code{identify_primary_eeps}.}
  \item{n_secondary_between}{Number of secondary EEPs to insert \emph{between} successive primaries; can be a single integer or a vector of length equal to the number of segments.}
  \item{metric_vars}{Variables used to compute the distance metric for interpolation in \code{build_eep_track}.}
  \item{metric_weights}{Optional weights for \code{metric_vars}.}
  \item{eep_track}{For \code{plot_eep}, the EEP-enhanced track returned by \code{build_eep_track}.}
  \item{xlab, ylab}{Axis labels for the HRD.}
  \item{show_secondary}{Logical; if \code{TRUE}, plot secondary EEPs as faint points.}
  \item{...}{Additional graphical parameters passed to \code{magicaxis::magplot} in \code{plot_eep}.}
}

\details{
\subsection{Primary EEP identification (\code{identify_primary_eeps})}{
EEPs are picked in the following order (indices recorded in a named integer vector):
\enumerate{
  \item \strong{PreMS}: first index where \code{log_center_T} exceeds \code{prems_logTc}; if the track already starts hotter, index 1.
  \item \strong{ZAMS}: (preferred) first step after PreMS where \eqn{L_\mathrm{H}/L_\mathrm{tot} >}\code{zams_lburn_frac} and \code{center_h1} has not dropped by more than \code{zams_dXc_max}. Fallback: first step where \code{center_h1} begins to decrease.
  \item \strong{IAMS}: first index with \code{center_h1} \eqn{\le} \code{iams_Xc}.
  \item \strong{TAMS}: first index with \code{center_h1} \eqn{\le} \code{tams_Xc}; if never reached, the minimum of \code{center_h1}.
  \item \strong{RGBTip}: after TAMS and while central He has not dropped by more than \code{rgbtip_Yc}, pick the earlier of (i) maximum \code{log_L} or (ii) minimum \code{log_Teff}.
  \item \strong{ZACHeB}: after RGBTip and while central He has not dropped by more than \code{zacheb_Yc}, pick the minimum \code{log_center_T}.
  \item \strong{TACHeB}: first index after ZACHeB with \code{center_he4} \eqn{\le} \code{tacheb_Yc}; otherwise, the last available index.
  \item \strong{TPAGB}: (if \code{m_shell_H} and \code{m_shell_He} are available) first step after TACHeB with \code{center_he4} \eqn{<} \code{tpagb_Yc} and \code{m_shell_H - m_shell_He} \eqn{<} \code{tpagb_delta_M}.
  \item \strong{CBurn}: (massive stars; if \code{center_c12} is available) first step after post-RGB/He phases where \code{center_c12} \eqn{\le} \code{cburn_Xc_X}.
  \item \strong{PostAGB}: (if \code{star_mass} and \code{he_core_mass} are available) first step where envelope fraction \eqn{(M_\mathrm{env}/M_\star) \le} \code{postagb_env_frac}, with \eqn{M_\mathrm{env} = M_\star - M_\mathrm{core}}.
  \item \strong{WDCS}: (if \code{center_gamma} is available) first step after PostAGB where \eqn{\Gamma_c \le} \code{wdcs_Gamma_max}.
}
Any of these may return \code{NA_integer_} if not identifiable with the available inputs.
}

\subsection{Metric distance (\code{metric_distance})}{
Computes a cumulative arc-length-like metric along indices \code{ind} in the space spanned by \code{vars} (defaults to \code{logL} and \code{logTeff}). The metric is
\deqn{D_k = \sum_{i=2}^k \sqrt{\sum_j w_j (x_{i,j} - x_{i-1,j})^2}},
with \eqn{w_j} provided by \code{weights} (defaults to 1). Returned vector length equals \code{length(ind)}; the first element is 0.
}

\subsection{Building an EEP track (\code{build_eep_track})}{
For each segment between successive primary EEPs (and the last primary to the end of the track), \code{build_eep_track}
interpolates \code{n_secondary_between} interior points using a cumulative metric based on \code{metric_vars}/\code{metric_weights}.
Numeric columns are interpolated (cubic spline via \code{\link[stats]{spline}}); non-numeric columns are carried from the nearest lower index.
The returned data frame contains the original columns plus:
\itemize{
  \item \code{EEP}: sequential integer index along the EEP-enhanced track
  \item \code{EEP_phase}: character label of the segment's starting primary EEP (e.g., \code{"ZAMS"})
  \item \code{EEP_is_primary}: logical flag; \code{TRUE} on the first row of each segment (and the last row of the last segment is reset to \code{FALSE} to avoid duplication)
  \item \code{label}: integer coarse phase code (see \emph{Phase code mapping})
}
}

\subsection{Plotting (\code{plot_eep})}{
Draws an HRD with the original track (grey), the EEP-enhanced track (blue), primary EEPs (open circles), and optional secondary EEPs (faint dots). Labels for primary EEPs are added next to their points. The function uses \code{magicaxis::magplot} for base plotting.
}
}

\section{Phase code mapping}{
The integer column \code{label} in the EEP track maps EEP phases to coarse MIST-like phases as follows:
\describe{  
  \item{\code{-1}}{PreMS}  
  \item{\code{0}}{MS (ZAMS/IAMS/TAMS)}  
  \item{\code{2}}{RGB (entire RGB; \code{RGB} or \code{RGBTip})}  
  \item{\code{3}}{CHeB (from ZACHeB through TACHeB)}  
  \item{\code{4}}{eAGB (if flagged as EAGB/EBGB)}  
  \item{\code{5}}{TP-AGB}  
  \item{\code{6}}{PostAGB}  
  \item{\code{7}}{CSPN (central star of PN; if present)}  
  \item{\code{8}}{WDCS}  
  \item{\code{9}}{WR (if present)}
}
}

\value{
\describe{
  \item{\code{identify_primary_eeps}}{A named integer vector of indices (one per primary EEP), with names \code{PreMS}, \code{ZAMS}, \code{IAMS}, \code{TAMS}, \code{RGBTip}, \code{ZACHeB}, \code{TACHeB}, \code{TPAGB}, \code{CBurn}, \code{PostAGB}, \code{WDCS}. Entries may be \code{NA_integer_}.}
  \item{\code{metric_distance}}{A numeric vector of length \code{length(ind)} giving the cumulative metric distance (first element 0).}
  \item{\code{build_eep_track}}{A \code{data.frame} of interpolated points with columns of \code{track} plus \code{EEP}, \code{EEP_phase}, \code{EEP_is_primary}, and \code{label}.}
  \item{\code{plot_eep}}{Invisibly returns \code{NULL}; called for its side effect of drawing on the active graphics device.}
}
}

\note{
\itemize{
  \item The default thresholds in \code{params} follow practical adaptations of EEP ideas (see Dotter 2016 for EEP methodology) but are not guaranteed to match any particular pipeline; tune as needed.
  \item \code{plot_eep} requires \pkg{magicaxis} for \code{magplot} (or provide a compatible \code{magplot} in your namespace).
  \item Optional columns listed in \code{cols} may be omitted (set to \code{NULL} or leave unmapped). Related EEPs are then skipped.
}
}

\references{
Aaron Dotter, 2016, ApJS, 222, 8D \cr
MESA instrument papers: Paxton et al. (2011, 2013, 2015, 2018, 2019), \emph{ApJS}. \cr 
MIST project documentation and isochrone column descriptions.
}

\seealso{\code{\link[stats]{spline}}, \pkg{magicaxis} (\code{magplot}).}

\examples{
## Minimal synthetic example (toy track; random walk in HRD)
track = read.table(system.file("extdata", '02000M.track',
  package="ProGeny"))

## Pick primaries (using defaults and available columns)
prim = identify_primary_eeps(track)

## Build EEP track with 20 secondaries between primaries
eep = build_eep_track(track, prim, n_secondary_between = 20)

## Plot HR diagram (notice glitches near ZAM)
plot_eep(track, eep)

## Build EEP track with 100 secondaries between primaries
eep = build_eep_track(track, prim, n_secondary_between = 100)

## Plot HR diagram (looks much cleaner)
plot_eep(track, eep)

## Classic log_Teff versus log_g plot:
plot_eep(track, eep, cols = c('log_Teff', 'log_g'))
}
