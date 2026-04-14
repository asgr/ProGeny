\name{EEP_tools}
\alias{EEP_tools}
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

The algorithms follow the Fortran \code{Iso} package (\url{https://github.com/aarondotter/iso})
as closely as possible, including the 5-dimensional weighted distance metric
and the primary EEP identification rules from \code{eep.f90} and \code{iso_eep_support.f90}.
}

\section{Shared Inputs}{
All functions expect a data frame \code{track} with at least the following columns (or mapped equivalents):
\itemize{
  \item \code{star_age}: stellar age in years (numeric)
  \item \code{log_L}: \eqn{log_{10}} of luminosity in solar units (numeric)
  \item \code{log_Teff}: \eqn{log_{10}} of effective temperature in Kelvin (numeric)
  \item \code{log_center_T}: \eqn{log_{10}} of central temperature (numeric)
  \item \code{center_h1}: central H mass fraction (numeric)
  \item \code{center_he4}: central He mass fraction (numeric)
}
Optional columns, if present (or if their names are provided in \code{cols}), enable additional phase detections:
\itemize{
  \item \code{log_g}: surface gravity (for ZAMS via max-logg method)
  \item \code{log_LH}: \eqn{log_{10}} H-burning luminosity
  \item \code{log_LHe}: \eqn{log_{10}} He-burning luminosity (for ZACHeB He-flash peak anchor)
  \item \code{log_center_Rho}: \eqn{log_{10}} central density (for distance metric)
  \item \code{center_gamma}: central Coulomb coupling parameter \eqn{\Gamma_c} (for WDCS detection)
  \item \code{star_mass}: total stellar mass (for PostAGB detection)
  \item \code{he_core_mass}: He-core mass (for TPAGB intershell calculation)
  \item \code{c_core_mass}: C-core mass (for TPAGB, PostAGB core CO fraction)
  \item \code{o_core_mass}: O-core mass (for TPAGB, PostAGB core CO fraction)
  \item \code{center_c12}: central carbon mass fraction (for core C-burning detection)
}
The \code{cols} argument maps the semantic names used by these functions to the actual column names in your \code{track}.
The default values in \code{cols} correspond to standard MESA history column names used by MIST.
}

\usage{
identify_primary_eeps(
  track,
  cols = list(
    star_age = "star_age", logL = "log_L", logTeff = "log_Teff",
    logTc = "log_center_T", logRhoc = "log_center_Rho",
    Xc = "center_h1", Yc = "center_he4",
    logg = "log_g", logLH = "log_LH", logLHe = "log_LHe",
    Gamma_c = "center_gamma", M_star = "star_mass",
    M_He_core = "he_core_mass", M_C_core = "c_core_mass",
    M_O_core = "o_core_mass", Xc_C = "center_c12"
  ),
  params = list(
    prems_logTc = 5.0, prems_logTc_adj = 0.01,
    zams_dXc_max = 0.001, iams_Xc = 0.35, tams_Xc = 1e-12,
    tams_max_age = 1.5e10,
    rgbtip_Yc = 0.01, zacheb_Yc = 0.03, tacheb_Yc = 1e-4,
    tpagb_Yc = 1e-6, tpagb_He_shell = 0.1,
    cburn_XY_lim = 1e-8, cburn_C_lim = 1e-4,
    postagb_core_frac = 0.8, wdcs_Gamma_min = 100
  )
)

metric_distance(
  track, ind,
  Teff_scale = 2.0, logL_scale = 0.125, Rhoc_scale = 1.0,
  Tc_scale = 1.0, age_scale = 0.05,
  weight_center_rho_T_by_Xc = TRUE,
  cols = list(star_age  = "star_age", logL = "log_L",
    logTeff = "log_Teff", logTc = "log_center_T",
    logRhoc = "log_center_Rho", Xc = "center_h1"
  )
)

build_eep_track(
  track, primary_idx, n_secondary_between = NULL,
  Teff_scale = 2.0, logL_scale = 0.125, Rhoc_scale = 1.0,
  Tc_scale = 1.0, age_scale = 0.05,
  weight_center_rho_T_by_Xc = TRUE,
  cols = list(star_age  = "star_age", logL = "log_L",
    logTeff = "log_Teff", logTc = "log_center_T",
    logRhoc = "log_center_Rho", Xc = "center_h1"
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
              to column names present in \code{track}. Entries may be set to \code{NULL} to indicate absence.
              Default names follow MESA history column conventions used by MIST.}
  \item{params}{For \code{identify_primary_eeps}, a named \code{list} of thresholds controlling phase picks.
  All defaults match the Fortran \code{Iso} \code{iso} package (\code{eep.f90} and \code{iso_eep_support.f90}):
  \describe{
    \item{\code{prems_logTc}}{\eqn{log_{10}} central temperature threshold for \strong{PreMS} start (default 5.0).}
    \item{\code{prems_logTc_adj}}{Adjustment added to \code{logTc[1]} when the track already starts hotter than
      \code{prems_logTc} (default 0.01, matching Fortran \code{Iso} \code{PreMS_Tc} high-mass fallback).}
    \item{\code{zams_dXc_max}}{Maximum allowed central-H depletion used to find the Xc-drop point
      (default 0.001, matching Fortran \code{Iso} \code{Xmin = Xmax - 1.0d-3}).
      ZAMS is then placed at max(\code{log_g}) from index 1 to the Xc-drop point.}
    \item{\code{iams_Xc}}{Central H abundance for \strong{IAMS} (default 0.35, matching Fortran \code{Iso} \code{TAMS(t,3.5d-1,...)}).}
    \item{\code{tams_Xc}}{Central H abundance for \strong{TAMS} (default 1e-12, matching Fortran \code{Iso} \code{TAMS(t,1d-12,...)}).}
    \item{\code{tams_max_age}}{Maximum age (years) for low-mass TAMS fallback (default 1.5e10 = 15 Gyr).}
    \item{\code{rgbtip_Yc}}{Allowed drop in central He prior to choosing \strong{RGBTip} (default 0.01).}
    \item{\code{zacheb_Yc}}{Allowed drop in central He for \strong{ZACHeB} search (default 0.03).
      ZACHeB uses a two-phase search: first find the He-luminosity peak (\code{log_LHe}),
      then find min \code{log_center_T} after the peak.}
    \item{\code{tacheb_Yc}}{Central He threshold for \strong{TACHeB} (default 1e-4).
      Returns \code{NA} if the threshold is never reached (matching Fortran \code{Iso} returning 0).}
    \item{\code{tpagb_Yc}}{Central He threshold for TP-AGB detection (default 1e-6).}
    \item{\code{tpagb_He_shell}}{He intershell mass threshold: \code{he_core_mass - c_core_mass < 0.1} (default 0.1, matching Fortran \code{Iso} \code{HeShellMin}).}
    \item{\code{cburn_XY_lim}}{H and He exhaustion threshold for C burning (default 1e-8, matching Fortran \code{Iso} \code{limit_XY}).}
    \item{\code{cburn_C_lim}}{Central C depletion threshold for C burning (default 1e-4, matching Fortran \code{Iso} \code{center_carbon_limit}).}
    \item{\code{postagb_core_frac}}{CO-core mass fraction threshold for PostAGB: \code{c_core_mass/star_mass > 0.8} (default 0.8, matching Fortran \code{Iso} \code{core_mass_frac_limit}). A Tc guard requires \code{logTc} at the PostAGB search starting point (after TPAGB or CBurn) to exceed \code{logTc} at the track end.}
    \item{\code{wdcs_Gamma_min}}{Minimum \eqn{\Gamma_c} for entry onto WD cooling sequence (default 100).
      Looks for \code{center_gamma >= wdcs_Gamma_min}, matching Fortran \code{Iso} \code{>= center_gamma_limit}.}
  }}
  \item{ind}{For \code{metric_distance}, integer indices defining the segment along \code{track}.}
  \item{Teff_scale}{Weight for \code{log_Teff} in the distance metric (default 2.0).}
  \item{logL_scale}{Weight for \code{log_L} in the distance metric (default 0.125).}
  \item{Rhoc_scale}{Weight for \code{log_center_Rho} in the distance metric (default 1.0). Xc-weighted when \code{weight_center_rho_T_by_Xc} is \code{TRUE}.}
  \item{Tc_scale}{Weight for \code{log_center_T} in the distance metric (default 1.0). Xc-weighted when \code{weight_center_rho_T_by_Xc} is \code{TRUE}.}
  \item{age_scale}{Weight for \eqn{log_{10}}(age) in the distance metric (default 0.05).}
  \item{weight_center_rho_T_by_Xc}{Logical (default \code{TRUE}). If \code{TRUE}, contributions from
    \code{log_center_Rho} and \code{log_center_T} are weighted by \code{max(0, Xc/max(Xc))} at each step,
    matching the Fortran \code{Iso} \code{weight_center_rho_T_by_Xc} option. This effectively down-weights interior
    structure changes after core-H exhaustion.}
  \item{primary_idx}{For \code{build_eep_track}, a named integer vector of primary EEP indices, as returned by \code{identify_primary_eeps}.}
  \item{n_secondary_between}{\code{NULL} for Fortran \code{Iso}-style defaults (variable per interval, e.g. 200 between PreMS-ZAMS, 150 between ZAMS-IAMS, etc.),
    a single integer (applied uniformly), or a vector of length equal to the number of intervals between consecutive valid primaries.}
  \item{eep_track}{For \code{plot_eep}, the EEP-enhanced track returned by \code{build_eep_track}.}
  \item{xlab, ylab}{Axis labels for the HRD.}
  \item{show_secondary}{Logical; if \code{TRUE}, plot secondary EEPs as faint points.}
  \item{...}{Additional graphical parameters passed to \code{magicaxis::magplot} in \code{plot_eep}.}
}

\details{
\subsection{Primary EEP identification (\code{identify_primary_eeps})}{
EEPs are picked in the following order (indices recorded in a named integer vector).
All algorithms match the Fortran \code{Iso} \code{iso} package (\code{eep.f90}):
\enumerate{
  \item \strong{PreMS}: first index where \code{log_center_T} exceeds \code{prems_logTc}.
    If the track already starts hotter, uses \code{logTc[1] + prems_logTc_adj}.
  \item \strong{ZAMS}: (Fortran \code{Iso} ZAMS3) find where central H has dropped by \code{zams_dXc_max} from the initial value;
    then pick the index of maximum \code{log_g} from the start to that point.
    Fallback if \code{log_g} is not available: use the Xc-drop point directly (ZAMS1).
  \item \strong{IAMS}: first index with \code{center_h1} \eqn{<} \code{iams_Xc} (default 0.35).
  \item \strong{TAMS}: first index with \code{center_h1} \eqn{<} \code{tams_Xc}; for very-low-mass stars (\eqn{\le 0.5\,M_\odot}), accepts the last track point if its age exceeds \code{tams_max_age}.
  \item \strong{RGBTip}: after TAMS and while central He has not dropped by more than \code{rgbtip_Yc} from its value at TAMS+1, pick the earlier of (i) maximum \code{log_L} or (ii) minimum \code{log_Teff}.
  \item \strong{ZACHeB}: (Fortran \code{Iso} ZAHB) two-phase search after RGBTip:
    (a) find the He-luminosity peak (\code{log_LHe}, if available) while Yc constraint holds;
    (b) find minimum \code{log_center_T} after the peak. If \code{log_LHe} is not available, phase (a) is skipped.
  \item \strong{TACHeB}: (Fortran \code{Iso} TAHB) first index after ZACHeB with \code{center_he4} \eqn{<} \code{tacheb_Yc}.
    Returns \code{NA} if the threshold is never reached (matching Fortran \code{Iso} returning 0).
  \item \strong{TPAGB}: (low-mass branch) first step after TACHeB with \code{center_he4} \eqn{<} \code{tpagb_Yc} and intershell mass \code{he_core_mass - c_core_mass} \eqn{<} \code{tpagb_He_shell}.
  \item \strong{CBurn}: (high-mass branch; attempted only if TPAGB was not found) first step after TACHeB where \code{center_h1} \eqn{<} \code{cburn_XY_lim}, \code{center_he4} \eqn{<} \code{cburn_XY_lim}, and \code{center_c12} \eqn{<} \code{cburn_C_lim}.
  \item \strong{PostAGB}: first step where \code{co_core_mass / star_mass >} \code{postagb_core_frac}, with a Tc guard requiring that \code{logTc} at the starting point exceeds \code{logTc} at the track end.
  \item \strong{WDCS}: first step after PostAGB where \code{center_gamma >=} \code{wdcs_Gamma_min}.
}
If any EEP is not identifiable (threshold not met, required column missing, etc.), it and all subsequent EEPs return \code{NA_integer_}.
}

\subsection{Metric distance (\code{metric_distance})}{
Computes a cumulative arc-length-like metric along indices \code{ind}, matching the Fortran \code{Iso} \code{distance_along_track} subroutine. The metric is:
\deqn{D_k = \sum_{i=2}^k \sqrt{w_\mathrm{Teff}\,\Delta\log T_\mathrm{eff}^2 + w_L\,\Delta\log L^2 + w(X_c)\,w_\rho\,\Delta\log\rho_c^2 + w(X_c)\,w_T\,\Delta\log T_c^2 + w_\mathrm{age}\,\Delta\log(\mathrm{age})^2}},
where \eqn{w(X_c) = \max(0, X_c / \max(X_c))} when \code{weight_center_rho_T_by_Xc} is \code{TRUE} (default).
Columns not present in the track are silently skipped (their contribution is 0).
}

\subsection{Building an EEP track (\code{build_eep_track})}{
For each interval between successive primary EEPs, \code{build_eep_track}
interpolates \code{n_secondary_between} interior points using the cumulative metric.
Matching the Fortran \code{Iso} behaviour, the track is only constructed between the first and last
valid primary EEPs --- it is not extended past the last primary.

Numeric columns are interpolated using monotone Hermite cubic splines
(\code{method = "monoH.FC"} in \code{\link[stats]{spline}}), which approximates the
Fortran \code{Iso} \code{interp_4pt_pm} (Steffen's piecewise monotonic method). Both are
monotonicity-preserving cubic Hermite interpolants that produce similar results,
though minor numerical differences may exist due to different derivative estimation
formulas (Fritsch-Carlson vs Steffen).

The returned data frame contains the original columns plus:
\itemize{
  \item \code{EEP}: sequential integer index along the EEP-enhanced track
  \item \code{EEP_phase}: character label of the segment's starting primary EEP (e.g., \code{"ZAMS"})
  \item \code{EEP_is_primary}: logical flag; \code{TRUE} on the first row of each segment and on the last row of the last segment
  \item \code{label}: integer coarse phase code (see \emph{Phase code mapping})
}
When \code{n_secondary_between} is \code{NULL}, the function uses Fortran \code{Iso}-style default
counts per interval: PreMS\eqn{\to}ZAMS: 200, ZAMS\eqn{\to}IAMS: 150,
IAMS\eqn{\to}TAMS: 100, TAMS\eqn{\to}RGBTip: 150, RGBTip\eqn{\to}ZACHeB: 25,
ZACHeB\eqn{\to}TACHeB: 75, TACHeB\eqn{\to}TPAGB: 100, TPAGB\eqn{\to}PostAGB: 600,
PostAGB\eqn{\to}WDCS: 300.
}

\subsection{Plotting (\code{plot_eep})}{
Draws an HRD with the original track (grey), the EEP-enhanced track (blue), primary EEPs (open circles), and optional secondary EEPs (faint dots). Labels for primary EEPs are added next to their points. The function uses \code{magicaxis::magplot} for base plotting.
}
}

\section{Phase code mapping}{
The integer column \code{label} in the EEP track maps EEP phases to coarse MIST-like phases as follows:
\itemize{
  \item{\code{-1}: PreMS}
  \item{\code{0}: MS (ZAMS/IAMS/TAMS)}
  \item{\code{2}: RGB (entire RGB; \code{RGB} or \code{RGBTip})}
  \item{\code{3}: CHeB (from ZACHeB through TACHeB)}
  \item{\code{4}: eAGB (if flagged as EAGB/EBGB)}
  \item{\code{5}: TPAGB}
  \item{\code{6}: PostAGB}
  \item{\code{7}: CSPN (central star of PN; if present)}
  \item{\code{8}: WDCS}
  \item{\code{9}: WR (if present)}
}

Note different mass ranges may or may not experience certain phases, e.g.:

Low-mass (Mini [Msol] < 2): \cr
PreMS → ZAMS → IAMS → TAMS → RGBTip → ZACHeB → TACHeB → TPAGB → PostAGB → WDCS \cr

Intermediate-mass (2 < Mini [Msol] < 8 ): \cr
PreMS → ZAMS → IAMS → TAMS → ZACHeB → TACHeB → TPAGB → WDCS \cr

Super-AGB (8 < Mini [Msol] < 10 Msol): \cr
PreMS → ZAMS → IAMS → TAMS → ZACHeB → TACHeB → CBurn → WDCS \cr

Massive-mass (Mini [Msol] > 10): \cr
PreMS → ZAMS → IAMS → TAMS → ZACHeB → TACHeB → CBurn → Core collapse\cr

Note in particular the RGBTip only occurs for lower mass stars, and this phase is entirely missed in larger stars.
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
  \item All default thresholds and algorithm choices match the Fortran \code{Iso} package (\url{https://github.com/aarondotter/iso}), specifically \code{eep.f90} and \code{iso_eep_support.f90}. See Dotter (2016) for the underlying methodology.
  \item The spline method \code{"monoH.FC"} (Fritsch-Carlson monotone Hermite) closely approximates the Fortran \code{Iso} \code{interp_4pt_pm} (Steffen's piecewise monotonic method). Both preserve monotonicity.
  \item \code{plot_eep} requires \pkg{magicaxis} for \code{magplot} (or provide a compatible \code{magplot} in your namespace).
  \item Optional columns listed in \code{cols} may be omitted (set to \code{NULL} or leave unmapped). Related EEPs or distance-metric components are then skipped.
  \item The MESA column names in the default \code{cols} follow MIST conventions (\url{https://github.com/aarondotter/iso/blob/master/README.md}).
}
}

\references{
Aaron Dotter, 2016, ApJS, 222, 8D \cr
Fortran \code{Iso} package: \url{https://github.com/aarondotter/iso} \cr
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

## Plot HR diagram
plot_eep(track, eep)

## Build EEP track with Fortran \code{Iso}-style default secondary counts
eep = build_eep_track(track, prim)

## Plot HR diagram
plot_eep(track, eep)

## Classic log_Teff versus log_g plot:
plot_eep(track, eep, cols = c('log_Teff', 'log_g'))
}
