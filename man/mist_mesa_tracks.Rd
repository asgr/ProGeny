\name{mist_mesa_tracks}
\alias{mist_mesa_tracks}
\alias{track000p1}
\alias{track001p0}
\alias{track020p0}
\alias{track300p0}
\title{Bundled MIST v1.2 MESA EEP Stellar Tracks}
\description{
Example MIST v1.2 EEP tracks from MESA, stored in
\code{inst/extdata} as EEP tables for four initial stellar masses:
\code{0.1}, \code{1}, \code{20}, and \code{300}~\eqn{M_\odot}.
}
\format{
Each file is a tabular dataset with columns:
\describe{
  \item{\code{star_age}}{Age (years).}
  \item{\code{star_mass}}{Current stellar mass (\eqn{M_\odot}).}
  \item{\code{he_core_mass}}{Helium core mass (\eqn{M_\odot}).}
  \item{\code{c_core_mass}}{Carbon core mass (\eqn{M_\odot}).}
  \item{\code{o_core_mass}}{Oxygen core mass (\eqn{M_\odot}).}
  \item{\code{log_L}}{Log10 luminosity relative to solar (\eqn{\log_{10}(L/L_\odot)}).}
  \item{\code{log_LH}}{Log10 hydrogen-burning luminosity relative to solar.}
  \item{\code{log_LHe}}{Log10 helium-burning luminosity relative to solar.}
  \item{\code{log_Teff}}{Log10 effective temperature (K).}
  \item{\code{log_g}}{Log10 surface gravity (cgs).}
  \item{\code{log_center_T}}{Log10 central temperature (K).}
  \item{\code{log_center_Rho}}{Log10 central density (cgs).}
  \item{\code{center_gamma}}{Central Coulomb coupling parameter \eqn{\Gamma}.}
  \item{\code{center_h1}}{Central mass fraction of \eqn{^1\mathrm{H}}.}
  \item{\code{center_he4}}{Central mass fraction of \eqn{^4\mathrm{He}}.}
  \item{\code{center_c12}}{Central mass fraction of \eqn{^{12}\mathrm{C}}.}
  \item{\code{phase}}{Integer evolutionary phase flag (MIST phase coding).}
}
}
\details{
The bundled files are:
\itemize{
  \item \code{track000p1}: 0.1~\eqn{M_\odot}
  \item \code{track001p0}: 1~\eqn{M_\odot}
  \item \code{track020p0}: 20~\eqn{M_\odot}
  \item \code{track300p0}: 300~\eqn{M_\odot}
}

These are plain-text, whitespace-delimited tables with a header row.
Column naming and definitions follow the MIST README tables documentation.
}
\section{Phase code mapping}{
The integer column \code{label} in the EEP track maps EEP phases to coarse MIST-like phases as follows:
\itemize{
  \item{\code{-1}: Pre-Main Sequence (PreMS)}
  \item{\code{0}: Main Sequence (ZAMS, IAMS)}
  \item{\code{2}: Red Giant Branch (RGB, TAMS)}
  \item{\code{3}: Core Helium Burning (RGBTip, ZACHeB)}
  \item{\code{4}: Early AGB (TACHeB; if present)}
  \item{\code{5}: TP-AGB (TPAGB; if present)}
  \item{\code{6}: Post AGB (PostAGB; if present)}
  \item{\code{7}: Central Star Planetary Nebula (CSPN; if present)}
  \item{\code{8}: White Dwarf Central Star (WDCS; if present)}
  \item{\code{9}: Wolf-Rayet (WR; if present)}
}

See the MIST README tables for full phase-code definitions.
}
\source{
MIST (MESA Isochrones and Stellar Tracks), v1.2 tables and documentation:
\url{https://mist.science/README_tables.pdf}
}
\examples{
f1 <- system.file("extdata", "track001p0.eep", package = "ProGeny")
if (nzchar(f1)) track1 <- read.table(f1, header = TRUE)

f20 <- system.file("extdata", "track020p0.eep", package = "ProGeny")
if (nzchar(f20)) track20 <- read.table(f20, header = TRUE)
}
\keyword{datasets}
