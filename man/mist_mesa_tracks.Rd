\name{mist_mesa_tracks}
\alias{mist_mesa_tracks}
\alias{track000p1.eep}
\alias{track001p0.eep}
\alias{track020p0.eep}
\alias{track300p0.eep}
\title{Bundled MIST v1.2 MESA EEP Stellar Tracks}
\description{
Example MIST v1.2 evolutionary tracks from MESA, stored in
\code{inst/extdata} as EEP tables for four initial stellar masses:
\code{0.1}, \code{1}, \code{20}, and \code{300}~\eqn{M_\odot}.
}
\details{
The bundled files are:
\itemize{
  \item \code{track000p1.eep}: 0.1~\eqn{M_\odot}
  \item \code{track001p0.eep}: 1~\eqn{M_\odot}
  \item \code{track020p0.eep}: 20~\eqn{M_\odot}
  \item \code{track300p0.eep}: 300~\eqn{M_\odot}
}

These are plain-text, whitespace-delimited tables with a header row.
Column naming and definitions follow the MIST README tables documentation.
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

Phase codes present in these bundled tracks are:
\code{-1}, \code{0}, \code{2}, \code{3}, \code{4}, \code{5}, \code{6}, and \code{9}.
See the MIST README tables for full phase-code definitions.
}
\source{
MIST (MESA Isochrones and Stellar Tracks), v1.2 tables and documentation:
\url{https://mist.science/README_tables.pdf}
}
\examples{
f1 <- system.file("extdata", "track001p0.eep", package = "ProGeny")
if (f1 == "") stop("track001p0.eep not found in extdata")
track1 <- read.table(f1, header = TRUE)

f20 <- system.file("extdata", "track020p0.eep", package = "ProGeny")
if (f20 == "") stop("track020p0.eep not found in extdata")
track20 <- read.table(f20, header = TRUE)
}
\keyword{datasets}
