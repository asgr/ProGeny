\name{progenyIsochrone}
\alias{progenyIsochrone}
\alias{progenyTrackInterp}
\alias{progenyFindMass}
\alias{progenyExtendIso}
\title{
Isochrone Creation Tools
}
\description{
Functions to convert MESA/Iso EEP tracks into MIST like isochrone ready for use in \code{ProGeny}.
}
\usage{
progenyTrackInterp(tracklist, target, make_iso = TRUE, logAge_lim = c(5, 10.3),
  logAge_bin = 0.05, iso_type = "approx", cores = 8, logZ_use = 0, ...)

progenyFindMass(tracklist, logAge_lim = c(5,10.3), logAge_bin = 0.05, logZ_use = 0)

progenyExtendIso(Iso_base, Iso_extend, label = NA, logA = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracklist}{
A full track list output from MESA that has already been passed through Aaron Dotter's "Iso" software, creating tags for key phases of stellar evolution. The original files will have the ending "eep", and probably look something like "17500M.track.eep" for a 175 solar mass initial mass star etc. Basically, all the tracks should be individually loaded and concatenated together into one large data.table. Only a few columns are needed for \code{ProGeny}: logZ (log10(Z/Zol) from MESA/Iso); logAge (log10('star_age') from MESA/Iso); Mini (initial 'star_mass' from MESA/Iso); Mass (current 'star_mass' from MESA/Iso); Lum (10^'log_L' from MESA/Iso); Teff (10^'log_Teff' from MESA/Iso); logG ('log_g' from MESA/Iso); label ('phase' from Iso).
}
  \item{target}{
The target masses to create interpolated tracks/isochrones for. If this is a vector input the same masses will be used for all ages, else it can be a data.table input with columns Mini, logAge_lo and logAge_hi. In this cases Mini should be in units of linear Msol; logAge_lo and logAge_hi should be log10(Age/Yr). If the input is a data.table with these columns then the output tracks use only the mass ranges specified for the relevant age limits. The idea is you target the appropriate mass resolution for each age using \code{progenyFindMass} first, where the output of that becomes the input for \code{progenyTrackInterp}.
}
  \item{make_iso}{
Logical; should the tracks be interpolated onto isochrones (TRUE, default) or the full resolution EEPs (FALSE).
}
  \item{logAge_lim}{
Numeric vector; the lower and upper age limits for the isochrones in log10(Age/Yr) units. Default is as used for MIST.
}
  \item{logAge_bin}{
Numeric scalar; the isochrone binning to use in dex. Default is as used for MIST.
}
  \item{iso_type}{
Character scalar; type of isochrone calculation to use. Default ("approx") uses the \code{\link{approxfun}} code, the other option ("integral") uses the \code{ProSpect} \code{specReBin} function to integrate inside the bin. The latter is more accurate and both take a similar amount of time.
}
  \item{cores}{
Integer scalar; the number of cores to run on.
}
  \item{logZ_use}{
Numeric scalar; the value of \option{logZ} to use. Isochrones can only be computed meaningfully for a specific metallicity, i.e. these functions need to be run for one value of metallicity at a time.
}
  \item{Iso_base}{
Data.table; isochrone that you want to be extended. This should be missing post-AGB stellar evolution (i.e. remnants) since they are what will be added via \option{Iso_extend}.
}
  \item{Iso_extend}{
Data.table; isochrone that you want to extend with. This should be possessing post-AGB stellar evolution (i.e. remnants) since they are what will be added to \option{Iso_base}.
}
  \item{label}{
Scalar; label to use for the extended isochrones. Options are a value (relevant for \option{Iso_base}); 'base_get', which will find the current largest \option{label} value in \option{Iso_base} and add 1 to it; 'extend_get', which will use whatever is currently in \option{Iso_extend}. Default is NA, which will just set the extended isochrone \option{label} to NA in the output, which is useful for identifying what entries have been added.
}
  \item{logA}{
Scalar; logA to use for the extended isochrones. Options are a value (relevant for \option{Iso_base}); 'base_get', which will find the current largest \option{logA} value in \option{Iso_base} and use that; 'extend_get', which will use whatever is currently in \option{Iso_extend}. Default is NA, which will just set the extended isochrone \option{label} to NA in the output, which is useful for identifying what entries have been added. Note this is only relevant if one or both of \option{Iso_base} and/or \option{Iso_extend} contain a \option{logA} column.
}
  \item{\dots}{
Other arguments to pass to \code{\link{approxfun}} and \code{specReBin} depending on the setting of \option{iso_type}.
}
}
\details{
Generally once you have loaded all your track EEP files into the required \option{tracklist} data.table format, you would pass this into \code{progenyFindMass} for you target isochrone grid, and then pass the output of this into the \option{target} argument of \code{progenyTrackInterp} with the isochrone arguments matching. The output of this should then be an isochrone ready to use with other \code{ProGeny} functions for a single metallicity (i.e. you will need to loop through all the \option{logZ} values you want for the final isochrone).

Note there is currently no option of interpolating between metallicity values. Aaron Dotter's Iso also only interpolates within a specific metallicity. A path to doing this manually (for now) would be to run \code{progenyFindMass} for the logZ values either side of the one you want to create; combine them (\code{\link{rbind}}) and compute the unique rows with \code{\link{unique}}; create two target isochrones using the same combined input \option{target} argument (so we know the resulting rows will match 1-1); combine the two output isochrones with linear logZ weights (i.e. if you have logZ 0 and 1 and want 0.5 the weights would be 0.5:0.5 respectively, but if you want 0.2 the weights would be 0.8,0.2 respectively).

\code{progenyExtendIso} will extend \option{Iso_base} using \option{Iso_extend}. The assumption is the provided \option{Iso_base} is missing evolved post-AGB stars and the current end points of stellar evolution will be found; then using \option{Iso_extend} an effort is made to extend these where possible. \option{Iso_extend} must therefore contain post-AGB stars (e.g. MIST or PARSEC). The best matching \option{logZ} is used (with no interpolation), so there may be issues if \option{Iso_base} has a very different \option{logZ} coverage than \option{Iso_extend}. The problem is most serious if it extends outside the range (lower or higher).
}
\value{
\code{progenyTrackInterp} outputs an isochrone data.table ready to use in the rest of \code{ProGeny} (usually the \option{Iso} argument). \cr

\code{progenyFindMass} creates the \option{target} input ready to use in \code{progenyTrackInterp} for a target isochrone. \cr

\code{progenyExtendIso} outputs an isochrone data.table ready to use in the rest of \code{ProGeny} (usually the \option{Iso} argument).
}
\references{
Aaron Dotter, 2016, ApJS, 222, 8D
}
\author{
Aaron Robotham (not Dotter!)
}
\seealso{
\code{\link{progenyMakeSSP}}
}
\examples{
\dontrun{
# Assume you have loaded all your EEP track files into a single data.table
# called 'tracklist', with columns: logZ, star_age, star_mass (initial),
# star_mass (current), log_L, log_Teff, log_g, and phase.
# Then rename/create the required ProGeny columns:
# logZ, logAge (log10(star_age)), Mini (initial mass), Mass (current mass),
# Lum (10^log_L), Teff (10^log_Teff), logG (log_g), label (phase).

# Find the optimal mass resolution for each age bin at solar metallicity:
target = progenyFindMass(tracklist, logAge_lim = c(5, 10.3), logAge_bin = 0.05,
  logZ_use = 0)

# Interpolate tracks onto a regular isochrone grid for solar metallicity:
Iso_solar = progenyTrackInterp(tracklist, target = target, make_iso = TRUE,
  logAge_lim = c(5, 10.3), logAge_bin = 0.05, logZ_use = 0)

# Repeat for each metallicity in your tracklist and then row-bind the results:
# Iso = rbind(Iso_solar, Iso_metal_poor, ...)

# If you have a base isochrone missing post-AGB phases, extend it:
# Iso_full = progenyExtendIso(Iso_base = Iso, Iso_extend = Iso_PARSEC)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ ~isochrone }
% \concept{ ~cpt2 }
% Use only one concept per line.
