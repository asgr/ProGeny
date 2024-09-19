\name{progenyInterpGrid}
\alias{progenyInterpGrid}
\alias{progenyInterpGrid_All}
\alias{progenyInterpBest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate Isochrone Grids onto Stellar Atmosphere Grids
}
\description{
Fast interpolation between target isochrone stellar paramaters and spectra from stellar atmospheres. Interpolation can use any combination of \option{Teff} (Kelvin), \option{logG} (log10(SurfaceGravity / [cm/s^2])) and \option{logZ} (log10(Z / 0.02).
}
\usage{
progenyInterpGrid(loc, info, rescale, radius = 2, weight_pow = 2, k = 8)

progenyInterpGrid_All(Iso, Spec_combine, radius = 2, weight_pow = 2, k = 8)

progenyInterpBest(Iso, Interp_combine, do_base = TRUE, do_extend = TRUE, do_hot = TRUE,
  do_AGB = TRUE, do_white = TRUE, b2e = 1.5, label_AGB = 7:8, label_white = 9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc}{
Data.frame; locations of \option{Teff} (Kelvin), \option{logG} (log10(SurfaceGravity / [cm/s^2])) and \option{logZ} (log10(Z / 0.02), which mean a value of 0 is solar metallicity) locations to interpolate.
}
  \item{info}{
Data.frame; information for each row of spectra. Can contain columns \option{Teff}, \option{logG} and \option{logZ} as appropriate for the data.
}
  \item{Iso}{
Data.frame / data.table; the input table of Isochrones. This must have a particular format (see ParsecIso in \code{\link{progenyIsoDownload}}).
}
  \item{Spec_combine}{
List; the required spectral atmospheres. This should probably be loaded via \code{\link{progenyAtmosLoad}}, where the inputs for that come via \code{\link{progenyAtmosDownload}}.
}
  \item{rescale}{
Data.table; rescaling to be applied to the \option{info} columns. Can contain columns \option{Teff}, \option{logG} and \option{logZ} as appropriate for the data in \option{info}. The relevant names \option{info} columns are divided by the values of \option{rescale}, so a \option{rescale} value for column \option{Teff} of 2 means all values are halfed. By this mechanism we bring all grid values onto a common search radius scale (\option{radius}).
}
  \item{radius}{
Numeric scalar; the post \option{rescale} search radius to use when interpolating \option{info} column values.
}
  \item{weight_pow}{
Numeric scalar; the interpolation weighting to apply. 1/distance^\option{weight_pow} becomes the per value weighting, where the sum of all \option{rescale} grid values within the specified \option{radius} are normalised to 1 using this scheme. A larger weight mean we bias out interpolated solution even more to closer grid values. A value between 1-2 (inclusive) is probably about optimal.
}
  \item{k}{
Integer scalar; the maximum number of parameter nearest neighbours to compute. Setting this lower is an additional method for limiting interpolation grid positions (beyond the \option{radius}). The default of 8 reflects the vertices in a 3D cube, so is reasonable for a 3D search. For a 2D search \option{k} = 4 might be more appropriate.
}
  \item{Interp_combine}{
List; the output of \code{progenyInterpGrid_All}.
}
  \item{do_base}{
Logical; should the 'base' spectra be considered when trying to determine the best interpolated spectra solution?
}
  \item{do_extend}{
Logical; should the 'extend' spectra be considered when trying to determine the best interpolated spectra solution?
}
  \item{do_hot}{
Logical; should the 'hot' spectra be considered when trying to determine the best interpolated spectra solution?
}
  \item{do_AGB}{
Logical; should the 'AGB' spectra be considered when trying to determine the best interpolated spectra solution?
}
  \item{do_white}{
Logical; should the 'white' dwarf spectra be considered when trying to determine the best interpolated spectra solution?
}
  \item{b2e}{
Numeric scalar; how much closer in re-scaled grid space does the 'extend' spectral atmosphere need to be before it is used instead of the 'base'. Values larger than \option{b2e} = 1 mean we are biased to preferring the 'base' solution even if it is not as close.
}
  \item{label_AGB}{
Integer vector; the \option{Iso} column \option{label} value that corresponds specifically to AGB periods (where the AGB spectra is then preferred). For MIST the setting should probably be 4:5, for PARSEC it should probably be 7:8, and for BaSTI it should probably be 3:5. If left as the default NULL then it will try to match the AGBs purely based of Teff and logG (ignoring any specific isochrone labels).
}
  \item{label_white}{
Integer vector; the \option{Iso} column \option{label} value that corresponds specifically to white dwarf periods (where the white dwarf spectra is then preferred). For MIST the setting should probably be 6:8 (note 9 is WR [Wolf-Rayet], which is a very different and very early phase of star formation, so avoid this!), for PARSEC it should probably be 9, and for BaSTI it should probably be 6. If left as the default NULL then it will try to match the AGBs purely based of Teff and logG (ignoring any specific isochrone labels).
}
}
\details{
TBD
}
\note{
In MIST the isochrones labels are provided as: -1 = PreMS, 0 = MS, 2 = RGB, 3 = CHeB, 4 = eAGB, 5 = TP-AGB, 6-8 = PostAGB, 9 = WR.

In PARSEC the isochrones labels are provided as: 0 = PreMS, 1 = MS, 2 = SubGB, 3 = RGB, 4-6 = CHeB, 7 = eAGB, 8 = TP-AGB, 9 = PostAGB.

In BaSTI the isochrones labels are provided as: -1 = PreMS, 0 = MS, 1 = RGB, 2 = CHeB, 3 = eAGB, 4-5 = AGB, 6 = PostAGB

In general:
PMS: pre main sequence
MS: main sequence
WR: Wolf-Rayet (generally only a few million years old)
SubGB: subgiant branch, or Hertzsprung gap for more intermediate+massive stars
RGB: red giant branch, or the quick stage of red giant for intermediate+massive stars
CHEB: core He-burning for low mass stars, or the very initial stage of CHeB for intermediate+massive stars
eAGB: the early asymptotic giant branch, or a quick stage of red giant for massive stars
TP-AGB: the thermally pulsing asymptotic giant branch
Post-AGB: post asymptotic giant branch
}
\value{
\code{progenyInterpGrid} returns a list output that contain the various outputs of \code{\link{nn2}}, with an additional weights list item that contains our bespokely scaled weights that will be used when interpolating an Isochrone vs Atmosphere solution.

\code{progenyInterpGrid_All} returns a list output from \code{progenyInterpGrid} for each spectral atmosphere present in \option{Spec_combine}.

\code{progenyInterpBest} adds a modifed \option{label} column to \option{Iso}, otherwise the format is the same as the \option{Iso} input.
}
\references{
Nothing yet
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{progenyMakeSSP}}
}
\examples{
## Nothing yet
}

