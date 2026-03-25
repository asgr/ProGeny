.onLoad<-function(libname, pkgname){
  packageStartupMessage(
'Isochrones: Download with progenyIsoDownload. Load with read.fst
Atmospheres: Download with progenyAtmosDownload. Load with progenyAtmosLoad
Interpolate grids with progenyInterpGrid_All
Update the best match between Isochrones and Atmospheres with progenyInterpBest
Generate your SSP using the input above with progenyMakeSSP, e.g. minimally:\n
Iso  = read.fst(\'path/to/Iso.fst\', as.data.table=TRUE)
Spec_combine = progenyAtmosLoad(\'path/to/atmos/\')
Interp_combine = progenyInterpGrid_All(Iso=Iso, Spec_combine=Spec_combine)
Iso = progenyInterpBest(Iso=Iso, Interp_combine=Interp_combine)
progenyMakeSSP(Iso=Iso, Spec_combine=Spec_combine, Interp_combine=Interp_combine)'
)
}

utils::globalVariables(
  c("..col_keep", "..xsel", "..zsel", "IMFint", "Isochrone", "Lum", "Mass", "Mini",
    "Teff", "V1", "V2", "best", "i", "label", "logA", "logAge", "logAge_hi", "logAge_lo",
    "logG", "logG_diff", "logT_diff", "logZ", "logZ_diff")
)
