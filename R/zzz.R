.onLoad<-function(libname, pkgname){
  packageStartupMessage(
'Isochrones: Download with progenyIsoDownload. Load with read.fst
Atmospheres: Download with progenyAtmosDownload. Load with atmos_load
Interpolate grids with progenyInterpGrid_All
Update the best match between Isochrones and Atmospheres with progenyInterpBest
Generate your SSP using the input above with progenyMakeSSP, minimally:\n
Iso  = read.fst(\'path/to/Iso.fst\', as.data.table=TRUE)
Spec_combine = progenyAtmosLoad(\'path/to/atmos/\')
Interp_combine = progenyInterpGrid_All(Iso=Iso, Spec_combine=Spec_combine)
Iso = progenyInterpBest(Iso=Iso, Interp_combine=Interp_combine)
progenyMakeSSP(Iso=Iso, Spec_combine=Spec_combine, Interp_combine=Interp_combine)'
)
}
