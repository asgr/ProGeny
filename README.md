# ProGeny

**ProGeny** is a stellar population library (SPL) that generates highly flexible simple/single stellar populations (SSPs) that can be directly loaded into the spectral energy distribution (SED) generating and fitting code **ProSpect**.

Details will be added over time, but the base package can be installed from this repo. Below we show a very simple (and fully default) example of how you can construct an SSP.

## Installation

You will need some other CRAN packages to install **ProGeny**, but these should install automatically. To install ProGeny itseld you should run:

```r
install.packages("remotes") # if you don't already have remotes
remotes::install_github("asgr/ProGeny")
```

That should be pretty painless since it doesn't require any compilation.

## Using ProGeny

### First download the desired isochrones and stellar atmospheres

```r
progenyIsoDownload()
progenyAtmosDownload()
```

Where below we assume the isochrone being used has been saved at path/to/Iso.fst, and the stellar atmospheres are in path/to/atmos/.

Load isochrone:

```r
install.packages("fst") # if you don't already have fst
Iso  = fst::read.fst('path/to/Iso.fst', as.data.table=TRUE)
```

Load stellar atmospheres:

```r
Spec_combine = progenyAtmosLoad('path/to/atmos/')
```

Interpolate stellar atmosheres onto isochrone grid:

```r
Interp_combine = progenyInterpGrid_All(Iso=Iso, Spec_combine=Spec_combine)
```

Select best stellar atmosphere library:

```r
Iso = progenyInterpBest(Iso=Iso, Interp_combine=Interp_combine)
```

Note above you might want to set label_AGB and label_white to restrict how AGB and white dwarf atmospheres are added.

Generate a full ProSpect ready SSP:

```r
SSP = progenyMakeSSP(Iso=Iso, IMFfunc=IMF_Chabrier, Spec_combine=Spec_combine,
  Interp_combine=Interp_combine)
```

Check this with the ProSpect function:

```r
library(ProSpect)
speclib_check(SSP)
```
