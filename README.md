# ProGeny: Stellar Population Library

**ProGeny** is a stellar population library (SPL) that generates highly flexible simple/single stellar populations (SSPs) that can be directly loaded into the spectral energy distribution (SED) generating and fitting code **ProSpect**.

Details will be added over time, but the base package can be installed from this repo. Below we show a very simple (and fully default) example of how you can construct an SSP. For more complex vignettes see (rpubs.com/asgr/)[https://rpubs.com/asgr/].

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

We can save it to use later:

```r
Rfits_write(SSP, 'path/to/SSP.fits', flatten = T)
```

We can load this into **ProSpect** easily:

```r
SSP2 = speclib_FITSload('path/to/SSP.fits')
```

Run it through **ProSpect** with the default SFH (constant):

```r
SSP2_pro = ProSpectSED(speclib=SSP2)
```

Plot it!

```r
plot(SSP2_pro, ylim=c(1e3,1e7))
```

We can compare that to the classic BC03 (high res) SSP (one line for convenience):

```r
data(BC03hr)
plot(ProSpectSED(speclib=BC03hr), ylim=c(1e3,1e7))
```

That is just the start, see (rpubs.com/asgr/)[https://rpubs.com/asgr/] for longer vignettes!
