# ProGeny: Stellar Population Library

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/ProGeny/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

**ProGeny** is a stellar population library (SPL) and software package that generates highly flexible simple/single stellar populations (SSPs) that can be directly loaded into the spectral energy distribution (SED) generating and fitting code **ProSpect**. The software is written in **R** and is released under a permissive LGPL-3 license.

If you just want to create SSPs with ProGeny using a simple UI over the web, then please try [progeny.datacentral.org.au](https://progeny.datacentral.org.au). This contains a good fraction of the core functionality, but lacks the ultra high resolution spectral atmospheres, so less well suited to higher resolution (over R~5000 in optical) work.

**ProGeny** is described in detail in the core software paper (ProGeny I: Robotham & Bellstedt, RASTI, 4, 19) and the associated application paper (ProGeny II: Bellstedt & Robotham, MNRAS, 540, 2703), and these papers should both be cited in any work making use of the package, interactive tool (either locally or through the web interface), or its outputs (the static FITS SSPs). 

Details will be added over time, but the base package can be installed from this repo. Below we show a very simple (and fully default) example of how you can construct an SSP. For more complex vignettes see [rpubs.com/asgr](https://rpubs.com/asgr).

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

Interpolate stellar atmospheres onto isochrone grid:

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

Check this with the **ProSpect** function:

```r
remotes::install_github("asgr/ProSpect") #Install ProSpect if needed
library(ProSpect)
speclib_check(SSP)
```

We can save it to use later:

```r
remotes::install_github("asgr/Rfits") #Install Rfits if needed
library(Rfits)
Rfits_write(SSP, 'path/to/SSP.fits', flatten = T)
```

We can load this back into the **ProSpect** format:

```r
SSP2 = speclib_FITSload('path/to/SSP.fits')
```

Now we run it through **ProSpect** with the default SFH (constant) to make our first SED:

```r
SSP_pro = ProSpectSED(speclib=SSP2) # Obviously ProSpectSED(speclib=SSP) also works
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

That is just the start, see [rpubs.com/asgr](https://rpubs.com/asgr) for longer vignettes!
