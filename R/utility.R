.blackbody_simp = function(wave, Temp, norm=1){
  return(norm*cosplanckLawRadWave(lambda=wave/1e10, Temp=Temp)/cosplanckSBLawRad_sr(Temp)/1e10)
}

.binlims = function(input, log=FALSE, max_width=NULL){
  # max_width: maximum bin width.
  #   - if log=FALSE: additive width (same units as input)
  #   - if log=TRUE : multiplicative width in dex, i.e. cap=0.1 means x10^0.1
  # max_width can be useful to catch pathological gaps in isochrones, but this can leave gaps.

  N_in = length(input)

  if(log){
    input = log10(input)
  }

  in_diff = diff(input)

  bins = input[1:(N_in - 1L)] + in_diff/2
  bins = c(input[1] - in_diff[1]/2, bins, input[N_in] + in_diff[N_in - 1L]/2)

  if(log){
    bins = 10^bins
  }

  bin_lo = bins[1:N_in]
  bin_hi = bins[2:(N_in + 1L)]

  if(!is.null(max_width)){
    if (!is.numeric(max_width) || length(max_width) != 1L || is.na(max_width)) stop("max_width must be a single non-NA number")
    if (max_width < 0) stop("max_width must be >= 0")

    if(log){
      bin_lim = bin_lo*10^max_width
    }else{
      bin_lim = bin_lo + max_width
    }

    bin_hi = pmin(bin_hi, bin_lim)
  }

  return(list(lo = bin_lo, hi = bin_hi))
}

progenyIsoFormat = function(Iso){
  setDT(Iso)
  setkeyv(Iso, c('logZ', 'logAge', 'Mini'))
  return(Iso)
}

#Scalar interpolation
.interp_quick = function(x, params, log=FALSE){
  if(length(x) > 1){stop('x must be scalar!')}
  if(x < min(params)){
    return(c(ID_lo=1, ID_hi=1, wt_lo=1, wt_hi=0))
  }
  if(x > max(params)){
    return(c(ID_lo=length(params), ID_hi=length(params), wt_lo=0, wt_hi=1))
  }
  if(log){
    params = log(params)
    x = log(x)
  }
  interp = approx(params, 1:length(params), xout=x)$y
  IDlo = floor(interp)
  IDhi = ceiling(interp)
  return(c(ID_lo=IDlo, ID_hi=IDhi, wt_lo=1-(interp-IDlo), wt_hi=interp-IDlo))
}
