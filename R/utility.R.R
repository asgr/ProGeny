.blackbody_simp = function(wave, Temp, norm=1){
  return(norm*cosplanckLawRadWave(lambda=wave/1e10, Temp=Temp)/cosplanckSBLawRad_sr(Temp)/1e10)
}

.binlims = function(input, log=FALSE){
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

  return(list(lo = bin_lo, hi = bin_hi))
}

progenyIsoFormat = function(Iso){
  setDT(Iso)
  setkeyv(Iso, c('logZ', 'logAge', 'Mini'))
  return(Iso)
}
