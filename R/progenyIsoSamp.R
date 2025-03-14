progenyIsoSamp = function(Iso, IMFfunc, ..., Nsamp=1e4, Lum_weight_func=NULL,
  logAge_weight_func=NULL, logZ_weight_func=NULL, label_use = seq(-1,9), seed=666){
  set.seed(seed)
  setDT(Iso)

  #To be safe to not alter the original Isochrones
  Iso_temp = copy(Iso)
  setkeyv(Iso_temp, c('logZ', 'logAge', 'Mini'))

  logZ = logAge = Mini = Mass = NULL

  Iso_temp[,IMFint := 0]

  message("Integrating the IMF")
  #The below should work more generically (regardless of isochrone)
  Iso_temp[Lum > 1e-6,IMFint := progenyUpdateIMF(Iso_temp[Lum > 1e-6,], IMFfunc, ...)$IMFint]

  message("Applying weights")
  if(!is.null(Lum_weight_func)){
    Iso_temp[,IMFint:=IMFint*Lum_weight_func(Lum)]
  }

  if(!is.null(logAge_weight_func)){
    Iso_temp[,IMFint:=IMFint*logAge_weight_func(logAge)]
  }

  if(!is.null(logZ_weight_func)){
    Iso_temp[,IMFint:=IMFint*logZ_weight_func(logZ)]
  }

  message("Creating inversion function")

  Iso_temp[IMFint == 0,IMFint := 1e-300]
  Iso_temp = Iso_temp[label %in% label_use,]

  temp_cumsum = cumsum(Iso_temp$IMFint)
  temp_cumsum = temp_cumsum/max(temp_cumsum)

  sel_dup = which(duplicated(temp_cumsum)) - 1L
  if(length(sel_dup) > 0){
    temp_cumsum = temp_cumsum[-sel_dup]
    Iso_temp = Iso_temp[-sel_dup,]
  }

  temp_map = approxfun(temp_cumsum, 1:length(temp_cumsum), yleft = 1, yright = length(temp_cumsum))

  message("Sampling")
  sel = ceiling(temp_map(runif(Nsamp)))

  return(Iso_temp[sel,])
}
