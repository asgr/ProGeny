progenyIso2Spec = function(logAge=8.4, logZ=0, Iso, IMFint, Spec_combine, Interp_combine, interp=FALSE){
  setDT(Iso)

  Lum = ID = wt = wt_sum = Teff = best = NULL

  #something here - change allZ to ParsecIso
  wave_grid = Spec_combine$base$wave
  spec_stack = rep(0, length(wave_grid))

  use = sort(unique(Iso$best))
  use = use[use > 0]

  IsoAge = sort(unique(Iso$logAge))
  IsoZ = sort(unique(Iso$logZ))

  if(interp){
    if(logAge < min(IsoAge) | logAge > max(IsoAge)){
      stop('Requested logAge not in range of isochrone logAge values!')
    }

    if(logZ < min(IsoZ) | logZ > max(IsoZ)){
      stop('Requested logZ not in range of isochrone logZ values!')
    }

    temp_Age = interp_quick(logAge, IsoAge)
    logAge_loc = as.integer(temp_Age[1:2])
    if(logAge_loc[1] == logAge_loc[2]){
      logAge_steps = IsoAge[logAge_loc[1]]
      logAge_wt = 1
    }else{
      #logAge_diff = IsoAge[logAge_loc[2]] - IsoAge[logAge_loc[1]]
      #logAge nearer to IsoAge[logAge_loc[1]] will have more weight (first element larger) since distance to IsoAge[logAge_loc[2]] is larger
      #logAge_wt = c((IsoAge[logAge_loc[2]] - logAge)/logAge_diff, (logAge - IsoAge[logAge_loc[1]])/logAge_diff)
      logAge_steps = IsoAge[logAge_loc]
      logAge_wt = temp_Age[3:4]
    }

    temp_Z = interp_quick(logZ, IsoZ)
    logZ_loc = as.integer(temp_Z[1:2])
    #logZ_loc = c(max(which(IsoZ <= logZ)), min(which(IsoZ >= logZ)))
    if(logZ_loc[1] == logZ_loc[2]){
      logZ_steps = IsoZ[logZ_loc[1]]
      logZ_wt = 1
    }else{
      #logZ_diff = IsoZ[logZ_loc[2]] - IsoZ[logZ_loc[1]]
      #logAge nearer to IsoZ[logZ_loc[1]] will have more weight (first element larger) since distance to IsoZ[logZ_loc[2]] is larger
      #logZ_wt = c((IsoZ[logZ_loc[2]] - logZ)/logZ_diff, (logZ - IsoZ[logZ_loc[1]])/logZ_diff)
      logZ_steps = IsoZ[logZ_loc]
      logZ_wt = temp_Z[3:4]
    }
  }else{
    logAge_steps = IsoAge[which.min(abs(IsoAge - logAge))]
    logZ_steps = IsoZ[which.min(abs(IsoZ - logZ))]
    logAge_wt = 1
    logZ_wt = 1
    if(logAge_steps != logAge){
      message('Using ', logAge_steps, ' as nearest logAge to the requested ', logAge)
    }
    if(logZ_steps != logZ){
      message('Using ', logZ_steps, ' as nearest logZ to the requested ', logZ)
    }
  }

  for(i in use){
    temp_stack = {}
    for(logAge_j in seq_along(logAge_steps)){
      logAge_step = logAge_steps[logAge_j]
      for(logZ_k in seq_along(logZ_steps)){
        logZ_step = logZ_steps[logZ_k]
        interp_wt = logAge_wt[logAge_j]*logZ_wt[logZ_k]
        subset = Iso[logZ==logZ_step & logAge==logAge_step & best==i, which=TRUE] #faster!
        if(length(subset) > 0){
          temp_DT = data.table(ID=as.vector(Interp_combine[[i]]$nn.idx[subset,]), wt=as.vector(Interp_combine[[i]]$weights[subset,])*Iso[subset,Lum]*IMFint[subset])
          temp_DT = temp_DT[ID > 0,]
          pre_stack = temp_DT[,list(wt_sum=sum(wt)),by=ID]
          pre_stack[,wt_sum:=wt_sum*interp_wt]
          temp_stack = rbind(temp_stack, pre_stack)
        }
      }
    }
    if(!is.null(temp_stack)){
      temp_stack = temp_stack[,list(wt_sum=sum(wt_sum)),by=ID]
      #spec_stack = spec_stack + as.numeric(t(temp_stack$V1) %*% Spec_combine[[i]]$spec[temp_stack$ID,])*interp_wt
      spec_stack = spec_stack + as.numeric(crossprod(temp_stack$wt_sum, Spec_combine[[i]]$spec[temp_stack$ID,]))
    }
  }

  subset = Iso[logZ==logZ_step & logAge==logAge_step & best==0L, which=TRUE] #faster!

  if(length(subset > 0)){
    for(sel in subset){
      spec_stack = spec_stack + .blackbody_simp(wave=wave_grid, Temp=Iso[sel,Teff], norm=Iso[sel,Lum]*IMFint[sel])*interp_wt
    }
  }

  return(invisible(spec_stack))
}
