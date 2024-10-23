progenyIso2Spec = function(logAge=8.4, logZ=0, Iso, IMFint, Spec_combine, Interp_combine){
  setDT(Iso)

  Lum = ID = wt = Teff = best = NULL

  #something here - change allZ to ParsecIso
  wave_grid = Spec_combine$base$wave
  spec_stack = rep(0, length(wave_grid))

  use = sort(unique(Iso$best))
  use = use[use > 0]

  logAge_step = logAge
  logZ_step = logZ

  for(i in use){
    #subset = which(Iso$logZ == logZ & Iso$logAge == logAge & Iso$best == i)

    subset = Iso[logZ==logZ_step & logAge==logAge_step & best==i, which=TRUE] #faster!

    temp_DT = data.table(ID=as.vector(Interp_combine[[i]]$nn.idx[subset,]), wt=as.vector(Interp_combine[[i]]$weights[subset,])*Iso[subset,Lum]*IMFint[subset])
    temp_DT = temp_DT[ID > 0,]
    temp_stack = temp_DT[,sum(wt),by=ID]

    # temp_stack_vec = rep(0, dim(Spec_combine[[i]]$info)[1])
    # temp_stack_vec[temp_stack$ID] = temp_stack$V1
    #
    # spec_stack = spec_stack + as.numeric(t(temp_stack_vec) %*% Spec_combine[[i]]$spec)

    spec_stack = spec_stack + as.numeric(t(temp_stack$V1) %*% Spec_combine[[i]]$spec[temp_stack$ID,])
  }

  #subset = which(Iso$logZ == logZ & Iso$logAge == logAge & Iso$best == 0)
  subset = Iso[logZ==logZ_step & logAge==logAge_step & best==0L, which=TRUE] #faster!

  if(length(subset > 0)){
    for(sel in subset){
      spec_stack = spec_stack + .blackbody_simp(wave=wave_grid, Temp=Iso[sel,Teff], norm=Iso[sel,Lum]*IMFint[sel])
    }
  }

  return(invisible(spec_stack))
}
