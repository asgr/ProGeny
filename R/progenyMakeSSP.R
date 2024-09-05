progenyMakeSSP = function(Iso, IMFfunc, ..., Spec_combine,
                          Interp_combine, Zsol=0.02, cores=8, Labels = list(
                            Zlab = "Metallicity", Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
                            Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)")){
  #currently takes about a minute or two to generate a target SSP on 8 cores
  logZ_steps = unique(Iso$logZ)
  logAge_steps = unique(Iso$logAge)

  setDT(Iso)

  Iso_temp = copy(Iso)
  setkeyv(Iso_temp, c('logZ', 'logAge', 'Mini'))

  logZ = logAge = Mini = Mass = NULL

  Iso_temp[,IMFint := progenyUpdateIMF(Iso_temp, IMFfunc, ...)$IMFint]
  #Iso_temp[,ID:= 1:dim(Iso_temp)[1]]

  cores = min(cores, length(logZ_steps), detectCores())
  registerDoParallel(cores=cores)
  #
  logAge_step = NULL
  logZ_step = NULL
  #
  message('Generating spec grids:')

  Zspec = foreach(logZ_step = logZ_steps)%dopar%{
    message('  ',logZ_step)
    output = foreach(logAge_step = logAge_steps)%do%{
      #message('    ',logAge_step)
      progenyIso2Spec(logAge_step, logZ_step, Iso=Iso_temp, IMFint=Iso_temp$IMFint, Spec_combine, Interp_combine)
    }
    return(do.call(rbind, output))
  }

  # Iso_stack = Iso_temp[,list(ID=list(ID)), by=list(logAge,logZ)]
  #
  # Zspec = foreach(i = 1:dim(Iso_stack))%dopar%{
  #   #for(i in 1:dim(Iso_stack)[1]){
  #   .progenyIso2SpecSub(subset=unlist(Iso_stack[i,ID]), Iso=Iso_temp, Iso_temp$IMFint, Spec_combine, Interp_combine)
  # }
  #need to check this a bit more carefully!

  message('Generating evo grids:')

  Zevo = foreach(logZ_step = logZ_steps)%dopar%{
    temp_cut = Iso_temp[logZ == logZ_step & logAge >= 8,data.table(Mini,Mass/Mini)[which.max(Mini),],by=logAge]
    temp_func = approxfun(temp_cut[,Mini], temp_cut[,V2], rule=2)

    missing_func = function(lo_lim){
        integrate(function(x){IMFfunc(x, massmult=TRUE, ...)*temp_func(x)}, lower=lo_lim, upper=Inf)$value
    }

    message('  ',logZ_step)
    temp_out = foreach(logAge_step = logAge_steps, .combine='rbind')%do%{
      mass_lim = Iso_temp[logZ == logZ_step & logAge == logAge_step, max(Mini)]
      SMrem_miss = missing_func(mass_lim)
      SMstar = Iso_temp[logZ == logZ_step & logAge == logAge_step, sum(Mass*IMFint)] + SMrem_miss
      SMstar[SMstar > 1] = 1
      SMgas = 1 - SMstar
      SMtot = 1
      SFR = 0
      SMrem = 0
      return(data.frame(SMstar=SMstar, SMgas=SMgas, SMtot=SMtot, SFR=SFR, SMrem=SMrem))
    }
    temp_out[1,'SFR'] = 1
    return(temp_out)
  }

  Age_lims = .binlims(10^logAge_steps, log=T)
  AgeBins = c(Age_lims$lo[1], Age_lims$hi)
  AgeWeights = diff(AgeBins)

  SSP = list(
    Z = Zsol * 10^logZ_steps,
    Age = 10^logAge_steps,
    AgeBins = AgeBins,
    AgeWeights = AgeWeights,
    Wave = Spec_combine$base$wave,
    Labels = Labels,
    Zspec = Zspec,
    Zevo = Zevo
  )

  return(SSP)
}
