progenyMakeSSP = function(Iso, IMFfunc = IMF_Chabrier, ..., rem_frac = 'get', Spec_combine,
                          Interp_combine, Zsol=0.02, cores=8, Labels = list(
                            Zlab = "Metallicity", Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
                            Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)")){
  #currently takes about a minute or two to generate a target SSP on 8 cores
  logZ_steps = unique(Iso$logZ)
  logAge_steps = unique(Iso$logAge)

  setDT(Iso)

  #To be safe to not alter the original Isochrones
  Iso_temp = copy(Iso)
  setkeyv(Iso_temp, c('logZ', 'logAge', 'Mini'))

  logZ = logAge = Mini = Mass = NULL

  Iso_temp[,IMFint := 0]

  #The below should work more generically (regardless of isochrone)
  Iso_temp[Lum > 1e-6,IMFint := progenyUpdateIMF(Iso_temp[Lum > 1e-6,], IMFfunc, ...)$IMFint]

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

  message('Generating evo grids:')

  Zevo = foreach(logZ_step = logZ_steps)%dopar%{
    if(rem_frac == 'get'){
      temp_cut = Iso_temp[logZ == logZ_step & logAge >= 8,data.table(Mini,Mass/Mini)[which.max(Mini),],by=logAge]
      temp_func = approxfun(temp_cut[,Mini], temp_cut[,V2], rule=2)
    }else{
      temp_func = function(x){rem_frac}
    }

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
      SMrem = SMrem_miss #might be missing some for some isochrones though (this is literally remnants we stop tracking, hence SM_rem_miss)
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

  return(invisible(SSP))
}
