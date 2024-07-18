progenyMakeSSP = function(Iso, IMFfunc, masslow = 0.1, massmax = 100, ..., Spec_combine,
                          Interp_combine, Zsol=0.02, cores=8, Labels = list(
  Zlab = "Metallicity", Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
  Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)")){
  #currently takes about a minute or two to generate a target SSP on 8 cores
  logZ_steps = unique(Iso$logZ)
  logAge_steps = unique(Iso$logAge)

  setDT(Iso)

  Iso_temp = copy(Iso)

  logZ = logAge = Mini = Mass = NULL

  Iso_temp[,IMFint:= progenyUpdateIMF(Iso_temp, IMFfunc, masslow=masslow, massmax=massmax, ...)$IMFint]

  cores = min(cores, length(logZ_steps))
  registerDoParallel(cores=cores)

  logAge_step = NULL
  logZ_step = NULL

  message('Generating spec grids:')

  Zspec = foreach(logZ_step = logZ_steps)%dopar%{
    message('  ',logZ_step)
    foreach(logAge_step = logAge_steps, .combine='rbind')%do%{
      #message('    ',logAge_step)
      progenyIso2Spec(logAge_step, logZ_step, Iso=Iso_temp, Iso_temp$IMFint, Spec_combine, Interp_combine)
    }
  }

  #need to check this a bit more carefully!

  message('Generating evo grids:')

  Zevo = foreach(logZ_step = logZ_steps)%dopar%{
    message('  ',logZ_step)
    foreach(logAge_step = logAge_steps, .combine='rbind')%do%{
      message('    ',logAge_step)
      SMstar = Iso_temp[logZ == logZ_step & logAge == logAge_step, sum(Mass*IMFint)]
      SMstar[SMstar > 1] = 1
      SMgas = 1 - SMstar
      SMtot = 1
      SFR = rep(0, length(SMstar))
      SFR[1] = 1
      SMrem = 0
      return(data.frame(SMstar=SMstar, SMgas=SMgas, SMtot=SMtot, SFR=SFR, SMrem=SMrem))
    }
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
