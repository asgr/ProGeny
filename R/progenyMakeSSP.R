progenyMakeSSP = function(Iso, IMFfunc = IMF_Chabrier, ..., rem_frac = 'get', Spec_combine,
                          Interp_combine, logAge_steps = NULL, logZ_steps = NULL, Zsol=0.02,
                          cores=8, Labels = list(Zlab = "Metallicity", Agelab = "Time since ZAM / Yrs",
                          Wavelab = "Wavelength / Ang", Lumlab = "Lsun / Ang (for 1 Msun SF)",
                          LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)")){
  #currently takes about a minute or two to generate a target SSP on 8 cores

  interp = FALSE

  if(is.null(logAge_steps)){
    logAge_steps = sort(unique(Iso$logAge))
    interp = TRUE
    message('Interpolating onto logAge grid different to isochrone!')
  }else{
    Iso_logAge_steps = sort(unique(Iso$logAge))
  }

  if(is.null(logZ_steps)){
    logZ_steps = sort(unique(Iso$logZ))
    interp = TRUE
    message('Interpolating onto logZ grid different to isochrone!')
  }else{
    Iso_logZ_steps = sort(unique(Iso$logZ))
  }

  setDT(Iso)

  #To be safe to not alter the original Isochrones
  Iso_temp = copy(Iso)
  setkeyv(Iso_temp, c('logZ', 'logAge', 'Mini'))

  logZ = logAge = Mini = Mass = NULL

  Iso_temp[,IMFint := 0]

  #The below should work more generically (regardless of isochrone)
  Iso_temp[Lum > 1e-6, IMFint := progenyUpdateIMF(Iso_temp[Lum > 1e-6,], IMFfunc, ...)$IMFint]

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
      progenyIso2Spec(logAge_step, logZ_step, Iso=Iso_temp, IMFint=Iso_temp$IMFint, Spec_combine, Interp_combine, interp=interp)
    }
    return(do.call(rbind, output))
  }

  message('Generating evo grids:')

  if(interp==FALSE){
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
  }else{
    #Need to really test this! 08/05/2025
    Zevo = foreach(logZ_step = logZ_steps)%dopar%{
      if(rem_frac == 'get'){
        temp_Z = interp_quick(logZ_step, Iso_logZ_steps)
        temp_cut_Zlo = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_lo']] & logAge >= 8,data.table(Mini,Mass/Mini)[which.max(Mini),],by=logAge]
        temp_cut_Zhi = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_hi']] & logAge >= 8,data.table(Mini,Mass/Mini)[which.max(Mini),],by=logAge]
        temp_func_Zlo = approxfun(temp_cut_Zlo[,Mini], temp_cut_Zlo[,V2], rule=2)*temp_Z['wt_lo']
        temp_func_Zhi = approxfun(temp_cut_Zhi[,Mini], temp_cut_Zhi[,V2], rule=2)*temp_Z['wt_hi']
      }else{
        temp_func_Zlo = function(x){rem_frac}
        temp_func_lo = function(x){rem_frac}
      }

      missing_func = function(lo_lim){
        integrate(function(x){IMFfunc(x, massmult=TRUE, ...)*temp_func(x)}, lower=lo_lim, upper=Inf)$value
      }

      message('  ',logZ_step)
      temp_out = foreach(logAge_step = logAge_steps, .combine='rbind')%do%{
        temp_Age = interp_quick(logAge_step, Iso_logAge_steps)
        mass_lim_Zlo_Alo = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_lo']] & logAge == Iso_logAge_steps[temp_Age['ID_lo']], max(Mini)]
        mass_lim_Zlo_Ahi = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_lo']] & logAge == Iso_logAge_steps[temp_Age['ID_hi']], max(Mini)]
        mass_lim_Zhi_Alo = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_hi']] & logAge == Iso_logAge_steps[temp_Age['ID_lo']], max(Mini)]
        mass_lim_Zhi_Ahi = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_hi']] & logAge == Iso_logAge_steps[temp_Age['ID_hi']], max(Mini)]
        SMrem_miss_Zlo = temp_func_Zlo(mass_lim_Zlo_Alo)*temp_Age['wt_lo'] + temp_func_Zlo(mass_lim_Zlo_Ahi)*temp_Age['wt_hi']
        SMrem_miss_Zhi = temp_func_Zhi(mass_lim_Zhi_Alo)*temp_Age['wt_lo'] + temp_func_Zhi(mass_lim_Zhi_Ahi)*temp_Age['wt_hi']
        SMrem_miss = SMrem_miss_Zlo + SMrem_miss_Zhi

        SMstar_Zlo_Alo = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_lo']] & logAge == Iso_logAge_steps[temp_Age['ID_lo']], sum(Mass*IMFint)] + SMrem_miss_Zlo
        SMstar_Zlo_Ahi = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_lo']] & logAge == Iso_logAge_steps[temp_Age['ID_hi']], sum(Mass*IMFint)] + SMrem_miss_Zlo
        SMstar_Zhi_Alo = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_hi']] & logAge == Iso_logAge_steps[temp_Age['ID_lo']], sum(Mass*IMFint)] + SMrem_miss_Zhi
        SMstar_Zhi_Ahi = Iso_temp[logZ == Iso_logZ_steps[temp_Z['ID_hi']] & logAge == Iso_logAge_steps[temp_Age['ID_hi']], sum(Mass*IMFint)] + SMrem_miss_Zhi
        SMstar = SMstar_Zlo_Alo + SMstar_Zlo_Ahi + SMstar_Zhi_Alo + SMstar_Zhi_Ahi

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
