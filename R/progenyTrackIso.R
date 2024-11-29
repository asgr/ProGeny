progenyTrackInterp = function(tracklist, mass_target, logAge_vec=NULL,
                              iso_type='approx', cores=8, logZ=0, ...){
  #track_masses = as.numeric(names(tracklist))
  track_masses = unique(tracklist$Mini)

  if(is.data.frame(mass_target)){
    logAge_lo = mass_target$logAge_lo
    logAge_hi = mass_target$logAge_hi
    mass_target = mass_target$Mini
  }else{
    logAge_lo = NULL
    logAge_hi = NULL
  }

  track_bins = as.integer(cut(mass_target, breaks=track_masses, right=FALSE))
  track_bins[mass_target == max(track_masses)] = length(track_masses) - 1L

  sel_track_bins = unique(track_bins)

  cores = min(cores, length(weights), detectCores())
  registerDoParallel(cores=cores)

  bin = NULL
  combine_out = foreach(bin = sel_track_bins, .combine=c)%dopar%{
    #mass_loc_start = which(names(tracklist) == mass_start)
    #track_start = tracklist[[bin]]
    track_start = tracklist[Mini == track_masses[bin],]
    mass_start = log10(track_masses[bin])

    #mass_loc_end = which(names(tracklist) == mass_end)
    #track_end = tracklist[[bin + 1L]]
    track_end = tracklist[Mini == track_masses[bin + 1L],]
    mass_end = log10(track_masses[bin + 1L])

    mass_target_loc = log10(mass_target)
    mass_target_loc = mass_target_loc[(mass_target_loc >= mass_start & mass_target_loc < mass_end) | (bin == max(sel_track_bins) & mass_target_loc == max(log10(mass_target)))] #The RHS of the OR is to catch masses at the extreme massive end that exactly match the max input tracklist (these can still be computed)

    weights = (mass_target_loc - mass_start)/(mass_end - mass_start)

    track_start_labels_tab = table(track_start$label)
    track_end_labels_tab = table(track_end$label)

    all_phases = sort(unique(names(track_start_labels_tab), names(track_end_labels_tab)))

    i = NULL
    max_res_phases = foreach(i = all_phases, .combine='c')%do%{
      max(track_start_labels_tab[names(track_start_labels_tab) == i], track_end_labels_tab[names(track_end_labels_tab) == i])
    }

    #here we calculate the best mapping of phases, making sure the all end up having the same length in terms of indexing

    remap_start = {}
    current_max = 0
    for(i in seq_along(all_phases)){
      if(all_phases[i] %in% names(track_start_labels_tab)){
        current_phase_map = seq(1,track_start_labels_tab[names(track_start_labels_tab) == all_phases[i]], length = max_res_phases[i])
      }else{
        current_phase_map = rep(NA, max_res_phases[i])
      }

      remap_start = c(remap_start,current_phase_map + current_max)
      current_max = max(remap_start, na.rm=TRUE)
    }

    remap_end = {}
    current_max = 0
    for(i in seq_along(all_phases)){
      if(all_phases[i] %in% names(track_end_labels_tab)){
        current_phase_map = seq(1,track_end_labels_tab[names(track_end_labels_tab) == all_phases[i]], length = max_res_phases[i])
      }else{
        current_phase_map = rep(NA, max_res_phases[i])
      }

      remap_end = c(remap_end,current_phase_map + current_max)
      current_max = max(remap_end, na.rm=TRUE)
    }

    tag_vec_start = 1:dim(track_start)[1]
    tag_vec_end = 1:dim(track_end)[1]


    #logAge2tag_start = approxfun(track_start$logAge, tag_vec_start, yleft=NA, yright=NA)
    #logAge2tag_end = approxfun(track_end$logAge, tag_vec_end, yleft=NA, yright=NA)

    start_func_logAge = approxfun(tag_vec_start, track_start$logAge, yleft=NA, yright=NA)
    end_func_logAge = approxfun(tag_vec_end, track_end$logAge, yleft=NA, yright=NA)

    start_func_logMass = approxfun(tag_vec_start, log10(track_start$Mass), yleft=NA, yright=NA)
    end_func_logMass = approxfun(tag_vec_end, log10(track_end$Mass), yleft=NA, yright=NA)

    start_func_logLum = approxfun(tag_vec_start, log10(track_start$Lum), yleft=NA, yright=NA)
    end_func_logLum = approxfun(tag_vec_end, log10(track_end$Lum), yleft=NA, yright=NA)

    start_func_logTeff = approxfun(tag_vec_start, log10(track_start$Teff), yleft=NA, yright=NA)
    end_func_logTeff = approxfun(tag_vec_end, log10(track_end$Teff), yleft=NA, yright=NA)

    start_func_logG = approxfun(tag_vec_start, track_start$logG, yleft=NA, yright=NA)
    end_func_logG = approxfun(tag_vec_end, track_end$logG, yleft=NA, yright=NA)

    start_func_label = approxfun(tag_vec_start, track_start$label, yleft=NA, yright=NA)
    end_func_label = approxfun(tag_vec_end, track_end$label, yleft=NA, yright=NA)

    #Still working on this: map ages to mean tag, then mean tag to properties.

    tag_vec_out = seq_along(remap_start)

    i = NULL
    temp_loop = foreach(i = seq_along(mass_target_loc))%do%{
      weight = weights[i]

      temp_mat_logAge = cbind(start_func_logAge(remap_start[tag_vec_out])*(1 - weight),
                              end_func_logAge(remap_end[tag_vec_out])*weight)
      temp_vec_logAge = rowSums(temp_mat_logAge, na.rm = FALSE)

      temp_mat_logLum = cbind(start_func_logLum(remap_start[tag_vec_out])*(1 - weight),
                              end_func_logLum(remap_end[tag_vec_out])*weight)
      temp_vec_logLum = rowSums(temp_mat_logLum, na.rm = FALSE)

      temp_mat_logMass = cbind(start_func_logMass(remap_start[tag_vec_out])*(1 - weight),
                               end_func_logMass(remap_end[tag_vec_out])*weight)
      temp_vec_logMass = rowSums(temp_mat_logMass, na.rm = FALSE)

      temp_mat_logTeff = cbind(start_func_logTeff(remap_start[tag_vec_out])*(1 - weight),
                               end_func_logTeff(remap_end[tag_vec_out])*weight)
      temp_vec_logTeff = rowSums(temp_mat_logTeff, na.rm = FALSE)

      temp_mat_logG = cbind(start_func_logG(remap_start[tag_vec_out])*(1 - weight),
                            end_func_logG(remap_end[tag_vec_out])*weight)
      temp_vec_logG = rowSums(temp_mat_logG, na.rm = FALSE)

      temp_mat_label = cbind(start_func_label(remap_start[tag_vec_out])*(1 - weight),
                             end_func_label(remap_end[tag_vec_out])*weight)
      temp_vec_label = rowSums(temp_mat_label, na.rm = FALSE)

      if(!is.null(logAge_vec)){
        if(iso_type == 'approx'){
          temp_vec_logMass = approxfun(temp_vec_logAge, temp_vec_logMass, ties = "ordered", ...)(logAge_vec)
          temp_vec_logLum = approxfun(temp_vec_logAge, temp_vec_logLum, ties = "ordered", ...)(logAge_vec)
          temp_vec_logTeff = approxfun(temp_vec_logAge, temp_vec_logTeff, ties = "ordered", ...)(logAge_vec)
          temp_vec_logG = approxfun(temp_vec_logAge, temp_vec_logG, ties = "ordered", ...)(logAge_vec)
          temp_vec_label = approxfun(temp_vec_logAge, temp_vec_label, ties = "ordered", ...)(logAge_vec)
        }else if(iso_type == 'integral' | iso_type == 'integrate'){
          temp_vec_logMass = ProSpect::specReBin(temp_vec_logAge, temp_vec_logMass, logAge_vec, ...)$flux
          temp_vec_logLum = ProSpect::specReBin(temp_vec_logAge, temp_vec_logLum, logAge_vec, ...)$flux
          temp_vec_logTeff = ProSpect::specReBin(temp_vec_logAge, temp_vec_logTeff, logAge_vec, ...)$flux
          temp_vec_logG = ProSpect::specReBin(temp_vec_logAge, temp_vec_logG, logAge_vec, ...)$flux
          temp_vec_label = ProSpect::specReBin(temp_vec_logAge, temp_vec_label, logAge_vec, ...)$flux
        }

        temp_vec_logAge = logAge_vec
      }

      temp_out = data.table(logAge = temp_vec_logAge,
                       Mini = mass_target[i],
                       Mass = 10^temp_vec_logMass,
                       Lum = 10^temp_vec_logLum,
                       Teff = 10^temp_vec_logTeff,
                       logG = temp_vec_logG,
                       label = ceiling(temp_vec_label)
      )
      return(temp_out)
    }
    return(temp_loop)
  }

  if(length(combine_out) == 1){
    return(combine_out[[1]])
  }else{
    return(rbindlist(combine_out))
  }
}

#This will extract the locations of important features
#Next step is then to log interpolate to find the implied mass they would occur at given a logAge_vec
progenyLabelFeatures = function(tracklist, logAge_vec, feature_sample=seq(-0.01,0.01,len=5)){
  label_vec = unique(tracklist$label)

  i = NULL
  feature_mass = foreach(i = label_vec, .combine=c)%do%{
    temp = tracklist[label == i, list(logAge=max(logAge)), by=Mini]
    temp_func = approxfun(temp[,list(logAge, log10(Mini))], yleft=NA, yright=NA)
    temp_mass = 10^temp_func(logAge_vec)
    return(temp_mass[!is.na(temp_mass)])
  }

  feature_mass = sort(c(unique(tracklist$Mini), as.numeric(outer(feature_mass, feature_sample, '+'))))

  return(feature_mass)
}

progenyGradFeatures = function(tracklist, logAge_vec, feature_sample=seq(-0.01,0.01,len=5),
                               deriv=1, n=1e4, keep_orig_mass=TRUE){
  mass_vec = unique(tracklist$Mini)

  i = NULL
  Mini = NULL
  feature_out = foreach(i = 1:(length(mass_vec) - 1L))%do%{
  #for(i in 1:(length(mass_vec) - 1L)){
    temp_func = splinefun(tracklist[Mini==mass_vec[i], list(10^logAge, log10(Lum))])
    temp_logAge = log10(rootSolve::uniroot.all(function(x){temp_func(x,deriv)}, c(0,10^max(tracklist[Mini==mass_vec[i],logAge])), n=n))
    age_bins = as.integer(cut(temp_logAge, breaks=logAge_vec, right=FALSE))
    age_bins = age_bins[!is.na(age_bins)]
    sel_age_bins = unique(age_bins)

    #Here we both log interpole track upwards and downwards. The logic is there will be contributing features from both less and more massive stars at an interpolated value, e.g. a 13 Msol star will be (roughly) 50% a 10 Msol and 50% 20 Msol star, so need to map features from both directions to not miss important things.
    temp_mass_lo = {}
    temp_mass_hi = {}
    for(bin in sel_age_bins){
      #For fraction above the current mass track the next more massive track
      temp_logAge_loc = temp_logAge[(temp_logAge >= logAge_vec[bin] & temp_logAge < logAge_vec[bin + 1L])]
      temp_logAge_frac = (temp_logAge_loc - logAge_vec[bin])/(logAge_vec[bin + 1L] - logAge_vec[bin])
      temp_mass_lo = c(temp_mass_lo, mass_vec[i] + temp_logAge_frac*(mass_vec[i + 1L] - mass_vec[i]))

      if(bin > sel_age_bins[1]){
        #For fraction below the current mass track to next less massive track
        temp_logAge_loc = temp_logAge[(temp_logAge >= logAge_vec[bin - 1L] & temp_logAge < logAge_vec[bin]) | (bin == max(sel_age_bins) & temp_logAge == max(logAge_vec))]
        temp_logAge_frac = (logAge_vec[bin] - temp_logAge_loc)/(logAge_vec[bin] - logAge_vec[bin - 1L])
        temp_mass_hi = c(temp_mass_hi, mass_vec[i + 1L] - temp_logAge_frac*(mass_vec[i + 1L] - mass_vec[i]))
      }
    }

    temp_mass_lo = as.numeric(outer(temp_mass_lo, feature_sample, '+'))
    temp_mass_hi = as.numeric(outer(temp_mass_hi, feature_sample, '+'))

    temp_DF_lo = data.table(Mini=temp_mass_lo, logAge_lo=logAge_vec[bin], logAge_hi=logAge_vec[bin + 1L])
    temp_DF_hi = data.table(Mini=temp_mass_hi, logAge_lo=logAge_vec[bin - 1L], logAge_hi=logAge_vec[bin])

    return(rbind(temp_DF_lo, temp_DF_hi))
  }

  feature_out = rbindlist(feature_out)
  #this seems to be doing something sensible now...

  feature_out = feature_out[Mini >= min(mass_vec) & Mini <= max(mass_vec),]

  if(keep_orig_mass){
    logAge_lo = logAge_vec[1:(length(logAge_vec) - 1L)]
    logAge_hi = logAge_vec[2:length(logAge_vec)]
    Mini = rep(mass_vec, each=length(logAge_lo))
    feature_out_add = data.table(Mini=Mini, logAge_lo=logAge_lo, logAge_hi=logAge_hi)
    feature_out = rbind(feature_out, feature_out_add)
  }

  setkey(feature_out, logAge_lo, Mini)

  return(feature_out)
}
