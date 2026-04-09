progenyTrackInterp = function(tracklist, target,
                              make_iso = TRUE, logAge_lim = c(5,10.3), logAge_bin=0.05,
                              iso_type='approx', cores=8, logZ_use=0, ...){
  tracklist = copy(tracklist[logZ == logZ_use,])
  target_type = 'Mini' #no longer thinking of Z interp, at least not for now.
  track_vals = unique(tracklist[[target_type]])
  logAge_vec = seq(logAge_lim[1], logAge_lim[2], by=logAge_bin)

  if(iso_type == 'integral' | iso_type == 'integrate'){
    if(!requireNamespace("ProSpect", quietly = TRUE)){
      stop('The ProSpect package is needed for this function to work with selected iso_type. Please install it from GitHub asgr/ProSpect', call. = FALSE)
    }
  }

  if(target_type == 'Mini'){
    target_log = TRUE
  }else if(target_type == 'logZ'){
    target_log = FALSE
  }else{
    stop('target type must be one of Mini or logZ')
  }

  if(is.data.frame(target)){
    target_orig = copy(target)
    target = unique(unlist(target[,1]))
  }else{
    target_orig = NULL
  }

  track_bins = findInterval(x=target, vec=track_vals, all.inside=TRUE)
  track_bins[target == max(track_vals)] = length(track_vals) - 1L

  sel_track_bins = sort(unique(track_bins))

  cores = min(cores, length(sel_track_bins), parallel::detectCores())

  registerDoParallel(cores=cores)

  bin = NULL
  combine_out = foreach(bin = sel_track_bins, .combine=c)%dopar%{
    track_start = tracklist[get(target_type) == track_vals[bin],]
    track_end = tracklist[get(target_type) == track_vals[bin + 1L],]

    target_start = track_vals[bin]
    target_end = track_vals[bin + 1L]
    target_loc = target

    #The RHS of the OR is to catch masses at the extreme massive end that exactly match the max input tracklist (these can still be computed)
    sel = (target_loc >= target_start & target_loc < target_end) | (bin == max(sel_track_bins) & target_loc == max(track_vals))
    target_loc = target_loc[sel]

    if(target_log){
      target_start = log10(target_start)
      target_end = log10(target_end)
      target_loc = log10(target_loc)
    }

    weights = (target_loc - target_start)/(target_end - target_start)

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

      remap_start = c(remap_start, current_phase_map + current_max)
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


    #logAge2tag_start = approxfun(track_start$logAge, tag_vec_start, yleft=NA, yright=NA, na.rm=TRUE)
    #logAge2tag_end = approxfun(track_end$logAge, tag_vec_end, yleft=NA, yright=NA, na.rm=TRUE)

    start_func_logAge = approxfun(tag_vec_start, track_start$logAge, yleft=NA, yright=NA, na.rm=TRUE)
    end_func_logAge = approxfun(tag_vec_end, track_end$logAge, yleft=NA, yright=NA, na.rm=TRUE)

    start_func_logMini = approxfun(tag_vec_start, log10(track_start$Mini), yleft=NA, yright=NA, na.rm=TRUE)
    end_func_logMini = approxfun(tag_vec_end, log10(track_end$Mini), yleft=NA, yright=NA, na.rm=TRUE)

    start_func_logMass = approxfun(tag_vec_start, log10(track_start$Mass), yleft=NA, yright=NA, na.rm=TRUE)
    end_func_logMass = approxfun(tag_vec_end, log10(track_end$Mass), yleft=NA, yright=NA, na.rm=TRUE)

    start_func_logLum = approxfun(tag_vec_start, log10(track_start$Lum), yleft=NA, yright=NA, na.rm=TRUE)
    end_func_logLum = approxfun(tag_vec_end, log10(track_end$Lum), yleft=NA, yright=NA, na.rm=TRUE)

    start_func_logTeff = approxfun(tag_vec_start, log10(track_start$Teff), yleft=NA, yright=NA, na.rm=TRUE)
    end_func_logTeff = approxfun(tag_vec_end, log10(track_end$Teff), yleft=NA, yright=NA, na.rm=TRUE)

    start_func_logG = approxfun(tag_vec_start, track_start$logG, yleft=NA, yright=NA, na.rm=TRUE)
    end_func_logG = approxfun(tag_vec_end, track_end$logG, yleft=NA, yright=NA, na.rm=TRUE)

    start_func_label = approxfun(tag_vec_start, track_start$label, yleft=NA, yright=NA, na.rm=TRUE)
    end_func_label = approxfun(tag_vec_end, track_end$label, yleft=NA, yright=NA, na.rm=TRUE)

    #Still working on this: map ages to mean tag, then mean tag to properties.

    tag_vec_out = seq_along(remap_start)

    # Pre-compute approxfun evaluations once, outside the per-mass loop
    start_logAge_vals  = start_func_logAge(remap_start[tag_vec_out])
    end_logAge_vals    = end_func_logAge(remap_end[tag_vec_out])
    start_logMini_vals = start_func_logMini(remap_start[tag_vec_out])
    end_logMini_vals   = end_func_logMini(remap_end[tag_vec_out])
    start_logMass_vals = start_func_logMass(remap_start[tag_vec_out])
    end_logMass_vals   = end_func_logMass(remap_end[tag_vec_out])
    start_logLum_vals  = start_func_logLum(remap_start[tag_vec_out])
    end_logLum_vals    = end_func_logLum(remap_end[tag_vec_out])
    start_logTeff_vals = start_func_logTeff(remap_start[tag_vec_out])
    end_logTeff_vals   = end_func_logTeff(remap_end[tag_vec_out])
    start_logG_vals    = start_func_logG(remap_start[tag_vec_out])
    end_logG_vals      = end_func_logG(remap_end[tag_vec_out])
    start_label_vals   = start_func_label(remap_start[tag_vec_out])
    end_label_vals     = end_func_label(remap_end[tag_vec_out])

    i = NULL
    temp_loop = foreach(i = seq_along(target_loc))%do%{
      weight = weights[i]

      temp_vec_logAge  = start_logAge_vals  * (1 - weight) + end_logAge_vals  * weight
      temp_vec_logMini = start_logMini_vals * (1 - weight) + end_logMini_vals * weight
      temp_vec_logMass = start_logMass_vals * (1 - weight) + end_logMass_vals * weight
      temp_vec_logLum  = start_logLum_vals  * (1 - weight) + end_logLum_vals  * weight
      temp_vec_logTeff = start_logTeff_vals * (1 - weight) + end_logTeff_vals * weight
      temp_vec_logG    = start_logG_vals    * (1 - weight) + end_logG_vals    * weight
      temp_vec_label   = start_label_vals   * (1 - weight) + end_label_vals   * weight

      if(make_iso){
        if(iso_type == 'approx'){
          temp_vec_logMini = approx(temp_vec_logAge, temp_vec_logMini, xout=logAge_vec, ties = "ordered", na.rm=TRUE, ...)$y
          temp_vec_logMass = approx(temp_vec_logAge, temp_vec_logMass, xout=logAge_vec, ties = "ordered", na.rm=TRUE, ...)$y
          temp_vec_logLum  = approx(temp_vec_logAge, temp_vec_logLum,  xout=logAge_vec, ties = "ordered", na.rm=TRUE, ...)$y
          temp_vec_logTeff = approx(temp_vec_logAge, temp_vec_logTeff, xout=logAge_vec, ties = "ordered", na.rm=TRUE, ...)$y
          temp_vec_logG    = approx(temp_vec_logAge, temp_vec_logG,    xout=logAge_vec, ties = "ordered", na.rm=TRUE, ...)$y
          temp_vec_label   = approx(temp_vec_logAge, temp_vec_label,   xout=logAge_vec, ties = "ordered", na.rm=TRUE, ...)$y
        }else if(iso_type == 'integral' | iso_type == 'integrate'){
          temp_vec_Age = 10^temp_vec_logAge
          Age_vec = 10^logAge_vec
          temp_vec_logMass = log10(ProSpect::specReBin(temp_vec_Age, 10^temp_vec_logMass, Age_vec, ...)$flux)
          temp_vec_logMini = log10(ProSpect::specReBin(temp_vec_Age, 10^temp_vec_logMini, Age_vec, ...)$flux)
          temp_vec_logLum = log10(ProSpect::specReBin(temp_vec_Age, 10^temp_vec_logLum, Age_vec, ...)$flux)
          temp_vec_logTeff = log10(ProSpect::specReBin(temp_vec_Age, 10^temp_vec_logTeff, Age_vec, ...)$flux)
          temp_vec_logG = log10(ProSpect::specReBin(temp_vec_Age, 10^temp_vec_logG, Age_vec, ...)$flux)
          temp_vec_label = ProSpect::specReBin(temp_vec_Age, temp_vec_label, Age_vec, ...)$flux
        }

        temp_vec_logAge = logAge_vec
      }

      if(target_type == 'Mini'){
        temp_out = data.table(
          logZ = unique(tracklist$logZ),
          logAge = temp_vec_logAge,
          Mini = target[sel][i],
          Mass = 10^temp_vec_logMass,
          Lum = 10^temp_vec_logLum,
          Teff = 10^temp_vec_logTeff,
          logG = temp_vec_logG,
          label = round(temp_vec_label)
        )
      }else if(target_type == 'logZ'){
        temp_out = data.table(
          logZ = target[sel][i],
          logAge = temp_vec_logAge,
          Mini = 10^temp_vec_logMini,
          Mass = 10^temp_vec_logMass,
          Lum = 10^temp_vec_logLum,
          Teff = 10^temp_vec_logTeff,
          logG = temp_vec_logG,
          label = round(temp_vec_label)
        )
      }

      return(temp_out)
    }
    return(temp_loop)
  }

  if(length(combine_out) == 1){
    combine_out = combine_out[[1]]
  }else{
    combine_out = rbindlist(combine_out)
  }

  setkey(combine_out, logZ, logAge, Mini)

  logAge = logAge_lo = logAge_hi = NULL

  if(!is.null(target_orig)){
    use_age = unique(target_orig[,list(logAge_lo, logAge_hi)])

    combine_out = foreach(i = 1:dim(use_age)[1])%do%{
      logAge_temp = (use_age[i,logAge_lo] + use_age[i,logAge_hi])/2
      combine_temp = combine_out[(get(target_type) %in% target_orig[logAge_lo == logAge_temp - logAge_bin/2, get(target_type)]) & (logAge == logAge_temp),]
      return(combine_temp)
    }

    combine_out = rbindlist(combine_out)
  }

  Mini = Mass = Lum = Teff = logG = NULL

  combine_out = combine_out[is.finite(Mini) &
                    is.finite(Mass) &
                    is.finite(Lum) &
                    is.finite(Teff) &
                    is.finite(logG)
                  , ]

  return(combine_out)
}

progenyFindMass = function(tracklist, logAge_lim = c(5,10.3), logAge_bin=0.05, logZ_use=0, label_use=-1:9){

  if(!requireNamespace("akima", quietly = TRUE)){
    stop('The akima package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }

  label = label_val = logAge = Mini = eep_i = NULL

  logAge_vec = seq(logAge_lim[1], logAge_lim[2], by=logAge_bin)

  temp = tracklist[logZ==logZ_use, list(logAge, eep_i=1:.N), by=list(Mini, label)]

  #to compute the interpolated linear fit onto a regular isochrone grid

  #need to loop round all labels, because some phases have different numbers of EEPs, and we want to make sure we are interpolating within the same phase (e.g. MS to MS, not MS to RGB)
  label_use = label_use[label_use %in% unique(temp$label)]

  output_all = foreach(label_val = label_use)%do%{
    temp_sub = temp[label==label_val,]
    akima.si = akima::interp(temp_sub[,logAge], temp_sub[, eep_i], temp_sub[, log10(Mini)],
                             xo=logAge_vec, yo=min(temp_sub[, eep_i]):max(temp_sub[, eep_i]),
                             linear = TRUE, extrap = FALSE, duplicate = 'median')
    #clean the top part for bits that go out of domain
    clean = ProPane::propaneBin2D(temp_sub[,logAge], temp_sub[, eep_i], image=matrix(0,length(akima.si$x), length(akima.si$y)),
                                  xlim=range(akima.si$x), ylim=range(akima.si$y))
    ymax = apply(clean$z != 0, 1, function(x) {if (any(x)) max(which(x)) else NA_integer_})
    if(length(which(!is.na(ymax))) > 1){
      temp_func = approxfun(akima.si$x, ymax)
      temp_grid = expand.grid(akima.si$x, akima.si$y)
      akima.si$z[temp_grid[,2] > (temp_func(temp_grid[,1]) + 0.5)] = NA
    }

    #make sure all mass values are in domain
    lims = range(log10(temp_sub[, Mini]), na.rm = TRUE)
    akima.si$z[] = pmax(pmin(akima.si$z, lims[2]), lims[1])

    i = NULL

    output = {
      z_mat = akima.si$z
      row_idx_all = rep(seq_len(nrow(z_mat)), times = ncol(z_mat))
      z_vec = as.vector(z_mat)
      finite_sel = is.finite(z_vec)
      dt = data.table(row_idx = row_idx_all[finite_sel], Mini = 10^z_vec[finite_sel])
      dt = unique(dt)
      dt[, `:=`(
        logAge_lo = akima.si$x[row_idx] - logAge_bin/2,
        logAge_hi = akima.si$x[row_idx] + logAge_bin/2
      )]
      dt[, row_idx := NULL]
      dt
    }

    return(output)
  }
  target = rbindlist(output_all)
  setkey(target, logAge_lo, Mini)
  return(target)
}

progenyExtendIso = function(Iso_base, Iso_extend, label=NA, logA=NA){
  label_use = label
  logA_use = logA

  setDT(Iso_base)
  setDT(Iso_extend)

  Mini = logZ = logAge = NULL

  Mini_end = Iso_base[,max(Mini),by=list(logZ, logAge)]

  if(identical(label_use, 'base_get')){
    label_use = Iso_base[,max(label, na.rm=TRUE)] + 1
  }

  if(identical(logA_use, 'base_get')){
    logA_use = Iso_base[,max(logA, na.rm=TRUE)] + 1
  }

  i = V1 = NULL

  output = foreach(i = 1:dim(Mini_end)[1])%do%{
    best_logZ = Iso_extend[which.min(abs(logZ - Mini_end[i,logZ])), logZ]
    best_logAge = Iso_extend[logZ == best_logZ & Mini >= Mini_end[i,V1], max(logAge)]
    temp = Iso_extend[logZ == best_logZ & logAge == best_logAge & Mini >= Mini_end[i,V1],]
    temp[,logZ:= Mini_end[i,logZ]]
    temp[,logAge:= Mini_end[i,logAge]]

    if(!identical(label_use, 'extend_get')){
        temp[,label:=label_use]
    }

    return(temp)
  }

  output = rbindlist(output)
  output = rbind(Iso_base, output, fill=TRUE)

  if('logA' %in% colnames(output) & !is.na(logA_use)){
    if(!identical(logA_use, 'extend_get')){
      output[is.na(logA), logA:=logA_use]
    }
  }

  setkey(output, logZ, logAge, Mini)
  return(output)
}

# .progenyTagFeatures = function(tracklist, deriv=1, n=1e4, expand=-1:1){
#   tracklist_temp = copy(tracklist)
#
#   if(!is.null(tracklist_temp$logZ)){
#     logZ_vals = unique(tracklist_temp$logZ)
#   }else{
#     tracklist_temp$logZ = -99
#     logZ_vals = -99
#   }
#
#   tracklist_temp[,tag:=FALSE]
#
#   logZ_use = NULL
#   Mini = NULL
#
#   for(logZ_use in logZ_vals){
#     i = NULL
#     target = NULL
#     track_vals = unique(tracklist_temp[logZ == logZ_use,Mini])
#     for(i in seq_along(track_vals)){
#       suppressWarnings({
#         sel = tracklist_temp[Mini == track_vals[i] & logZ == logZ_use, which=TRUE]
#         interval = c(0, 10^max(tracklist_temp[Mini == track_vals[i],logAge]))
#
#         temp_func = splinefun(tracklist_temp[sel, list(10^logAge, Mass)])
#         Mass_root = rootSolve::uniroot.all(function(x){temp_func(x,deriv)}, interval, n=n)
#
#         temp_func = splinefun(tracklist_temp[sel, list(10^logAge, Lum)])
#         Lum_root = rootSolve::uniroot.all(function(x){temp_func(x,deriv)}, interval, n=n)
#
#         temp_func = splinefun(tracklist_temp[sel, list(10^logAge, Teff)])
#         Teff_root = rootSolve::uniroot.all(function(x){temp_func(x,deriv)}, interval, n=n)
#
#         temp_func = splinefun(tracklist_temp[sel, list(10^logAge, logG)])
#         logG_root = rootSolve::uniroot.all(function(x){temp_func(x,deriv)}, interval, n=n)
#
#         all_root = log10(sort(unique(c(Mass_root, Lum_root, Teff_root, logG_root))))
#
#         tag_temp = findInterval(all_root, tracklist_temp[sel, logAge], all.inside = T)
#         tag_temp = c(1, tag_temp, length(sel))
#         tag_temp = unique(as.integer(outer(tag_temp, expand, FUN='+')))
#         tracklist_temp[sel, tag:= (1:length(sel) %in% tag_temp)]
#       })
#     }
#   }
#   return(tracklist_temp)
# }
#
# .progenyTagWeights = function(tracklist, logAge_vec, logZ_vec=NULL){
#   tracklist_temp = copy(tracklist[tag == TRUE])
#
#   age_bin = findInterval(tracklist_temp$logAge, logAge_vec, all.inside=TRUE)
#   tracklist_temp = tracklist_temp[which(age_bin > 0),]
#   age_bin = age_bin[age_bin > 0]
#   age_bin_unique = sort(unique(age_bin))
#
#   if(!is.null(logZ_vec)){
#     Z_bin = findInterval(tracklist_temp$logZ, logZ_vec, all.inside=TRUE)
#     tracklist_temp = tracklist_temp[which(Z_bin > 0),]
#     Z_bin = Z_bin[Z_bin > 0]
#     Z_bin_unique = sort(unique(Z_bin))
#   }
#
#   Mini_vec = sort(unique(tracklist$Mini))
#   Mini_match = match(tracklist_temp$Mini, Mini_vec)
#   Mini_match_lo = Mini_match
#   Mini_match_lo[Mini_match_lo == 1L] = 2L
#   Mini_match_hi = Mini_match
#   Mini_match_hi[Mini_match_hi == length(Mini_vec)] = length(Mini_vec) - 1L
#
#   tracklist_temp[,Mini_prev:=Mini_vec[Mini_match_lo - 1L]]
#   tracklist_temp[,Mini_next:=Mini_vec[Mini_match_hi + 1L]]
#
#   output = {}
#
#   for(i in age_bin_unique){
#     if(is.null(logZ_vec)){
#       #interpolating over just logAge
#       temp = tracklist_temp[which(age_bin == i),list(Mini, logAge, Mini_prev, Mini_next)]
#       logAge_weight = (temp$logAge - logAge_vec[i]) / (logAge_vec[i + 1L] - logAge_vec[i])
#       logAge_weight[logAge_weight < 0] = 0
#       logAge_weight[logAge_weight > 1] = 1
#       temp_Mini = 10^unlist(temp[,log10(Mini) + (log10(Mini_next) - log10(Mini))*logAge_weight])
#       output = rbind(output, data.table(Mini = temp_Mini,
#                                         Mini_lo = temp$Mini,
#                                         Mini_hi = temp$Mini_next,
#                                         logAge = temp$logAge,
#                                         logAge_weight = logAge_weight,
#                                         logAge_lo = logAge_vec[i],
#                                         logAge_hi = logAge_vec[i + 1L])
#                      )
#     }else{
#       for(j in Z_bin_unique){
#         #if we need to also interpolate over logZ
#         #need to check this weight a bit more carefully...
#         temp = tracklist_temp[which(age_bin == i & Z_bin == j),list(Mini, logAge, Mini_prev, Mini_next)]
#         logAge_weight = (temp$logAge - logAge_vec[i]) / (logAge_vec[i + 1L] - logAge_vec[i])
#         logAge_weight[logAge_weight < 0] = 0
#         logAge_weight[logAge_weight > 1] = 1
#         logZ_weight = (temp$logZ - logZ_vec[j]) / (logZ_vec[j + 1L] - logZ_vec[j])
#         logZ_weight[logZ_weight < 0] = 0
#         logZ_weight[logZ_weight > 1] = 1
#         temp_Mini = 10^unlist(temp[,log10(Mini) + (log10(Mini_next) - log10(Mini))*logAge_weight])
#         output = rbind(output, data.table(Mini = temp_Mini,
#                                           Mini_lo = temp$Mini,
#                                           Mini_hi = temp$Mini_next,
#                                           logAge = temp$logAge,
#                                           logAge_weight = logAge_weight,
#                                           logAge_lo = logAge_vec[i],
#                                           logAge_hi = logAge_vec[i + 1L],
#                                           logZ = temp$logZ,
#                                           logZ_weight = logZ_weight,
#                                           logZ_lo = logZ_vec[j],
#                                           logZ_hi = logZ_vec[j + 1L]
#                                           )
#         )
#       }
#     }
#   }
#   return(output)
# }
#
# .progenyGradFeatures = function(tracklist, logAge_vec,
#                                feature_sample=seq(-0.01,0.01,len=5),
#                                deriv=1, n=1e4, keep_orig_mass=TRUE, cores=8){
#   tracklist_temp = copy(tracklist)
#   target_type = 'Mini'
#   target_log = TRUE
#
#   if(!is.null(tracklist_temp$logZ)){
#     logZ_vals = unique(tracklist_temp$logZ)
#   }else{
#     tracklist_temp$logZ = -99
#     logZ_vals = -99
#   }
#
#   cores = min(cores, length(logZ_vals), detectCores())
#   registerDoParallel(cores=cores)
#
#   logZ_use = NULL
#   Mini = NULL
#   final_out = foreach(logZ_use = logZ_vals)%dopar%{
#     i = NULL
#     target = NULL
#     track_vals = unique(tracklist_temp[logZ == logZ_use,Mini])
#     feature_out = foreach(i = 1:(length(track_vals) - 1L))%do%{
#     #for(i in 1:(length(track_vals) - 1L)){
#       temp_func = splinefun(tracklist_temp[Mini == track_vals[i] & logZ == logZ_use, list(10^logAge, log10(Lum))])
#       temp_logAge = log10(rootSolve::uniroot.all(function(x){temp_func(x,deriv)}, c(0,10^max(tracklist_temp[Mini == track_vals[i],logAge])), n=n))
#       age_bins = as.integer(cut(temp_logAge, breaks=logAge_vec, right=FALSE))
#       age_bins = age_bins[!is.na(age_bins)]
#       sel_age_bins = unique(age_bins)
#
#       #Here we both log interpolate track upwards and downwards. The logic is there will be contributing features from both less and more massive stars at an interpolated value, e.g. a 13 Msol star will be (roughly) 50% a 10 Msol and 50% 20 Msol star, so need to map features from both directions to not miss important things.
#       temp_target_lo = {}
#       temp_target_hi = {}
#       for(bin in sel_age_bins){
#         #For fraction above the current mass track the next more massive track
#         temp_logAge_loc = temp_logAge[(temp_logAge >= logAge_vec[bin] & temp_logAge < logAge_vec[bin + 1L])]
#         temp_logAge_frac = (temp_logAge_loc - logAge_vec[bin])/(logAge_vec[bin + 1L] - logAge_vec[bin])
#         if(target_log){
#           temp_target_lo = c(temp_target_lo, 10^(log10(track_vals[i]) + temp_logAge_frac*(log10(track_vals[i + 1L]) - log10(track_vals[i]))))
#         }else{
#           temp_target_lo = c(temp_target_lo, track_vals[i] + temp_logAge_frac*(track_vals[i + 1L] - track_vals[i]))
#         }
#
#         if(bin > sel_age_bins[1]){
#           #For fraction below the current mass track to next less massive track
#           temp_logAge_loc = temp_logAge[(temp_logAge >= logAge_vec[bin - 1L] & temp_logAge < logAge_vec[bin]) | (bin == max(sel_age_bins) & temp_logAge == max(logAge_vec))]
#           temp_logAge_frac = (logAge_vec[bin] - temp_logAge_loc)/(logAge_vec[bin] - logAge_vec[bin - 1L])
#           if(target_log){
#             temp_target_hi = c(temp_target_hi, 10^(log10(track_vals[i + 1L]) - temp_logAge_frac*(log10(track_vals[i + 1L]) - log10(track_vals[i]))))
#           }else{
#             temp_target_hi = c(temp_target_hi, track_vals[i + 1L] - temp_logAge_frac*(track_vals[i + 1L] - track_vals[i]))
#           }
#         }
#       }
#
#       temp_target_lo = as.numeric(outer(temp_target_lo, feature_sample, '+'))
#       temp_target_hi = as.numeric(outer(temp_target_hi, feature_sample, '+'))
#
#       temp_DF_lo = data.table(Mini=temp_target_lo, logAge_lo=logAge_vec[bin], logAge_hi=logAge_vec[bin + 1L])
#       temp_DF_hi = data.table(Mini=temp_target_hi, logAge_lo=logAge_vec[bin - 1L], logAge_hi=logAge_vec[bin])
#
#       return(rbind(temp_DF_lo, temp_DF_hi))
#     }
#
#     feature_out = rbindlist(feature_out)
#     #this seems to be doing something sensible now...
#
#     feature_out = feature_out[Mini >= min(track_vals) & Mini <= max(track_vals),]
#
#     if(keep_orig_mass){
#       logAge_lo = logAge_vec[1:(length(logAge_vec) - 1L)]
#       logAge_hi = logAge_vec[2:length(logAge_vec)]
#       target = rep(track_vals, each=length(logAge_lo))
#       feature_out_add = data.table(Mini=target, logAge_lo=logAge_lo, logAge_hi=logAge_hi)
#       feature_out = rbind(feature_out, feature_out_add)
#     }
#     return(feature_out)
#   }
#
#   final_out = rbindlist(final_out)
#   final_out= unique(final_out)
#   setkey(final_out, logAge_lo, Mini)
#   return(final_out)
# }
#
# #This will extract the locations of important features
# #Next step is then to log interpolate to find the implied mass they would occur at given a logAge_vec
# .progenyLabelFeatures = function(tracklist, logAge_vec, feature_sample=seq(-0.01,0.01,len=5)){
#   label_vec = unique(tracklist$label)
#
#   i = NULL
#   feature_mass = foreach(i = label_vec, .combine=c)%do%{
#     temp = tracklist[label == i, list(logAge=max(logAge)), by=Mini]
#     temp_func = approxfun(temp[,list(logAge, log10(Mini))], yleft=NA, yright=NA)
#     temp_mass = 10^temp_func(logAge_vec)
#     return(temp_mass[!is.na(temp_mass)])
#   }
#
#   feature_mass = sort(c(unique(tracklist$Mini), as.numeric(outer(feature_mass, feature_sample, '+'))))
#
#   return(feature_mass)
# }
