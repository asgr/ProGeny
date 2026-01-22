progenyInterpGrid = function(loc, info, rescale, radius=2, weight_pow=2, k=8){
  setDT(loc)
  setDT(info)

  col_keep = Teff = logG = logZ = NULL

  if(rescale$Teff != -999 & 'Teff' %in% colnames(loc)){
    has_Teff = TRUE
    col_keep = c(col_keep, 'Teff')
  }else{
    has_Teff = FALSE
  }

  if(rescale$logG != -999 & 'logG' %in% colnames(loc)){
    has_logG = TRUE
    col_keep = c(col_keep, 'logG')
  }else{
    has_logG = FALSE
  }

  if(rescale$logZ != -999 & 'logZ' %in% colnames(loc)){
    has_logZ = TRUE
    col_keep = c(col_keep, 'logZ')
  }else{
    has_logZ = FALSE
  }

  if(!is.null(rescale$logA)){
    if(rescale$logA != -999 & 'logA' %in% colnames(loc)){
      has_logA = TRUE
      col_keep = c(col_keep, 'logA')
    }else{
      has_logA = FALSE
    }
  }else{
    has_logA = FALSE
  }

  loc_local = copy(loc[,..col_keep])
  info_local = copy(info[,..col_keep])

  if(has_Teff){
    loc_local[,Teff:=Teff/rescale$Teff]
    info_local[,Teff:=Teff/rescale$Teff]
  }

  if(has_logG){
    loc_local[,logG:=logG/rescale$logG]
    info_local[,logG:=logG/rescale$logG]
  }

  if(has_logZ){
    loc_local[,logZ:=logZ/rescale$logZ]
    info_local[,logZ:=logZ/rescale$logZ]
  }

  if(has_logA){
    loc_local[,logA:=logA/rescale$logA]
    info_local[,logA:=logA/rescale$logA]
  }

  temp = nn2(info_local, loc_local, k=k, radius=radius, searchtype = 'radius')
  temp$weights = 1/(temp$nn.dists^weight_pow)
  temp$weights = temp$weights / rowSums(temp$weights)
  temp$weights[is.na(temp$weights)] = 1
  temp$nn.dists[temp$nn.idx == 0L] = Inf
  temp$weights[temp$nn.idx == 0L] = 0
  return(temp)
}

progenyInterpGrid_All = function(Iso, Spec_combine, radius=2, weight_pow=2, k=8){
  setDT(Iso)

  dobase = TRUE #always need base
  doextend = !is.null(Spec_combine$extend)
  dohot = !is.null(Spec_combine$hot)
  doAGB = !is.null(Spec_combine$AGB)
  dowhite = !is.null(Spec_combine$white)
  doWR = !is.null(Spec_combine$WR)

  #cores = max(dobase + doextend + dohot + doAGB + dowhite + doWR, cores)

  #with(plan(multisession), local=TRUE)

  message('Interpolating base spectra')
  Interp_base = progenyInterpGrid(Iso, Spec_combine$base$info, Spec_combine$base$rescale, radius=radius, weight_pow=weight_pow)

  if(doextend){
    message('Interpolating extend spectra')
    Interp_extend = progenyInterpGrid(Iso, Spec_combine$extend$info, Spec_combine$extend$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_extend = NULL
  }

  if(dohot){
    message('Interpolating hot spectra')
    Interp_hot = progenyInterpGrid(Iso, Spec_combine$hot$info, Spec_combine$hot$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_hot = NULL
  }

  if(doAGB){
    message('Interpolating AGB spectra')
    Interp_AGB = progenyInterpGrid(Iso, Spec_combine$AGB$info, Spec_combine$AGB$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_AGB = NULL
  }

  if(dowhite){
    message('Interpolating white spectra')
    Interp_white = progenyInterpGrid(Iso, Spec_combine$white$info, Spec_combine$white$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_white = NULL
  }

  if(doWR){
    message('Interpolating WR spectra')
    Interp_WR = progenyInterpGrid(Iso, Spec_combine$WR$info, Spec_combine$WR$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_WR = NULL
  }

  # Interp_base = value(Interp_base)
  #
  # if(doextend){
  #   Interp_extend = value(Interp_extend)
  # }
  #
  # if(dohot){
  #   Interp_hot = value(Interp_hot)
  # }
  #
  # if(doAGB){
  #   Interp_AGB = value(Interp_AGB)
  # }
  #
  # if(dowhite){
  #   Interp_white = value(Interp_white)
  # }
  #
  # if(doWR){
  #   Interp_WR = value(Interp_WR)
  # }

  Interp_combine = list(
    base = Interp_base,
    extend = Interp_extend,
    hot = Interp_hot,
    AGB = Interp_AGB,
    white = Interp_white,
    WR = Interp_WR
  )

  return(invisible(Interp_combine))
}

progenyInterpBest = function(Iso, Interp_combine, do_base=TRUE, do_extend=TRUE, b2e=1.5,
                             do_hot=TRUE, do_AGB=TRUE, do_white=TRUE, do_WR=TRUE,
                             prefer_hot=FALSE, prefer_AGB=FALSE, prefer_white=FALSE, prefer_WR=FALSE,
                             label_AGB=NULL, label_white=NULL, label_WR=NULL){
  setDT(Iso)

  #To be safe to not alter the original Isochrones
  Iso_temp = copy(Iso)

  best_spec = rep(0L, dim(Iso_temp)[1])

  if(do_base){
    best_spec[best_spec == 0L & Interp_combine$base$nn.idx[,1] > 0 & (Interp_combine$base$nn.dists[,1] < Interp_combine$extend$nn.dists[,1]*b2e)] = 1L
  }

  if(do_extend & !is.null(Interp_combine$extend)){
    best_spec[best_spec == 0L & Interp_combine$extend$nn.idx[,1] > 0] = 2L
  }

  if(do_hot & !is.null(Interp_combine$hot)){
    if(prefer_hot){
      best_spec[Interp_combine$hot$nn.idx[,1] > 0] = 3L
    }else{
      best_spec[best_spec == 0L & Interp_combine$hot$nn.idx[,1] > 0] = 3L
    }
  }

  if(do_AGB & !is.null(Interp_combine$AGB)){
    if(prefer_AGB){
      if(is.null(label_AGB)){
        best_spec[Interp_combine$AGB$nn.idx[,1] > 0] = 4L
      }else{
        best_spec[Interp_combine$AGB$nn.idx[,1] > 0 & Iso_temp$label %in% label_AGB] = 4L
      }
    }else{
      if(is.null(label_AGB)){
        best_spec[best_spec == 0L & Interp_combine$AGB$nn.idx[,1] > 0] = 4L
      }else{
        best_spec[best_spec == 0L & Interp_combine$AGB$nn.idx[,1] > 0 & Iso_temp$label %in% label_AGB] = 4L
      }
    }
  }

  if(do_white & !is.null(Interp_combine$white)){
    if(prefer_white){
      if(is.null(label_white)){
        best_spec[Interp_combine$white$nn.idx[,1] > 0] = 5L
      }else{
        best_spec[Interp_combine$white$nn.idx[,1] > 0 & Iso_temp$label %in% label_white] = 5L
      }
    }else{
      if(is.null(label_white)){
        best_spec[best_spec == 0L & Interp_combine$white$nn.idx[,1] > 0] = 5L
      }else{
        best_spec[best_spec == 0L & Interp_combine$white$nn.idx[,1] > 0 & Iso_temp$label %in% label_white] = 5L
      }
    }
  }

  if(do_WR & !is.null(Interp_combine$WR)){
    if(prefer_WR){
      if(is.null(label_WR)){
        best_spec[Interp_combine$WR$nn.idx[,1] > 0] = 6L
      }else{
        best_spec[Interp_combine$WR$nn.idx[,1] > 0 & Iso_temp$label %in% label_WR] = 6L
      }
    }else{
      if(is.null(label_WR)){
        best_spec[best_spec == 0L & Interp_combine$WR$nn.idx[,1] > 0] = 6L
      }else{
        best_spec[best_spec == 0L & Interp_combine$WR$nn.idx[,1] > 0 & Iso_temp$label %in% label_WR] = 6L
      }
    }
  }

  best = NULL
  return(invisible(Iso_temp[,best:=best_spec]))
}

progenyInterpStat = function(Iso, Spec_combine, Interp_combine){
  Iso_use = copy(Iso)

  Iso_len = dim(Iso_use)[1]
  mat_dim = dim(Interp_combine[[1]]$nn.idx)
  best_ID = Iso_use$best
  best_ID[best_ID == 0L] = NA

  do1 = !is.null(Spec_combine[[1]])
  do2 = !is.null(Spec_combine[[2]])
  do3 = !is.null(Spec_combine[[3]])
  do4 = !is.null(Spec_combine[[4]])
  do5 = !is.null(Spec_combine[[5]])
  do6 = !is.null(Spec_combine[[6]])

  if(do1){
    doT1 = !is.null(Spec_combine[[1]]$info$Teff)
    doG1 = !is.null(Spec_combine[[1]]$info$logG)
    doZ1 = !is.null(Spec_combine[[1]]$info$logZ)
  }else{
    doT1 = FALSE
    doG1 = FALSE
    doZ1 = FALSE
  }

  if(do2){
    doT2 = !is.null(Spec_combine[[2]]$info$Teff)
    doG2 = !is.null(Spec_combine[[2]]$info$logG)
    doZ2 = !is.null(Spec_combine[[2]]$info$logZ)
  }else{
    doT2 = FALSE
    doG2 = FALSE
    doZ2 = FALSE
  }

  if(do3){
    doT3 = !is.null(Spec_combine[[3]]$info$Teff)
    doG3 = !is.null(Spec_combine[[3]]$info$logG)
    doZ3 = !is.null(Spec_combine[[3]]$info$logZ)
  }else{
    doT3 = FALSE
    doG3 = FALSE
    doZ3 = FALSE
  }

  if(do4){
    doT4 = !is.null(Spec_combine[[4]]$info$Teff)
    doG4 = !is.null(Spec_combine[[4]]$info$logG)
    doZ4 = !is.null(Spec_combine[[4]]$info$logZ)
  }else{
    doT4 = FALSE
    doG4 = FALSE
    doZ4 = FALSE
  }

  if(do5){
    doT5 = !is.null(Spec_combine[[5]]$info$Teff)
    doG5 = !is.null(Spec_combine[[5]]$info$logG)
    doZ5 = !is.null(Spec_combine[[5]]$info$logZ)
  }else{
    doT5 = FALSE
    doG5 = FALSE
    doZ5 = FALSE
  }

  if(do6){
    doT6 = !is.null(Spec_combine[[6]]$info$Teff)
    doG6 = !is.null(Spec_combine[[6]]$info$logG)
    doZ6 = !is.null(Spec_combine[[6]]$info$logZ)
  }else{
    doT6 = FALSE
    doG6 = FALSE
    doZ6 = FALSE
  }

  stats = matrix(NA, Iso_len, 3)

  if(do1){
    # ID1 = as.integer(Interp_combine[[1]]$nn.idx)
    # ID1[ID1 == 0L] = NA
    #
    # if(doT1){
    #   logT1 = Spec_combine[[1]]$info[ID1,log10(Teff)] - log10(Iso_use$Teff)
    # }else{
    #   logT1 = rep(NA, Iso_len)
    # }
    #
    # if(doG1){
    #   logG1 = Spec_combine[[1]]$info[ID1,logG] - Iso_use$logG
    # }else{
    #   logG1 = rep(NA, Iso_len)
    # }
    #
    # if(doZ1){
    #   logZ1 = Spec_combine[[1]]$info[ID1,logZ] - Iso_use$logZ
    # }else{
    #   logZ1 = rep(NA, Iso_len)
    # }
    sel1 = Iso_use[best == 1L,, which = TRUE]
    SpecID1 = as.integer(Interp_combine[[1]]$nn.idx[sel1,1])

    if(doT1){
      stats[sel1,1] = Spec_combine[[1]]$info[SpecID1,log10(Teff)] - log10(Iso_use[sel1,Teff])
    }

    if(doG1){
      stats[sel1,2] = Spec_combine[[1]]$info[SpecID1,logG] - Iso_use[sel1,logG]
    }

    if(doZ1){
      stats[sel1,3] = Spec_combine[[1]]$info[SpecID1,logZ] - Iso_use[sel1,logZ]
    }
  }

  if(do2){
    sel2 = Iso_use[best == 2L,, which = TRUE]
    SpecID2 = as.integer(Interp_combine[[2]]$nn.idx[sel2,1])

    if(doT2){
      stats[sel2,1] = Spec_combine[[2]]$info[SpecID2,log10(Teff)] - log10(Iso_use[sel2,Teff])
    }

    if(doG2){
      stats[sel2,2] = Spec_combine[[2]]$info[SpecID2,logG] - Iso_use[sel2,logG]
    }

    if(doZ2){
      stats[sel2,3] = Spec_combine[[2]]$info[SpecID2,logZ] - Iso_use[sel2,logZ]
    }
  }

  if(do3){
    sel3 = Iso_use[best == 3L,, which = TRUE]
    SpecID3 = as.integer(Interp_combine[[3]]$nn.idx[sel3,1])

    if(doT3){
      stats[sel3,1] = Spec_combine[[3]]$info[SpecID3,log10(Teff)] - log10(Iso_use[sel3,Teff])
    }

    if(doG3){
      stats[sel3,2] = Spec_combine[[3]]$info[SpecID3,logG] - Iso_use[sel3,logG]
    }

    if(doZ3){
      stats[sel3,3] = Spec_combine[[3]]$info[SpecID3,logZ] - Iso_use[sel3,logZ]
    }
  }

  if(do4){
    sel4 = Iso_use[best == 4L,, which = TRUE]
    SpecID4 = as.integer(Interp_combine[[4]]$nn.idx[sel4,1])

    if(doT4){
      stats[sel4,1] = Spec_combine[[4]]$info[SpecID4,log10(Teff)] - log10(Iso_use[sel4,Teff])
    }

    if(doG4){
      stats[sel4,2] = Spec_combine[[4]]$info[SpecID4,logG] - Iso_use[sel4,logG]
    }

    if(doZ4){
      stats[sel4,3] = Spec_combine[[4]]$info[SpecID4,logZ] - Iso_use[sel4,logZ]
    }
  }

  if(do5){
    sel5 = Iso_use[best == 5L,, which = TRUE]
    SpecID5 = as.integer(Interp_combine[[5]]$nn.idx[sel5,1])

    if(doT5){
      stats[sel5,1] = Spec_combine[[5]]$info[SpecID5,log10(Teff)] - log10(Iso_use[sel5,Teff])
    }

    if(doG5){
      stats[sel5,2] = Spec_combine[[5]]$info[SpecID5,logG] - Iso_use[sel5,logG]
    }

    if(doZ5){
      stats[sel5,3] = Spec_combine[[5]]$info[SpecID5,logZ] - Iso_use[sel5,logZ]
    }
  }

  if(do6){
    sel6 = Iso_use[best == 6L,, which = TRUE]
    SpecID6 = as.integer(Interp_combine[[6]]$nn.idx[sel6,1])

    if(doT6){
      stats[sel6,1] = Spec_combine[[6]]$info[SpecID6,log10(Teff)] - log10(Iso_use[sel6,Teff])
    }

    if(doG6){
      stats[sel6,2] = Spec_combine[[6]]$info[SpecID6,logG] - Iso_use[sel6,logG]
    }

    if(doZ6){
      stats[sel6,3] = Spec_combine[[6]]$info[SpecID6,logZ] - Iso_use[sel6,logZ]
    }
  }

  #logT = array(c(logT1, logT2, logT3, logT4, logT5, logT6), dim = c(mat_dim, 6))
  #logG = array(c(logG1, logG2, logG3, logG4, logG5, logG6), dim = c(mat_dim, 6))
  #logZ = array(c(logZ1, logZ2, logZ3, logZ4, logZ5, logZ6), dim = c(mat_dim, 6))

  # minT = cbind(
  #   rowMins(abs(logT[,,1]), na.rm=TRUE),
  #   rowMins(abs(logT[,,2]), na.rm=TRUE),
  #   rowMins(abs(logT[,,3]), na.rm=TRUE),
  #   rowMins(abs(logT[,,4]), na.rm=TRUE),
  #   rowMins(abs(logT[,,5]), na.rm=TRUE),
  #   rowMins(abs(logT[,,6]), na.rm=TRUE)
  # )
  #
  # minG = cbind(
  #   rowMins(abs(logG[,,1]), na.rm=TRUE),
  #   rowMins(abs(logG[,,2]), na.rm=TRUE),
  #   rowMins(abs(logG[,,3]), na.rm=TRUE),
  #   rowMins(abs(logG[,,4]), na.rm=TRUE),
  #   rowMins(abs(logG[,,5]), na.rm=TRUE),
  #   rowMins(abs(logG[,,6]), na.rm=TRUE)
  # )
  #
  # minZ = cbind(
  #   rowMins(abs(logZ[,,1]), na.rm=TRUE),
  #   rowMins(abs(logZ[,,2]), na.rm=TRUE),
  #   rowMins(abs(logZ[,,3]), na.rm=TRUE),
  #   rowMins(abs(logZ[,,4]), na.rm=TRUE),
  #   rowMins(abs(logZ[,,5]), na.rm=TRUE),
  #   rowMins(abs(logZ[,,6]), na.rm=TRUE)
  # )
  #
  # meanT = cbind(
  #   rowMeans(logT[,,1], na.rm=TRUE),
  #   rowMeans(logT[,,2], na.rm=TRUE),
  #   rowMeans(logT[,,3], na.rm=TRUE),
  #   rowMeans(logT[,,4], na.rm=TRUE),
  #   rowMeans(logT[,,5], na.rm=TRUE),
  #   rowMeans(logT[,,6], na.rm=TRUE)
  # )
  #
  # meanG = cbind(
  #   rowMeans(logG[,,1], na.rm=TRUE),
  #   rowMeans(logG[,,2], na.rm=TRUE),
  #   rowMeans(logG[,,3], na.rm=TRUE),
  #   rowMeans(logG[,,4], na.rm=TRUE),
  #   rowMeans(logG[,,5], na.rm=TRUE),
  #   rowMeans(logG[,,6], na.rm=TRUE)
  # )
  #
  # meanZ = cbind(
  #   rowMeans(logZ[,,1], na.rm=TRUE),
  #   rowMeans(logZ[,,2], na.rm=TRUE),
  #   rowMeans(logZ[,,3], na.rm=TRUE),
  #   rowMeans(logZ[,,4], na.rm=TRUE),
  #   rowMeans(logZ[,,5], na.rm=TRUE),
  #   rowMeans(logZ[,,6], na.rm=TRUE)
  # )
  #
  # logT_min = minT[cbind(1:Iso_len, best_ID)]
  # logT_min[!is.finite(logT_min)] = NA
  # logG_min = minG[cbind(1:Iso_len, best_ID)]
  # logG_min[!is.finite(logG_min)] = NA
  # logZ_min = minZ[cbind(1:Iso_len, best_ID)]
  # logZ_min[!is.finite(logZ_min)] = NA
  #
  # logT_mean = meanT[cbind(1:Iso_len, best_ID)]
  # logT_mean[!is.finite(logT_mean)] = NA
  # logG_mean = meanG[cbind(1:Iso_len, best_ID)]
  # logG_mean[!is.finite(logG_mean)] = NA
  # logZ_mean = meanZ[cbind(1:Iso_len, best_ID)]
  # logZ_mean[!is.finite(logZ_mean)] = NA

  # stats = data.frame(
  #   logT_min = logT_min,
  #   logG_min = logG_min,
  #   logZ_min = logZ_min,
  #   logT_mean = logT_mean,
  #   logG_mean = logG_mean,
  #   logZ_mean = logZ_mean
  # )

  colnames(stats) = c('logT_diff', 'logG_diff', 'logZ_diff')

  Iso_use = cbind(Iso_use[,1:9], as.data.table(stats))

  logT_diff_med = median(Iso_use$logT_diff, na.rm=TRUE)
  logG_diff_med = median(Iso_use$logG_diff, na.rm=TRUE)
  logZ_diff_med = median(Iso_use$logZ_diff, na.rm=TRUE)

  logT_diff_sd = sd(Iso_use$logT_diff, na.rm=TRUE)
  logG_diff_sd = sd(Iso_use$logG_diff, na.rm=TRUE)
  logZ_diff_sd = sd(Iso_use$logZ_diff, na.rm=TRUE)

  stat_out = data.table(stat = c('logT_diff_med', 'logG_diff_med', 'logZ_diff_med', 'logT_diff_sd', 'logG_diff_sd', 'logZ_diff_sd'),
             val = c(logT_diff_med, logG_diff_med, logZ_diff_med, logT_diff_sd, logG_diff_sd, logZ_diff_sd)
            )

  return(list(Iso = Iso_use, stat = stat_out))
}
