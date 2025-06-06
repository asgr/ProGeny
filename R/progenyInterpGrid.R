progenyInterpGrid = function(loc, info, rescale, radius=2, weight_pow=2, k=8){
  setDT(loc)
  setDT(info)

  col_keep = Teff = logG = logZ = NULL

  if(rescale$Teff != -999){
    has_Teff = TRUE
    col_keep = c(col_keep, 'Teff')
  }else{
    has_Teff = FALSE
  }

  if(rescale$logG != -999){
    has_logG = TRUE
    col_keep = c(col_keep, 'logG')
  }else{
    has_logG = FALSE
  }

  if(rescale$logZ != -999){
    has_logZ = TRUE
    col_keep = c(col_keep, 'logZ')
  }else{
    has_logZ = FALSE
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
