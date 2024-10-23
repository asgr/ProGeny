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

  # if(length(rescale) == 3){
  #   loc = loc[,list(Teff/rescale[1], logG/rescale[2], logZ/rescale[3])]
  #   info = info[,list(Teff/rescale[1], logG/rescale[2], logZ/rescale[3])]
  # }else if(length(rescale) == 2){
  #   loc = loc[,list(Teff/rescale[1], logG/rescale[2])]
  #   info = info[,list(Teff/rescale[1], logG/rescale[2])]
  # }else if(length(rescale) == 1){
  #   loc = loc[,list(Teff/rescale[1])]
  #   info = info[,list(Teff/rescale[1])]
  # }

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

  Interp_base = progenyInterpGrid(Iso, Spec_combine$base$info, Spec_combine$base$rescale, radius=radius, weight_pow=weight_pow)

  if(!is.null(Spec_combine$extend)){
    Interp_extend = progenyInterpGrid(Iso, Spec_combine$extend$info, Spec_combine$extend$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_extend = NULL
  }

  if(!is.null(Spec_combine$hot)){
    Interp_hot = progenyInterpGrid(Iso, Spec_combine$hot$info, Spec_combine$hot$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_hot = NULL
  }

  if(!is.null(Spec_combine$AGB)){
    Interp_AGB = progenyInterpGrid(Iso, Spec_combine$AGB$info, Spec_combine$AGB$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_AGB = NULL
  }

  if(!is.null(Spec_combine$white)){
    Interp_white = progenyInterpGrid(Iso, Spec_combine$white$info, Spec_combine$white$rescale, radius=radius, weight_pow=weight_pow, k=k)
  }else{
    Interp_white = NULL
  }

  Interp_combine = list(
    base = Interp_base,
    extend = Interp_extend,
    hot = Interp_hot,
    AGB = Interp_AGB,
    white = Interp_white
  )

  return(invisible(Interp_combine))
}

progenyInterpBest = function(Iso, Interp_combine, do_base=TRUE, do_extend=TRUE, do_hot=TRUE,
                             do_AGB=TRUE, do_white=TRUE, b2e=1.5, label_AGB=NULL, label_white=NULL){
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
    best_spec[best_spec == 0L & Interp_combine$hot$nn.idx[,1] > 0] = 3L
  }
  if(do_AGB & !is.null(Interp_combine$AGB)){
    if(is.null(label_AGB)){
      best_spec[best_spec == 0L & Interp_combine$AGB$nn.idx[,1] > 0] = 4L
    }else{
      best_spec[best_spec == 0L & Interp_combine$AGB$nn.idx[,1] > 0 & Iso_temp$label %in% label_AGB] = 4L
    }
  }
  if(do_white & !is.null(Interp_combine$white)){
    if(is.null(label_white)){
      best_spec[best_spec == 0L & Interp_combine$white$nn.idx[,1] > 0] = 5L
    }else{
      best_spec[best_spec == 0L & Interp_combine$white$nn.idx[,1] > 0 & Iso_temp$label %in% label_white] = 5L
    }
  }

  best = NULL
  return(invisible(Iso_temp[,best:=best_spec]))
}
