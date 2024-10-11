progenyIsoDownload = function(URL = 'https://tinyurl.com/progeny-isochrone/', ...) {
  # if (Iso == 'avail') {
  #   url.show(paste0(URL, 'avail.txt?raw=1'))
  # } else{
  #   Iso = paste0(Iso, '.fst')
  #
  #   download.file(paste0(URL, Iso, '?raw=1'),
  #                 destfile = paste0(destpath, Iso))
  #   return(paste0(destpath, Iso))
  # }
  browseURL(URL, ...)
}

progenyAtmosDownload = function(URL = 'https://tinyurl.com/progeny-atmos/', ...) {
  # if (atmos == 'avail') {
  #   url.show(paste0(URL, 'avail.txt?raw=1'))
  # } else{
  #   atmos = paste0(atmos, '.fits')
  #
  #   download.file(paste0(URL, atmos, '?raw=1'),
  #                 destfile = paste0(destpath, atmos))
  #   return(paste0(destpath, atmos))
  # }
  browseURL(URL, ...)
}

progenyAtmosLoad = function(destpath = '',
                      base = 'combine_C3K_conroy',
                      extend = 'combine_PHOENIX_allard',
                      hot = 'combine_TMAP_werner',
                      AGB = 'combine_AGB_lancon',
                      white = NULL,
                      wavegrid = NULL,
                      ...
                      ){
  if(!is.null(wavegrid)){
    if(!requireNamespace("ProSpect", quietly = TRUE)){
      stop('The ProSpect package is needed for spectral re-binning. Please install from GitHub asgr/ProSpect.', call. = FALSE)
    }
  }
  
  if(is.character(base)){
    base = Rfits_read(paste0(destpath,'/',base,'.fits'), pointer=F, header=F)
    
    if(anyNA(base$spec)){
      base$spec[is.na(base$spec)] = 0
    }
    
    base$spec[base$spec < 0] = 0
    
    if(!is.null(wavegrid)){
      spec_new = matrix(0, dim(base$spec)[1], length(wavegrid))
      
      for(i in 1:dim(spec_new)[1]){
        spec_new[i,] = ProSpect::specReBin(wave=base$wave, flux=base$spec[i,], wavegrid=wavegrid, ...)$flux
      }
      
      base$spec = spec_new
      base$wave = wavegrid
    }
    
  }else{
    stop('base is required!')
  }

  if(is.character(extend)){
    extend = Rfits_read(paste0(destpath,'/',extend,'.fits'), pointer=F, header=F)
    
    if(anyNA(extend$spec)){
      extend$spec[is.na(extend$spec)] = 0
    }
    
    extend$spec[extend$spec < 0] = 0
    
    if(!is.null(wavegrid)){
      spec_new = matrix(0, dim(extend$spec)[1], length(wavegrid))
      
      for(i in 1:dim(spec_new)[1]){
        spec_new[i,] = ProSpect::specReBin(wave=extend$wave, flux=extend$spec[i,], wavegrid=wavegrid, ...)$flux
      }
      
      extend$spec = spec_new
      extend$wave = wavegrid
    }
    
  }else{
    extend = NULL
  }

  if(is.character(hot)){
    hot = Rfits_read(paste0(destpath,'/',hot,'.fits'), pointer=F, header=F)
    
    if(anyNA(hot$spec)){
      hot$spec[is.na(hot$spec)] = 0
    }
    
    hot$spec[hot$spec < 0] = 0
    
    if(!is.null(wavegrid)){
      spec_new = matrix(0, dim(hot$spec)[1], length(wavegrid))
      
      for(i in 1:dim(spec_new)[1]){
        spec_new[i,] = ProSpect::specReBin(wave=hot$wave, flux=hot$spec[i,], wavegrid=wavegrid, ...)$flux
      }
      
      hot$spec = spec_new
      hot$wave = wavegrid
    }
    
  }else{
    hot = NULL
  }

  if(is.character(AGB)){
    AGB = Rfits_read(paste0(destpath,'/',AGB,'.fits'), pointer=F, header=F)
    
    if(anyNA(AGB$spec)){
      AGB$spec[is.na(AGB$spec)] = 0
    }
    
    AGB$spec[AGB$spec < 0] = 0
    
    if(!is.null(wavegrid)){
      spec_new = matrix(0, dim(AGB$spec)[1], length(wavegrid))
      
      for(i in 1:dim(spec_new)[1]){
        spec_new[i,] = ProSpect::specReBin(wave=AGB$wave, flux=AGB$spec[i,], wavegrid=wavegrid, ...)$flux
      }
      
      AGB$spec = spec_new
      AGB$wave = wavegrid
    }
    
  }else{
    AGB = NULL
  }

  if(is.character(white)){
    white = Rfits_read(paste0(destpath,'/',white,'.fits'), pointer=F, header=F)
    
    if(anyNA(white$spec)){
      white$spec[is.na(white$spec)] = 0
    }
    
    white$spec[white$spec < 0] = 0
    
    if(!is.null(wavegrid)){
      spec_new = matrix(0, dim(white$spec)[1], length(wavegrid))
      
      for(i in 1:dim(spec_new)[1]){
        spec_new[i,] = ProSpect::specReBin(wave=white$wave, flux=white$spec[i,], wavegrid=wavegrid, ...)$flux
      }
      
      white$spec = spec_new
      white$wave = wavegrid
    }
    
  }else{
    white = NULL
  }
  
  wavegrid = base$wave
  
  if(length(wavegrid) != length(base$wave)){
    stop('base wave does not much wavegrid!')
  }
  
  if(length(wavegrid) != dim(base$spec)[2]){
    stop('base spec matrix does not much wavegrid!')
  }
  
  if(!is.null(extend)){
    if(length(wavegrid) != length(extend$wave)){
      stop('extend wave does not much wavegrid!')
    }
    
    if(length(wavegrid) != dim(extend$spec)[2]){
      stop('extend spec matrix does not much wavegrid!')
    }
  }
  
  if(!is.null(hot)){
    if(length(wavegrid) != length(hot$wave)){
      stop('hot wave does not much wavegrid!')
    }
    
    if(length(wavegrid) != dim(hot$spec)[2]){
      stop('hot spec matrix does not much wavegrid!')
    }
  }
  
  if(!is.null(AGB)){
    if(length(wavegrid) != length(AGB$wave)){
      stop('AGB wave does not much wavegrid!')
    }
    
    if(length(wavegrid) != dim(AGB$spec)[2]){
      stop('AGB spec matrix does not much wavegrid!')
    }
  }
  
  if(!is.null(white)){
    if(length(wavegrid) != length(white$wave)){
      stop('white wave does not much wavegrid!')
    }
    
    if(length(wavegrid) != dim(white$spec)[2]){
      stop('white spec matrix does not much wavegrid!')
    }
  }

  Spec_combine = list(
    base = base,
    extend = extend,
    hot = hot,
    AGB = AGB,
    white = white
  )
  return(Spec_combine)
}
