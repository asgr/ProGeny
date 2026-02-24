progenyIsoDownload = function(URL = 'https://drive.google.com/drive/folders/1SnUUDgXGOiRZgcUZUKnWpG0WN0Lzg1Zs?usp=sharing', ...) {
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

progenyAtmosDownload = function(URL = 'https://drive.google.com/drive/folders/1oYS1JBQzP56cvTXuvDaw8b1ckwIZSnzw?usp=sharing', ...) {
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
                      hot = 'combine_OB_PoWR',
                      AGB = 'combine_AGB_lancon',
                      white = 'combine_TMAP_werner',
                      WR = 'combine_WNE_PoWR',
                      wavegrid = NULL,
                      cores = 8,
                      pointer = FALSE,
                      ...
                      ){
  if(!is.null(wavegrid)){
    if(pointer){
      stop('wavegrid rebinning requires the full spec matrix in memory; use pointer=FALSE to rebin spectra.')
    }
    if(!requireNamespace("ProSpect", quietly = TRUE)){
      stop('The ProSpect package is needed for spectral re-binning. Please install from GitHub asgr/ProSpect.', call. = FALSE)
    }
    if(cores > 1){
      registerDoParallel(cores = cores)
    }
  }

  if(is.character(base)){
    base_file = paste0(destpath,'/',base,'.fits')
    if(file.exists(base_file)){
      base = Rfits_read(base_file, pointer=pointer, header=F)
    }else{
      stop('No file at: ', base_file)
    }

    if(pointer){
      if(inherits(base$wave, 'Rfits_pointer')){
        base$wave = Rfits_read_image(base$wave$filename, ext=base$wave$ext, header=FALSE)
      }
    }else{
      if(anyNA(base$spec)){
        base$spec[is.na(base$spec)] = 0
      }

      base$spec[base$spec < 0] = 0
    }

    if(!is.null(wavegrid)){
      if(wavegrid[1] == 'get'){
        wavegrid = base$wave
      }else{
        Nrow = dim(base$spec)[1]
        Nwave = length(wavegrid)
        if(cores == 1){
          spec_new = matrix(0, Nrow, Nwave)
          for(i in 1:Nrow){
            spec_new[i,] = ProSpect::specReBin(wave=base$wave, flux=base$spec[i,], wavegrid=wavegrid, ...)$flux
          }
        }else{
          spec_new = foreach(i = 1:Nrow)%dopar%{
            ProSpect::specReBin(wave=base$wave, flux=base$spec[i,], wavegrid=wavegrid, ...)$flux
          }
          spec_new = do.call(rbind, spec_new)
        }

        base$spec = spec_new
        base$wave = wavegrid
      }
    }

  }else{
    stop('base is required!')
  }

  if(is.character(extend)){
    extend_file = paste0(destpath,'/',extend,'.fits')
    if(file.exists(extend_file)){
      extend = Rfits_read(extend_file, pointer=pointer, header=F)
    }else{
      stop('No file at: ', extend_file)
    }

    if(pointer){
      if(inherits(extend$wave, 'Rfits_pointer')){
        extend$wave = Rfits_read_image(extend$wave$filename, ext=extend$wave$ext, header=FALSE)
      }
    }else{
      if(anyNA(extend$spec)){
        extend$spec[is.na(extend$spec)] = 0
      }

      extend$spec[extend$spec < 0] = 0
    }

    if(!is.null(wavegrid)){
      Nrow = dim(extend$spec)[1]
      if(cores == 1){
        spec_new = matrix(0, Nrow, Nwave)
        for(i in 1:Nrow){
          spec_new[i,] = ProSpect::specReBin(wave=extend$wave, flux=extend$spec[i,], wavegrid=wavegrid, ...)$flux
        }
      }else{
        spec_new = foreach(i = 1:Nrow)%dopar%{
          ProSpect::specReBin(wave=extend$wave, flux=extend$spec[i,], wavegrid=wavegrid, ...)$flux
        }
        spec_new = do.call(rbind, spec_new)
      }

      extend$spec = spec_new
      extend$wave = wavegrid
    }

  }else{
    extend = NULL
  }

  if(is.character(hot)){
    hot_file = paste0(destpath,'/',hot,'.fits')
    if(file.exists(hot_file)){
      hot = Rfits_read(hot_file, pointer=pointer, header=F)
    }else{
      stop('No file at: ', hot_file)
    }

    if(pointer){
      if(inherits(hot$wave, 'Rfits_pointer')){
        hot$wave = Rfits_read_image(hot$wave$filename, ext=hot$wave$ext, header=FALSE)
      }
    }else{
      if(anyNA(hot$spec)){
        hot$spec[is.na(hot$spec)] = 0
      }

      hot$spec[hot$spec < 0] = 0
    }

    if(!is.null(wavegrid)){
      Nrow = dim(hot$spec)[1]
      if(cores == 1){
        spec_new = matrix(0, Nrow, Nwave)
        for(i in 1:Nrow){
          spec_new[i,] = ProSpect::specReBin(wave=hot$wave, flux=hot$spec[i,], wavegrid=wavegrid, ...)$flux
        }
      }else{
        spec_new = foreach(i = 1:Nrow)%dopar%{
          ProSpect::specReBin(wave=hot$wave, flux=hot$spec[i,], wavegrid=wavegrid, ...)$flux
        }
        spec_new = do.call(rbind, spec_new)
      }

      hot$spec = spec_new
      hot$wave = wavegrid
    }

  }else{
    hot = NULL
  }

  if(is.character(AGB)){
    AGB_file = paste0(destpath,'/',AGB,'.fits')
    if(file.exists(AGB_file)){
      AGB = Rfits_read(AGB_file, pointer=pointer, header=F)
    }else{
      stop('No file at: ', AGB_file)
    }

    if(pointer){
      if(inherits(AGB$wave, 'Rfits_pointer')){
        AGB$wave = Rfits_read_image(AGB$wave$filename, ext=AGB$wave$ext, header=FALSE)
      }
    }else{
      if(anyNA(AGB$spec)){
        AGB$spec[is.na(AGB$spec)] = 0
      }

      AGB$spec[AGB$spec < 0] = 0
    }

    if(!is.null(wavegrid)){
      Nrow = dim(AGB$spec)[1]
      if(cores == 1){
        spec_new = matrix(0, Nrow, Nwave)
        for(i in 1:Nrow){
          spec_new[i,] = ProSpect::specReBin(wave=AGB$wave, flux=AGB$spec[i,], wavegrid=wavegrid, ...)$flux
        }
      }else{
        spec_new = foreach(i = 1:Nrow)%dopar%{
          ProSpect::specReBin(wave=AGB$wave, flux=AGB$spec[i,], wavegrid=wavegrid, ...)$flux
        }
        spec_new = do.call(rbind, spec_new)
      }

      AGB$spec = spec_new
      AGB$wave = wavegrid
    }

  }else{
    AGB = NULL
  }

  if(is.character(white)){
    white_file = paste0(destpath,'/',white,'.fits')
    if(file.exists(white_file)){
      white = Rfits_read(white_file, pointer=pointer, header=F)
    }else{
      stop('No file at: ', white_file)
    }

    if(pointer){
      if(inherits(white$wave, 'Rfits_pointer')){
        white$wave = Rfits_read_image(white$wave$filename, ext=white$wave$ext, header=FALSE)
      }
    }else{
      if(anyNA(white$spec)){
        white$spec[is.na(white$spec)] = 0
      }

      white$spec[white$spec < 0] = 0
    }

    if(!is.null(wavegrid)){
      Nrow = dim(white$spec)[1]
      if(cores == 1){
        spec_new = matrix(0, Nrow, Nwave)
        for(i in 1:Nrow){
          spec_new[i,] = ProSpect::specReBin(wave=white$wave, flux=white$spec[i,], wavegrid=wavegrid, ...)$flux
        }
      }else{
        spec_new = foreach(i = 1:Nrow)%dopar%{
          ProSpect::specReBin(wave=white$wave, flux=white$spec[i,], wavegrid=wavegrid, ...)$flux
        }
        spec_new = do.call(rbind, spec_new)
      }

      white$spec = spec_new
      white$wave = wavegrid
    }

  }else{
    white = NULL
  }

  if(is.character(WR)){
    WR_file = paste0(destpath,'/',WR,'.fits')
    if(file.exists(WR_file)){
      WR = Rfits_read(WR_file, pointer=pointer, header=F)
    }else{
      stop('No file at: ', WR_file)
    }

    if(pointer){
      if(inherits(WR$wave, 'Rfits_pointer')){
        WR$wave = Rfits_read_image(WR$wave$filename, ext=WR$wave$ext, header=FALSE)
      }
    }else{
      if(anyNA(WR$spec)){
        WR$spec[is.na(WR$spec)] = 0
      }

      WR$spec[WR$spec < 0] = 0
    }

    if(!is.null(wavegrid)){
      Nrow = dim(WR$spec)[1]
      if(cores == 1){
        spec_new = matrix(0, Nrow, Nwave)
        for(i in 1:Nrow){
          spec_new[i,] = ProSpect::specReBin(wave=WR$wave, flux=WR$spec[i,], wavegrid=wavegrid, ...)$flux
        }
      }else{
        spec_new = foreach(i = 1:Nrow)%dopar%{
          ProSpect::specReBin(wave=WR$wave, flux=WR$spec[i,], wavegrid=wavegrid, ...)$flux
        }
        spec_new = do.call(rbind, spec_new)
      }

      WR$spec = spec_new
      WR$wave = wavegrid
    }

  }else{
    WR = NULL
  }

  wavegrid = base$wave

  if(length(wavegrid) != length(base$wave)){
    stop('base wave does not match wavegrid!')
  }

  if(length(wavegrid) != dim(base$spec)[2]){
    stop('base spec matrix does not match wavegrid!')
  }

  if(!is.null(extend)){
    if(length(wavegrid) != length(extend$wave)){
      stop('extend wave does not match wavegrid!')
    }

    if(length(wavegrid) != dim(extend$spec)[2]){
      stop('extend spec matrix does not match wavegrid!')
    }
  }

  if(!is.null(hot)){
    if(length(wavegrid) != length(hot$wave)){
      stop('hot wave does not match wavegrid!')
    }

    if(length(wavegrid) != dim(hot$spec)[2]){
      stop('hot spec matrix does not match wavegrid!')
    }
  }

  if(!is.null(AGB)){
    if(length(wavegrid) != length(AGB$wave)){
      stop('AGB wave does not match wavegrid!')
    }

    if(length(wavegrid) != dim(AGB$spec)[2]){
      stop('AGB spec matrix does not match wavegrid!')
    }
  }

  if(!is.null(white)){
    if(length(wavegrid) != length(white$wave)){
      stop('white wave does not match wavegrid!')
    }

    if(length(wavegrid) != dim(white$spec)[2]){
      stop('white spec matrix does not match wavegrid!')
    }
  }

  Spec_combine = list(
    base = base,
    extend = extend,
    hot = hot,
    AGB = AGB,
    white = white,
    WR = WR
  )
  return(invisible(Spec_combine))
}
