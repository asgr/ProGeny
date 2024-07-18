progenyIsoDownload = function(Iso = 'avail',
                          destpath = '',
                          URL = 'https://tinyurl.com/progeny_isochrone/') {
  if (Iso == 'avail') {
    url.show(paste0(URL, 'avail.txt?raw=1'))
  } else{
    Iso = paste0(Iso, '.fst')

    download.file(paste0(URL, Iso, '?raw=1'),
                  destfile = paste0(destpath, Iso))
    return(paste0(destpath, Iso))
  }
}

progenyAtmosDownload = function(atmos = 'avail',
                            destpath = '',
                            URL = 'https://tinyurl.com/progeny_atmos/') {
  if (atmos == 'avail') {
    url.show(paste0(URL, 'avail.txt?raw=1'))
  } else{
    atmos = paste0(atmos, '.fits')

    download.file(paste0(URL, atmos, '?raw=1'),
                  destfile = paste0(destpath, atmos))
    return(paste0(destpath, atmos))
  }
}

progenyAtmosLoad = function(destpath= '',
                      base = 'combine_PHOENIX_husser',
                      extend = 'combine_PHOENIX_allard',
                      hot = 'combine_TMAP_werner',
                      AGB = 'combine_AGB',
                      white = NULL){
  if(is.character(base)){
    base = Rfits_read(paste0(destpath,'/',base,'.fits'), pointer=F, header=F)
    if(anyNA(base$spec)){
      base$spec[is.na(base$spec)] = 0
    }
    base$spec[base$spec < 0] = 0
  }else{
    base = NULL
  }

  if(is.character(extend)){
    extend = Rfits_read(paste0(destpath,'/',extend,'.fits'), pointer=F, header=F)
    if(anyNA(extend$spec)){
      extend$spec[is.na(extend$spec)] = 0
    }
    extend$spec[extend$spec < 0] = 0
  }else{
    extend = NULL
  }

  if(is.character(hot)){
    hot = Rfits_read(paste0(destpath,'/',hot,'.fits'), pointer=F, header=F)
    if(anyNA(hot$spec)){
      hot$spec[is.na(hot$spec)] = 0
    }
    hot$spec[hot$spec < 0] = 0
  }else{
    hot = NULL
  }

  if(is.character(AGB)){
    AGB = Rfits_read(paste0(destpath,'/',AGB,'.fits'), pointer=F, header=F)
    if(anyNA(AGB$spec)){
      AGB$spec[is.na(AGB$spec)] = 0
    }
    AGB$spec[AGB$spec < 0] = 0
  }else{
    AGB = NULL
  }

  if(is.character(white)){
    white = Rfits_read(paste0(destpath,'/',white,'.fits'), pointer=F, header=F)
    if(anyNA(white$spec)){
      white$spec[is.na(white$spec)] = 0
    }
    white$spec[white$spec < 0] = 0
  }else{
    white = NULL
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
