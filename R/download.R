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
                      base = 'combine_PHOENIX_alpha+0.0',
                      extend = 'combine_basel_cut',
                      hot = 'combine_hot',
                      AGB = 'combine_AGB',
                      white = 'combine_white'){
  Spec_combine = list(
    base = Rfits_read(paste0(destpath,'/',base,'.fits'), pointer=F, header=F),
    extend = Rfits_read(paste0(destpath,'/',extend,'.fits'), pointer=F, header=F),
    hot = Rfits_read(paste0(destpath,'/',hot,'.fits'), pointer=F, header=F),
    AGB = Rfits_read(paste0(destpath,'/',AGB,'.fits'), pointer=F, header=F),
    white = Rfits_read(paste0(destpath,'/',white,'.fits'), pointer=F, header=F)
  )
  return(Spec_combine)
}
