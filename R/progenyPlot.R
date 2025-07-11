progenyAtmosPlot = function(Spec_combine, add=FALSE,
                            do_base=TRUE, do_extend=TRUE, do_hot=TRUE, do_AGB=TRUE, do_white=TRUE, do_WR=TRUE,
                            col = c('black', 'grey', 'blue', 'red', 'darkgreen', 'orange'),
                            pch=16, cex=1, xlim=c(2e3,4e5), ylim=c(-1.5,9.5), log='x', ...){

  legvec = {}
  if(length(pch) == 1){pch = rep(pch, 6)}
  if(length(cex) == 1){cex = rep(cex, 6)}

  colvec = col[c(do_base, do_extend, do_hot, do_AGB, do_white, do_WR)]
  pchvec = pch[c(do_base, do_extend, do_hot, do_AGB, do_white, do_WR)]
  cexvec = cex[c(do_base, do_extend, do_hot, do_AGB, do_white, do_WR)]

  if(add == FALSE){
    magplot(NA, NA, xlim=xlim, ylim=ylim, log=log, xlab='Teff / K', ylab=expression(log(g / cm.s^{-2})))
  }

  if(do_WR & !is.null(Spec_combine$WR)){
    points(Spec_combine$WR$info[,list(Teff, logG)], col=col[6], pch=pch[6], cex=cex[6], ...)
    legvec = c(legvec, 'WR')
  }

  if(do_white & !is.null(Spec_combine$white)){
    points(Spec_combine$white$info[,list(Teff, logG)], col=col[5], pch=pch, cex=cex, ...)
    legvec = c(legvec, 'white')
  }

  if(do_AGB & !is.null(Spec_combine$AGB)){
    points(Spec_combine$AGB$info[,list(Teff, -1)], col=col[4], pch=pch, cex=cex, ...)
    legvec = c(legvec, 'AGB')
  }

  if(do_hot & !is.null(Spec_combine$hot)){
    points(Spec_combine$hot$info[,list(Teff, logG)], col=col[3], pch=pch, cex=cex, ...)
    legvec = c(legvec, 'hot')
  }

  if(do_extend & !is.null(Spec_combine$extend)){
    points(Spec_combine$extend$info[,list(Teff, logG + logZ*0.05)], pch=pch, col=col[2], cex=cex, ...)
    legvec = c(legvec, 'extend')
  }

  if(do_base & !is.null(Spec_combine$base)){
    points(Spec_combine$base$info[,list(Teff, logG + logZ*0.05)], pch=pch, col=col[1], cex=cex, ...)
    legvec = c(legvec, 'base')
  }

  legend('topleft', legend=rev(legvec), col=colvec, pch=pchvec, cex=1)
}

progenyIsoPlot = function(Iso, add=FALSE, Nsamp=1e4, seed=666, zsel='Lum', zunit='Lsol',
                          zstretch = 'log', zlegend = NULL,
                          col = hcl.colors(101, 'geyser'), pch=16, cex=0.5,
                          xlim=c(2e3,4e5), ylim=c(-1.5,9.5), zlim=NULL, log='x', draw_regions=FALSE, ...){

  if(add == FALSE){
    magplot(NA, NA, xlim=xlim, ylim=ylim, log=log, xlab='Teff / K', ylab=expression(log(g / cm.s^{-2})))
  }

  if(Nsamp == 'all'){
    Iso_use = copy(Iso)
  }else{
    set.seed(seed)
    Nsamp = min(Nsamp, dim(Iso)[1])
    Iso_use = copy(Iso[sample(dim(Iso)[1], Nsamp),])
  }

  if(is.null(zlim)){
    zlim = range(unlist(Iso_use[,..zsel]), na.rm=TRUE)
  }

  locut = zlim[1]
  hicut = zlim[2]

  ParmOff(plot.xy, .args = list(xy=Iso_use[,list(Teff, logG)], type='p',
                                col=col[magmap(unlist(Iso_use[,..zsel]), range=c(1,length(col)),
                                               stretch=zstretch, type='num', locut=locut, hicut=hicut)$map],
                                pch=pch, cex=cex), .pass_dots=FALSE, ...)

  if(is.null(zlegend)){
    ParmOff(magbar, .args = list(position='bottomright', range=zlim, col=col, log=(zstretch=='log'),
                                 title = paste(zsel, zunit, sep=' / ')), .pass_dots=FALSE, ...)
  }else{
    if(!isFALSE(zlegend)){
      legend('bottomright', legend=zlegend, col=col, pch=16, ...)
    }
  }

  if(draw_regions){
    polygon(c(2e3,9e3,5e4,2e3), c(2,2,5.5,5.5))
    text(2.5e3, 4.3, 'MS')
    lines(c(2e3,2.3e4), c(4,4), lty=3)
    text(2.8e3, c(3,2.3), c('RGB', 'CHeB'))
    polygon(c(4e3, 2e3, 3e3, 7e3), c(2, -1.55, -1.55, 2))
    text(2.5e3, -1.8, 'AGB')
    polygon(c(4.9e3, 1e5, 1.3e4, 1.3e4, 4e5, 4e5, 3e3), c(0.5, 6.5, 7, 10, 8.5, 6, -1.55))
    text(3e4, 6.5, 'Post-AGB')
  }
}

progenySampPlot = function(wave, z=0, v=NULL, h=NULL, add=FALSE, xlim='auto',
                           ylim='auto', xlab='Wavelength / Ang', ylab='Sampling Resolution', ...){
  wave_grid = (1+z)*(wave[-1] + wave[-length(wave)])/2
  samp_res = wave[-1]/diff(wave)
  if(add){
    lines(wave_grid, samp_res, ...)
  }else{
    if(ylim[1] == 'auto'){
      ylim = c(0, max(samp_res, na.rm=TRUE))
    }

    if(xlim[1] == 'auto'){
      #don't really need an xlim argument, but makes things seem more consistent
      xlim = range(wave_grid, na.rm=TRUE)
    }

    magplot(wave_grid, samp_res,
            log='x', type='l', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  }
  abline(v=v, lty=3, col='red')
  abline(h=h, lty=3, col='blue')
  return(invisible(data.frame(wave=wave_grid, res=samp_res)))
}

progenyIsoBestPlot = function(Iso, xsel='logZ', col = c('black', 'grey', 'blue', 'red', 'darkgreen', 'orange', 'purple'), log='', Nsamp=1e4, seed=666){

  if(tolower(Nsamp) == 'all'){
    sel = TRUE
  }else{
    set.seed(seed)
    sel = sample(dim(Iso)[1], Nsamp)
  }
  
  if(xsel == 'Teff'){
    xlab = 'Isochrone Teff / K'
  }else if(xsel == 'logG'){
    xlab = expression('Isochrone '*log(g / cm.s^{-2}))
  }else if(xsel == 'logZ'){
    xlab = 'Isochrone log(Z / Zsol)'
  }else{
    xlab = paste(Isochrone, xsel)
  }
  ylab = 'Atmos Type Frac'

  magplot(NA, NA, xlim=range(Iso[,..xsel], na.rm = TRUE), ylim=c(0,1), xlab=xlab, ylab=ylab, log=log)

  for(i in 1:6){
    lines(magrun(Iso[sel,list(.SD, best == i), .SDcols=xsel], type = 'mean', log=log, equalN=FALSE, bins=10), col=col[i], lwd=3)
  }
  lines(magrun(Iso[sel,list(.SD, best == 0), .SDcols=xsel], type = 'mean', log=log, equalN=FALSE, bins=10), col=col[7], lwd=3)
  legend('bottomleft', legend=c('base','extend','hot','AGB','white','WR', 'BB'), col=col, lty=1, lwd=3)
}

progenyAtmosErrorPlot = function(Iso, xsel='logZ', log='', Nsamp=1e4, seed=666){

  if(tolower(Nsamp) == 'all'){
    sel = TRUE
  }else{
    set.seed(seed)
    sel = sample(dim(Iso)[1], Nsamp)
  }
  
  if(xsel == 'Teff'){
    xlab = 'Isochrone Teff / K'
  }else if(xsel == 'logG'){
    xlab = expression('Isochrone '*log(g / cm.s^{-2}))
  }else if(xsel == 'logZ'){
    xlab = 'Isochrone log(Z / Zsol)'
  }else{
    xlab = paste(Isochrone, xsel)
  }
    ylab = expression('Median '*Delta)

  magplot(NA, NA, xlim=range(Iso[,..xsel], na.rm = TRUE), ylim=c(-0.3,0.3), xlab=xlab, ylab=ylab, log=log)

  temp_run = magrun(Iso[sel,list(.SD, logZ_diff), .SDcols=xsel], type = 'median', log=log, equalN=FALSE, bins=10)
  lines(temp_run, col='blue', lwd=3)
  lines(temp_run$x, temp_run$yquan[,1], col='blue', lwd=1, lty=2)
  lines(temp_run$x, temp_run$yquan[,2], col='blue', lwd=1, lty=2)

  temp_run = magrun(Iso[sel,list(.SD, logG_diff), .SDcols=xsel], type = 'median', log=log, equalN=FALSE, bins=10)
  lines(temp_run, col='darkgreen', lwd=3)
  lines(temp_run$x, temp_run$yquan[,1], col='darkgreen', lwd=1, lty=2)
  lines(temp_run$x, temp_run$yquan[,2], col='darkgreen', lwd=1, lty=2)

  temp_run = magrun(Iso[sel,list(.SD, logT_diff), .SDcols=xsel], type = 'median', log=log, equalN=FALSE, bins=10)
  lines(temp_run, col='red', lwd=3)
  lines(temp_run$x, temp_run$yquan[,1], col='red', lwd=1, lty=2)
  lines(temp_run$x, temp_run$yquan[,2], col='red', lwd=1, lty=2)

  legend('bottomleft', legend=c('logZ', 'logG', 'logT'), col=c('blue', 'darkgreen', 'red'), lty=1, lwd=3)
}
