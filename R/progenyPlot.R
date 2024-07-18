progenyAtmosPlot = function(Spec_combine, add=FALSE,
                            do_base=TRUE, do_extend=TRUE, do_hot=TRUE, do_AGB=TRUE, do_white=TRUE,
                            col = c('black', 'grey', 'blue', 'red', 'darkgreen'),
                            pch=16, cex=1, xlim=c(2e3,4e5), ylim=c(-1,9), log='x', ...){

  legvec = {}
  colvec = col[c(do_base, do_extend, do_AGB, do_hot, do_white)]

  if(add == FALSE){
    magplot(NA, NA, xlim=xlim, ylim=ylim, log=log, xlab='Teff / K', ylab=expression(log(g / cm.s^{-2})))
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

  legend('topleft', legend=rev(legvec), col=colvec, pch=1)
}

progenyIsoPlot = function(Iso, add=FALSE, sampleN=1e4, zsel='Lum', zunit='Lsol', zstretch = 'log',
                          col = hcl.colors(1000, 'geyser'), pch=16, cex=0.5,
                          xlim=c(2e3,4e5), ylim=c(-1,9), log='x', ...){

  if(add == FALSE){
    magplot(NA, NA, xlim=xlim, ylim=ylim, log=log, xlab='Teff / K', ylab=expression(log(g / cm.s^{-2})))
  }

  if(sampleN == 'all'){
    Iso_use = copy(Iso)
  }else{
    sampleN = min(sampleN, dim(Iso)[1])
    Iso_use = copy(Iso[sample(dim(Iso)[1], sampleN),])
  }

  zlim = range(unlist(Iso_use[,..zsel]))

  points(Iso_use[,list(Teff, logG)],
         col=col[magmap(unlist(Iso_use[,..zsel]), range=c(1,1e3), stretch=zstretch)$map],
         pch=pch, cex=cex, ...)
  magbar('bottomright', range=zlim, col=col, log=(zstretch=='log'), title = paste(zsel,zunit,sep=' / '))
}
