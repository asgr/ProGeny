progenyUpdateIMF = function(Iso, IMFfunc, max_width=NULL, ...){
  setDT(Iso)

  order_check = order(Iso$logZ, Iso$logAge, Iso$Mini)

  if(any(diff(order_check) != 1L)){
    reorder = order(order_check) #this will reorder the results correctly
    reorder_flag = TRUE
  }else{
    reorder_flag = FALSE
  }

  i = j = k = NULL


  split_lib = split(Iso, by=c('logZ', 'logAge'), flatten=FALSE)

  #If we are using an evolving form of the IMF then we need to compute
  output = foreach(i = 1:length(split_lib), .combine = 'rbind') %do% {
    foreach(j = 1:length(split_lib[[i]]), .combine = 'rbind') %do% {
      #when set to log, using 0.1/0.2 means the Mini can only change by 26/58% between bins
      sub_bins = data.frame(.binlims(split_lib[[i]][[j]]$Mini, log = TRUE, max_width = max_width))
      if ('Age' %in% names(formals(IMFfunc))){
        sub_bins$IMFint = (sub_bins$hi - sub_bins$lo) * IMFfunc(split_lib[[i]][[j]]$Mini,
                                                                Age = 10^split_lib[[i]][[j]]$logAge[1] / 1e9, ...)
      }else if('redshift' %in% names(formals(IMFfunc))){
        sub_bins$IMFint = (sub_bins$hi - sub_bins$lo) * IMFfunc(split_lib[[i]][[j]]$Mini,
                                                                redshift = split_lib[[i]][[j]]$redshift[1], ...)
      }else if('logZ' %in% names(formals(IMFfunc))){
        sub_bins$IMFint = (sub_bins$hi - sub_bins$lo) * IMFfunc(split_lib[[i]][[j]]$Mini,
                                                                logZ = split_lib[[i]][[j]]$logZ[1], ...)
      }else{
        #This should be most standard IMF functions in practice
        sub_bins$IMFint = (sub_bins$hi - sub_bins$lo) * IMFfunc(split_lib[[i]][[j]]$Mini, ...)
      }
      return(sub_bins)
    }
  }

  #Old, more accurate than we need (and slow)
  #this takes about 20 seconds
  #all_IMF_int = IMFfunc(all_bins, ...) #this will be the number count of stars for a total mass of 1 Msol of star formation

  if(reorder_flag){
    #This is to reorder the output in case our isochrones are not logZ, logAge ordered
    output = output[reorder,]
  }

  return(invisible(output))
}
