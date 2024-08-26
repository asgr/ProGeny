.binlims = function(input, log=FALSE){
  N_in = length(input)

  if(log){
    input = log10(input)
  }

  in_diff = diff(input)

  bins = input[1:(N_in - 1L)] + in_diff/2
  bins = c(input[1] - in_diff[1]/2, bins, input[N_in] + in_diff[N_in - 1L]/2)

  if(log){
    bins = 10^bins
  }

  bin_lo = bins[1:N_in]
  bin_hi = bins[2:(N_in + 1L)]

  return(list(lo = bin_lo, hi = bin_hi))
}

progenyUpdateIMF = function(Iso, IMFfunc, ...){
  setDT(Iso)

  i = j = k = NULL

  split_lib = split(Iso, by=c('logZ', 'logAge'), flatten=FALSE)

  if(any(c('Age', 'redshift') %in% names(formals(IMFfunc)))){
    #If we are using an evolving form of the IMF then we need to compute
    output = foreach(i = 1:length(split_lib), .combine='rbind')%do%{
      foreach(j = 1:length(split_lib[[i]]), .combine='rbind')%do%{
        sub_bins = data.frame(.binlims(split_lib[[i]][[j]]$Mini, log=TRUE))
        if('Age' %in% names(formals(IMFfunc))){
          return(data.table(lo=sub_bins$lo, hi=sub_bins$hi, IMFint=(sub_bins$hi - sub_bins$lo)*IMFfunc(split_lib[[i]][[j]]$Mini, Age=10^split_lib[[i]][[j]]$logAge[1]/1e9, ...)))
        }else if('redshift' %in% names(formals(IMFfunc))){
          return(data.table(lo=sub_bins$lo, hi=sub_bins$hi, IMFint=(sub_bins$hi - sub_bins$lo)*IMFfunc(split_lib[[i]][[j]]$Mini, redshift=split_lib[[i]][[j]]$redshift[1], ...)))
        }
      }
    }
  }else{
    #about a second
    all_bins = foreach(i = 1:length(split_lib), .combine='rbind')%do%{
      foreach(j = 1:length(split_lib[[i]]), .combine='rbind')%do%{
        data.frame(.binlims(split_lib[[i]][[j]]$Mini, log=TRUE))
      }
    }
    #very fast
    output = data.table(lo=all_bins$lo, hi=all_bins$hi, IMFint=(all_bins$hi - all_bins$lo)*IMFfunc(Iso$Mini, ...))
  }
    #Old, more accurate than we need (and slow)
    #this takes about 20 seconds
    #all_IMF_int = IMFfunc(all_bins, ...) #this will be the number count of stars for a total mass of 1 Msol of star formation

  return(output)
}
