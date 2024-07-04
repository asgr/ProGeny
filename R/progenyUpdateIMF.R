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

progenyUpdateIMF = function(Iso, IMFfunc,  masslow=0.1, massmax=100, ...){
  setDT(Iso)

  i = j = k = NULL

  split_lib = split(Iso, by=c('logZ', 'logAge'), flatten=FALSE)

  #about a second
  all_bins = foreach(i = 1:length(split_lib), .combine='rbind')%do%{
    foreach(j = 1:length(split_lib[[i]]), .combine='rbind')%do%{
      data.frame(.binlims(split_lib[[i]][[j]]$Mini, log=TRUE))
    }
  }

    all_IMF_int = (all_bins$hi - all_bins$lo)*IMFfunc(Iso$Mini, masslow=masslow, massmax=massmax, ...)
    #Old, more accurate than we need (and slow)
    #this takes about 20 seconds
    #all_IMF_int = IMFfunc(all_bins, ...) #this will be the number count of stars for a total mass of 1 Msol of star formation

  #very fast
  return(data.table(lo=all_bins$lo, hi=all_bins$hi, IMFint=all_IMF_int))
}
