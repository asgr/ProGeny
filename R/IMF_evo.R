.interval = function(x, lim){
  ifelse(x < lim[1], lim[1], ifelse(x > lim[2], lim[2], x))
}

IMF_Kroupa_evo = function(mass, Age = 0, Age_lim = c(0,13.8),
                          alpha1_lim = c(0.3,0.3), alpha2_lim = c(1.3,1.3), alpha3_lim = c(2.3,2),
                          masslow_lim = c(0.01,0.01), mass1_lim = c(0.08,0.08), mass2_lim = c(0.5,0.5),
                          massmax_lim = c(150,150), Lookback_Age = 0, massform = 1,
                          massmult = FALSE, rel.tol = .Machine$double.eps^0.25, method = 'linear', ...){

  assertNumeric(Age_lim)
  ageLen = length(Age_lim)

  assertNumeric(alpha1_lim, len=ageLen)
  assertNumeric(alpha2_lim, len=ageLen)
  assertNumeric(alpha3_lim, len=ageLen)
  assertNumeric(masslow_lim, len=ageLen)
  assertNumeric(mass1_lim, len=ageLen)
  assertNumeric(mass2_lim, len=ageLen)
  assertNumeric(massmax_lim, len=ageLen)

  Age = .interval(Age + Lookback_Age, Age_lim)

  alpha1 = approxfun(Age_lim, alpha1_lim, method=method, ...)(Age)
  alpha2 = approxfun(Age_lim, alpha2_lim, method=method, ...)(Age)
  alpha3 = approxfun(Age_lim, alpha3_lim, method=method, ...)(Age)
  masslow = approxfun(Age_lim, masslow_lim, method=method, ...)(Age)
  mass1 = approxfun(Age_lim, mass1_lim, method=method, ...)(Age)
  mass2 = approxfun(Age_lim, mass2_lim, method=method, ...)(Age)
  massmax = approxfun(Age_lim, massmax_lim, method=method, ...)(Age)

  norm = integrate(.Kroupa_shape, lower=masslow, upper=massmax, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=TRUE, rel.tol=rel.tol)$value

  if(is.numeric(mass)){
    return(massform*.Kroupa_shape(mass=mass, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult)/norm)
  }else if(is.list(mass)){
    if(is.null(mass$lo) | is.null(mass$hi)){
      stop('lo and hi list elements must be present!')
    }
    if(length(mass$lo) != length(mass$hi)){
      stop('lo and hi must be the same length!')
    }
    if(any(mass$lo > mass$hi)){
      stop('some lo > hi')
    }
    i = NULL
    tempout = rep(0,length(mass$lo))
    for(i in 1:length(mass$lo)){
      tempout[i] = integrate(.Kroupa_shape, lower=mass$lo[i], upper=mass$hi[i], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult, rel.tol=rel.tol)$value
    }
    return(massform*tempout/norm)
  }
}

IMF_Kroupa_Zevo = function(mass, logZ = 0, logZ_lim = c(-4,0),
                          alpha1_lim = c(0.3,0.3), alpha2_lim = c(1.3,1.3), alpha3_lim = c(2,2.3),
                          masslow_lim = c(0.01,0.01), mass1_lim = c(0.08,0.08), mass2_lim = c(0.5,0.5),
                          massmax_lim = c(150,150), massform = 1,
                          massmult = FALSE, rel.tol = .Machine$double.eps^0.25, method = 'linear', ...){

  logZ = .interval(logZ, logZ_lim)

  alpha1 = approxfun(logZ_lim, alpha1_lim, method=method, ...)(logZ)
  alpha2 = approxfun(logZ_lim, alpha2_lim, method=method, ...)(logZ)
  alpha3 = approxfun(logZ_lim, alpha3_lim, method=method, ...)(logZ)
  masslow = approxfun(logZ_lim, masslow_lim, method=method, ...)(logZ)
  mass1 = approxfun(logZ_lim, mass1_lim, method=method, ...)(logZ)
  mass2 = approxfun(logZ_lim, mass2_lim, method=method, ...)(logZ)
  massmax = approxfun(logZ_lim, massmax_lim, method=method, ...)(logZ)

  norm = integrate(.Kroupa_shape, lower=masslow, upper=massmax, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=TRUE, rel.tol=rel.tol)$value

  if(is.numeric(mass)){
    return(massform*.Kroupa_shape(mass=mass, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult)/norm)
  }else if(is.list(mass)){
    if(is.null(mass$lo) | is.null(mass$hi)){
      stop('lo and hi list elements must be present!')
    }
    if(length(mass$lo) != length(mass$hi)){
      stop('lo and hi must be the same length!')
    }
    if(any(mass$lo > mass$hi)){
      stop('some lo > hi')
    }
    i = NULL
    tempout = rep(0,length(mass$lo))
    for(i in 1:length(mass$lo)){
      tempout[i] = integrate(.Kroupa_shape, lower=mass$lo[i], upper=mass$hi[i], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult, rel.tol=rel.tol)$value
    }
    return(massform*tempout/norm)
  }
}

IMF_Lacey_evo = function(mass, Age = 0, Age_lim = c(0, 10), alpha1_lim = c(0.4, 2),
              alpha2_lim = c(0.4, 2), alpha3_lim = c(2.3, 2), masslow_lim = c(0.1, 0.1),
              mass1_lim = c(0.1, 0.1), mass2_lim = c(1, 1), massmax_lim = c(100, 100),
              Lookback_Age = 0, massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25,
              method = 'constant', ...){
  return(IMF_Kroupa_evo(
    mass = mass,
    Age = Age,
    Age_lim = Age_lim,
    alpha1_lim = alpha1_lim,
    alpha2_lim = alpha2_lim,
    alpha3_lim = alpha3_lim,
    masslow_lim = masslow_lim,
    mass1_lim = mass1_lim,
    mass2_lim = mass2_lim,
    massmax_lim = massmax_lim,
    Lookback_Age = Lookback_Age,
    massform = massform,
    massmult = massmult,
    rel.tol = rel.tol,
    method = method,
    ...
  ))
}

IMF_Lacey_Zevo = function(mass, logZ = 0, logZ_lim = c(-4,0), alpha1_lim = c(2, 0.4),
                         alpha2_lim = c(2, 0.4), alpha3_lim = c(2, 2.3), masslow_lim = c(0.1, 0.1),
                         mass1_lim = c(0.1, 0.1), mass2_lim = c(1, 1), massmax_lim = c(100, 100),
                         massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25,
                         method = 'constant', ...){
  return(IMF_Kroupa_Zevo(
    mass = mass,
    logZ = logZ,
    logZ_lim = logZ_lim,
    alpha1_lim = alpha1_lim,
    alpha2_lim = alpha2_lim,
    alpha3_lim = alpha3_lim,
    masslow_lim = masslow_lim,
    mass1_lim = mass1_lim,
    mass2_lim = mass2_lim,
    massmax_lim = massmax_lim,
    massform = massform,
    massmult = massmult,
    rel.tol = rel.tol,
    method = method,
    ...
  ))
}
