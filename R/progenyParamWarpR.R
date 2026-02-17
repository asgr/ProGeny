progenyParamWarp = function(x_src, y_src, x_tar, y_tar, smooth=FALSE, check_order=FALSE,
                            open.end=TRUE, open.begin=FALSE, ...) {

  if(!requireNamespace("dtw", quietly = TRUE)){
    stop('The dtw package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }

  if (is.matrix(x_src) || is.data.frame(x_src)) {
    if (ncol(x_src) == 1) {
      x_src = unlist(x_src[, 1])
    }else if (ncol(x_src) >= 2) {
      y_src = as.matrix(x_src[, 2:ncol(x_src)])
      x_src = unlist(x_src[, 1])
    }
  }

  if(is.data.frame(y_src)){
    y_src = as.matrix(y_src)
  }

  if (is.matrix(x_tar) || is.data.frame(x_tar)) {
    if (ncol(x_tar) == 1) {
      x_tar = unlist(x_tar[, 1])
    }else if (ncol(x_tar) >= 2) {
      y_tar = as.matrix(x_tar[, 2:ncol(x_tar)])
      x_tar = unlist(x_tar[, 1])
    }
  }

  if(is.data.frame(y_tar)){
    y_tar = as.matrix(y_tar)
  }

  if(check_order){
    # Sort inputs
    ord_s = order(x_src)
    x_src = x_src[ord_s]
    if(is.matrix(y_src)){
      y_src = y_src[ord_s,]
      stopifnot(length(x_src) == dim(y_src)[1])
    }else{
      y_src = y_src[ord_s]
      stopifnot(length(x_src) == length(y_src))
    }

    ord_t = order(x_tar)
    x_tar = x_tar[ord_t]
    if(is.matrix(y_tar)){
      y_tar = y_tar[ord_t,]
      stopifnot(length(x_tar) == dim(y_tar)[1])
      stopifnot(dim(y_src)[2] == dim(y_tar)[2])
    }else{
      y_tar = y_tar[ord_t]
      stopifnot(length(x_tar) == length(y_tar))
    }
  }

  #x_eval = sort(c(x_src, x_tar))

  if(smooth){
    if(is.matrix(y_src)){
      k_src = matrix(NA, dim(y_src)[1], dim(y_src)[2])
      for(i in 1:dim(y_src)[2]){
        k_src[,i] = smooth.spline(y_src[,i])$y
        #k_src_obj = smooth.spline(x_src, y_src[,i])
        #k_src[,i] = predict(k_src_obj, x_eval)$y
      }
    }else{
      k_src = smooth.spline(y_src)$y
      #k_src_obj = smooth.spline(x_src, y_src)
      #k_src = predict(k_src_obj, x_eval)$y
    }

    if(is.matrix(y_tar)){
      k_tar = matrix(NA, dim(y_tar)[1], dim(y_tar)[2])
      for(i in 1:dim(y_tar)[2]){
        k_tar[,i] = smooth.spline(y_tar[,i])$y
        #k_tar_obj = smooth.spline(x_tar, y_tar[,i])
        #k_tar[,i] = predict(k_tar_obj, x_eval)$y
      }
    }else{
      k_tar = smooth.spline(y_tar)$y
      #k_tar_obj = smooth.spline(x_tar, y_tar)
      #k_tar = predict(k_tar_obj, x_eval)$y
    }
  }else{
    k_src = y_src
    k_tar = y_tar
  }

  #to help pin the final value better
  # k_src_N = length(k_src)
  # k_tar_N = length(k_tar)
  # k_src_end = k_src[k_src_N]
  # k_tar_end = k_tar[k_tar_N]
  # if(k_src_end < k_tar_end){
  #   k_src = c(k_src, k_tar[max(which(k_tar < k_src_end)):k_tar_N])
  # }else{
  #   k_src = c(k_src, k_tar[max(which(k_tar > k_src_end)):k_tar_N])
  # }

  # DTW alignment of curvature sequences
  dtw_alignment = dtw::dtw(k_src, k_tar, open.end=open.end, open.begin=open.begin, ...)

  idx_src = dtw_alignment$index1
  idx_tar = dtw_alignment$index2

  suppressWarnings({
    # Build a continuous warp function: x_tar -> x_src
    warp_tar2src = function(x, wt=1){
      suppressWarnings({
        temp = approxfun(x = x_tar[idx_tar], y = x_src[idx_src], rule = 1)
        #temp = smooth.spline(x = x_tar[idx_tar], y = x_src[idx_src])
      })
      temp(x)*wt + x*(1 - wt)
      #predict(temp,x)$y*wt + x*(1 - wt)
    }

    # Build a continuous warp function: x_src x  -> x_tar
    warp_src2tar = function(x, wt=1){
      suppressWarnings({
        temp = approxfun(x = x_src[idx_src], y = x_tar[idx_tar], rule = 1)
        #temp = smooth.spline(x = x_src[idx_src], y = x_tar[idx_tar])
      })
      temp(x)*wt + x*(1 - wt)
      #predict(temp,x)$y*wt + x*(1 - wt)
    }

  })

  list(
    src = as.data.frame(cbind(x_src, y_src)),
    tar = as.data.frame(cbind(x_tar, y_tar)),
    warp_src2tar    = warp_src2tar,        # function: x_src -> x_tar
    warp_tar2src    = warp_tar2src,        # function: x_tar -> x_src
    dtw_alignment   = dtw_alignment
  )
}

progenyWarpInterp = function(x_src, y_src, x_tar, y_tar, ParamWarp_out, wt=0){

  if(missing(x_src)){
    x_src = ParamWarp_out$src
  }

  if(missing(x_tar)){
    x_tar = ParamWarp_out$tar
  }

  if (is.matrix(x_src) || is.data.frame(x_src)) {
    if (ncol(x_src) == 1) {
      x_src = unlist(x_src[, 1])
    }else if (ncol(x_src) >= 2) {
      src_colnames = colnames(x_src)
      y_src = as.matrix(x_src[, 2:ncol(x_src)])
      x_src = unlist(x_src[, 1])
    }
  }

  if(is.data.frame(y_src)){
    y_src = as.matrix(y_src)
  }

  if (is.matrix(x_tar) || is.data.frame(x_tar)) {
    if (ncol(x_tar) == 1) {
      x_tar = unlist(x_tar[, 1])
    }else if (ncol(x_tar) >= 2) {
      y_tar = as.matrix(x_tar[, 2:ncol(x_tar)])
      x_tar = unlist(x_tar[, 1])
    }
  }

  if(is.data.frame(y_tar)){
    y_tar = as.matrix(y_tar)
  }

  if(is.matrix(y_src)){
    i = NULL
    suppressWarnings({
      output = foreach(i = 1:dim(y_src)[2], .combine='cbind')%do%{
        y_tar[,i]*wt + approx(ParamWarp_out$warp_src2tar(x_src, wt=1), y_src[,i], x_tar)$y*(1 - wt)
      }

      output = cbind(ParamWarp_out$warp_tar2src(x_tar, wt = 1 - wt), output)
      colnames(output) = src_colnames
      output = as.data.frame(output)
    })
  }else{
    suppressWarnings({
      output = data.frame(
        x = ParamWarp_out$warp_tar2src(x_tar, wt = 1 - wt),
        y = y_tar*wt + approx(ParamWarp_out$warp_src2tar(x_src, wt=1), y_src, x_tar)$y*(1 - wt)
      )
    })
  }
  return(output)
}
