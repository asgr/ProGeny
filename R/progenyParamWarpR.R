progenyParamWarp = function(x_src, y_src, x_tar, y_tar, check_order = FALSE, ...) {

  if(!requireNamespace("dtw", quietly = TRUE)){
    stop('The dtw package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }

  if (is.matrix(x_src) || is.data.frame(x_src)) {
    if (ncol(x_src) == 1) {
      x_src = unlist(x_src[, 1])
    }else if (ncol(x_src) == 2) {
      y_src = unlist(x_src[, 2])
      x_src = unlist(x_src[, 1])
    }
  }

  if (is.matrix(x_tar) || is.data.frame(x_tar)) {
    if (ncol(x_tar) == 1) {
      x_tar = unlist(x_tar[, 1])
    }else if (ncol(x_tar) == 2) {
      y_tar = unlist(x_tar[, 2])
      x_tar = unlist(x_tar[, 1])
    }
  }

  stopifnot(length(x_src) == length(y_src),
            length(x_tar) == length(y_tar))

  if(check_order){
    # Sort inputs
    ord_s = order(x_src)
    x_src = x_src[ord_s]
    y_src = y_src[ord_s]

    ord_t = order(x_tar)
    x_tar = x_tar[ord_t]
    y_tar = y_tar[ord_t]
  }

  x_eval = sort(c(x_src, x_tar))

  k_src_obj = smooth.spline(x_src, y_src)
  k_src = predict(k_src_obj, x_eval)$y

  k_tar_obj = smooth.spline(x_tar, y_tar)
  k_tar = predict(k_tar_obj, x_eval)$y

  #to help pin the final value better
  k_src_N = length(k_src)
  k_tar_N = length(k_tar)
  k_src_end = k_src[k_src_N]
  k_tar_end = k_tar[k_tar_N]
  if(k_src_end < k_tar_end){
    k_src = c(k_src, k_tar[max(which(k_tar < k_src_end)):k_tar_N])
  }else{
    k_src = c(k_src, k_tar[max(which(k_tar > k_src_end)):k_tar_N])
  }

  # DTW alignment of curvature sequences
  dtw_alignment = dtw::dtw(k_src, k_tar, ...)

  idx_src = dtw_alignment$index1
  idx_tar = dtw_alignment$index2

  suppressWarnings({
    # Build a continuous warp function: x_tar -> x_src
    warp_tar2src = function(x, wt=1){
      suppressWarnings({
        temp = approxfun(x = x_eval[idx_tar], y = x_eval[idx_src], rule = 1)
      })
      temp(x)*wt + x*(1 - wt)
    }

    # Build a continuous warp function: x_src x  -> x_tar
    warp_src2tar = function(x, wt=1){
      suppressWarnings({
        temp = approxfun(x = x_eval[idx_src], y = x_eval[idx_tar], rule = 1)
      })
      temp(x)*wt + x*(1 - wt)
    }

  })

  list(
    src = data.frame(x = x_src, y = y_src),
    tar = data.frame(x = x_tar, y = y_tar),
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
    }else if (ncol(x_src) == 2) {
      y_src = unlist(x_src[, 2])
      x_src = unlist(x_src[, 1])
    }
  }

  if (is.matrix(x_tar) || is.data.frame(x_tar)) {
    if (ncol(x_tar) == 1) {
      x_tar = unlist(x_tar[, 1])
    }else if (ncol(x_tar) == 2) {
      y_tar = unlist(x_tar[, 2])
      x_tar = unlist(x_tar[, 1])
    }
  }

  suppressWarnings({
    output = data.frame(
      x = ParamWarp_out$warp_tar2src(x_tar, wt = 1 - wt),
      y = y_tar*wt + approx(ParamWarp_out$warp_src2tar(x_src, wt=1), y_src, x_tar)$y*(1 - wt)
    )
  })
  return(output)
}
