progenyParamWarp = function(x_src, y_src, x_tar, y_tar, check_order = FALSE, ...) {

  if(!requireNamespace("dtw", quietly = TRUE)){
    stop('The dtw package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
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

  # DTW alignment of curvature sequences
  dtw_alignment = dtw::dtw(k_src, k_tar, ...)

  idx_src = dtw_alignment$index1
  idx_tar = dtw_alignment$index2

  suppressWarnings({
    # Build a continuous warp function: x_tar -> x_src
    warp_tar2src = approxfun(
      x = x_eval[idx_tar],
      y = x_eval[idx_src],
      rule = 1   # set to NA at boundaries
    )

    # Build a continuous warp function: x_src x  -> x_tar
    warp_src2tar = approxfun(
      x = x_eval[idx_src],
      y = x_eval[idx_tar],
      rule = 1   # set to NA at boundaries
    )

  })

  list(
    warp_src2tar    = warp_src2tar,        # function: x_src -> x_tar
    warp_tar2src    = warp_tar2src,        # function: x_tar -> x_src
    dtw_alignment   = dtw_alignment
  )
}
