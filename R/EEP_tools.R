## Column mapping – adjust this to your actual MIST output
# default_cols = list(
#   star_age = "star_age",
#   logL     = "log_L",
#   logTeff  = "log_Teff",
#   logTc    = "log_center_T",
#   Xc       = "center_h1",
#   Yc       = "center_he4",
#   logLH      = "log_LH",      # optional
#   Gamma_c  = "center_gamma",     # optional
#   M_star   = "star_mass",   # optional
#   M_Hshell = "m_shell_H",   # optional
#   M_Heshell= "m_shell_He",   # optional
#   Xc_C     = "center_c12",   # optional (for C burning)
#   M_core  = "he_core_mass"   # optional (for PostAGB)
# )

## Simple helper
.first_or_na = function(idx) if (length(idx) == 0) NA_integer_ else idx[1]

identify_primary_eeps = function(track,
                                  cols = list(
                                    star_age = "star_age",
                                    logL     = "log_L",
                                    logTeff  = "log_Teff",
                                    logTc    = "log_center_T",
                                    Xc       = "center_h1",
                                    Yc       = "center_he4",
                                    logLH      = "log_LH",      # optional
                                    Gamma_c  = "center_gamma",     # optional
                                    M_star   = "star_mass",   # optional
                                    M_Hshell = "m_shell_H",   # optional
                                    M_Heshell= "m_shell_He",   # optional
                                    Xc_C     = "center_c12",   # optional (for C burning)
                                    M_core  = "he_core_mass"   # optional (for PostAGB)
                                  ),
                                  params = list(
                                    prems_logTc      = 5.0,
                                    zams_lburn_frac  = 0.999,
                                    zams_dXc_max     = 0.0015,
                                    iams_Xc          = 0.3,
                                    tams_Xc          = 1e-12,
                                    rgbtip_Yc        = 0.01,
                                    zacheb_Yc        = 0.03,
                                    tacheb_Yc        = 1e-4,
                                    tpagb_Yc         = 1e-6,
                                    tpagb_delta_M    = 0.1,
                                    cburn_Xc_X       = 1e-4,
                                    postagb_env_frac = 0.2,
                                    wdcs_Gamma_max   = 100
                                  )) {
  n = nrow(track)
  get = function(name) track[[cols[[name]]]]

  logTc = get("logTc")
  logL  = get("logL")
  logTeff = get("logTeff")
  Xc   = get("Xc")
  Yc   = get("Yc")

  ## Optional columns (may be NULL if not present)
  has_col = function(name) !is.null(cols[[name]]) && cols[[name]] %in% names(track)
  logLH     = if (has_col("logLH"))  get("logLH") else NULL
  Gamma_c = if (has_col("Gamma_c")) get("Gamma_c") else NULL
  M_star  = if (has_col("M_star")) get("M_star") else NULL
  M_Hshell= if (has_col("M_Hshell")) get("M_Hshell") else NULL
  M_Heshell= if (has_col("M_Heshell")) get("M_Heshell") else NULL
  Xc_C    = if (has_col("Xc_C")) get("Xc_C") else NULL
  M_core    = if (has_col("M_core")) get("M_core") else NULL
  if(!is.null(M_core)){
    M_env = M_star - M_core
  }else{
    M_env = NULL
  }

  ## Result container
  eep_idx = c(
    PreMS   = NA_integer_,
    ZAMS    = NA_integer_,
    IAMS    = NA_integer_,
    TAMS    = NA_integer_,
    RGBTip  = NA_integer_,
    ZACHeB  = NA_integer_,
    TACHeB  = NA_integer_,
    TPAGB   = NA_integer_,
    CBurn   = NA_integer_,
    PostAGB = NA_integer_,
    WDCS    = NA_integer_
  )

  ## Below is all based on Dotter 2016 paper (rather than Iso code)
  ## 1. Pre-MS (PreMS): first point where log Tc > threshold
  i_prems = .first_or_na(which(logTc > params$prems_logTc))
  #if (is.na(i_prems)) i_prems = 1L   # if track starts already hotter, use first
  eep_idx["PreMS"] = i_prems

  ## 2. ZAMS
  Xc0 = Xc[1]
  if (!is.null(logLH)) {
    ## Using Dotter’s definition: logLH / L_tot > 0.999, before Xc has dropped by 0.0015
    Ltot  = 10^logL
    LH    = 10^logLH
    fracH = LH / Ltot
    cand  = which(
      seq_len(n) >= i_prems &
        fracH > params$zams_lburn_frac &
        Xc > (Xc0 - params$zams_dXc_max)
    )
    i_zams = .first_or_na(cand)
  } else if(is.na(eep_idx["PreMS"])){
    ## Fallback: pick first point after PreMS where Xc has begun to decrease slightly
    dXc   = c(0, diff(Xc))
    cand  = which(dXc < 0)
    i_zams = .first_or_na(cand)
    if (is.na(i_zams)) i_zams = 1L
  } else {
    ## Fallback: pick first point after PreMS where Xc has begun to decrease slightly
    dXc   = c(0, diff(Xc))
    cand  = which(seq_len(n) > i_prems & dXc < 0)
    i_zams = .first_or_na(cand)
    if (is.na(i_zams)) i_zams = i_prems + 1L
  }
  eep_idx["ZAMS"] = i_zams

  ## 3. IAMS: Xc ~ 0.3 (absolute)
  cand_iams = which(seq_len(n) >= i_zams & Xc <= params$iams_Xc)
  i_iams = .first_or_na(cand_iams)
  eep_idx["IAMS"] = i_iams

  ## 4. TAMS: Xc ~ 1e-12
  cand_tams = which(seq_len(n) >= i_iams & Xc <= params$tams_Xc)
  i_tams = .first_or_na(cand_tams)
  if (is.na(i_tams)) {
    ## If Xc never goes that low, use last point with minimum Xc
    i_tams = which.min(Xc)
  }
  eep_idx["TAMS"] = i_tams

  ## 5. RGB Tip: L max or Teff min after TAMS, before Yc drops much
  Yc_TAMS = Yc[i_tams]
  idx_range = seq.int(i_tams + 1L, n)
  if (length(idx_range) > 0) {
    ok = idx_range[Yc[idx_range] >= (Yc_TAMS - params$rgbtip_Yc)]
    if (length(ok) == 0) ok = idx_range

    logL_sub   = logL[ok]
    logTeff_sub= logTeff[ok]

    idx_Lmax   = ok[which.max(logL_sub)]
    idx_Teffmin= ok[which.min(logTeff_sub)]

    i_rgbtip = min(idx_Lmax, idx_Teffmin)  # whichever occurs first
    eep_idx["RGBTip"] = i_rgbtip
  }

  ## 6. ZACHeB: Tc minimum after RGBTip while Yc not yet much reduced
  if (!is.na(eep_idx["RGBTip"])) {
    i_rgbtip = eep_idx["RGBTip"]
    Yc_RGBTip = Yc[i_rgbtip]
    idx_range = seq.int(i_rgbtip + 1L, n)
    ok = idx_range[Yc[idx_range] > (Yc_RGBTip - params$zacheb_Yc)]
    if (length(ok) > 0) {
      logTc_sub = logTc[ok]
      i_zacheb = ok[which.min(logTc_sub)]
      eep_idx["ZACHeB"] = i_zacheb
    }
  }

  ## 7. TACHeB: Yc ~ 1e-4 (end of core He burning)
  if (!is.na(eep_idx["ZACHeB"])) {
    i_zacheb = eep_idx["ZACHeB"]
    idx_range = seq.int(i_zacheb + 1L, n)
    Yc_sub = Yc[idx_range]
    cand_tacheb = which(Yc_sub <= params$tacheb_Yc)
    if (length(cand_tacheb) > 0) {
      i_tacheb = idx_range[cand_tacheb[1]]
    } else {
      i_tacheb = tail(idx_range, 1L)
    }
    eep_idx["TACHeB"] = i_tacheb
  }

  ## 8a. TP-AGB: when (M_Hshell - M_Heshell) < 0.1 Msun after He burning (if data available) & Yc_sub < 1e-6
  if (!is.na(eep_idx["TACHeB"]) && !is.null(M_Hshell) && !is.null(M_Heshell)) {
    i_tacheb = eep_idx["TACHeB"]
    idx_range = seq.int(i_tacheb + 1L, n)
    Yc_sub = Yc[idx_range]
    dM = M_Hshell[idx_range] - M_Heshell[idx_range]
    cand_tp = which(Yc_sub < params$tpagb_Yc & dM < params$tpagb_delta_M)
    if (length(cand_tp) > 0) {
      eep_idx["TPAGB"] = idx_range[cand_tp[1]]
    }
  }

  ## 8b. CBurn: end of core C burning (Xc_C ~ 1e-4) – massive stars only
  if (!is.null(Xc_C) && !all(is.na(Xc_C)) & is.na(eep_idx["TPAGB"])) {
    ## Start after RGBTip / ZACHeB / TACHeB if available
    start_i = max(eep_idx[c("RGBTip", "ZACHeB", "TACHeB")], na.rm = TRUE)
    if (is.finite(start_i) && start_i < n) {
      idx_range = seq.int(start_i + 1L, n)
      XcC_sub = Xc_C[idx_range]
      cand_cb = which(XcC_sub <= params$cburn_Xc_X)
      if (length(cand_cb) > 0) {
        eep_idx["CBurn"] = idx_range[cand_cb[1]]
      }
    }
  }

  ## 9. Post-AGB: envelope mass < 0.2 M_star (if envelope mass can be inferred)
  ## This requires some combination of M_star and core mass; this is very
  ## model-specific
  if (!is.null(M_env) && !is.null(M_star)) {
    start_i = max(eep_idx[c("RGBTip", "ZACHeB", "TACHeB", "TPAGB", "CBurn")], na.rm = TRUE)
    if (is.finite(start_i) && start_i < n) {
      idx_range = seq.int(start_i + 1L, n)
      M_env_sub = M_env[idx_range]
      M_star_sub = M_star[idx_range]
      cand_pagb = which(M_env_sub/M_star_sub <= params$postagb_env_frac)
      if (length(cand_pagb) > 0) {
        eep_idx["PostAGB"] = idx_range[cand_pagb[1]]
      }
    }
  }

  ## 10. WDCS: Gamma_c <= 100 (if available)
  if (!is.null(Gamma_c) && !is.na(eep_idx["PostAGB"])) {
    i_post = eep_idx["PostAGB"]
    idx_range = seq.int(i_post + 1L, n)
    cand_wd = which(Gamma_c[idx_range] <= params$wdcs_Gamma_max)
    if (length(cand_wd) > 0) {
      eep_idx["WDCS"] = idx_range[cand_wd[1]]
    }
  }

  ## Return
  return(eep_idx)
}

## Compute cumulative metric distance along indices ind
## Note any log (dex) movement positively accumulates monotonically
## I.e. this captures any change, up or down, in the same sense
## Can trace other properties like log_center_T and log_center_Rho
## No limit (does not just have to be 2)
metric_distance = function(track,
                            ind,
                            vars = c("logL", "logTeff"),
                            weights = NULL,
                            cols = list(
                              star_age = "star_age",
                              logL     = "log_L",
                              logTeff  = "log_Teff",
                              logTc    = "log_center_T",
                              Xc       = "center_h1",
                              Yc       = "center_he4",
                              logLH      = "log_LH",      # optional
                              Gamma_c  = "center_gamma",     # optional
                              M_star   = "star_mass",   # optional
                              M_Hshell = "m_shell_H",   # optional
                              M_Heshell= "m_shell_He",   # optional
                              Xc_C     = "center_c12",   # optional (for C burning)
                              M_core  = "he_core_mass"   # optional (for PostAGB)
                            )) {
  if (is.null(weights)) weights = rep(1, length(vars))
  stopifnot(length(weights) == length(vars))

  #x = sapply(vars, function(v) track[[cols[[v]]]][ind])
  x = as.matrix(track[unlist(cols[vars])][ind,,drop=FALSE])
  # for (i in 2:n) {
  #   dx = x[i, ] - x[i - 1, ]
  #   D[i] = D[i - 1] + sqrt(sum(weights * dx^2))
  # }
  dx = diff(x)
  dx = rbind(rep(0, length(vars)), dx)
  metric = cumsum(sqrt(dx^2 %*% weights))
  return(metric)
}

## Interpolate all columns in 'track' along a segment defined by indices 'ind',
## placing 'n_secondary' EEPs between the endpoints (i.e. n_secondary+2 total).
.interpolate_segment_eep = function(track,
                                    ind,
                                    n_secondary,
                                    metric_vars = c("logL", "logTeff"),
                                    metric_weights = NULL,
                                    cols = list(
                                      star_age = "star_age",
                                      logL     = "log_L",
                                      logTeff  = "log_Teff",
                                      logTc    = "log_center_T",
                                      Xc       = "center_h1",
                                      Yc       = "center_he4",
                                      logLH      = "log_LH",      # optional
                                      Gamma_c  = "center_gamma",     # optional
                                      M_star   = "star_mass",   # optional
                                      M_Hshell = "m_shell_H",   # optional
                                      M_Heshell= "m_shell_He",   # optional
                                      Xc_C     = "center_c12",   # optional (for C burning)
                                      M_core  = "he_core_mass"   # optional (for PostAGB)
                                    )) {
  if (length(ind) < 2) message("Segment contain only 1 point!")

  dist_metric = metric_distance(track, ind, vars = metric_vars,
                       weights = metric_weights, cols = cols)
  
  n_points = n_secondary + 2
  D_total = dist_metric[length(dist_metric)]
  target_D = seq(0, D_total, length.out = n_points)
  
  seg = track[ind, , drop = FALSE]
  
  out = matrix(NA, n_points, dim(seg)[2])
  for(i in 1:dim(seg)[2]){
    if (is.numeric(seg[,i])) {
      out[,i] = spline(dist_metric, seg[,i], xout=target_D)$y
    } else {
      best_match = ceiling(spline(dist_metric, 1:D_total, xout=target_D))
      out[,i] = seg[best_match,i]
    }
  }
  
  return(out)
  # D_total = D[length(D)]
  # 
  # ## Number of points including boundaries
  # n_points = n_secondary + 2
  # target_D = seq(0, D_total, length.out = n_points)
  # 
  # seg = track[ind, , drop = FALSE]
  # 
  # ## For each target_D, find enclosing D interval and linearly interpolate
  # interp_row = function(tD) {
  #   if (tD <= min(D)) {
  #     return(seg[1, , drop = FALSE])
  #   }
  #   if (tD >= max(D)) {
  #     return(seg[nrow(seg), , drop = FALSE])
  #   }
  #   k = max(which(D <= tD))
  #   if (k == length(D)) return(seg[nrow(seg), , drop = FALSE])
  #   w = (tD - D[k]) / (D[k + 1] - D[k])
  #   ## linear interpolation for all numeric columns; factor/character copied from lower index
  #   out = seg[k, , drop = FALSE]
  #   for (col in names(seg)) {
  #     v = seg[[col]]
  #     if (is.numeric(v)) {
  #       out[[col]] = v[k] + w * (v[k + 1] - v[k])
  #     } else {
  #       out[[col]] = v[k]
  #     }
  #   }
  #   out
  # }
  # 
  # return(do.call(rbind, lapply(target_D, interp_row)))
}

build_eep_track = function(track,
                            primary_idx,
                            n_secondary_between = 50,
                            metric_vars = c("logL", "logTeff"),
                            metric_weights = NULL,
                            cols = list(
                              star_age = "star_age",
                              logL     = "log_L",
                              logTeff  = "log_Teff",
                              logTc    = "log_center_T",
                              Xc       = "center_h1",
                              Yc       = "center_he4",
                              logLH      = "log_LH",      # optional
                              Gamma_c  = "center_gamma",     # optional
                              M_star   = "star_mass",   # optional
                              M_Hshell = "m_shell_H",   # optional
                              M_Heshell= "m_shell_He",   # optional
                              Xc_C     = "center_c12",   # optional (for C burning)
                              M_core  = "he_core_mass"   # optional (for PostAGB)
                            )) {
  ## primary_idx: named integer vector from identify_primary_eeps()
  ## n_secondary_between: either single integer or vector of length (#segments)
  ##
  ## Returns:
  ##   track with interpolated EEPs + columns:
  ##     - EEP          : integer index along EEP track
  ##     - EEP_phase    : phase label (character/factor), inherited by secondaries
  ##     - EEP_is_primary: logical, TRUE for primary EEP rows

  # Keep only valid primary EEPs (non-NA) in chronological order
  valid     = !is.na(primary_idx)
  prim_names = names(primary_idx)[valid]
  prim_i     = as.integer(primary_idx[valid])

  o = order(prim_i)
  prim_names = prim_names[o]
  prim_i     = prim_i[o]

  n_seg = length(prim_i)
  if (n_seg < 1) stop("Not enough primary EEPs to build segments.")

  if (length(n_secondary_between) == 1L) {
    n_secondary_between = rep(n_secondary_between, n_seg)
  }
  stopifnot(length(n_secondary_between) == n_seg)

  eep_list      = list()
  phase_list    = list()
  primary_flag  = list()

  for (s in seq_len(n_seg)) {
    i_start = prim_i[s]
    if(s < n_seg){
      i_end   = prim_i[s + 1L]
    }else{
      i_end = dim(track)[1]
    }

    ind = seq.int(i_start, i_end)

    seg_eep = .interpolate_segment_eep(
      track,
      ind            = ind,
      n_secondary    = n_secondary_between[s],
      metric_vars    = metric_vars,
      metric_weights = metric_weights,
      cols           = cols
    )

    n_rows_seg = nrow(seg_eep)

    # Phase label for the whole segment: you can use start, end, or something more
    # sophisticated. For simplicity, use the start primary's label.
    seg_phase = rep(prim_names[s], n_rows_seg)

    # Primary flags: first and last rows in the segment are primary
    primary_seg = rep(FALSE, n_rows_seg)
    primary_seg[1]           = TRUE
    primary_seg[n_rows_seg]  = TRUE


    # Drop the last row of all segments except the last to avoid duplicate primaries
    if (s < n_seg) {
      seg_eep     = seg_eep[-n_rows_seg, , drop = FALSE]
      seg_phase   = seg_phase[-n_rows_seg]
      primary_seg = primary_seg[-n_rows_seg]
    }else{
      primary_seg[n_rows_seg] = FALSE
    }

    eep_list[[s]]     = seg_eep
    phase_list[[s]]   = seg_phase
    primary_flag[[s]] = primary_seg
  }

  eep_track = do.call(rbind, eep_list)
  eep_track = as.data.frame(eep_track)
  rownames(eep_track) = NULL
  colnames(eep_track) = colnames(track)

  EEP_phase    = unlist(phase_list, use.names = FALSE)
  EEP_is_primary = unlist(primary_flag, use.names = FALSE)

  eep_track$EEP            = seq_len(nrow(eep_track))
  #eep_track$EEP_phase      = factor(EEP_phase, levels = unique(prim_names))
  eep_track$EEP_phase      = EEP_phase
  eep_track$EEP_is_primary = EEP_is_primary
  
  eep_track$label <- NA_integer_
  
  # Pre–Main Sequence
  eep_track$label[eep_track$EEP_phase %in% c("PreMS")] = -1L  # MIST phase = -1
  
  # Main Sequence (ZAMS, IAMS, TAMS all map to MS)
  eep_track$label[eep_track$EEP_phase %in% c("ZAMS", "IAMS", "TAMS")] = 0L
  
  # Red Giant Branch (entire RGB, including the tip)
  # Currently no specific RGB flag
  eep_track$label[eep_track$EEP_phase %in% c("RGB", "RGBTip")] = 2L
  
  # Core He-burning (from ZACHeB through TACHeB)
  eep_track$label[eep_track$EEP_phase %in% c("ZACHeB", "CHeB", "TACHeB")] = 3L
  
  # Early AGB (after TACHeB; EEP names vary: EAGB/EBGB)
  # Currently no specific EAGB flag
  eep_track$label[eep_track$EEP_phase %in% c("EAGB", "EBGB")] = 4L
  
  # Thermally Pulsing AGB
  eep_track$label[eep_track$EEP_phase %in% c("TPAGB")] = 5L
  
  # Post-AGB bucket (6–8):
  # PostAGB
  eep_track$label[eep_track$EEP_phase %in% "PostAGB"] = 6L
  # Central Star Planetary Nebula
  # Currently no specific CSPN flag
  eep_track$label[eep_track$EEP_phase %in% "CSPN"]   = 7L
  # White Dwarf Central Star
  eep_track$label[eep_track$EEP_phase %in% "WDCS"]   = 8L
  
  # Wolf–Rayet phase (if present)
  # Currently no specific WR flag
  eep_track$label[eep_track$EEP_phase %in% c("WR")] <- 9L

  return(eep_track)
}

plot_eep = function(track,
                         eep_track,
                         cols = list(
                           x = "log_Teff",
                           y = "log_L"
                         ),
                         xlab = cols$x, ylab = cols$y,
                         show_secondary = TRUE,
                         ...) {
  # Extract originals
  x_orig   = track[[cols$x]]
  y_orig= track[[cols$y]]

  # EEP-based
  x_eep = eep_track[[cols$x]]
  y_eep    = eep_track[[cols$y]]

  # Prep primary/secondary split
  prim_mask    = eep_track$EEP_is_primary

  # Base HRD plot (Teff decreases to the right)
  magplot(x_orig, y_orig,
       type = "l", lwd = 1, col = "grey50",
       xlab = xlab, ylab = ylab, ...
       )

  # EEP-based continuous track
  lines(x_eep, y_eep,
        col = "blue", lwd = 2)

  # Secondary EEPs (small points)
  if (show_secondary) {
    points(x_eep[!prim_mask],
           y_eep[!prim_mask],
           pch = 20, cex = 0.2, col = adjustcolor("blue", alpha.f = 0.4))
  }

  # Primary EEPs (larger, coloured points)
  points(x_eep[prim_mask],
         y_eep[prim_mask],
         pch = 1, bg = "red", col = "black", cex = 1.2)

  # Add text labels for primary EEPs
  labs = as.character(eep_track$EEP_phase[prim_mask])
  # Remove "secondary" from level set if present
  valid = !is.na(labs)

  text(x_eep[prim_mask][valid],
       y_eep[prim_mask][valid],
       labels = labs[valid],
       pos = 4, cex = 0.7, col = "red")

  legend("bottomleft",
         legend = c("Original track", "EEP-based track", "Primary EEPs", "Secondary EEPs"),
         col    = c("grey50", "blue", "black", adjustcolor("blue", alpha.f = 0.4)),
         pt.bg  = c(NA, NA, "red", adjustcolor("blue", alpha.f = 0.4)),
         pch    = c(NA, NA, 21, 20),
         lty    = c(1, 1, NA, NA),
         lwd    = c(1, 2, NA, NA),
         bty    = "n", cex = 0.8)
}
