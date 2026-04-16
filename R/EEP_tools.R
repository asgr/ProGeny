## Simple helpers
.first_or_na = function(idx) if (length(idx) == 0) NA_integer_ else idx[1]

## Default secondary EEP counts per interval, matching Fortran Iso input.eep
## (value = number of secondary EEPs in the interval starting at the named primary)
.default_n_secondary = c(
  PreMS   = 200L,
  ZAMS    = 150L,
  IAMS    = 100L,
  TAMS    = 150L,
  RGBTip  = 25L,
  ZACHeB  = 75L,
  TACHeB  = 100L,
  TPAGB   = 600L,
  CBurn   = 100L,
  PostAGB = 300L
)

identify_primary_eeps = function(track,
                                  cols = list(
                                    star_age  = "star_age",
                                    logL      = "log_L",
                                    logTeff   = "log_Teff",
                                    logTc     = "log_center_T",
                                    logRhoc   = "log_center_Rho",  # for distance metric
                                    Xc        = "center_h1",
                                    Yc        = "center_he4",
                                    logg      = "log_g",           # optional (for ZAMS)
                                    logLH     = "log_LH",          # optional
                                    logLHe    = "log_LHe",         # optional (for ZACHeB)
                                    Gamma_c   = "center_gamma",    # optional
                                    M_star    = "star_mass",       # optional
                                    M_He_core = "he_core_mass",    # optional (for TPAGB)
                                    M_C_core  = "c_core_mass",     # optional (for TPAGB, PostAGB)
                                    M_O_core  = "o_core_mass",     # optional (for TPAGB, PostAGB)
                                    Xc_C      = "center_c12"       # optional (for C burning)
                                  ),
                                  params = list(
                                    prems_logTc       = 5.0,
                                    prems_logTc_adj   = 0.01,     # high-mass Tc adjustment
                                    zams_dXc_max      = 0.001,    # Fortran Iso: Xmax - 1.0d-3
                                    iams_Xc           = 0.35,     # Fortran Iso: 3.5d-1
                                    tams_Xc           = 1e-12,    # Fortran Iso: 1d-12
                                    tams_max_age      = 1.5e10,   # 15 Gyr low-mass fallback
                                    rgbtip_Yc         = 0.01,     # Fortran Iso: 1d-2
                                    zacheb_Yc         = 0.03,     # Fortran Iso: 3d-2
                                    tacheb_Yc         = 1e-4,     # Fortran Iso: 1d-4
                                    tpagb_Yc          = 1e-6,     # Fortran Iso: 1d-6
                                    tpagb_He_shell    = 0.1,      # Fortran Iso: HeShellMin = 1d-1
                                    cburn_XY_lim      = 1e-8,     # Fortran Iso: limit_XY = 1d-8
                                    cburn_C_lim       = 1e-4,     # Fortran Iso: center_carbon_limit
                                    postagb_core_frac = 0.8,      # Fortran Iso: core_mass_frac_limit
                                    wdcs_Gamma_min    = 100       # Fortran Iso: center_gamma_limit
                                  )) {
  nrow = nrow(track)

  .get = function(name, req = TRUE) {
    if (!is.null(cols[[name]]) && cols[[name]] %in% names(track)) {
      return(track[[cols[[name]]]])
    } else {
      if (req) {
        stop('Column ', name, ' is missing and is required!')
      } else {
        return(NULL)
      }
    }
  }

  ## Required columns
  logTc   = .get("logTc")
  logL    = .get("logL")
  logTeff = .get("logTeff")
  Xc      = .get("Xc")
  Yc      = .get("Yc")
  ## Optional columns (may be NULL if not present)
  logg      = .get("logg", req = FALSE)
  logLH     = .get("logLH", req = FALSE)
  logLHe    = .get("logLHe", req = FALSE)
  Gamma_c   = .get("Gamma_c", req = FALSE)
  M_star    = .get("M_star", req = FALSE)
  M_He_core = .get("M_He_core", req = FALSE)
  M_C_core = .get("M_C_core", req = FALSE)
  M_O_core = .get("M_O_core", req = FALSE)
  Xc_C      = .get("Xc_C", req = FALSE)
  star_age  = .get("star_age", req = FALSE)

  M_CO_core = M_C_core
  if(!is.null(M_O_core)){
    M_CO_core = M_CO_core + M_O_core
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

  ## 1. Pre-MS (PreMS): first point where log Tc > threshold
  ## Fortran Iso: if track starts already hotter, use logTc[1] + adjustment
  my_logTc = params$prems_logTc
  if (logTc[1] > my_logTc) {
    my_logTc = logTc[1] + params$prems_logTc_adj
  }
  i_prems = .first_or_na(which(logTc > my_logTc))
  eep_idx["PreMS"] = i_prems
  if (is.na(i_prems) || i_prems == nrow) return(eep_idx)

  ## 2. ZAMS: Fortran Iso algorithm -- find where Xc drops by zams_dXc_max, then
  ## pick max(logg) from 1 to that point.
  Xc0  = Xc[1]
  Xmin = Xc0 - params$zams_dXc_max

  ## Find first point (starting from PreMS) where Xc has dropped below Xmin
  i_start = max(1L, i_prems, na.rm=TRUE)
  i_Xcdrop = which(Xc <= Xmin & seq_along(Xc) >= i_start)
  if(length(i_Xcdrop) > 0){
    i_Xcdrop = min(i_Xcdrop)
  }else{
    i_Xcdrop = NA_integer_
  }

  if (is.na(i_Xcdrop)) i_Xcdrop = nrow

  ## ZAMS = max(logg) from 1 to i_Xcdrop (matching Fortran Iso ZAMS3)
  if (!is.null(logg)) {
    i_zams = which.max(logg[1:i_Xcdrop])
  } else {
    ## Fallback: use the Xc-drop point itself (ZAMS1)
    i_zams = i_Xcdrop
  }
  eep_idx["ZAMS"] = i_zams
  if (i_zams == 0L || i_zams == nrow) return(eep_idx)

  ## 3. IAMS: Xc < 0.35 (Fortran Iso: TAMS(t,3.5d-1,...))
  cand_iams = which(seq_len(nrow) >= (i_zams + 1L) & Xc < params$iams_Xc)
  i_iams = .first_or_na(cand_iams)
  eep_idx["IAMS"] = i_iams
  if (is.na(i_iams) || i_iams == nrow) return(eep_idx)

  ## 4. TAMS: Xc < 1e-12 (Fortran Iso: TAMS(t,1d-12,...))
  cand_tams = which(seq_len(nrow) >= (i_iams + 1L) & Xc < params$tams_Xc)
  i_tams = .first_or_na(cand_tams)
  if (is.na(i_tams)) {
    ## Fortran Iso fallback: if age > max_age for very low mass, accept last point
    if (!is.null(star_age) && !is.null(M_star) &&
        M_star[1] <= 0.5 && star_age[nrow] >= params$tams_max_age) {
      i_tams = nrow
    }
  }
  eep_idx["TAMS"] = i_tams
  if (is.na(i_tams) || i_tams == nrow) return(eep_idx)

  ## 5. RGB Tip: max logL or min logTeff after TAMS, while Yc has not
  ## dropped by more than rgbtip_Yc from the starting value.
  ## Fortran Iso: Ymin = Yc(guess) - 0.01; result = min(idx_Lmax, idx_Teffmin)
  Yc_start = Yc[i_tams + 1L]
  Ymin_rgb = Yc_start - params$rgbtip_Yc
  idx_range = seq.int(i_tams + 1L, nrow)
  if (length(idx_range) > 0) {
    ok = idx_range[Yc[idx_range] > Ymin_rgb]
    if (length(ok) > 0) {
      idx_Lmax    = ok[which.max(logL[ok])]
      idx_Teffmin = ok[which.min(logTeff[ok])]
      i_rgbtip = min(idx_Lmax, idx_Teffmin)
      eep_idx["RGBTip"] = i_rgbtip
    }
  }
  if (is.na(eep_idx["RGBTip"]) || eep_idx["RGBTip"] == nrow) return(eep_idx)
  i_rgbtip = eep_idx["RGBTip"]

  ## 6. ZACHeB (ZAHB): Fortran Iso two-phase search:
  ##   Phase 1 -- find max logLHe while Yc > Ymin (He flash peak)
  ##   Phase 2 -- find min logTc after phase-1 peak while Yc > Ymin
  Yc_start_zahb = Yc[i_rgbtip + 1L]
  Ymin_zahb = Yc_start_zahb - params$zacheb_Yc
  idx_range = seq.int(i_rgbtip + 1L, nrow)
  ok_zahb = idx_range[Yc[idx_range] > Ymin_zahb]

  if (length(ok_zahb) > 0) {
    ## Phase 1: find He luminosity peak (if logLHe available)
    if (!is.null(logLHe)) {
      i_LHe_peak = ok_zahb[which.max(logLHe[ok_zahb])]
    } else {
      i_LHe_peak = ok_zahb[1]  # fallback: start of range
    }
    ## Phase 2: find min logTc from the peak onward, still under Yc constraint
    ok_phase2 = ok_zahb[ok_zahb >= i_LHe_peak]
    if (length(ok_phase2) > 0) {
      i_zacheb = ok_phase2[which.min(logTc[ok_phase2])]
      ## Fortran Iso guard: if result equals starting guess, treat as not found
      # Not sure we actually want this?
      if (i_zacheb == (i_rgbtip + 1L)) {
        i_zacheb = i_rgbtip
        eep_idx["RGBTip"] = NA #This doesn't occur for more massive stars

      }

      eep_idx["ZACHeB"] = i_zacheb
    }
  }
  if (is.na(eep_idx["ZACHeB"]) || eep_idx["ZACHeB"] == nrow) return(eep_idx)
  i_zacheb = eep_idx["ZACHeB"]

  ## 7. TACHeB (TAHB): Yc < 1e-4  (Fortran Iso returns 0 if not found)
  idx_range = seq.int(i_zacheb + 1L, nrow)
  if (length(idx_range) > 0) {
    Yc_sub = Yc[idx_range]
    cand_tacheb = which(Yc_sub < params$tacheb_Yc)
    if (length(cand_tacheb) > 0) {
      eep_idx["TACHeB"] = idx_range[cand_tacheb[1]]
    }
    ## Fortran Iso: no fallback -- if not found, TACHeB stays NA
  }
  if (is.na(eep_idx["TACHeB"]) || eep_idx["TACHeB"] == nrow) return(eep_idx)
  i_tacheb = eep_idx["TACHeB"]

  ## 8a. TPAGB: Yc < 1e-6 AND (he_core_mass - co_core_mass) < 0.1
  ## Fortran Iso: HeShell = i_He_Core - i_CO_Core
  if (!is.null(M_He_core) && !is.null(M_CO_core)) {
    idx_range = seq.int(i_tacheb + 1L, nrow)
    if (length(idx_range) > 0) {
      Yc_sub    = Yc[idx_range]
      HeShell   = M_He_core[idx_range] - M_CO_core[idx_range]
      cand_tp   = which(Yc_sub < params$tpagb_Yc & HeShell < params$tpagb_He_shell)
      if (length(cand_tp) > 0) {
        eep_idx["TPAGB"] = idx_range[cand_tp[1]]
      }
    }
  }

  ## 8b. CBurn: Xc < 1e-8 AND Yc < 1e-8 AND Cc < 1e-4 (massive stars only)
  ## Fortran Iso: limit_XY=1d-8, center_carbon_limit=1d-4
  ## Only attempted if TPAGB was not found (mutually exclusive branches)
  if (!is.null(Xc_C) && !all(is.na(Xc_C)) && is.na(eep_idx["TPAGB"])) {
    start_i = i_tacheb
    if (is.finite(start_i) && start_i < nrow) {
      idx_range = seq.int(start_i + 1L, nrow)
      cand_cb = which(
        Xc[idx_range]   < params$cburn_XY_lim &
        Yc[idx_range]   < params$cburn_XY_lim &
        Xc_C[idx_range] < params$cburn_C_lim
      )
      if (length(cand_cb) > 0) {
        eep_idx["CBurn"] = idx_range[cand_cb[1]]
      }
    }
  }

  ## 9. PostAGB: co_core_mass / star_mass > 0.8
  ## Fortran Iso: core_mass_frac = i_co_core / i_mass, with Tc guard
  if (!is.null(M_CO_core) && !is.null(M_star)) {
    start_i = max(eep_idx[c("TPAGB", "CBurn")], na.rm=TRUE)
    if (is.finite(start_i) && start_i < nrow) {
      ## Fortran Iso Tc guard: only search if Tc is decreasing (has TP-AGB)
      Tc_now = logTc[start_i]
      Tc_end = logTc[nrow]
      if (Tc_now > Tc_end) {
        idx_range = seq.int(start_i + 1L, nrow)
        core_frac = M_CO_core[idx_range] / M_star[idx_range]
        cand_pagb = which(core_frac > params$postagb_core_frac)
        if (length(cand_pagb) > 0) {
          eep_idx["PostAGB"] = idx_range[cand_pagb[1]]
        }
      }
    }
  }

  ## 10. WDCS: center_gamma >= Gamma_min (Fortran Iso: >= center_gamma_limit)
  if (!is.null(Gamma_c) && !is.na(eep_idx["PostAGB"])) {
    i_post = eep_idx["PostAGB"]
    idx_range = seq.int(i_post + 1L, nrow)
    if (length(idx_range) > 0) {
      cand_wd = which(Gamma_c[idx_range] >= params$wdcs_Gamma_min)
      if (length(cand_wd) > 0) {
        eep_idx["WDCS"] = idx_range[cand_wd[1]]
      }
    }
  }

  ## Return
  return(eep_idx)
}

## Compute cumulative metric distance along indices ind
## Matches the Fortran Iso distance_along_track subroutine:
##   5D weighted Euclidean distance with Xc-dependent weighting
##   for log_center_Rho and log_center_T.
## Columns that are absent from the track are silently skipped
## (their contribution is zero).
metric_distance = function(track,
                            ind,
                            Teff_scale = 2.0,
                            logL_scale = 0.125,
                            Rhoc_scale = 1.0,
                            Tc_scale   = 1.0,
                            age_scale  = 0.05,
                            weight_center_rho_T_by_Xc = TRUE,
                            cols = list(
                              star_age  = "star_age",
                              logL      = "log_L",
                              logTeff   = "log_Teff",
                              logTc     = "log_center_T",
                              logRhoc   = "log_center_Rho",
                              Xc        = "center_h1"
                            )) {
  nn = length(ind)
  if (nn < 2) return(0)

  ## Helper to safely extract a column vector for the given indices
  .get = function(name, ind) {
    if (!is.null(cols[[name]]) && cols[[name]] %in% names(track)) {
      return(track[[cols[[name]]]][ind])
    } else {
      return(NULL)
    }
  }

  logTeff_v = .get("logTeff", ind)
  logL_v    = .get("logL", ind)
  logRhoc_v = .get("logRhoc", ind)
  logTc_v   = .get("logTc", ind)
  star_age_v = .get("star_age", ind)
  Xc_v      = .get("Xc", ind)

  ## Xc-dependent weighting (Fortran Iso: weight = max(0, Xc/max(Xc)))
  if (weight_center_rho_T_by_Xc && !is.null(Xc_v)) {
    max_Xc = max(Xc_v, na.rm = TRUE)
    if (max_Xc <= 0) max_Xc = 1.0
  } else {
    max_Xc = 1.0
  }

  ## Compute log10(age), protecting against non-positive ages
  log_age = NULL
  if (!is.null(star_age_v)) {
    log_age = log10(pmax(star_age_v, 1.0))
  }


  ## Indices 2:nn correspond to diffs
  idx = 2:nn

  ## Weight w[j]
  if (weight_center_rho_T_by_Xc && !is.null(Xc_v)) {
    w = pmax(0, Xc_v[idx] / max_Xc)
  } else {
    w = rep(1.0, length(idx))
  }

  ## Initialise tmp
  tmp = numeric(length(idx))

  ## Add contributions conditionally
  if (!is.null(logTeff_v)) {
    d = diff(logTeff_v)
    tmp = tmp + Teff_scale * d^2
  }

  if (!is.null(logL_v)) {
    d = diff(logL_v)
    tmp = tmp + logL_scale * d^2
  }

  if (!is.null(logRhoc_v)) {
    d = diff(logRhoc_v)
    tmp = tmp + w * Rhoc_scale * d^2
  }

  if (!is.null(logTc_v)) {
    d = diff(logTc_v)
    tmp = tmp + w * Tc_scale * d^2
  }

  if (!is.null(log_age)) {
    d = diff(log_age)
    tmp = tmp + age_scale * d^2
  }

  ## Cumulative distance
  dist_metric = c(0, cumsum(sqrt(tmp)))

  return(dist_metric)
}

## Interpolate all columns in 'track' along a segment defined by indices 'ind',
## placing 'n_secondary' EEPs between the endpoints (i.e. n_secondary+2 total).
## Uses monotone Hermite splinefun (monoH.FC), matching the spirit of the Fortran Iso
## interp_4pt_pm (Steffen's piecewise monotonic method).
## NOTE: monoH.FC (Fritsch-Carlson) and Steffen's method are both
## monotonicity-preserving cubic Hermite interpolants; results are very close.
.interpolate_segment_eep = function(track,
                                    ind,
                                    n_secondary,
                                    Teff_scale = 2.0,
                                    logL_scale = 0.125,
                                    Rhoc_scale = 1.0,
                                    Tc_scale   = 1.0,
                                    age_scale  = 0.05,
                                    weight_center_rho_T_by_Xc = TRUE,
                                    cols = list(
                                      star_age  = "star_age",
                                      logL      = "log_L",
                                      logTeff   = "log_Teff",
                                      logTc     = "log_center_T",
                                      logRhoc   = "log_center_Rho",
                                      Xc        = "center_h1"
                                    )) {
  if (length(ind) < 2) message("Segment contains only 1 point!")

  dist_metric = metric_distance(track, ind,
                                Teff_scale = Teff_scale,
                                logL_scale = logL_scale,
                                Rhoc_scale = Rhoc_scale,
                                Tc_scale   = Tc_scale,
                                age_scale  = age_scale,
                                weight_center_rho_T_by_Xc = weight_center_rho_T_by_Xc,
                                cols = cols)

  n_points = n_secondary + 2
  D_total  = tail(dist_metric, 1)

  ## Handle degenerate case where distance is zero
  if (D_total <= 0) {
    target_D = seq(0, 1, length.out = n_points)
    dist_metric = seq(0, 1, length.out = length(ind))
  } else {
    target_D = seq(0, D_total, length.out = n_points)
  }

  ## Remove any duplicate distance values to ensure strictly increasing x for splinefun
  dup = duplicated(dist_metric)
  if (any(dup)) {
    keep = !dup
    dist_metric = dist_metric[keep]
    ind = ind[keep]
  }

  seg = track[ind, , drop = FALSE]

  out = matrix(NA_real_, n_points, ncol(seg))
  for (i in seq_len(ncol(seg))) {
    if (is.numeric(seg[, i])) {
      if (length(dist_metric) >= 2) {
        out[, i] = splinefun(dist_metric, seg[, i], method = "monoH.FC")(target_D)
      } else {
        out[, i] = rep(seg[1, i], n_points)
      }
    } else {
      ## Non-numeric: carry forward from nearest lower distance point
      idx = findInterval(target_D, dist_metric, all.inside = TRUE)
      out[, i] = seg[idx, i]
    }
  }

  return(out)
}

build_eep_track = function(track,
                            primary_idx,
                            n_secondary_between = NULL,
                            Teff_scale = 2.0,
                            logL_scale = 0.125,
                            Rhoc_scale = 1.0,
                            Tc_scale   = 1.0,
                            age_scale  = 0.05,
                            weight_center_rho_T_by_Xc = TRUE,
                           cols = list(
                             star_age  = "star_age",
                             logL      = "log_L",
                             logTeff   = "log_Teff",
                             logTc     = "log_center_T",
                             logRhoc   = "log_center_Rho",
                             Xc        = "center_h1"
                           )) {
  ## primary_idx: named integer vector from identify_primary_eeps()
  ## n_secondary_between: NULL for Fortran Iso defaults, single integer,
  ##   or vector of length (number_of_intervals)
  ##
  ## Returns:
  ##   data.frame with interpolated EEPs + columns:
  ##     - EEP          : integer index along EEP track
  ##     - EEP_phase    : phase label (character), inherited by secondaries
  ##     - EEP_is_primary: logical, TRUE for primary EEP rows
  ##     - label        : integer coarse phase code

  # Keep only valid primary EEPs (non-NA) in chronological order
  valid      = !is.na(primary_idx)
  prim_names = names(primary_idx)[valid]
  prim_i     = as.integer(primary_idx[valid])

  o = order(prim_i)
  prim_names = prim_names[o]
  prim_i     = prim_i[o]

  n_prim = length(prim_i)
  if (n_prim == 0) stop("Need at least 1 primary EEP to build segments.")

  ## Number of intervals between consecutive primaries (matching Fortran Iso)
  #n_intervals = n_prim - 1
  n_intervals = n_prim

  ## Determine secondary EEP counts per interval
  if (is.null(n_secondary_between)) {
    ## Use Fortran Iso-style defaults keyed by primary name
    n_secondary_between = integer(n_intervals)
    for (s in seq_len(n_intervals)) {
      nm = prim_names[s]
      if (nm %in% names(.default_n_secondary)) {
        n_secondary_between[s] = .default_n_secondary[[nm]]
      } else {
        n_secondary_between[s] = 50L  # safe fallback
      }
    }
  } else if (length(n_secondary_between) == 1L) {
    n_secondary_between = rep(as.integer(n_secondary_between), n_intervals)
  }
  stopifnot(length(n_secondary_between) == n_intervals)

  eep_list      = list()
  phase_list    = list()
  primary_flag  = list()

  for (s in 1:n_intervals) {
    i_start = prim_i[s]

    if(s < n_intervals){
      i_end   = prim_i[s + 1L]
    }else{
      i_end = dim(track)[1]
    }

    ind     = seq.int(i_start, i_end)

    seg_eep = .interpolate_segment_eep(
      track,
      ind         = ind,
      n_secondary = n_secondary_between[s],
      Teff_scale  = Teff_scale,
      logL_scale  = logL_scale,
      Rhoc_scale  = Rhoc_scale,
      Tc_scale    = Tc_scale,
      age_scale   = age_scale,
      weight_center_rho_T_by_Xc = weight_center_rho_T_by_Xc,
      cols        = cols
    )

    n_rows_seg = nrow(seg_eep)

    # Phase label: use the start primary's label
    seg_phase   = rep(prim_names[s], n_rows_seg)

    # Primary flags: first rows in the segment are primary
    primary_seg            = rep(FALSE, n_rows_seg)
    primary_seg[1]         = TRUE

    # Drop the last row of all segments except the last to avoid
    # duplicate primaries at interval boundaries
    if (s < n_intervals) {
      seg_eep     = seg_eep[-n_rows_seg, , drop = FALSE]
      seg_phase   = seg_phase[-n_rows_seg]
      primary_seg = primary_seg[-n_rows_seg]
    }

    eep_list[[s]]     = seg_eep
    phase_list[[s]]   = seg_phase
    primary_flag[[s]] = primary_seg
  }

  eep_track = do.call(rbind, eep_list)
  eep_track = as.data.frame(eep_track)
  rownames(eep_track) = NULL
  colnames(eep_track) = colnames(track)

  EEP_phase      = unlist(phase_list, use.names = FALSE)
  EEP_is_primary = unlist(primary_flag, use.names = FALSE)

  eep_track$EEP            = seq_len(nrow(eep_track))
  eep_track$EEP_phase      = EEP_phase
  eep_track$EEP_is_primary = EEP_is_primary

  eep_track$label = NA_integer_

  # Pre-Main Sequence
  eep_track$label[eep_track$EEP_phase %in% c("PreMS")] = -1L  # MIST phase = -1

  # Main Sequence (ZAMS, IAMS )
  eep_track$label[eep_track$EEP_phase %in% c("ZAMS", "IAMS")] = 0L

  # Red Giant Branch (TAMS, RGB)
  # Note TAMS ends up asigned to 2 in original Fortran Iso (which seems odd to me)
  eep_track$label[eep_track$EEP_phase %in% c("TAMS", "RGB")] = 2L

  # Core He-burning (ZACHeB, CHeB)
  # Note RGBTip ends up asigned to 3 in original Fortran Iso (which seems odd to me)
  eep_track$label[eep_track$EEP_phase %in% c("RGBTip", "ZACHeB", "CHeB")] = 3L

  # Early AGB (EEP names vary: EAGB / EBGB)
  # Note TACHeB ends up asigned to 4 in original Fortran Iso (which seems odd to me)
  eep_track$label[eep_track$EEP_phase %in% c("TACHeB", "EAGB", "EBGB")] = 4L

  # Thermally Pulsing AGB (TPAGB)
  eep_track$label[eep_track$EEP_phase %in% c("TPAGB")] = 5L

  # Post-AGB bucket (6-8):
  # PostAGB
  eep_track$label[eep_track$EEP_phase %in% "PostAGB"] = 6L
  # Central Star Planetary Nebula
  eep_track$label[eep_track$EEP_phase %in% "CSPN"]   = 7L
  # White Dwarf Central Star
  eep_track$label[eep_track$EEP_phase %in% "WDCS"]   = 8L

  # Wolf-Rayet phase (if present, currently do not flag this in ProGeny)
  eep_track$label[eep_track$EEP_phase %in% c("WR")] = 9L

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
  if(is.vector(cols)){
    cols = as.list(cols)
  }

  if(!'x' %in% names(cols)){
    names(cols)[1] = 'x'
  }

  if(!'y' %in% names(cols)){
    names(cols)[2] = 'y'
  }

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
       type = "l", lwd = 2, col = "grey50",
       xlab = xlab, ylab = ylab, ...
       )

  # EEP-based continuous track
  lines(x_eep, y_eep,
        col = "lightblue", lwd = 1)

  # Secondary EEPs (small points)
  if (show_secondary) {
    points(x_eep[!prim_mask],
           y_eep[!prim_mask],
           pch = 20, cex = 0.5, col = "blue")
  }

  # Primary EEPs (larger, coloured points)
  points(x_eep[prim_mask],
         y_eep[prim_mask],
         pch = 21, bg = "red", col = "black", cex = 1.2)

  # Add text labels for primary EEPs
  labs = as.character(eep_track$EEP_phase[prim_mask])
  # Remove "secondary" from level set if present
  valid = !is.na(labs)

  text(x_eep[prim_mask][valid],
       y_eep[prim_mask][valid],
       labels = labs[valid],
       pos = 4, cex = 0.7, col = "red")

  legend("bottomright",
         legend = c("Original track", "EEP-based track", "Primary EEPs", "Secondary EEPs"),
         col    = c("grey50", "blue", "black", "blue"),
         pt.bg  = c(NA, NA, "red", "blue"),
         pch    = c(NA, NA, 21, 20),
         lty    = c(1, 1, NA, NA),
         lwd    = c(2, 1, NA, NA),
         bty    = "n", cex = 0.8)
}
