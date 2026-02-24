progenyIso2Spec = function(logAge=8.4, logZ=0, logA=NULL, Iso, IMFint, Spec_combine,
                           Interp_combine, interp=FALSE){
  setDT(Iso)

  Lum = ID = wt = wt_sum = Teff = best = NULL

  #something here - change allZ to ParsecIso
  wave_grid = Spec_combine$base$wave
  spec_stack = rep(0, length(wave_grid))

  use = sort(unique(Iso$best))
  use = use[use > 0]

  IsoAge = sort(unique(Iso$logAge))
  IsoZ = sort(unique(Iso$logZ))
  if(!is.null(logA)){
    if(!'logA' %in% colnames(Iso)){
      stop('No logA column in Iso!')
    }

    # if(!'logA' %in% colnames(Spec_combine[[1]]$info)){
    #   stop('No logA column in base Spec_combine!')
    # }
    IsoA = sort(unique(Iso$logA))
  }

  if(interp){
    if(logAge < min(IsoAge) | logAge > max(IsoAge)){
      stop('Requested logAge not in range of isochrone logAge values!')
    }

    if(logZ < min(IsoZ) | logZ > max(IsoZ)){
      stop('Requested logZ not in range of isochrone logZ values!')
    }

    if(!is.null(logA)){
      if(logA < min(IsoA) | logA > max(IsoA)){
        stop('Requested logA not in range of isochrone logA values!')
      }
    }

    temp_Age = .interp_quick(logAge, IsoAge)
    logAge_loc = as.integer(temp_Age[1:2])
    if(logAge_loc[1] == logAge_loc[2]){
      logAge_steps = IsoAge[logAge_loc[1]]
      logAge_wt = 1
    }else{
      logAge_steps = IsoAge[logAge_loc]
      logAge_wt = temp_Age[3:4]
    }

    temp_Z = .interp_quick(logZ, IsoZ)
    logZ_loc = as.integer(temp_Z[1:2])
    if(logZ_loc[1] == logZ_loc[2]){
      logZ_steps = IsoZ[logZ_loc[1]]
      logZ_wt = 1
    }else{
      logZ_steps = IsoZ[logZ_loc]
      logZ_wt = temp_Z[3:4]
    }

    if(!is.null(logA)){
      temp_A = .interp_quick(logA, IsoA)
      logA_loc = as.integer(temp_A[1:2])
      if(logA_loc[1] == logA_loc[2]){
        logA_steps = IsoA[logA_loc[1]]
        logA_wt = 1
      }else{
        logA_steps = IsoA[logA_loc]
        logA_wt = temp_A[3:4]
      }
    }
  }else{
    logAge_steps = IsoAge[which.min(abs(IsoAge - logAge))]
    logAge_wt = 1
    if(logAge_steps != logAge){
      message('Using ', logAge_steps, ' as nearest logAge to the requested ', logAge)
    }

    logZ_steps = IsoZ[which.min(abs(IsoZ - logZ))]
    logZ_wt = 1
    if(logZ_steps != logZ){
      message('Using ', logZ_steps, ' as nearest logZ to the requested ', logZ)
    }

    if(!is.null(logA)){
      logA_steps = IsoA[which.min(abs(IsoA - logA))]
      logA_wt = 1
      if(logA_steps != logA){
        message('Using ', logA_steps, ' as nearest logA to the requested ', logA)
      }
    }
  }

  bb_spec = numeric(length(wave_grid))

  for(i in use){
    temp_stack_list = list()
    for(logAge_j in seq_along(logAge_steps)){
      logAge_step = logAge_steps[logAge_j]
      for(logZ_k in seq_along(logZ_steps)){
        logZ_step = logZ_steps[logZ_k]
        interp_wt = logAge_wt[logAge_j]*logZ_wt[logZ_k]
        if(is.null(logA)){
          subset = Iso[logZ==logZ_step & logAge==logAge_step & best==i, which=TRUE] #faster!
          if(length(subset) > 0){
            temp_DT = data.table(
              ID=as.vector(Interp_combine[[i]]$nn.idx[subset,]),
              wt=as.vector(Interp_combine[[i]]$weights[subset,])*Iso[subset,Lum]*IMFint[subset]
            )
            temp_DT = temp_DT[ID > 0,]
            pre_stack = temp_DT[,list(wt_sum=sum(wt)),by=ID]
            pre_stack[,wt_sum:=wt_sum*interp_wt]
            temp_stack_list[[length(temp_stack_list) + 1L]] = pre_stack
          }
        }else{
#          if('logA' %in% colnames(Spec_combine[[i]]$info)){
            for(logA_m in seq_along(logA_steps)){
              logA_step = logA_steps[logA_m]
              interp_wt = logAge_wt[logAge_j]*logZ_wt[logZ_k]*logA_wt[logA_m]
              subset = Iso[logZ==logZ_step & logAge==logAge_step & logA==logA_step & best==i, which=TRUE] #faster!
              if(length(subset) > 0){
                temp_DT = data.table(
                  ID=as.vector(Interp_combine[[i]]$nn.idx[subset,]),
                  wt=as.vector(Interp_combine[[i]]$weights[subset,])*Iso[subset,Lum]*IMFint[subset]
                )
                temp_DT = temp_DT[ID > 0,]
                pre_stack = temp_DT[,list(wt_sum=sum(wt)),by=ID]
                pre_stack[,wt_sum:=wt_sum*interp_wt]
                temp_stack_list[[length(temp_stack_list) + 1L]] = pre_stack
              }
            }
#          }
        }
      }
    }
    if(length(temp_stack_list) > 0L){
      temp_stack = rbindlist(temp_stack_list)
      temp_stack = temp_stack[,list(wt_sum=sum(wt_sum)),by=ID]
      #This version didn't seem to help, so removed to reduce complication installing ProGeny
      #.vec_add_cpp(spec_stack, crossprod(temp_stack$wt_sum, Spec_combine[[i]]$spec[temp_stack$ID,]))
      spec_mat = Spec_combine[[i]]$spec
      if(inherits(spec_mat, 'Rfits_pointer')){
        id_lo = min(temp_stack$ID)
        id_hi = max(temp_stack$ID)
        spec_block = spec_mat[id_lo:id_hi, ]  # reads only the needed row range from disk
        # relative indexing is safe: all IDs are in [id_lo, id_hi] by construction
        spec_stack = spec_stack + as.numeric(crossprod(temp_stack$wt_sum, spec_block[temp_stack$ID - id_lo + 1L, ]))
      }else{
        spec_stack = spec_stack + as.numeric(crossprod(temp_stack$wt_sum, spec_mat[temp_stack$ID,]))
      }
    }
  }

  #below now works properly with interp=TRUE. This didn't matter for Robotham 2025 since we always generated SSPs on isochrone logAge and logZ values.
  if(is.null(logA)){
    bb_check = Iso[logZ %in% logZ_steps & logAge %in% logAge_steps & best==0L, which=TRUE]
  }else{
    bb_check = Iso[logZ %in% logZ_steps & logAge %in% logAge_steps & logA %in% logA_steps & best==0L, which=TRUE]
  }

  if(length(bb_check) > 0){
    for(logAge_j in seq_along(logAge_steps)){
      logAge_step = logAge_steps[logAge_j]
      for(logZ_k in seq_along(logZ_steps)){
        logZ_step = logZ_steps[logZ_k]
        if(is.null(logA)){
          interp_wt = logAge_wt[logAge_j]*logZ_wt[logZ_k]
          bb_subset = Iso[logZ==logZ_step & logAge==logAge_step & best==0L, which=TRUE] #faster!
          if(length(bb_subset > 0)){
            for(sel in bb_subset){
              #This version didn't seem to help, so removed to reduce complication installing ProGeny
              #.vec_add_cpp(spec_stack, .blackbody_simp(wave=wave_grid, Temp=Iso[sel,Teff], norm=Iso[sel,Lum]*IMFint[sel])*interp_wt)
              spec_stack = spec_stack + .blackbody_simp(wave=wave_grid, Temp=Iso[sel,Teff], norm=Iso[sel,Lum]*IMFint[sel])*interp_wt
            }
          }
        }else{
          for(logA_m in seq_along(logA_steps)){
            logA_step = logA_steps[logA_m]
            interp_wt = logAge_wt[logAge_j]*logZ_wt[logZ_k]*logA_wt[logA_m]
            bb_subset = Iso[logZ==logZ_step & logAge==logAge_step & logA==logA_step & best==i, which=TRUE] #faster!
            if(length(bb_subset > 0)){
              for(sel in bb_subset){
                #This version didn't seem to help, so removed to reduce complication installing ProGeny
                #.vec_add_cpp(spec_stack, .blackbody_simp(wave=wave_grid, Temp=Iso[sel,Teff], norm=Iso[sel,Lum]*IMFint[sel])*interp_wt)
                spec_stack = spec_stack + .blackbody_simp(wave=wave_grid, Temp=Iso[sel,Teff], norm=Iso[sel,Lum]*IMFint[sel])*interp_wt
              }
            }
          }
        }
      }
    }
  }

  return(invisible(spec_stack))
}

# progenySSP2Spec = function(speclib, logAge=8.4, logZ=0, Zsol=0.02){
#   interp = FALSE
#
#   SSP_logAge_steps = log10(speclib$Age)
#   if(min(logAge) < min(SSP_logAge_steps)){
#     stop('logAge minimum is less than minimum logAge present in the provided SSP!')
#   }
#   if(max(logAge) > max(SSP_logAge_steps)){
#     stop('logAge maximum is more than maximum logAge present in the provided SSP!')
#   }
#   interp = TRUE
#
#   SSP_logZ_steps = log10(speclib$Z/Zsol)
#   if(min(logZ) < min(SSP_logZ_steps)){
#     stop('logAge minimum is less than minimum logAge present in the provided SSP!')
#   }
#   if(max(logAge) > max(Iso_logAge_steps)){
#     stop('logAge maximum is more than maximum logAge present in the provided SSP!')
#   }
#   interp = TRUE
# }
