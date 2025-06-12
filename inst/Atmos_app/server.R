server = function(input, output, session) {
  ASGR_atmos_path = '/Volumes/Macintosh HD/Users/aaron/Google Drive/My Drive/ProGeny_atmos'
  if(file.exists(ASGR_atmos_path)){
    volumes = c(ASGR = ASGR_atmos_path, wd = getwd(), Home = '~/', shinyFiles::getVolumes()())
  }else{
    volumes = c(wd = getwd(), Home = '~/', shinyFiles::getVolumes()())
  }

  shinyFiles::shinyDirChoose(input, id = 'destpath', roots = volumes)

  atmos_result = reactiveVal(NULL)
  wave_grid_result = reactiveVal(NULL)
  iso_result = reactiveVal(NULL)
  interp_all_result = reactiveVal(NULL)
  SSP_result = reactiveVal(NULL)
  iso_info_result = reactiveVal(NULL)
  atmos_info_result = reactiveVal(NULL)
  interp_info_result = reactiveVal(NULL)
  #IMF_info_result = reactiveVal(NULL)

  observeEvent(input$destpath,{
    output$selectedPath <- renderText({
      req(input$destpath)
      paste(
        'Path to use: \n',
        shinyFiles::parseDirPath(volumes, input$destpath)
      )
    })
  })

  #output$status <- renderText("Waiting...")

  observeEvent(input$iso_file, {
    req(input$iso_file)
    shinybusy::show_spinner()
    tryCatch({
      iso_out = fst::read_fst(input$iso_file$datapath, as.data.table = TRUE)
      iso_result(iso_out)
      output$iso_summary <- renderPrint(summary(iso_out))
      output$iso_status <- renderText("Isochrone file loaded successfully!")
      output$plot_iso = renderPlot({
        progenyIsoPlot(iso_result())
      })
    }, error = function(e) {
      output$iso_status <- renderText(paste("Error loading isochrone file:", e$message))
    })
    shinybusy::hide_spinner()

    if(grepl('Mist', input$iso_file$name)){
      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_AGB',
        choices = -1:10,
        selected = c(4,5)
      )

      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_white',
        choices = -1:10,
        selected = c(6,7,8)
      )

      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_WR',
        choices = -1:10,
        selected = 9
      )

      Iso_info = c(Iso_type = 'MIST')
    }else if(grepl('Parsec', input$iso_file$name)){
      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_AGB',
        choices = -1:10,
        selected = c(7,8)
      )

      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_white',
        choices = -1:10,
        selected = 9
      )

      Iso_info = c(Iso_type = 'Parsec')
    }else if(grepl('Basti', input$iso_file$name)){
      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_AGB',
        choices = -1:10,
        selected = c(3,4,5)
      )

      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_white',
        choices = -1:10,
        selected = 6
      )

      Iso_info = c(Iso_type = 'BaSTI')
    }else if(grepl('Padova', input$iso_file$name)){
      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_AGB',
        choices = -1:10,
        selected = 5
      )

      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_white',
        choices = -1:10,
        selected = 6
      )

      updateSelectInput(
        session = getDefaultReactiveDomain(),
        inputId = 'label_WR',
        choices = -1:10,
        selected = 9
      )
      Iso_info = c(Iso_type = 'Padova')
    }
    iso_info_result(Iso_info)
  })

  # observeEvent(input$load_atmos, {
  #   output$atmos_status <- renderText("Loading atmosphere data, please be patient (~20 seconds)...")
  # }, priority = 10)

  observeEvent(input$load_atmos, {
    shinybusy::show_spinner()
    input_base = input$base
    input_extend = input$extend
    input_hot = input$hot
    input_AGB = input$AGB
    input_white = input$white
    input_WR = input$WR

    if(input_extend == 'None'){
      input_extend = NULL
    }

    if(input_hot == 'None'){
      input_hot = NULL
    }

    if(input_AGB == 'None'){
      input_AGB = NULL
    }

    if(input_white == 'None'){
      input_white = NULL
    }

    if(input_WR == 'None'){
      input_WR = NULL
    }

    destpath = shinyFiles::parseDirPath(volumes, input$destpath)

    if(!file.exists(paste0(destpath,'/',input_base,'.fits'))){
      output$atmos_status <- renderText(paste("No base atmosphere at:", paste0(destpath,'/',input_base,'.fits')))
      return()
    }

    if(!is.null(input_extend)){
      if(!file.exists(paste0(destpath,'/',input_extend,'.fits'))){
        output$atmos_status <- renderText(paste("No extend atmosphere at:", paste0(destpath,'/',input_extend,'.fits')))
        return()
      }
    }

    if(!is.null(input_hot)){
      if(!file.exists(paste0(destpath,'/',input_hot,'.fits'))){
        output$atmos_status <- renderText(paste("No hot atmosphere at:", paste0(destpath,'/',input$hot,'.fits')))
        return()
      }
    }

    if(!is.null(input_AGB)){
      if(!file.exists(paste0(destpath,'/',input_AGB,'.fits'))){
        output$atmos_status <- renderText(paste("No AGB atmosphere at:", paste0(destpath,'/',input_AGB,'.fits')))
        return()
      }
    }

    if(!is.null(input_white)){
      if(!file.exists(paste0(destpath,'/',input_white,'.fits'))){
        output$atmos_status <- renderText(paste("No white atmosphere at:", paste0(destpath,'/',input_white,'.fits')))
        return()
      }
    }

    if(!is.null(input_WR)){
      if(!file.exists(paste0(destpath,'/',input_WR,'.fits'))){
        output$atmos_status <- renderText(paste("No WR atmosphere at:", paste0(destpath,'/',input_WR,'.fits')))
        return()
      }
    }

    #shinyjs::delay(100, {
    tryCatch({
      atmos_out = progenyAtmosLoad(
        destpath = destpath,
        base = input_base,
        extend = input_extend,
        hot = input_hot,
        AGB = input_AGB,
        white = input_white,
        WR = input_WR,
        wavegrid = wave_grid_result(),
        cores = input$atmos_cores
      )
      atmos_result(atmos_out)

      output$atmos_status <- renderText("Atmosphere data loaded successfully!")

      output$plot_atmos = renderPlot({
        if(is.null(iso_result())){
          progenyAtmosPlot(atmos_out, add=FALSE)
        }else{
          progenyIsoPlot(iso_result())
          progenyAtmosPlot(atmos_out, add=TRUE)
        }
      })

      # output$loaded_wave_samp <- renderText({
      #   'Wavelength grid info below:'
      # })

      output$summary_wave_samp <- renderPrint({
        temp = summary(atmos_out$base$wave)
        temp_diff = summary(diff(atmos_out$base$wave))
        data.frame(Stat=names(temp), 'Wave.Ang'=as.numeric(temp), 'Bin.Ang'=as.numeric(temp_diff))
      })

      output$plot_wave_samp <- renderPlot({
        progenySampPlot(atmos_out$base$wave)
      })
    }, error = function(e) {
      output$atmos_status <- renderText(paste("Error:", e$message))
    })
    shinybusy::hide_spinner()

    atmos_info = c(
      atmos_base = input_base,
      atmos_extend = input_extend,
      atmos_hot = input_hot,
      atmos_AGB = input_AGB,
      atmos_white = input_white,
      atmos_WR = input_WR
    )

    atmos_info_result(atmos_info)
  })

  observeEvent(input$wave_file, {
    req(input$wave_file)
    tryCatch({
      wave_grid = as.numeric(read.table(input$wave_file$datapath, header = FALSE)[,1])
      wave_grid_result(wave_grid)
    }, error = function(e) {
      output$wave_status <- renderText(paste("Error loading wavegrid file:", e$message))
    })
  })

  # output$iso_tab <- renderUI({
  #   req(wave_grid_result())
  #   if(is.null(wave_grid_result())){
  #     tabPanel("Isochrone",
  #              verbatimTextOutput("iso_status"),
  #              plotOutput("plot_iso", height = "600px"),
  #              verbatimTextOutput("iso_summary")
  #     )
  #   } else {
  #     tabPanel(HTML("Isochrone <span style='color:green;'>&#10003;</span>"),
  #              verbatimTextOutput("iso_status"),
  #              plotOutput("plot_iso", height = "600px"),
  #              verbatimTextOutput("iso_summary")
  #     )
  #   }
  # })

  observeEvent(input$run_interp, {
    shinybusy::show_spinner()
    tryCatch({
      interp_all_out = progenyInterpGrid_All(
        Iso = iso_result(),
        Spec_combine = atmos_result(),
        radius = input$radius,
        weight_pow = input$weight_pow,
        k = input$k
      )

      interp_all_result(interp_all_out)

      new_iso = progenyInterpBest(
        Iso = iso_result(),
        Interp_combine = interp_all_out,
        do_hot = input$do_hot,
        do_AGB = input$do_AGB,
        do_white = input$do_white,
        do_WR = input$do_WR,
        prefer_hot = input$prefer_hot,
        prefer_AGB = input$prefer_AGB,
        prefer_white = input$prefer_white,
        prefer_WR = input$prefer_WR,
        label_AGB = input$label_AGB,
        label_white = input$label_white,
        label_WR = input$label_WR
      )

      iso_result(new_iso)
      output$interp_status <- renderText("Interpolation successfully run!")
    }, error = function(e) {
      output$interp_status <- renderText(paste("Error:", e$message))
    })
    shinybusy::hide_spinner()

    output$plot_iso_best_logZ = renderPlot({progenyIsoBestPlot(new_iso, xsel = 'logZ')})
    output$plot_iso_best_logG = renderPlot({progenyIsoBestPlot(new_iso, xsel = 'logG')})
    output$plot_iso_best_logT = renderPlot({progenyIsoBestPlot(new_iso, xsel = 'Teff', log='x')})

    interp_stat = progenyInterpStat(new_iso, atmos_result(), interp_all_out)

    output$interp_summary = renderPrint({interp_stat$stat})

    output$plot_atmos_error_logZ = renderPlot({progenyAtmosErrorPlot(interp_stat$Iso, xsel = 'logZ')})
    output$plot_atmos_error_logG = renderPlot({progenyAtmosErrorPlot(interp_stat$Iso, xsel = 'logG')})
    output$plot_atmos_error_logT = renderPlot({progenyAtmosErrorPlot(interp_stat$Iso, xsel = 'Teff', log='x')})

    interp_info = c(
      interp_radius = input$radius,
      interp_weight_pow = input$weight_pow,
      interp_k = input$k,
      interp_do_hot = input$do_hot,
      interp_do_AGB = input$do_AGB,
      interp_do_white = input$do_white,
      interp_do_WR = input$do_WR,
      interp_prefer_hot = input$prefer_hot,
      interp_prefer_AGB = input$prefer_AGB,
      interp_prefer_white = input$prefer_white,
      interp_prefer_WR = input$prefer_WR,
      interp_lab_AGB = paste(input$label_AGB, collapse=' / '),
      interp_lab_white =  paste(input$label_white, collapse=' / '),
      interp_lab_WR =  paste(input$label_WR, collapse=' / ')
    )

    interp_info_result(interp_info)
  })

  color_bool <- function(value) {
    if (value) {
      '<span style="color:green;">TRUE</span>'
    } else {
      '<span style="color:red;">FALSE</span>'
    }
  }

  output$used_imf <- renderUI({
    Iso_exits = !is.null(iso_result())
    Atmos_exist = !is.null(atmos_result())
    Interp_exist = !is.null(interp_all_result())
    SSP_exist = !is.null(SSP_result())

    can_run_SSP = Iso_exits & Atmos_exist & Interp_exist

    HTML(paste0('Isochrone loaded: ', color_bool(Iso_exits), '<br>',
          'Isochrone for SSP: ', input$iso_file$name, '<br><br>',
          'Atmospheres for SSP loaded: ', color_bool(Atmos_exist), '<br>',
          'See (Atmospheres) tab for parameter details.', '<br><br>',
          'Interpolated grid generated: ', color_bool(Interp_exist), '<br>',
          'See (Interpolate) tab for parameter details.', '<br><br>',
          'IMF for SSP: ', input$imf, '<br>',
          'See (IMF) tab for parameter details.', '<br><br>',
          'Ready to [Make SSP]: ', color_bool(can_run_SSP), '<br>',
          'Ready to [Check SSP]: ', color_bool(SSP_exist), '<br>',
          'Ready to [Return SSP]: ', color_bool(SSP_exist), '<br>',
          'Ready to [Download SSP]: ', color_bool(SSP_exist)
          ))
  })

  output$plot_imf <- renderPlot({
    req(input$imf)

    if (input$imf == 'IMF_Chabrier') {
      req(input$chab_masslow, input$chab_massmax,
          input$chab_alpha, input$chab_a, input$chab_b)

      magicaxis::magcurve(IMF_Chabrier(x,
                            alpha = input$chab_alpha,
                            a = input$chab_a,
                            b = input$chab_b,
                            masslow = input$chab_masslow,
                            massmax = input$chab_massmax),
               input$chab_masslow, input$chab_massmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2)
    } else if (input$imf == 'IMF_Kroupa') {
      req(input$kroupa_masslow, input$kroupa_massmax,
          input$kroupa_alpha1, input$kroupa_alpha2, input$kroupa_alpha3,
          input$kroupa_mass1, input$kroupa_mass2)

      magicaxis::magcurve(IMF_Kroupa(x,
                          alpha1 = input$kroupa_alpha1,
                          alpha2 = input$kroupa_alpha2,
                          alpha3 = input$kroupa_alpha3,
                          mass1 = input$kroupa_mass1,
                          mass2 = input$kroupa_mass2,
                          masslow = input$kroupa_masslow,
                          massmax = input$kroupa_massmax),
               input$kroupa_masslow, input$kroupa_massmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2)
    } else if (input$imf == 'IMF_Salpeter') {
      req(input$salp_masslow, input$salp_massmax,
          input$salp_alpha)

      magicaxis::magcurve(IMF_Salpeter(x,
                            alpha = input$salp_alpha,
                            masslow = input$salp_masslow,
                            massmax = input$salp_massmax),
               input$salp_masslow, input$salp_massmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2)
    } else if (input$imf == 'IMF_Kroupa_evo') {
      req(input$kroupa_masslow_lim, input$kroupa_massmax_lim, input$kroupa_Age_lim,
          input$kroupa_alpha1_lim, input$kroupa_alpha2_lim, input$kroupa_alpha3_lim,
          input$kroupa_mass1_lim, input$kroupa_mass2_lim,
          input$kroupa_Lookback_Age,
          !is.null(input$kroupa_alpha1_lim_Rv),
          !is.null(input$kroupa_alpha2_lim_Rv),
          !is.null(input$kroupa_alpha3_lim_Rv),
          !is.null(input$kroupa_masslow_lim_Rv),
          !is.null(input$kroupa_massmax_lim_Rv)
      )

      alpha1_lim = if(input$kroupa_alpha1_lim_Rv){rev(input$kroupa_alpha1_lim)}else{input$kroupa_alpha1_lim}
      alpha2_lim = if(input$kroupa_alpha2_lim_Rv){rev(input$kroupa_alpha2_lim)}else{input$kroupa_alpha2_lim}
      alpha3_lim = if(input$kroupa_alpha3_lim_Rv){rev(input$kroupa_alpha3_lim)}else{input$kroupa_alpha3_lim}
      masslow_lim = if(input$kroupa_masslow_lim_Rv){rev(input$kroupa_masslow_lim)}else{input$kroupa_masslow_lim}
      massmax_lim = if(input$kroupa_massmax_lim_Rv){rev(input$kroupa_massmax_lim)}else{input$kroupa_massmax_lim}
      xmin = min(masslow_lim)
      xmax = max(massmax_lim)

      temp_func = function(x, Age=0){
        IMF_Kroupa_evo(x,
                       Age = Age,
                       Age_lim = input$kroupa_Age_lim,
                       alpha1_lim = alpha1_lim,
                       alpha2_lim = alpha2_lim,
                       alpha3_lim = alpha3_lim,
                       mass1 = input$kroupa_mass1_lim,
                       mass2 = input$kroupa_mass2_lim,
                       masslow = masslow_lim,
                       massmax = massmax_lim,
                       Lookback_Age = input$kroupa_Lookback_Age)
      }
      magicaxis::magcurve(temp_func(x, Age=0), xmin, xmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2, xlim=c(xmin, xmax))
      curve(temp_func(x, Age=6.9), xmin, xmax, log='xy', lty=2, add=TRUE, lwd=2)
      curve(temp_func(x, Age=13.8), xmin, xmax, log='xy', lty=3, add=TRUE, lwd=2)
      legend('bottomleft', legend=c('0 Gyrs', '6.9 Gyrs', '13.8 Gyrs'), lty=1:3, lwd=2)
    } else if (input$imf == 'IMF_Kroupa_Zevo') {
      req(input$kroupa_masslow_lim, input$kroupa_massmax_lim, input$kroupa_logZ_lim,
          input$kroupa_alpha1_lim, input$kroupa_alpha2_lim, input$kroupa_alpha3_lim,
          input$kroupa_mass1_lim, input$kroupa_mass2_lim,
          !is.null(input$kroupa_alpha1_lim_Rv),
          !is.null(input$kroupa_alpha2_lim_Rv),
          !is.null(input$kroupa_alpha3_lim_Rv),
          !is.null(input$kroupa_masslow_lim_Rv),
          !is.null(input$kroupa_massmax_lim_Rv)
        )

      alpha1_lim = if(input$kroupa_alpha1_lim_Rv){rev(input$kroupa_alpha1_lim)}else{input$kroupa_alpha1_lim}
      alpha2_lim = if(input$kroupa_alpha2_lim_Rv){rev(input$kroupa_alpha2_lim)}else{input$kroupa_alpha2_lim}
      alpha3_lim = if(input$kroupa_alpha3_lim_Rv){rev(input$kroupa_alpha3_lim)}else{input$kroupa_alpha3_lim}
      masslow_lim = if(input$kroupa_masslow_lim_Rv){rev(input$kroupa_masslow_lim)}else{input$kroupa_masslow_lim}
      massmax_lim = if(input$kroupa_massmax_lim_Rv){rev(input$kroupa_massmax_lim)}else{input$kroupa_massmax_lim}
      xmin = min(masslow_lim)
      xmax = max(massmax_lim)

      temp_func = function(x, logZ=0){
        IMF_Kroupa_Zevo(x,
                        logZ = logZ,
                        logZ_lim = input$kroupa_logZ_lim,
                        alpha1_lim = alpha1_lim,
                        alpha2_lim = alpha2_lim,
                        alpha3_lim = alpha3_lim,
                        mass1 = input$kroupa_mass1_lim,
                        mass2 = input$kroupa_mass2_lim,
                        masslow = masslow_lim,
                        massmax = massmax_lim
        )
      }
      magicaxis::magcurve(temp_func(x, logZ=0), xmin, xmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2, xlim=c(xmin, xmax))
      curve(temp_func(x, logZ=-2), xmin, xmax, log='xy', lty=2, add=TRUE, lwd=2)
      curve(temp_func(x, logZ=-4), xmin, xmax, log='xy', lty=3, add=TRUE, lwd=2)
      legend('bottomleft', legend=c('logZ: 0', 'logZ: -2', 'logZ: -4'), lty=1:3, lwd=2)
    } else if (input$imf == 'IMF_Lacey_evo') {
      req(input$lacey_masslow_lim, input$lacey_massmax_lim, input$lacey_Age_lim,
          input$lacey_alpha1_lim, input$lacey_alpha2_lim, input$lacey_alpha3_lim,
          input$lacey_mass1_lim, input$lacey_mass2_lim,
          input$lacey_Lookback_Age,
          !is.null(input$lacey_alpha1_lim_Rv),
          !is.null(input$lacey_alpha2_lim_Rv),
          !is.null(input$lacey_alpha3_lim_Rv),
          !is.null(input$lacey_masslow_lim_Rv),
          !is.null(input$lacey_massmax_lim_Rv)
      )

      alpha1_lim = if(input$lacey_alpha1_lim_Rv){rev(input$lacey_alpha1_lim)}else{input$lacey_alpha1_lim}
      alpha2_lim = if(input$lacey_alpha2_lim_Rv){rev(input$lacey_alpha2_lim)}else{input$lacey_alpha2_lim}
      alpha3_lim = if(input$lacey_alpha3_lim_Rv){rev(input$lacey_alpha3_lim)}else{input$lacey_alpha3_lim}
      masslow_lim = if(input$lacey_masslow_lim_Rv){rev(input$lacey_masslow_lim)}else{input$lacey_masslow_lim}
      massmax_lim = if(input$lacey_massmax_lim_Rv){rev(input$lacey_massmax_lim)}else{input$lacey_massmax_lim}
      xmin = min(masslow_lim)
      xmax = max(massmax_lim)

      temp_func = function(x, Age=0){
        IMF_Lacey_evo(x,
                      Age = Age,
                      Age_lim = input$lacey_Age_lim,
                      alpha1_lim = alpha1_lim,
                      alpha2_lim = alpha2_lim,
                      alpha3_lim = alpha3_lim,
                      mass1 = input$lacey_mass1_lim,
                      mass2 = input$lacey_mass2_lim,
                      masslow = masslow_lim,
                      massmax = massmax_lim,
                      Lookback_Age = input$lacey_Lookback_Age)
      }
      magicaxis::magcurve(temp_func(x, Age=0), xmin, xmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2, xlim=c(xmin, xmax))
      curve(temp_func(x, Age=6.9), xmin, xmax, log='xy', lty=2, add=TRUE, lwd=2)
      curve(temp_func(x, Age=13.8), xmin, xmax, log='xy', lty=3, add=TRUE, lwd=2)
      legend('bottomleft', legend=c('0 Gyrs', '6.9 Gyrs', '13.8 Gyrs'), lty=1:3, lwd=2, title='Lookback')
    } else if (input$imf == 'IMF_Lacey_Zevo') {

      req(input$lacey_masslow_lim, input$lacey_massmax_lim, input$lacey_logZ_lim,
          input$lacey_alpha1_lim, input$lacey_alpha2_lim, input$lacey_alpha3_lim,
          input$lacey_mass1_lim, input$lacey_mass2_lim,
          !is.null(input$lacey_alpha1_lim_Rv),
          !is.null(input$lacey_alpha2_lim_Rv),
          !is.null(input$lacey_alpha3_lim_Rv),
          !is.null(input$lacey_masslow_lim_Rv),
          !is.null(input$lacey_massmax_lim_Rv)
      )

      alpha1_lim = if(input$lacey_alpha1_lim_Rv){rev(input$lacey_alpha1_lim)}else{input$lacey_alpha1_lim}
      alpha2_lim = if(input$lacey_alpha2_lim_Rv){rev(input$lacey_alpha2_lim)}else{input$lacey_alpha2_lim}
      alpha3_lim = if(input$lacey_alpha3_lim_Rv){rev(input$lacey_alpha3_lim)}else{input$lacey_alpha3_lim}
      masslow_lim = if(input$lacey_masslow_lim_Rv){rev(input$lacey_masslow_lim)}else{input$lacey_masslow_lim}
      massmax_lim = if(input$lacey_massmax_lim_Rv){rev(input$lacey_massmax_lim)}else{input$lacey_massmax_lim}
      xmin = min(masslow_lim)
      xmax = max(massmax_lim)

      temp_func = function(x, logZ=0){
        IMF_Lacey_Zevo(x,
                       logZ = logZ,
                       logZ_lim = input$lacey_logZ_lim,
                       alpha1_lim = input$lacey_alpha1_lim,
                       alpha2_lim = input$lacey_alpha2_lim,
                       alpha3_lim = input$lacey_alpha3_lim,
                       mass1 = input$lacey_mass1_lim,
                       mass2 = input$lacey_mass2_lim,
                       masslow = masslow_lim,
                       massmax = massmax_lim
        )
      }
      magicaxis::magcurve(temp_func(x, logZ=0), xmin, xmax, log = 'xy',
               xlab = 'Star Mass / Msol', ylab = 'dN/dM (1 Msol)', lwd=2, xlim=c(xmin, xmax))
      curve(temp_func(x, logZ=-2), xmin, xmax, log='xy', lty=2, add=TRUE, lwd=2)
      curve(temp_func(x, logZ=-4), xmin, xmax, log='xy', lty=3, add=TRUE, lwd=2)
      legend('bottomleft', legend=c('logZ: 0', 'logZ: -2', 'logZ: -4'), lty=1:3, lwd=2)
    }
    legend('topright', legend = c('User IMF', 'Ref Chab [0.1 - 150 Msol]'), col=c('black', 'darkgreen'), lty=c(1,3), lwd=2)
    curve(IMF_Chabrier(x, masslow=0.1, massmax=150), 0.1, 150, add=TRUE, col='darkgreen', lty=3, lwd=2)
  })

  # observeEvent(input$make_ssp, {
  #   output$SSP_status <- renderUI(HTML('<span style="color:purple;">SSP running, please be patient (~one minute)...</span>'))
  # }, priority = 10)

  observeEvent(input$make_ssp, {
    shinybusy::show_spinner()
    tryCatch({
      iso_data <- iso_result()
      atmos_data <- atmos_result()
      interp_data <- interp_all_result()
      IMF_function <- match.fun(input$imf)

      if(is.null(iso_data)){
        message('Missing Isochrone!')
      }

      if(is.null(atmos_data)){
        message('Missing Atmospheres!')
      }

      if(is.null(interp_data)){
        message('Missing Interpolation Grid!')
      }

      if(is.null(IMF_function)){
        message('Missing IMF!')
      }

      req(iso_data, atmos_data, interp_data, IMF_function)

      if(input$imf == 'IMF_Chabrier'){
        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          alpha = input$chab_alpha,
          a = input$chab_a,
          b = input$chab_b,
          masslow = input$chab_masslow,
          massmax = input$chab_massmax,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Chabrier',
          IMF_alpha = paste(input$chab_alpha, collapse = ' - '),
          IMF_a = paste(input$chab_a, collapse = ' - '),
          IMF_b = paste(input$chab_b, collapse = ' - '),
          IMF_masslow = paste(input$chab_masslow, collapse = ' - '),
          IMF_massmax = paste(input$chab_massmax, collapse = ' - ')
        )
      }else if(input$imf == 'IMF_Kroupa'){
        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          alpha1 = input$kroupa_alpha1,
          alpha2 = input$kroupa_alpha2,
          alpha3 = input$kroupa_alpha3,
          mass1 = input$kroupa_mass1,
          mass2 = input$kroupa_mass2,
          masslow = input$kroupa_masslow,
          massmax = input$kroupa_massmax,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Kroupa',
          IMF_alpha1 = paste(input$kroupa_alpha1, collapse = ' - '),
          IMF_alpha2 = paste(input$kroupa_alpha2, collapse = ' - '),
          IMF_alpha3 = paste(input$kroupa_alpha3, collapse = ' - '),
          IMF_mass1 = paste(input$kroupa_mass1, collapse = ' - '),
          IMF_mass2 = paste(input$kroupa_mass2, collapse = ' - '),
          IMF_masslow = paste(input$kroupa_masslow, collapse = ' - '),
          IMF_massmax = paste(input$kroupa_massmax, collapse = ' - ')
        )
      }else if(input$imf == 'IMF_Salpeter'){
        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          alpha = input$salp_alpha,
          masslow = input$salp_masslow,
          massmax = input$salp_massmax,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Salpeter',
          IMF_alpha = paste(input$salp_alpha, collapse = ' - '),
          IMF_masslow = paste(input$salp_mass1, collapse = ' - '),
          IMF_massmax = paste(input$salp_massmax, collapse = ' - ')
        )
      }else if(input$imf == 'IMF_Kroupa_evo'){
        alpha1_lim = if(input$kroupa_alpha1_lim_Rv){rev(input$kroupa_alpha1_lim)}else{input$kroupa_alpha1_lim}
        alpha2_lim = if(input$kroupa_alpha2_lim_Rv){rev(input$kroupa_alpha2_lim)}else{input$kroupa_alpha2_lim}
        alpha3_lim = if(input$kroupa_alpha3_lim_Rv){rev(input$kroupa_alpha3_lim)}else{input$kroupa_alpha3_lim}
        masslow_lim = if(input$kroupa_masslow_lim_Rv){rev(input$kroupa_masslow_lim)}else{input$kroupa_masslow_lim}
        massmax_lim = if(input$kroupa_massmax_lim_Rv){rev(input$kroupa_massmax_lim)}else{input$kroupa_massmax_lim}

        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          Age_lim = input$kroupa_Age_lim,
          alpha1_lim = alpha1_lim,
          alpha2_lim = alpha2_lim,
          alpha3_lim = alpha3_lim,
          mass1_lim = input$kroupa_mass1_lim,
          mass2_lim = input$kroupa_mass2_lim,
          masslow_lim = masslow_lim,
          massmax_lim = massmax_lim,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Kroupa_evo',
          IMF_age_lim = paste(input$kroupa_Age_lim, collapse = ' - '),
          IMF_alpha1_lim = paste(alpha1_lim, collapse = ' - '),
          IMF_alpha2_lim = paste(alpha2_lim, collapse = ' - '),
          IMF_alpha3_lim = paste(alpha3_lim, collapse = ' - '),
          IMF_mass1_lim = paste(input$kroupa_mass1_lim, collapse = ' - '),
          IMF_mass2_lim = paste(input$kroupa_mass2_lim, collapse = ' - '),
          IMF_masslow_lim = paste(masslow_lim, collapse = ' - '),
          IMF_massmax_lim = paste(massmax_lim, collapse = ' - ')
        )
      }else if(input$imf == 'IMF_Kroupa_Zevo'){
        alpha1_lim = if(input$kroupa_alpha1_lim_Rv){rev(input$kroupa_alpha1_lim)}else{input$kroupa_alpha1_lim}
        alpha2_lim = if(input$kroupa_alpha2_lim_Rv){rev(input$kroupa_alpha2_lim)}else{input$kroupa_alpha2_lim}
        alpha3_lim = if(input$kroupa_alpha3_lim_Rv){rev(input$kroupa_alpha3_lim)}else{input$kroupa_alpha3_lim}
        masslow_lim = if(input$kroupa_masslow_lim_Rv){rev(input$kroupa_masslow_lim)}else{input$kroupa_masslow_lim}
        massmax_lim = if(input$kroupa_massmax_lim_Rv){rev(input$kroupa_massmax_lim)}else{input$kroupa_massmax_lim}

        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          logZ_lim = input$kroupa_logZ_lim,
          alpha1_lim = alpha1_lim,
          alpha2_lim = alpha2_lim,
          alpha3_lim = alpha3_lim,
          mass1_lim = input$kroupa_mass1_lim,
          mass2_lim = input$kroupa_mass2_lim,
          masslow_lim = masslow_lim,
          massmax_lim = massmax_lim,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Kroupa_Zevo',
          IMF_Z_lim = paste(input$kroupa_logZ_lim, collapse = ' - '),
          IMF_alpha1_lim = paste(alpha1_lim, collapse = ' - '),
          IMF_alpha2_lim = paste(alpha2_lim, collapse = ' - '),
          IMF_alpha3_lim = paste(alpha3_lim, collapse = ' - '),
          IMF_mass1_lim = paste(input$kroupa_mass1_lim, collapse = ' - '),
          IMF_mass2_lim = paste(input$kroupa_mass2_lim, collapse = ' - '),
          IMF_masslow_lim = paste(masslow_lim, collapse = ' - '),
          IMF_massmax_lim = paste(massmax_lim, collapse = ' - ')
        )
      }else if(input$imf == 'IMF_Lacey_evo'){
        alpha1_lim = if(input$lacey_alpha1_lim_Rv){rev(input$lacey_alpha1_lim)}else{input$lacey_alpha1_lim}
        alpha2_lim = if(input$lacey_alpha2_lim_Rv){rev(input$lacey_alpha2_lim)}else{input$lacey_alpha2_lim}
        alpha3_lim = if(input$lacey_alpha3_lim_Rv){rev(input$lacey_alpha3_lim)}else{input$lacey_alpha3_lim}
        masslow_lim = if(input$lacey_masslow_lim_Rv){rev(input$lacey_masslow_lim)}else{input$lacey_masslow_lim}
        massmax_lim = if(input$lacey_massmax_lim_Rv){rev(input$lacey_massmax_lim)}else{input$lacey_massmax_lim}

        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          Age_lim = input$lacey_Age_lim,
          alpha1_lim = alpha1_lim,
          alpha2_lim = alpha2_lim,
          alpha3_lim = alpha3_lim,
          mass1_lim = input$lacey_mass1_lim,
          mass2_lim = input$lacey_mass2_lim,
          masslow_lim = masslow_lim,
          massmax_lim = massmax_lim,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Lacey_evo',
          IMF_age_lim = paste(input$lacey_Age_lim, collapse = ' - '),
          IMF_alpha1_lim = paste(alpha1_lim, collapse = ' - '),
          IMF_alpha2_lim = paste(alpha2_lim, collapse = ' - '),
          IMF_alpha3_lim = paste(alpha3_lim, collapse = ' - '),
          IMF_mass1_lim = paste(input$lacey_mass1_lim, collapse = ' - '),
          IMF_mass2_lim = paste(input$lacey_mass2_lim, collapse = ' - '),
          IMF_masslow_lim = paste(masslow_lim, collapse = ' - '),
          IMF_massmax_lim = paste(massmax_lim, collapse = ' - ')
        )
      }else if(input$imf == 'IMF_Lacey_Zevo'){
        alpha1_lim = if(input$lacey_alpha1_lim_Rv){rev(input$lacey_alpha1_lim)}else{input$lacey_alpha1_lim}
        alpha2_lim = if(input$lacey_alpha2_lim_Rv){rev(input$lacey_alpha2_lim)}else{input$lacey_alpha2_lim}
        alpha3_lim = if(input$lacey_alpha3_lim_Rv){rev(input$lacey_alpha3_lim)}else{input$lacey_alpha3_lim}
        masslow_lim = if(input$lacey_masslow_lim_Rv){rev(input$lacey_masslow_lim)}else{input$lacey_masslow_lim}
        massmax_lim = if(input$lacey_massmax_lim_Rv){rev(input$lacey_massmax_lim)}else{input$lacey_massmax_lim}

        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          logZ_lim = input$lacey_logZ_lim,
          alpha1_lim = alpha1_lim,
          alpha2_lim = alpha2_lim,
          alpha3_lim = alpha3_lim,
          mass1_lim = input$lacey_mass1_lim,
          mass2_lim = input$lacey_mass2_lim,
          masslow_lim = masslow_lim,
          massmax_lim = massmax_lim,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )

        IMF_info = c(
          IMF_type = 'Lacey_Zevo',
          IMF_Z_lim = paste(input$lacey_logZ_lim, collapse = ' - '),
          IMF_alpha1_lim = paste(alpha1_lim, collapse = ' - '),
          IMF_alpha2_lim = paste(alpha2_lim, collapse = ' - '),
          IMF_alpha3_lim = paste(alpha3_lim, collapse = ' - '),
          IMF_mass1_lim = paste(input$lacey_mass1_lim, collapse = ' - '),
          IMF_mass2_lim = paste(input$lacey_mass2_lim, collapse = ' - '),
          IMF_masslow_lim = paste(masslow_lim, collapse = ' - '),
          IMF_massmax_lim = paste(massmax_lim, collapse = ' - ')
        )
      }

      SSP_out$PG_info = data.frame(
        name = c(
          names(iso_info_result()),
          names(atmos_info_result()),
          names(interp_info_result()),
          names(IMF_info)
        ),
        value = as.character(c(
          iso_info_result(),
          atmos_info_result(),
          interp_info_result(),
          IMF_info
        )),
        row.names = NULL
      )

      output$PG_info = renderPrint({SSP_out$PG_info})

      SSP_result(SSP_out)
      output$SSP_status <- renderUI(HTML('<span style="color:green;">SSP successfully run!</span>'))
    }, error = function(e) {
      output$SSP_status <- renderUI(HTML(paste('<span style="color:red;">Error:', e$message,'</span>')))
    })
    shinybusy::hide_spinner()
  })

  output$dynamic_imf <- renderUI({
    req(input$imf)
    if(input$imf == 'IMF_Chabrier'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput('chab_masslow', 'Low Mass IMF Cutoff', min=0.01, max=0.2, value=0.1, step=0.01),
        sliderInput('chab_massmax', 'High Mass IMF Cutoff', min=10, max=500, value=150, step=10),
        sliderInput('chab_alpha', 'alpha', min=0, max=4, value=2.3, step=0.1),
        sliderInput('chab_a', 'a', min=0, max=1, value=0.08, step=0.01),
        sliderInput('chab_b', 'b', min=0, max=1, value=0.69, step=0.01)
      )
    } else if(input$imf == 'IMF_Kroupa'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput('kroupa_masslow', 'Low Mass IMF Cutoff', min=0.01, max=0.2, value=0.1, step=0.01),
        sliderInput('kroupa_massmax', 'High Mass IMF Cutoff', min=10, max=500, value=150, step=10),
        sliderInput('kroupa_alpha1', 'alpha_1', min=0, max=4, value=0.3, step=0.1),
        sliderInput('kroupa_alpha2', 'alpha_2', min=0, max=4, value=1.3, step=0.1),
        sliderInput('kroupa_alpha3', 'alpha_3', min=0, max=4, value=2.3, step=0.1),
        sliderInput('kroupa_mass1', 'mass_1', min=0.01, max=0.2, value=0.08, step=0.01),
        sliderInput('kroupa_mass2', 'mass_2', min=0.2, max=2, value=0.5, step=0.1)
      )
    }else if(input$imf == 'IMF_Salpeter'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput('salp_masslow', 'Low Mass IMF Cutoff', min=0.01, max=0.2, value=0.1, step=0.01),
        sliderInput('salp_massmax', 'High Mass IMF Cutoff', min=10, max=500, value=150, step=10),
        sliderInput('salp_alpha', 'alpha', min=0, max=4, value=2.35, step=0.1)
      )
    }else if (input$imf == 'IMF_Kroupa_evo'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput("kroupa_Age_lim", "Age Lim (Gyr)", min = 0, max = 13.8, value = c(0,13.8), step=0.1),
        fluidRow(
          column(9, sliderInput("kroupa_masslow_lim", "Mass Low", min=0.01, max = 0.2, value = c(0.1, 0.1), step = 0.01)),
          column(3, checkboxInput('kroupa_masslow_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_massmax_lim", "Mass Max", min = 10, max = 500, value = c(150,150), step = 10)),
          column(3, checkboxInput('kroupa_massmax_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_alpha1_lim", "Alpha 1", min = 0, max = 4, value = c(0.3, 0.3), step=0.1)),
          column(3, checkboxInput('kroupa_alpha1_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_alpha2_lim", "Alpha 2", min = 0, max = 4, value = c(1.3, 1.3), step=0.1)),
          column(3, checkboxInput('kroupa_alpha2_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_alpha3_lim", "Alpha 3", min = 0, max = 4, value = c(2, 2.3), step=0.1)),
          column(3, checkboxInput('kroupa_alpha3_lim_Rv', 'Rv', value = TRUE))
        ),
        sliderInput("kroupa_mass1_lim", "Mass 1", min = 0.01, max = 0.2, value = c(0.08,0.08), step = 0.01),
        sliderInput("kroupa_mass2_lim", "Mass 2", min = 0.2, max = 2, value = c(0.5,0.5), step = 0.1),
        sliderInput("kroupa_Lookback_Age", "Lookback Age (Gyr)", min = 0, max = 13.8, value = 0, step=0.1)
      )
    }else if (input$imf == 'IMF_Kroupa_Zevo'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput("kroupa_logZ_lim", "log(Z/Zsol) Lim", min = -4, max = 1, value = c(-4,0), step=0.1),
        fluidRow(
          column(9, sliderInput("kroupa_masslow_lim", "Mass Low", min=0.01, max = 0.2, value = c(0.1, 0.1), step = 0.01)),
          column(3, checkboxInput('kroupa_masslow_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_massmax_lim", "Mass Max", min = 10, max = 500, value = c(150,150), step = 10)),
          column(3, checkboxInput('kroupa_massmax_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_alpha1_lim", "Alpha 1", min = 0, max = 4, value = c(0.3, 0.3), step=0.1)),
          column(3, checkboxInput('kroupa_alpha1_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_alpha2_lim", "Alpha 2", min = 0, max = 4, value = c(1.3, 1.3), step=0.1)),
          column(3, checkboxInput('kroupa_alpha2_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("kroupa_alpha3_lim", "Alpha 3", min = 0, max = 4, value = c(2, 2.3), step=0.1)),
          column(3, checkboxInput('kroupa_alpha3_lim_Rv', 'Rv', value = FALSE))
        ),
        sliderInput("kroupa_mass1_lim", "Mass 1", min = 0.01, max = 0.2, value = c(0.08,0.08), step = 0.01),
        sliderInput("kroupa_mass2_lim", "Mass 2", min = 0.2, max = 2, value = c(0.5,0.5), step = 0.1)
      )
    }else if (input$imf == 'IMF_Lacey_evo'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput("lacey_Age_lim", "Age Lim (Gyr)", min = 0, max = 13.8, value = c(0,13.8), step=0.1),
        fluidRow(
          column(9, sliderInput("lacey_masslow_lim", "Mass Low", min=0.01, max = 0.2, value = c(0.1, 0.1), step = 0.01)),
          column(3, checkboxInput('lacey_masslow_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_massmax_lim", "Mass Max", min = 10, max = 500, value = c(150,150), step = 10)),
          column(3, checkboxInput('lacey_massmax_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_alpha1_lim", "Alpha 1", min = 0, max = 4, value = c(0.4, 2), step=0.1)),
          column(3, checkboxInput('lacey_alpha1_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_alpha2_lim", "Alpha 2", min = 0, max = 4, value = c(0.4, 2), step=0.1)),
          column(3, checkboxInput('lacey_alpha2_lim_Rv', 'Rv', value = FALSE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_alpha3_lim", "Alpha 3", min = 0, max = 4, value = c(2, 2.3), step=0.1)),
          column(3, checkboxInput('lacey_alpha3_lim_Rv', 'Rv', value = TRUE))
        ),
        sliderInput("lacey_mass1_lim", "Mass 1", min = 0.01, max = 0.2, value = c(0.08,0.08), step = 0.01),
        sliderInput("lacey_mass2_lim", "Mass 2", min = 0.2, max = 2, value = c(0.5,0.5), step = 0.1),
        sliderInput("lacey_Lookback_Age", "Lookback Age (Gyr)", min = 0, max = 13.8, value = 0, step=0.1)

      )
    }else if (input$imf == 'IMF_Lacey_Zevo'){
      output$SSP_status <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP]</span>'))
      #output$SSP_status <- renderText("IMF updated, need to run [Make SSP]")
      output$SSP_check <- renderUI(HTML('<span style="color:orange;">IMF updated, need to run [Make SSP] and then [Check SSP]</span>'))
      #output$SSP_check <- renderText("IMF updated, need to run [Make SSP] and then [Check SSP]")
      tagList(
        sliderInput("lacey_logZ_lim", "log(Z/Zsol) Lim", min = -4, max = 1, value = c(-4,0), step=0.1),
        fluidRow(
          column(9, sliderInput("lacey_masslow_lim", "Mass Low", min=0.01, max = 0.2, value = c(0.1, 0.1), step = 0.01)),
          column(3, checkboxInput('lacey_masslow_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_massmax_lim", "Mass Max", min = 10, max = 500, value = c(150,150), step = 10)),
          column(3, checkboxInput('lacey_massmax_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_alpha1_lim", "Alpha 1", min = 0, max = 4, value = c(0.4, 2), step=0.1)),
          column(3, checkboxInput('lacey_alpha1_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_alpha2_lim", "Alpha 2", min = 0, max = 4, value = c(0.4, 2), step=0.1)),
          column(3, checkboxInput('lacey_alpha2_lim_Rv', 'Rv', value = TRUE))
        ),
        fluidRow(
          column(9, sliderInput("lacey_alpha3_lim", "Alpha 3", min = 0, max = 4, value = c(2, 2.3), step=0.1)),
          column(3, checkboxInput('lacey_alpha3_lim_Rv', 'Rv', value = FALSE))
        ),
        sliderInput("lacey_mass1_lim", "Mass 1", min = 0.01, max = 0.2, value = c(0.08,0.08), step = 0.01),
        sliderInput("lacey_mass2_lim", "Mass 2", min = 0.2, max = 2, value = c(0.5,0.5), step = 0.1)
      )
    }
  })


  observeEvent(input$check_ssp, {
    output$SSP_check <- renderUI({
      HTML(
        if(is.null(SSP_result())){
          '<span style="color:red;">No SSP has been generated - check failed!</span>'
        }else{
          messages <- character()
          withCallingHandlers(
            {
              ProSpect::speclib_check(SSP_result()[1:8])
            },
            message = function(m) {
              messages <<- c(messages, conditionMessage(m))
              invokeRestart("muffleMessage")
            }
          )
          paste(messages, collapse = "<br>")
        }
      )
    })
  })


  observeEvent(input$iso_done, {
    stopApp(iso_result())
  })

  observeEvent(input$atmos_done, {
    stopApp(atmos_result())
  })

  observeEvent(input$interp_done, {
    stopApp(interp_all_result())
  })

  observeEvent(input$return_ssp, {
    if(is.null(SSP_result())){
      output$SSP_return <- renderUI({
        renderUI(HTML('<span style="color:orange;">No SSP has been generated - nothing to return!</span>'))
      })
    }else{
      stopApp(SSP_result())
    }
  })


  output$download_ssp <- downloadHandler(
    filename = paste0("PG_SSP_", format(Sys.time(), "%y_%m_%d_%H_%M_%S"), ".fits"),
    content = function(file) {
      if(!is.null(SSP_result())){
        Rfits::Rfits_write_all(data = SSP_result(), filename = file, flatten = TRUE, compress = FALSE)
      }else{
        showNotification("SSP data is missing or invalid. Cannot generate FITS file.", type = "error")
        stop("SSP data is missing or invalid.")

      }
    },
    contentType = 'FITS'
  )

}
