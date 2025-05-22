ui = fluidPage(
  titlePanel("ProGeny SSP Generator"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tabs == 'Isochrones'",
        fileInput("iso_file", "Choose Isochrone File (.fst)", accept = ".fst"),
        br(), br(),
        actionButton("iso_done", "Return Isochrone"),
        verbatimTextOutput("iso_summary")
      ),

      conditionalPanel(
        condition = "input.tabs == 'Atmospheres'",
        shinyDirButton("destpath", "Atmos Path", "Select a folder"),
        selectInput("base", "Base", choices = c("C3K" = "combine_C3K_conroy", "Husser" = "combine_PHOENIX_husser", "Allard" = "combine_PHOENIX_Allard", "MILES" = "combine_MILES_vazdekis", "BASEL" = "combine_BASEL_wlbc")),
        selectInput("extend", "Extend", choices = c("Allard" = "combine_PHOENIX_Allard", "ATLAS9" = "combine_ATLAS9_castelli")),
        selectInput("hot", "Hot", choices = c("OB" = "combine_OB_PoWR", "hot" = "combine_hot")),
        selectInput("AGB", "AGB", choices = "combine_AGB_lancon"),
        selectInput("white", "White Dwarfs", choices = c("TMAP" = "combine_TMAP_werner", "white" = "combvine_white")),
        selectInput("WR", "Wolf-Rayet", choices = c("PoWR" = "combine_WNE_PoWR")),
        numericInput("atmos_cores", "Number of Cores", value = 8, min = 1, step=1),
        actionButton("load_atmos", "Load Atmospheres"),
        br(), br(),
        actionButton("atmos_done", "Return Atmospheres")
      ),

      conditionalPanel(
        condition = "input.tabs == 'Interpolate'",
        sliderInput('radius', 'Search Radius', min=1, max=4, value=2),
        sliderInput('weight_pow', 'Weight Power', min=1, max=4, value=2),
        sliderInput('k', 'Search k', min=4, max=16, value=8, step=1),
        checkboxInput('do_hot', 'Use Hot', value = TRUE),
        checkboxInput('do_AGB', 'Use AGB', value = TRUE),
        checkboxInput('do_white', 'Use White Dwarf', value = TRUE),
        checkboxInput('do_WR', 'Use Wolf-Rayet', value = TRUE),
        selectInput('label_AGB', 'AGB Phase', choices=-1:10, multiple=TRUE),
        selectInput('label_white', 'White Dwarf Phase', choices=-1:10, multiple=TRUE),
        selectInput('label_WR', 'AGB Wolf-Rayet', choices=-1:10, multiple=TRUE),
        actionButton("run_interp", "Run Interpolation"),
        br(), br(),
        actionButton("interp_done", "Return Interpolation Grids"),
      ),

      conditionalPanel(
        condition = "input.tabs == 'Make SSP'",
        selectInput("imf", "IMF", choices = c("Chabrier" = "IMF_Chabrier",
                                              "Kroupa" = "IMF_Kroupa",
                                              "Salpeter" = "IMF_Salpeter")),
        actionButton("update_imf", "Change IMF Type"),
        sliderInput('masslow', 'Low Mass IMF Cutoff', min=0.01, max=0.2, value=0.1, step=0.01),
        sliderInput('massmax', 'High Mass IMF Cutoff', min=10, max=500, value=150, step=10),

        tags$h4("IMF parameters:"),
        uiOutput("dynamic_imf"),

        numericInput("SSP_cores", "Number of Cores", value = 8, min = 1, step=1),
        actionButton("make_ssp", "Make SSP"),
        br(), br(),
        actionButton("check_ssp", "Check SSP"),
        br(), br(),
        actionButton("SSP_done", "Return SSP"),
      )
    ),

    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Isochrones",
                           verbatimTextOutput("iso_status"),
                           plotOutput("plot_iso")
                  ),
                  tabPanel("Atmospheres",
                           verbatimTextOutput("selectedPath"),
                           verbatimTextOutput("atmos_status"),
                           plotOutput("plot_atmos")
                  ),
                  tabPanel("Interpolate",
                           verbatimTextOutput("interp_status")
                  ),
                  tabPanel("Make SSP",
                           plotOutput("plot_imf"),
                           verbatimTextOutput("SSP_status"),
                           verbatimTextOutput("SSP_check")
                  )
      )
    )
  )
)

server = function(input, output, session) {
  ASGR_atmos_path = '/Volumes/Macintosh HD/Users/aaron/Google Drive/My Drive/ProGeny_atmos'
  if(file.exists(ASGR_atmos_path)){
    volumes = c(ASGR = ASGR_atmos_path, wd = getwd(), Home = '~/', getVolumes()())
  }else{
    volumes = c(wd = getwd(), Home = '~/', getVolumes()())
  }

  shinyDirChoose(input, id = 'destpath', roots = volumes)

  atmos_result = reactiveVal(NULL)
  iso_result = reactiveVal(NULL)
  interp_all_result = reactiveVal(NULL)
  SSP_result = reactiveVal(NULL)

  output$selectedPath <- renderText({
    req(input$destpath)
    parseDirPath(volumes, input$destpath)
  })

  output$status <- renderText("Waiting...")

  observeEvent(input$iso_file, {
    req(input$iso_file)
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

    print(input$iso_file$name)

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
    }
  })

  observeEvent(input$load_atmos, {
    output$atmos_status <- renderText("Loading atmosphere data!")
    #shinyjs::delay(100, {
      tryCatch({
        atmos_out = progenyAtmosLoad(
          destpath = shinyFiles::parseDirPath(volumes, input$destpath),
          base = input$base,
          extend = input$extend,
          hot = input$hot,
          AGB = input$AGB,
          white = input$white,
          WR = input$WR,
          wavegrid = NULL,
          cores = input$atmos_cores
        )
        atmos_result(atmos_out)
        output$atmos_status <- renderText("Atmosphere data loaded successfully!")

        output$plot_atmos = renderPlot({
          progenyIsoPlot(iso_result())
          progenyAtmosPlot(atmos_out, add=TRUE)
        })
      }, error = function(e) {
        output$atmos_status <- renderText(paste("Error:", e$message))
      })
    #})
  })

  observeEvent(input$run_interp, {
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
        label_AGB = input$label_AGB,
        label_white = input$label_white,
        label_WR = input$label_WR
      )

      iso_result(new_iso)
      output$interp_status <- renderText("Interpolation successfully run!")
    }, error = function(e) {
      output$interp_status <- renderText(paste("Error:", e$message))
    })
  })

  observeEvent(input$update_imf, {
    if(input$imf == 'IMF_Chabrier'){
      output$plot_imf = renderPlot({
        magcurve(IMF_Chabrier(x,
                              alpha = input$chab_alpha,
                              a = input$chab_a,
                              b = input$chab_b,
                              masslow = input$masslow,
                              massmax = input$massmax),
                 input$masslow, input$massmax, log='xy', xlab='Star Mass / Msol', ylab='dN/dM (1 Msol)',
                 )
      })
    }else if(input$imf == 'IMF_Kroupa'){
      output$plot_imf = renderPlot({
        magcurve(IMF_Kroupa(x,
                            alpha1 = input$kroupa_alpha1,
                            alpha2 = input$kroupa_alpha2,
                            alpha3 = input$kroupa_alpha3,
                            mass1 = input$kroupa_mass1,
                            mass2 = input$kroupa_mass2,
                            masslow = input$masslow,
                            massmax = input$massmax),

                 input$masslow, input$massmax, log='xy', xlab='Star Mass / Msol', ylab='dN/dM (1 Msol)',
                 )
       })
    }else if(input$imf == 'IMF_Salpeter'){
      output$plot_imf = renderPlot({
        magcurve(IMF_Salpeter(x,
                              alpha = input$salp_alpha,
                              masslow = input$masslow,
                              massmax = input$massmax),
                 input$masslow, input$massmax, log='xy', xlab='Star Mass / Msol', ylab='dN/dM (1 Msol)',
                 )
      })
    }
  }
  )

  observeEvent(input$make_ssp, {
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

      if(input$imf == 'IMF_Chabrier'){
        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          alpha = input$chab_alpha,
          a = input$chab_a,
          b = input$chab_b,
          masslow = input$masslow,
          massmax = input$massmax,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
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
          masslow = input$masslow,
          massmax = input$massmax,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )
      }else if(input$imf == 'IMF_Salpeter'){
        print(input$salp_alpha)
        print(input$masslow)
        print(input$massmax)
        SSP_out = progenyMakeSSP(
          Iso = iso_data,
          IMFfunc = IMF_function,
          alpha = input$salp_alpha,
          masslow = input$masslow,
          massmax = input$massmax,
          Spec_combine = atmos_data,
          Interp_combine = interp_data,
          cores = input$SSP_cores
        )
      }

      SSP_result(SSP_out)
      output$SSP_status <- renderText("SSP successfully run!")
    }, error = function(e) {
      output$SSP_status <- renderText(paste("Error:", e$message))
    })
  })

  output$dynamic_imf <- renderUI({
    req(input$imf)
    if(input$imf == 'IMF_Chabrier'){
      tagList(
        sliderInput('chab_alpha', 'alpha', min=1, max=3, value=2.3, step=0.1),
        sliderInput('chab_a', 'a', min=0, max=1, value=0.08, step=0.01),
        sliderInput('chab_b', 'b', min=0, max=1, value=0.69, step=0.01)
      )
    } else if(input$imf == 'IMF_Kroupa'){
      tagList(
        sliderInput('kroupa_alpha1', 'alpha_1', min=0, max=1, value=0.3, step=0.1),
        sliderInput('kroupa_alpha2', 'alpha_2', min=1, max=2, value=1.3, step=0.1),
        sliderInput('kroupa_alpha3', 'alpha_3', min=2, max=4, value=2.3, step=0.1),
        sliderInput('kroupa_mass1', 'mass_1', min=0, max=1, value=0.08, step=0.01),
        sliderInput('kroupa_mass2', 'mass_2', min=0, max=2, value=0.5, step=0.01)
      )
    }else if(input$imf == 'IMF_Salpeter'){
      tagList(
        sliderInput('salp_alpha', 'alpha', min=1, max=3, value=2.35, step=0.1)
      )
    }
  })


  # observeEvent(input$check_ssp, {
  #   output$SSP_check <- renderText({
  #     capture.output(ProSpect::speclib_check(SSP_result()))
  #   })
  # })


  observeEvent(input$check_ssp, {
    output$SSP_check <- renderText({
      messages <- character()
      withCallingHandlers(
        {
          ProSpect::speclib_check(SSP_result())
        },
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      paste(messages, collapse = "\n")
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

  observeEvent(input$SSP_done, {
    stopApp(SSP_result())
  })
}

shinyApp(ui=ui, server=server)
