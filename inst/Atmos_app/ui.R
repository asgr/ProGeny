ui = fluidPage(
  titlePanel("ProGeny SSP Generator"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tabs == 'Isochrone'",
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
        condition = "input.tabs == 'IMF'",
        selectInput("imf", "IMF", choices = c("Chabrier" = "IMF_Chabrier",
                                              "Kroupa" = "IMF_Kroupa",
                                              "Salpeter" = "IMF_Salpeter",
                                              "Kroupa_evo" = "IMF_Kroupa_evo",
                                              "Kroupa_Zevo" = "IMF_Kroupa_Zevo",
                                              "Lacey_evo" = "IMF_Lacey_evo",
                                              "Lacey_Zevo" = "IMF_Lacey_Zevo")),
        tags$h4("IMF parameters:"),
        uiOutput("dynamic_imf")
      ),

      conditionalPanel(
        condition = "input.tabs == 'Make SSP'",
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
                  tabPanel("Isochrone",
                           verbatimTextOutput("iso_status"),
                           plotOutput("plot_iso", height = "600px")
                  ),
                  tabPanel("Atmospheres",
                           verbatimTextOutput("selectedPath"),
                           verbatimTextOutput("atmos_status"),
                           plotOutput("plot_atmos", height = "600px")
                  ),
                  tabPanel("Interpolate",
                           verbatimTextOutput("interp_status")
                  ),
                  tabPanel("IMF",
                           plotOutput("plot_imf", height = "600px")
                  ),
                  tabPanel("Make SSP",
                           verbatimTextOutput("used_imf"),
                           verbatimTextOutput("SSP_status"),
                           verbatimTextOutput("SSP_check")
                  )
      )
    )
  )
)
