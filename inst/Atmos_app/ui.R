ui = fluidPage(
  titlePanel("Progeny Atmosphere Loader"),

  sidebarLayout(
    sidebarPanel(
      shinyDirButton("destpath", "Atmos Path", "Select a folder"),
      selectInput("base", "Base", choices = c("C3K" = "combine_C3K_conroy", "Husser" = "combine_PHOENIX_husser", "Allard" = "combine_PHOENIX_Allard", "MILES" = "combine_MILES_vazdekis", "BASEL = combine_BASEL_wlbc")),
      selectInput("extend", "Extend", choices = c("Allard" = "combine_PHOENIX_Allard", "ATLAS9" = "combine_ATLAS9_castelli")),
      selectInput("hot", "Hot", choices = c("OB" = "combine_OB_PoWR", "hot" = "combine_hot")),
      selectInput("AGB", "AGB", choices = "combine_AGB_lancon"),
      selectInput("white", "White Dwarfs", choices = c("TMAP" = "combine_TMAP_werner", "white" = "combvine_white")),
      selectInput("WR", "Wolf-Rayet", choices = c("PoWR" = "combine_WNE_PoWR")),
      numericInput("cores", "Number of Cores", value = 8, min = 1),
      actionButton("run", "Run progenyAtmosLoad"),
      br(), br(),
      actionButton("done", "Return Result")
    ),

    mainPanel(
      verbatimTextOutput("selectedPath"),
      verbatimTextOutput("status")
    )
  )
)
