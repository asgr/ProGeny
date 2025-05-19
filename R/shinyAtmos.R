shinyAtmos = function(){
    require(shiny)
    require(shinyFiles)

    # Define UI for application that draws a histogram

    runApp(list(
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
        ),

        server = function(input, output) {

            volumes = c(wd=getwd(), Home='~/', getVolumes()())

            shinyDirChoose(input, id='destpath', roots = volumes)

            atmos_result <- reactiveVal(NULL)


            output$selectedPath <- renderText({
                req(input$destpath)
                parseDirPath(volumes, input$destpath)
            })


            result = observeEvent(input$run, {
                output$status <- renderText("Running...")
                tryCatch({
                    atmos_out = progenyAtmosLoad(
                        destpath = parseDirPath(volumes,input$destpath),
                        base = input$base,
                        extend = input$extend,
                        hot = input$hot,
                        AGB = input$AGB,
                        white = input$white,
                        WR = input$WR,
                        wavegrid = NULL,
                        cores = input$cores
                    )
                    atmos_result(atmos_out)
                    output$status <- renderText("Done!")
                }, error = function(e) {
                    output$status <- renderText(paste("Error:", e$message))
                })
            })


            observeEvent(input$done, {
                stopApp(atmos_result())
            })
        }
    )
    )
}
