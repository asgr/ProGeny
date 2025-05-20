server = function(input, output){

  volumes = c(wd=getwd(), Home='~/', getVolumes()())

  shinyDirChoose(input, id='destpath', roots = volumes)

  atmos_result <- reactiveVal(NULL)


  output$selectedPath <- renderText({
    req(input$destpath)
    parseDirPath(volumes, input$destpath)
  })

  output$status <- renderText("Waiting...")

  result = observeEvent(input$run, {
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
