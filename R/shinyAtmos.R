shinyAtmos = function(){
    require(shiny)
    require(shinyFiles)

    appDir <- system.file("Atmos_app", package = "ProSpect")
    if (appDir == "") {
        stop("Could not find Atmos_app directory. Try re-installing ProGeny.", call. = FALSE)
    }

    return(runApp(appDir, display.mode = "normal"))
}
