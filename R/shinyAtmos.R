shinyAtmos = function(){
    require(shiny)
    require(shinyFiles)
    library(fst)

    options(shiny.maxRequestSize = 2000*1024^2)

    appDir <- system.file("Atmos_app", package = "ProGeny")
    if (!file.exists(appDir)) {
        stop("Could not find Atmos_app directory. Try re-installing ProGeny.", call. = FALSE)
    }

    return(runApp(appDir, display.mode = "normal"))
}
