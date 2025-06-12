runShinyProGeny = function(){
    if(!requireNamespace("shiny", quietly = TRUE)){
        stop("Need shiny package for runShinyProGeny to work! Get from CRAN.")
    }

    if(!requireNamespace("shinyFiles", quietly = TRUE)){
        stop("Need shinyFiles package for runShinyProGeny to work! Get from CRAN.")
    }

    if(!requireNamespace("fst", quietly = TRUE)){
        stop("Need fst package for runShinyProGeny to work! Get from CRAN.")
    }

    options(shiny.maxRequestSize = 2000*1024^2)

    appDir <- system.file("Atmos_app", package = "ProGeny")
    if (!file.exists(appDir)) {
        stop("Could not find Atmos_app directory. Try re-installing ProGeny.", call. = FALSE)
    }

    return(shiny::runApp(appDir, display.mode = "normal"))
}
