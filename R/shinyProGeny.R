runShinyProGeny = function(app_mode = 'user',
  iso_path = "~/Google Drive/My Drive/ProGeny_isochrone/",
  atmos_path = "~/Google Drive/My Drive/ProGeny_atmos/",
  cores_server = 4){
    if(!requireNamespace("ProSpect", quietly = TRUE)){
        stop("Need ProSpect package for runShinyProGeny to work! Get from GitHub asgr/ProSpect")
    }

    if(!requireNamespace("Rfits", quietly = TRUE)){
        stop("Need Rfits package for runShinyProGeny to work! Get from GitHub asgr/Rfits")
    }

    if(!requireNamespace("magicaxis", quietly = TRUE)){
        stop("Need magicaxis package for runShinyProGeny to work! Get from GitHub asgr/magicaxis")
    }

    if(!requireNamespace("shiny", quietly = TRUE)){
        stop("Need shiny package for runShinyProGeny to work! Get from CRAN.")
    }

    if(!requireNamespace("shinyFiles", quietly = TRUE)){
        stop("Need shinyFiles package for runShinyProGeny to work! Get from CRAN.")
    }

    if(!requireNamespace("shinybusy", quietly = TRUE)){
        stop("Need shinybusy package for runShinyProGeny to work! Get from CRAN.")
    }

    if(!requireNamespace("fst", quietly = TRUE)){
        stop("Need fst package for runShinyProGeny to work! Get from CRAN.")
    }

    if(!requireNamespace("fst", quietly = TRUE)){
        stop("Need fst package for runShinyProGeny to work! Get from CRAN.")
    }

    options(shiny.maxRequestSize = 2000*1024^2)
    options(ProGeny_app_mode = app_mode)

    if(app_mode == 'server'){
      options(ProGeny_iso_path = iso_path)
      options(ProGeny_atmos_path = atmos_path)
      options(ProGeny_cores = cores_server)
    }

    appDir <- system.file("Atmos_app", package = "ProGeny")
    if (!file.exists(appDir)) {
        stop("Could not find Atmos_app directory. Try re-installing ProGeny.", call. = FALSE)
    }

    return(shiny::runApp(appDir, display.mode = "normal"))
}
