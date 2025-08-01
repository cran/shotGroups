runGUI <-
function(app=c("analyze", "hitprob", "range", "angular"), ...) {
    app <- match.arg(toupper(app),
                     choices=c("ANALYZE", "HITPROB", "RANGE", "ANGULAR"),
                     several.ok=FALSE)

    appDir <- if(app == "ANALYZE") {
        system.file("shotGroups_AnalyzeGroups", package="shotGroups")
    } else if(app == "HITPROB") {
        system.file("shotGroups_HitProb",       package="shotGroups")
    } else if(app == "RANGE") {
        system.file("shotGroups_RangeStat",     package="shotGroups")
    } else if(app == "ANGULAR") {
        system.file("shotGroups_AngularSize",   package="shotGroups")
    } else {
    	NA_character_
    }

    if(is.na(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'shotGroups'.", call.=FALSE)
    }

    if(requireNamespace("shiny", quietly=TRUE)) {
        ## check if we have bs4Dash for newer GUI
        if(requireNamespace("bs4Dash", quietly=TRUE)) {
            ## breaking changes introduced in bs4Dash 2.0.0
            ## check which version is available
            bs4Dash_version <- packageVersion("bs4Dash")
            if(compareVersion("2.0.0", as.character(bs4Dash_version)) <= 0) {
                # appDir_legacy <- paste0(appDir, "_legacy")
                # shiny::runApp(appDir_legacy, ...)
                shiny::runApp(appDir, ...)
            } else {
                appDir_bs4Dash_old <- paste0(appDir, "_bs4Dash_05")
                shiny::runApp(appDir_bs4Dash_old, ...)
            }
        } else {
            warning("Package 'bs4Dash' not found - running legacy version")
            appDir_legacy <- paste0(appDir, "_legacy")
            shiny::runApp(appDir_legacy, ...)
        }
    } else {
        stop("Could not find package shiny - please install first")
    }
}
