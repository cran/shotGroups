fluidPage(
    fluidRow(
        box(# height = "600px",
            title = "About",
            width = 6,
            # status = NULL, 
            # closable = FALSE,
            # maximizable = FALSE, 
            collapsible = FALSE,
            p("The", a("shotGroups", href="https://CRAN.R-project.org/package=shotGroups"),
              "package for", a("R", href="https://www.r-project.org/"),
              "provides functions to read in, plot,
                      statistically describe, analyze, and compare shooting data with respect
                      to group shape, precision, and accuracy. This includes graphical methods,
                      descriptive statistics, and inference tests using standard, but also
                      non-parametric and robust statistical methods. The data can be imported
                      from files produced by",
              a("OnTarget PC and OnTarget TDS", href="https://ontargetshooting.com/"), ", ",
              a("TARAN", href="http://taran.ptosis.ch/"), ", ",
              a("ShotMarker e-target", href="https://autotrickler.com/pages/shotmarker"), ", ",
              a("Silver Mountain e-target", href="https://www.silvermountaintargets.com/"), ", ",
              a("SIUS e-target", href="https://www.sius.com/"), ", ",
              "or from custom data files in text format with a similar structure.
                      For further explanations and an example walkthrough, see the",
              a("package vignette",
                href="https://cloud.r-project.org/web/packages/shotGroups/vignettes/shotGroups.pdf"),
              ". Many statistical methods are described on",
              a("Ballistipedia", href="http://ballistipedia.com/"), "."),
            p("shotGroups and this web application are written by:", br(),
              "Daniel", HTML("Wollschl&auml;ger"),
              a("<dwoll@kuci.org>", href="mailto:dwoll@kuci.org"), br(),
              a("Source code shotGroups",
                href="https://github.com/dwoll/shotGroups/"), br(),
              a("Source code shiny app",
                href="https://github.com/dwoll/shotGroups/tree/master/inst/shotGroups_AnalyzeGroups")),
            
            h6("More shotGroups web applications"),
            p("Absolute", icon("resize-horizontal", lib="glyphicon"),
              "angular size conversions:",
              a("http://shiny.imbei.uni-mainz.de:3838/shotGroups_AngularSize/",
                href="shiny.imbei.uni-mainz.de:3838/shotGroups_AngularSize/"), br(),
              "Region", icon("resize-horizontal", lib="glyphicon"),
              "hit probability calculations:",
              a("http://shiny.imbei.uni-mainz.de:3838/shotGroups_HitProb/",
                href="http://shiny.imbei.uni-mainz.de:3838/shotGroups_HitProb/"), br(),
              "Estimate Rayleigh sigma from range statistics:",
              a("http://shiny.imbei.uni-mainz.de:3838/shotGroups_RangeStat/",
                href="http://shiny.imbei.uni-mainz.de:3838/shotGroups_RangeStat/"))
        ),
        box(title="References",
            width=6,
            collapsible = FALSE,
            p("This web application is built with R, shiny, and bs4Dash.
              The shotGroups package
              uses functionality provided by the R packages boot, coin, CompQuadForm,
              energy, mvoutlier, and robustbase:"),
            
            p("Canty A, & Ripley BD. (2025). boot: Bootstrap R (S-Plus) Functions.", br(),
              a("https://CRAN.R-project.org/package=boot",
                href="https://CRAN.R-project.org/package=boot")),
            p("Duchesne P & Lafaye de Micheaux P. (2010). Computing the distribution
                       of quadratic forms: Further comparisons between the Liu-Tang-Zhang
                       approximation and exact methods. Computational Statistics and Data
                       Analysis, 54 (4), 858-862.", br(),
              a("https://CRAN.R-project.org/package=CompQuadForm",
                href="https://CRAN.R-project.org/package=CompQuadForm")),
            p("Filzmoser P & Gschwandtner M. (2025). mvoutlier: Multivariate
                       outlier detection based on robust methods.", br(),
              a("https://CRAN.R-project.org/package=mvoutlier",
                href="https://CRAN.R-project.org/package=mvoutlier")),
            p("Granjon D. (2025). bs4Dash: A 'Bootstrap 4' Version of 'shinydashboard'.", br(),
              a("https://CRAN.R-project.org/package=bs4Dash",
                href="https://CRAN.R-project.org/package=bs4Dash")),
            p("Hothorn T, Hornik K, van de Wiel MA, Zeileis A. (2008).
                       Implementing a Class of Permutation Tests: The coin Package. Journal of
                       Statistical Software, 28 (8), 1-23.", br(),
              a("https://www.jstatsoft.org/v28/i08/",
                href="https://www.jstatsoft.org/v28/i08/"), br(),
              a("https://CRAN.R-project.org/package=coin",
                href="https://CRAN.R-project.org/package=coin")),
            p("R Core Team (2025). R: A language and environment for statistical computing.
                       R Foundation for Statistical Computing, Vienna, Austria.", br(),
              a("https://www.R-project.org/", br(), href="https://www.R-project.org/")),
            p("Rizzo ML & Szekely GJ. (2025). energy: E-statistics
                      (energy statistics).", br(),
              a("https://CRAN.R-project.org/package=energy",
                href="https://CRAN.R-project.org/package=energy")),
            p("Rousseeuw PJ, Croux C, Todorov V et al. (2025). robustbase: Basic Robust Statistics.", br(),
              a("https://CRAN.R-project.org/package=robustbase",
                href="https://CRAN.R-project.org/package=robustbase")),
            p("Chang W, Cheng J, Allaire JJ et al. (2025). shiny: Web Application Framework for R.", br(),
              a("https://CRAN.R-project.org/package=shiny", href="https://CRAN.R-project.org/package=shiny"), br(),
              a("https://shiny.rstudio.com/", href="https://shiny.rstudio.com/")),
            p("Wollschlaeger, D. (2025). shotGroups: Analyze shot group data.", br(),
              a("https://CRAN.R-project.org/package=shotGroups",
                href="https://CRAN.R-project.org/package=shotGroups"), br(),
              a("https://github.com/dwoll/shotGroups/",
                href="https://github.com/dwoll/shotGroups/"), "- development version")
        )
    )
)
