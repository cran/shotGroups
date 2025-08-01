fluidPage(
    fluidRow(
        box(title="Enter data",
            width=4,
            radioButtons("datIn", "",
                         list("Use built-in data"=1,
                              "Upload file"=2,
                              "Paste data"=3)),
            conditionalPanel(condition="input.datIn == '1'",
                             radioButtons("builtInData", h5("Built-in data:"),
                                          dataBuiltInInv, selected="2")),
            conditionalPanel(condition="(input.datIn == '2') || (input.datIn == '3')",
                             radioButtons("fileType", "File format:",
                                          list("OnTarget 1.*"=1,
                                               "OnTarget 2.*, 3.*"=2,
                                               "Silver Mountain"=3,
                                               "ShotMarker"=4,
                                               "SIUS"=5,
                                               "Other"=6), selected="2")),
            conditionalPanel(condition="input.datIn == '2'",
                             h5("Upload file: "),
                             fileInput("fileUpload", "Select file:", multiple=TRUE)),
            conditionalPanel(condition="input.datIn == '3'",
                             h5("Paste data:"),
                             textAreaInput("datPaste", label=NULL, rows=10, cols=30, resize="both")),
            actionButton("applyData", "Apply")
        ),
        box(title="Distance to target",
            width=8,
            uiOutput("unitDstXY"),
            h3("Information from imported file(s)"),
            uiOutput("fileInfo"),
            p("For details on how to read in data, see the documentation for",
              a("readDataOT1()",
                href="https://www.rdocumentation.org/packages/shotGroups/functions/readDataOT1"),
              ",",
              a("readDataOT2()",
                href="https://www.rdocumentation.org/packages/shotGroups/functions/readDataOT2"),
              ",",
              a("readDataSMT()",
                href="https://www.rdocumentation.org/packages/shotGroups/functions/readDataSMT"),
              ",",
              a("readDataShotMarker()",
                href="https://www.rdocumentation.org/packages/shotGroups/functions/readDataShotMarker"),
              ",",
              a("readDataMisc()",
                href="https://www.rdocumentation.org/packages/shotGroups/functions/readDataMisc"),
              "and the",
              a("shotGroups vignette",
                href="https://cran.rstudio.com/web/packages/shotGroups/vignettes/shotGroups.pdf"),
              "section 2.1")
        )
    )
)
