
library(assertthat)
library(magrittr)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(ggplot2)
library(fs)


options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# UI ----------------------------------------------------------------------

shinyUI(

    fixedPage(
        # shinythemes::themeSelector(),
        titlePanel("GEX"),
        theme = "maglab_theme.css",

        useShinyjs(),

        setBackgroundColor(
            color = c("#FFFFFF"),
            gradient = "linear",
            direction = "bottom"
        ),

        tags$head(
            tags$style(HTML("hr {border-top: 1px solid #A9A9A9;}")),
            tags$style(
                HTML(
                    ".shiny-notification {
                        position:fixed;
                        top: calc(20%);
                        left: calc(40%);
                    }"
                )
            )
        ),

        # Sidebar -----------------------------------------------------------------

        sidebarLayout(
            sidebarPanel(
                tabsetPanel(
                    id = "mainpanel",
                    type = "pills",
                    tabPanel(
                        "Make XICs",
                        br(),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            shinyFilesButton(
                                "rawfile",
                                "Select .raw file",
                                "Select a Thermo .raw file",
                                multiple = FALSE
                            )
                        ),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            shinyFilesButton(
                                "targetseqs",
                                "Select .csv file",
                                "Select a .csv file",
                                multiple = FALSE
                            )
                        ),
                        br(), br(),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            actionButton(
                                "GEXstart",
                                "Generate XICs"
                            )
                        ),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            downloadButton("downloadReport", label = "XIC Report")
                        ),
                        hr(),
                        # selectInput(
                        #     "plotchoice",
                        #     "Choose XIC to view",
                        #     choices = c("Generate XICs first")
                        # ),
                        selectInput(
                            "plotchoice",
                            "Choose XIC to view",
                            choices = c("Generate XICs first"),
                            selectize = FALSE,
                            size = 5
                        ),
                        br(),
                        uiOutput("spectrumchoice")
                        # selectInput(
                        #     "spectrumchoice",
                        #     "Choose spectrum to view",
                        #     choices = c("Select an XIC"),
                        #     selectize = FALSE,
                        #     size = 5
                        # )
                        # actionButton(
                        #     "updateplot",
                        #     "Update Plot"
                        # )
                    ),
                    tabPanel(
                        "Make Spectra"
                    ),
                    tabPanel(
                        "Settings",
                        br(),
                        selectInput(
                            "usedepleted",
                            "Use Depleted Isotopes?",
                            choices = c("Yes" = TRUE, "No" = FALSE)
                        ),
                        numericInput(
                            "sample_pforms",
                            "# of proteoforms to sample",
                            0,
                            min = 0,
                            max = 100,
                            step = 1
                        ),
                        numericRangeInput(
                            "massrange",
                            "Proteoform search mass range (Da)",
                            value = c(0,100000)
                        ),
                        numericRangeInput(
                            "charges",
                            "Charge state range",
                            value = c(1,50)
                        ),
                        numericRangeInput(
                            "mzrange",
                            "Spectrum search mass range (m/z)",
                            value = c(600,2000)
                        ),
                        numericInput(
                            "abundcutoff",
                            "Isotopic peak relative abundance cutoff (%)",
                            value = 5,
                            min = 1,
                            max = 100,
                            step = 1
                        ),
                        hr(),
                        h4("Enter names for target columns"),
                        br(),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            textInput(
                                "targetcolname",
                                "Name",
                                value = "UNIPROTKB"
                            )
                        ),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            textInput(
                                "targetseqname",
                                "Proteoform Seq.",
                                value = "ProteoformSequence"
                            )
                        ),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            textInput(
                                "ptmcolname",
                                "PTM",
                                value = "PTMname"
                            )
                        ),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            textInput(
                                "ptmaddname",
                                "Formula to add",
                                value = "FormulaToAdd"
                            )
                        ),
                        div(
                            style="display: inline-block;vertical-align:top; width: 150px;",
                            textInput(
                                "ptmsubname",
                                "Formula to subtract",
                                value = "FormulaToSubtract"
                            )
                        )
                    ),
                    tabPanel(
                        "About",
                        hr(),
                        includeMarkdown("about.md")
                    )
                )
            ),


            # Main Panel --------------------------------------------------------------

            mainPanel(
                htmlOutput("ULconfirm"),
                uiOutput(
                    "tdrep_fracassign"
                ),
                br(),
                textOutput("confirm"),
                br(),
                textOutput("error"),
                br(),
                plotOutput("outputPlot1"),
                br(), br(),
                plotOutput("outputPlot2")

            ),

            fluid = FALSE
        )


    )
)
