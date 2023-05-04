##########################################################
#################### GLOBAL FUNCTIONS ####################
##########################################################
source("module/helpers.R",encoding = "utf-8")

##########################################################
################## GLOBAL BUILT-IN DATASETS ##############
##########################################################

# These datasets are loaded in "global.R":
# BDmetadata_df
# BDtaxa_myco
# BDtaxa_micro
# BDtaxa_myco_phy
# BDtaxa_micro_phy
# Working datasets are initialized to built-in datasets:
# metadata_df <- BDmetadata_df
# taxa_myco <- BDtaxa_myco
# taxa_micro <- BDtaxa_micro
# taxa_myco_phy <- BDtaxa_myco_phy
# taxa_micro_phy <- BDtaxa_micro_phy



betaClient <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow( 

            # htmltools::HTML("<script>$(function(){ $('.btn').mouseup(function() { this.blur(); console.log('hi'); }) });</script>"),
            # The button class is "bttn-pill" (F12). Adds pop-back action:
            htmltools::HTML("<script>$(function(){ $('.bttn-pill').mouseup(function() { this.blur(); }) });</script>"),

            # Add datatable row selection style:
            # tags$style(HTML(names("css_row_select"))),

            ## Control box
            box(
                width=3,
# actionButton("browser", "browser"), #debug
# tags$script("$('#browser').hide();"), #debug
                # Use built-in files (can last up to 3m to produce the 6 plots)
                div(class = "buttonupload",
                  HTML('&nbsp;'), HTML('&nbsp;'),            
                  actionBttn(
                    inputId = ns("loadbd"),
                    label = "Use built-in dataset",
                    style = "pill", 
                    color = "success",
                    icon = icon("server"),
                    size = "md",
                  )),
                tags$style(".buttonupload .bttn-primary{ background-color: green; }"),                                                                            tags$style(".buttonupload .bttn-primary{ background-color: green; }"),
                tags$h4(HTML('&nbsp;'),HTML('&nbsp;'),"Or upload your own 2 files:"),
                tags$br(),
                tags$ul( 
                        # Use uploaded files
                        div(tags$li(strong("Metadata")," file:", style = "font-size:1.5em")),
                        tags$h5("(format: csv,txt,xls,xlsx; <5MB)"),
                        actionBttn(
                          inputId = ns("show1"),
                          label = "Built-in metadata",
                          style = "pill", 
                          color = "success",
                          size = "sm",
                          icon = icon("eye")
                        ),  
                        tags$br(),
                        tags$br(),
                        tags$head(tags$style(".progress-bar{background-color: green;}")), #Change file upload progress bar color
                        # tags$style(".fa-chart-column {color:#00a65a}"),  # Change icon color
                        fileInput(ns("file1"),NULL,
                                  multiple = FALSE,
                                  # buttonLabel = "Browse...",
                                  # color:#0072B2  # Color blue
                                  buttonLabel = tags$div(HTML('<i class="fa fa-cloud-upload" style = "color:#28b78d;"></i> Browse...')),
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),

                        div(tags$li(strong("Taxon abundance")," file:", style = "font-size:1.5em")),
                        tags$h5("(format: csv,txt,xls,xlsx; <5MB)"),
                        actionBttn(
                          inputId = ns("show2"),
                          label = "Built-in fungi taxa",
                          style = "pill", 
                          color = "success",
                          size = "sm",
                          icon = icon("eye")

                        ),
                        tags$br(),
                        tags$br(),
                        fileInput(ns("file2"),NULL,
                                  multiple = FALSE,
                                  buttonLabel = tags$div(HTML('<i class="fa fa-cloud-upload" style = "color:#28b78d;"></i> Browse...')),
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),

                        # div(tags$li(strong("Bacteria taxa")," file and update:", style = "font-size:1.5em")),
                        # tags$h5("(format: csv,txt,xls,xlsx; <5MB)"),
                        # actionBttn(
                        #   inputId = ns("show3"),
                        #   label = "Built-in bacteria taxa",
                        #   style = "pill", 
                        #   color = "success",
                        #   size = "sm",
                        #   icon = icon("eye")
                        # ),
                        # tags$br(),
                        # tags$br(),
                        # fileInput(ns("file3"),NULL,
                        #           multiple = FALSE,
                        #           buttonLabel = tags$div(HTML('<i class="fa fa-cloud-upload" style = "color:#28b78d;"></i> Browse...')),
                        #           accept = c("text/csv",
                        #                      "text/comma-separated-values,text/plain",
                        #                      ".csv"))                        
                ),
                tags$br(),
                tags$h4(HTML('&nbsp;'),"Distances are calculated using Bray-Curtis method."),
                tags$br(),
                # selectInput(inputId=ns('method'), label='Distance method:',  choices=c("unifrac","bray"), selected="bray"),                
                # selectInput(inputId=ns('weighted'), label='Weighted (unifrac only):',  choices=c(T,F), selected=T),                
                selectInput(inputId=ns('splitby'), label='Group by:',  choices=""),
                tags$br(),                                                           
                # Download plot.
                # HTML specs default screen resolution 96ppi and optimal for Win-PC 
                # Internet and Mac optimal resolution is 72ppi (for printing: 300dpi)
                # A4 (WxH): 8.3 x 11.7 in 
                # A3 : 11.7 x 16.5 in
                # A5 : 5.8 x 8.3 in
                # Book size: 6 x 9 in (octavo) (9 x 6 in in landscape mode)
                numericInput(ns("w"), label = "Figure Width (in) (def: A4 landscape)", value = 11.7),
                numericInput(ns("h"), label = "Figure Height (in) (def: A4 landscape)", value = 8.3),
                numericInput(ns("ppi"), label = "Figure Resolution (ppi)", value = 72), 
                dropdownButton(
                    downloadBttn(
                      outputId = ns("pdf"),
                      label="PDF figure",
                      style = "pill",
                      color = "success",
                      size='sm',
                      block=TRUE
                    ),
                    downloadBttn(
                      outputId = ns("png"),
                      label="PNG figure",
                      style = "pill",
                      color = "success",
                      size='sm',
                      block=TRUE
                    ),
                    downloadBttn(
                      outputId = ns("jpeg"),
                      label="JPEG figure",
                      style = "pill",
                      color = "success",
                      size='sm',
                      block=TRUE
                    ),
                    downloadBttn(
                      outputId = ns("tiff"),
                      label="TIFF figure",
                      style = "pill",
                      color = "success",
                      size='sm',
                      block=TRUE
                    ),
                    circle=TRUE,
                    label="Download Figure",
                    status="success",
                    icon=icon("download")
                )
            ),

            ## Plots box
            tags$style(HTML("#beta-pvalText {
                            background-color: white;
                            border: none;
                          }
                      ")),
            box (
                  # style = 'padding: 0px',
                  # background = "black",
                  width=9,  #width=10, #cant fit with control box
                  height="100%",                                                          
                  title = "Beta diversity analysis", 
                  status = "success", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                            plotOutput(ns("plotgraph3")) %>% withSpinner(color="green",size=0.5), 
                            htmlOutput(ns("pvalText"))
            )
        )
    )
}



betaServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {


            #Show built-in files if clicked:
            showBuiltinMetadata <- function(failed = FALSE) {
              modalDialog(
                easyClose=TRUE,
                size="l",
                footer = tagList(modalButton("Dismiss")),
                span(HTML(paste(tags$h1("Built-in metadata file"),tags$strong('(NOTE: File has to be comma "," delimited.)')))),
                tags$hr(),
                
                # ", server = FALSE)", speeds up table pop-up up to 6seconds
                withSpinner( tagList( renderDT({DT_style(BDmetadata_df)}, server = FALSE) ), color="green") 
                  # Or: 
                  # shiny::renderDataTable(
                  #           #metadata_df[1:5,1:5],
                  #           metadata_df , 
                  #           options = list(
                  #           scrollX = T,
                  #           autoWidth = F,
                  #           pageLength = 5
                  #           )
              )
            }

            showBuiltinFungiData <- function(failed = FALSE) {
              modalDialog(
                easyClose=TRUE,
                size="l",
                footer = tagList(modalButton("Dismiss")),
                span(HTML(paste(tags$h1("Built-in abundance fungi file"),tags$strong('(NOTE: File has to be tab delimited.)')))),
                tags$hr(),
                withSpinner( tagList( renderDT({DT_style(BDtaxa_myco)}, server = FALSE) ), color="green") 
              )
            }

          # React to show built-in dataset files.
          observeEvent(input$show1, {
            showModal(showBuiltinMetadata())
          })

          observeEvent(input$show2, {
            showModal(showBuiltinFungiData())
          })


# Start up empty render areas (canvas)
output$plotgraph3 <- renderPlot({ 

})

output$pvalText <- renderText({ 

})

# Render plot (update canvas)
observeEvent(vals$plotgraph3, { 
    printd("TRIGGER # vals$plotgraph3"); 
    output$plotgraph3 <- renderPlot({ vals$plotgraph3 }) 
})

# Render stats text output (update canvas)
observeEvent(vals$pvalue, { 
    printd("TRIGGER # vals$pvalue"); 
    output$pvalText <- renderText({ 

                printd("RUNNING observeEvent(vals$pvalue "); 
                print(vals$pvalue)

                if( !is.null(vals$pvalue) ) {
                    HTML( paste("<br/><br/><span style=\"font-family: 'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;font-size: 16px;\"> <p style='text-align:center;'>  Statistical significance (Adonis test). P-value:  ",
                        vals$pvalue,"</p></span>") )
                }
    })
})


# Check UI interaction and other reactive values that trigger the plot (plot input values)
observeEvent(c(input$splitby,vals$taxa_myco_phy), {
        if(!is.null(vals$taxa_myco_phy) & input$splitby != "") {
            printd("TRIGGER # ANY OF input$splitby,vals$taxa_myco_phy HAVE BEEN UPDATED")
            print(vals$taxa_myco_phy)
            printd("--")
            print(input$splitby)
            printd("#--#")
            if( !preloaded | (preloaded & input$splitby != "country" ) ) {
              if(input$splitby != "none") {
                printd("TRIGGER # splitby != 'none' # RECALCULATING - running 'plot_beta_diversity()' ")

                # Update UI before performing long calcs
                vals$plotgraph3 <<-NULL  # Trigger HTML spinner

                p=plot_beta_diversity(physeq =  isolate(vals$taxa_myco_phy), method = 'bray', weighted = T, variable = isolate(input$splitby))
                vals$plotgraph3 <<- p[["plot"]]
                vals$pvalue <<- p[["pvalue"]]
                printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_beta_diversity()' ")
                assign("beta.vals_plotgraph3", vals$plotgraph3, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                printd("--+")
                print(vals$pvalue)
                printd("#--+#")
              } else {
                printd("TRIGGER # splitby == 'none' #  RECALCULATING - running 'plot_beta_diversity()' ")
                
                # Update UI before performing long calcs
                vals$plotgraph3 <<-NULL  # Trigger HTML spinner

                p=plot_beta_diversity(physeq =  isolate(vals$taxa_myco_phy), method = 'bray', weighted = T)
                vals$plotgraph3 <<- p[["plot"]]+geom_point(color='#1b9e77')
                vals$pvalue <<- p[["pvalue"]]
                printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_beta_diversity()' ")
                assign("beta.vals_plotgraph3", vals$plotgraph3, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                printd("--+")
                print(vals$pvalue)
                printd("#--+#")
              }
            } else {
              if(preloaded & input$splitby == "country") {
                printd("TRIGGER # splitby == 'country' #  UPDATING PRELOAD PLOT ")
                
                # No need to clear UI, as we don't need to spend time ("spinner") re-calc a new plot
                vals$plotgraph3 <<- BDplot3
                vals$pvalue <<- BDpvalue

                printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_beta_diversity()' ")
                assign("beta.vals_plotgraph3", vals$plotgraph3, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                printd("--+")
                print(vals$pvalue)
                printd("#--+#")
              }
            }
        }
})



# Load Built-in Dataset
observeEvent(input$loadbd, {

        printd("TRIGGER # input$loadbd")

        metadata_df <<- BDmetadata_df
        taxa_myco <<- BDtaxa_myco

        vals$taxa_myco_phy <<- BDtaxa_myco_phy

        printd("TRIGGER # input$loadbd # WILL TRIGGER input$splitby")
        updateSelectInput(session, inputId="splitby", choices=c("none", colnames(metadata_df)), selected="country")

        printd("TRIGGER # input$loadbd # AFTER TRIGGERING input$splitby")
        print(input$splitby)  # Value still NOT UPDATED

        # vals$taxa_myco_phy <<- BDtaxa_myco_phy

        preloaded <<- T

        printd("LOADBD # BEFORE updating vals$plotgraph3")
        vals$plotgraph3 <<- BDplot3
        printd("LOADBD # AFTER updating vals$plotgraph3")
        vals$pvalue <<- BDpvalue
        printd("LOADBD # vals$pvalue")
        print(vals$pvalue) #debug

  })

# Load Built-in Dataset (2 files must be provided)
observeEvent(input$file2, {

        # In order to make plots, file1 (metadata), must be provided first.
        req(input$file1) 

        # Updates metadata_df, taxa_myco and taxa_myco_phy
        # "taxa_myco_phy" should trigger "observeEvent(c(input$method,input$splitby,vals$taxa_myco_phy)"
        # to update the plot
        readData(input,session) 

        printd("TRIGGER # input$loadbd # WILL TRIGGER input$splitby")
        updateSelectInput(session, inputId="splitby", choices=c("none", colnames(metadata_df)), selected=colnames(metadata_df)[1])
        printd("TRIGGER # input$loadbd # AFTER TRIGGERING input$splitby")
        print(input$splitby)  #debug. Value still NOT UPDATED, untill all observers(observeEvent) and outputs(render) are finished.

        preloaded <<- F
  })


    }
  )
}