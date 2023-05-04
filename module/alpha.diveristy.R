#######################################################
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


alphaClient <- function(id) {
    ns <- NS(id)
    tagList(

        fluidRow( 

            # htmltools::HTML("<script>$(function(){ $('.btn').mouseup(function() { this.blur(); console.log('hi'); }) });</script>"),
            # The button class is "bttn-pill" (F12). Adds pop-back action:
            htmltools::HTML("<script>$(function(){ $('.bttn-pill').mouseup(function() { this.blur(); }) });</script>"),

            # Add datatable row selection style:
            tags$style(HTML(names("css_row_select"))),

            # Control box
            box(

              shiny::tags$head(
                                tags$style(
                                  ".tabwid  table  caption { color: black;}
                                  .dropdown-menu { background-color: transparent; border-color: transparent }
                                  "
                                )
                              ),  # Set caption text in black in stats flextable
                width=3,

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
                selectInput(inputId=ns('method'), label='Alpha diversity index:',  choices=c("Shannon", "Chao1", "Obs"), selected="Shannon"),                
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
                                outputId = ns("downloadData"),
                                label="Diversity table",
                                style = "pill",
                                color = "success",
                                size='sm',
                                block=TRUE
                    ),
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
                ),
            ),

            ## Plots box
            box (
                  # style = 'padding: 0px',
                  # style="width: fit-content; overflow-x: auto",
                  sytle="overflow-x: auto",
                  # background = "black",
                  id=ns("render_box"),
                  width=9,  #width=10, #cant fit with control box
                  height="100%",                                                          
                  title = "Alpha diversity analysis", 
                  status = "success", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                            plotOutput(ns("plotgraph5")) %>% withSpinner(color="green",size=0.5), 
                            tags$br(),
                            tags$br(),
                            # DT::dataTableOutput(ns("wilcoxTable"))
                            #tableOutput(ns('wilcoxTable'))
                            uiOutput(ns("wilcoxTable")) %>% withSpinner(color="green",size=0.5), 
                            downloadBttn(
                                        outputId = ns("download"),
                                        label="Diversity table",
                                        style = "pill",
                                        color = "success",
                                        size='sm',
                                        # block=TRUE   # If TRUE, the button uses whole width
                            )
            ),

            # CSS options to fit content:
            # width: fit-content;
            # overflow-x: auto;
            tags$script(HTML("
                $(document).ready(function() {
                    
                    document.getElementById('alpha-render_box').style.overflowY = 'hidden';
                    document.getElementById('alpha-render_box').style.overflowX = 'auto';
                    //$('#alpha-render_box').parent().overflowX = 'auto';

                    // Get the content width
                    //var contentWidth = $('#alpha-wilcoxTable').width();

                    //alert(contentWidth);

                    // Add some margin
                    //var margin = 20;

                    // Set the box width
                    //if( contentWidth !=0 ) {
                      //$('#alpha-render_box').parent().width(contentWidth + margin);
                    //}

                });
            "))
        
        )
    )
}


alphaServer <- function(id) {
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
            span(HTML(paste(tags$h1("Built-in abundance fungi file"),tags$strong('(NOTE: File has to be comma "\\t" delimited.)')))),
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
output$plotgraph5 <- renderPlot({ 

})

output$wilcoxTable <- renderUI({ 

})


# Render plot (update canvas)
observeEvent(vals$plotgraph5, { 
    printd("TRIGGER # vals$plotgraph5"); 
    output$plotgraph5 <- renderPlot({ vals$plotgraph5 }) 
})

# Render stats table (update canvas)
observeEvent(vals$wtable, { 
    printd("TRIGGER # vals$wtable"); 
    output$wilcoxTable <- renderUI({ 

            inputsplitby <- isolate(input$splitby)
            # inputsplitby <- input$splitby
            printd("RUNNING observeEvent(vals$wtable # formatTable"); 
            print(head(vals$wtable))
            printd("--.")
            print(inputsplitby)
            printd("#--.#")

            if( !is.null(vals$wtable) & !is.null(inputsplitby) ) {
                if(inputsplitby != "") {
                    printd("FORMAT TABLE # ")
                    print(vals$wtable)
                    printd("--*")
                    print(inputsplitby)
                    printd("#--*#")

                    vals$wtable %>%
                        tibble::rownames_to_column(var = inputsplitby) %>%  # "rownames" cant be shown in flextable unless they are a column
                        flextable() %>%
                        theme_apa() %>%
                        colformat_double(big.mark = ",", digits = 3, na_str = "") %>%  # data.frame cols must be numeric
                        align(align = "center", part = "all") %>%
                        autofit() %>%
                        set_caption(caption = "Statistical significance (Wilcox test p-values)")  %>%
                        htmltools_value()
                }
            }

     }) 
})


# Check UI interaction and other reactive values that trigger the plot (plot input values)
observeEvent(c(input$method,input$splitby,vals$taxa_myco_phy), {
        if(!is.null(vals$taxa_myco_phy) & input$method != "" & input$splitby != "") {
            printd("TRIGGER # ANY OF input$method,input$splitby,vals$taxa_myco_phy HAVE BEEN UPDATED")
            print(vals$taxa_myco_phy)
            printd("--")
            print(input$method)
            printd("--")
            print(input$splitby)
            printd("#--#")
            if( !preloaded | (preloaded & (input$method != "Shannon" | input$splitby != "country" )) ) {
                printd("TRIGGER # RECALCULATING - running 'plot_alpha_diversity()' ")
                
                # Update UI before performing long calcs
                output$plotgraph5 <- renderPlot({  })    # Trigger HTML spinner
                output$wilcoxTable <- renderUI({ NULL }) # Clear canvas

                p=plot_alpha_diversity(physeq = isolate(vals$taxa_myco_phy), method = isolate(input$method),attribute = isolate(input$splitby))
                vals$plotgraph5 <<- p[["plot"]]
                vals$wtable <<- p[["wtable"]]
                diversity_df <<- p[["diversity"]]
                printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_alpha_diversity()' ")
                assign("alpha.vals_plotgraph5", vals$plotgraph5, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                printd("--+")
                print(vals$wtable)
                printd("--+")
                print(head(diversity_df))
                printd("#--+#")
                # output$plotgraph5 <- renderPlot({ vals$plotgraph5 })
                # output$wilcoxTable <- renderUI({ formatTable() })
            } else {
                if( preloaded & input$method == "Shannon" & input$splitby == "country" ) {
                    printd("TRIGGER # input$method == 'Shannon' & input$splitby == 'country' # UPDATING PRELOAD PLOT ")

                    vals$plotgraph5 <<- BDplot5
                    vals$wtable <<- BDwtable
                    diversity_df <<- BDdiversity

                    printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_alpha_diversity()' ")
                    assign("alpha.vals_plotgraph5", vals$plotgraph5, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                    printd("--+")
                    print(vals$wtable)
                    printd("--+")
                    print(head(diversity_df))
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

        printd("TRIGGER # input$method and input$loadbd # WILL TRIGGER input$splitby")
        updateSelectInput(session, inputId="splitby", choices=colnames(metadata_df), selected="country")
        updateSelectInput(session, inputId="method", selected="Shannon")
        printd("TRIGGER # input$loadbd # AFTER TRIGGERING input$splitby")
        print(input$method)   #debug. Value still NOT UPDATED
        print(input$splitby)  #debug

        # vals$taxa_myco_phy <<- BDtaxa_myco_phy

        preloaded <<- T
        
        printd("LOADBD # BEFORE updating vals$plotgraph5")
        vals$plotgraph5 <<- BDplot5
        printd("LOADBD # AFTER updating vals$plotgraph5")
        vals$wtable <<- BDwtable
        printd("LOADBD # vals$wtable")
        print(head(vals$wtable))
        diversity_df <<- BDdiversity
        printd("LOADBD # diversity_df")
        print(head(diversity_df))
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
        updateSelectInput(session, inputId="splitby", choices=colnames(metadata_df), selected=colnames(metadata_df)[1])
        printd("TRIGGER # input$loadbd # AFTER TRIGGERING input$splitby")
        print(input$splitby)  #debug. Value still NOT UPDATED, untill all observers(observeEvent) and outputs(render) are finished.


        preloaded <<- F

        # "observeEvent(c(input$method,input$splitby,vals$taxa_myco_phy)" must
        # update: vals$plotgraph5, vals$wtable, diversity_df
  })

    }
  )
}