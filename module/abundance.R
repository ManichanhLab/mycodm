##########################################################
#################### GLOBAL FUNCTIONS ####################
##########################################################
source("module/helpers.R",encoding = "utf-8")

##########################################################
################## GLOBAL BUILT-IN DATASETS ##############
##########################################################

abundanceClient <- function(id) {
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
                ),
                selectInput(inputId=ns('taxafill'), label='Taxa level:',  choices=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), selected="Phylum"),                
                # Select metadata column to split by: 
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
            box (
                  # style = 'padding: 0px',
                  # background = "black",
                  width=9,  #width=10, #cant fit with control box
                  height="100%",                                                          
                  #title = "Abundance analysis of taxa", 
                  title = "Microbial composition",
                  status = "success", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  # conditionalPanel(  # INITIAL loading (NULL graph). Keeps box space but plots are empty.
                  #       condition= c("input.loadbd == 0"), 
                  #       # Dummy plots to keep empty plot spaces.
                  #       plotOutput("plotgraph1"),
                  #       plotOutput("plotgraph2"),
                  #       plotOutput("plotgraph3"),
                  #       plotOutput("plotgraph4"),
                  #       plotOutput("plotgraph5"),
                  #       plotOutput("plotgraph6")
                  # ),
                      # style= "display: none;",
                      #splitLayout(cellArgs = list(style = "padding: 1px"),cellWidths = c("50%", "50%"), withSpinner(plotOutput("plotgraph1"),color="green",size=0.5),     plotOutput("plotgraph2") %>% withSpinner(color="green",size=0.5)), 
                            plotOutput(ns("plotgraph1")) %>% withSpinner(color="green",size=0.5),
                            plotOutput(ns("plotgraph1Legend")) %>% withSpinner(color="green",size=0.5)
            )
        )
    )
}



abundanceServer <- function(id) {
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
output$plotgraph1 <- renderPlot({ 

})

output$plotgraph1Legend <- renderPlot({ 

})


# Render plot (update canvas)
observeEvent(vals$plotgraph1, { 
    printd("TRIGGER # vals$plotgraph1"); 
    output$plotgraph1 <- renderPlot({ vals$plotgraph1 }) 
})

# Render plot (update canvas)
observeEvent(vals$plotgraph1Legend, { 
    printd("TRIGGER # vals$plotgraph1Legend"); 
    output$plotgraph1Legend <- renderPlot({ vals$plotgraph1Legend }) 
})


# Check UI interaction and other reactive values that trigger the plot (plot input values)
observeEvent(c(input$taxafill,input$splitby,vals$taxa_myco_phy), {
        if(!is.null(vals$taxa_myco_phy) & input$taxafill != "" & input$splitby != "") {
            printd("TRIGGER # ANY OF input$taxafill,input$splitby,vals$taxa_myco_phy HAVE BEEN UPDATED")
            print(vals$taxa_myco_phy)
            printd("--")
            print(input$taxafill)
            printd("--")
            print(input$splitby)
            printd("#--#")
            if( !preloaded | (preloaded & (input$taxafill != "Phylum" | input$splitby != "disease" )) ) {
                
                # Update UI before performing long calcs
                vals$plotgraph1 <<-NULL    # Trigger HTML spinner
                vals$plotgraph1Legend <<-NULL

                if(input$splitby != "none") {
                    printd("TRIGGER # splitby != 'none' # RECALCULATING - running 'plot_high_abundance()' ")                    
                    if(isolate(input$taxafill) != "Species") {
                        p=plot_high_abundance(physeq = isolate(vals$taxa_myco_phy), level=isolate(input$taxafill), title = "Fungal composition", legend = "bottom", facet=isolate(input$splitby))
                        vals$plotgraph1 <<- p
                    } else {
                        p=plot_high_abundance(physeq = isolate(vals$taxa_myco_phy), level=isolate(input$taxafill), title = "Fungal composition", legend = "bottom", n_taxa=25, facet=isolate(input$splitby))
                        legend <- get_legend(p) %>% as_ggplot()
                        vals$plotgraph1Legend <<- legend
                        vals$plotgraph1 <<- p + theme(legend.position = "none")
                    }

                } else {
                    printd("TRIGGER # splitby != 'none' # RECALCULATING - running 'plot_high_abundance()' ")
                    if(isolate(input$taxafill) != "Species") {
                        p=plot_high_abundance(physeq = isolate(vals$taxa_myco_phy), level=isolate(input$taxafill), title = "Fungal composition", legend = "bottom")
                        vals$plotgraph1 <<- p
                    } else {
                        p=plot_high_abundance(physeq = isolate(vals$taxa_myco_phy), level=isolate(input$taxafill), title = "Fungal composition", legend = "bottom", n_taxa=25)
                        legend <- get_legend(p) %>% as_ggplot()
                        vals$plotgraph1Legend <<- legend
                        vals$plotgraph1 <<- p + theme(legend.position = "none")
                    }
                }
                printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_high_abundance()' ")
                assign("higha.vals_plotgraph1", vals$plotgraph1, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                printd("#--+#")
            } else {
                if( preloaded & input$taxafill == "Phylum" & input$splitby == "disease" ) {
                    printd("TRIGGER # input$taxafill == 'Phylum' & input$splitby == 'disease' # UPDATING PRELOAD PLOT ")

                    # Clear Species-Legend
                    vals$plotgraph1Legend <<-NULL
                    
                    vals$plotgraph1 <<- BDplot1

                    printd("TRIGGER # AFTER RECALCULATING - AFTER running 'plot_high_abundance()' ")
                    assign("higha.vals_plotgraph1", vals$plotgraph1, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
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

        printd("TRIGGER # input$loadbd # WILL TRIGGER  input$taxafill and input$splitby")
        updateSelectInput(session, inputId="taxafill", selected="Phylum")
        updateSelectInput(session, inputId="splitby", choices=c("none", colnames(metadata_df)), selected="disease")

        printd("TRIGGER # input$loadbd # AFTER TRIGGERING input$splitby")
        print(input$taxafill)   #debug. Value still NOT UPDATED
        print(input$splitby)    #debug

        preloaded <<- T

        printd("LOADBD # BEFORE updating vals$plotgraph1")
        vals$plotgraph1 <<- BDplot1
        printd("LOADBD # AFTER updating vals$plotgraph1")
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
        updateSelectInput(session, inputId="splitby", choices=c("none",colnames(metadata_df)), selected=colnames(metadata_df)[1])
        printd("TRIGGER # input$loadbd # AFTER TRIGGERING input$splitby")
        print(input$splitby)  #debug. Value still NOT UPDATED, untill all observers(observeEvent) and outputs(render) are finished.

        preloaded <<- F
  })

# Download plots hanlders:

#PDF (up to 15s)
output$pdf <- downloadHandler(
    filename="summary.pdf",
    content = function(file){
      pdf(file,width=input$w,height=input$h,onefile = TRUE) #units: in
          print(vals$plotgraph1)
      dev.off()
    }
  )

#PNG (up to 40s)
output$png <- downloadHandler(
  filename="summary_png.zip",
  # filename = function() { paste(input$dataset, '.png', sep='') },
  # contentType= "image/png",
  content = function(file){
    # A temp file is created in AppData\Local\Temp\RtmpXXX
    # this file is renamed at the end of this function with "filename" and returned as a downloaded file

    # png(file,width=input$w,height=input$h,units="in",res=input$ppi)
    #   print( (vals$plotgraph1+vals$plotgraph2)/(vals$plotgraph3+vals$plotgraph4)/(vals$plotgraph5+vals$plotgraph6 ) )   #requires "library(patchwork)"
    # dev.off()
    # Saves each plot in a different .png file (named by numbering)
    # multi.page <- ggarrange(vals$plotgraph1, vals$plotgraph2, vals$plotgraph3, vals$plotgraph4,vals$plotgraph5,vals$plotgraph6,
    #                     nrow = 2, ncol = 1)
    # ggexport(multi.page, filename = file,res=input$ppi)
    list_plots <- list()
    list_plots[[1]] <- vals$plotgraph1
    # #device <- function(..., width, height) grDevices::png(..., width = input$w, height = input$h, res = input$ppi, units = "in")
    imgs_url <- sapply( 1:length(list_plots), function(i) {
              #file_n <- paste0(sub('\\.png$', '\\.', file) ,i,".png")
              file_n <- paste0(dirname(file),"/plot_" ,i,".png")
              # If not using "dirname(file)" (which is the temp folder), the files is saved in he app directory (where app.R is run)
              ggsave( 
                  filename =  file_n , 
                  plot = list_plots[[i]], 
                    width = input$w, height = input$h,units="in",dpi=input$ppi
                  # device = device
                ) #Save files in temp folder: AppData\Local\Temp\RtmpXXX
              return(file_n)
            }
        )
    # Set wd (currently set to the app directory) to the current temporary working directory to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    imgs <- image_read(imgs_url)
    # imgs <- c(imgs[1:2], image_blank(width = 3, height = 5), imgs[3])
    # https://www.imagemagick.org/Magick++/Geometry.html
    # A4 595x842 pixels (geometry must be given in px) 
    im<-image_montage(imgs, tile = '2x3',geometry = "842x595")
    file_combined <- "plot_combined.png"
    image_write(im, path = file_combined, format = "png")
    # Create the zip file
    zip(file,c(basename(imgs_url),file_combined))
  }
)

#JPEG
output$jpeg <- downloadHandler(
  filename="summary_jpeg.zip",
  content = function(file){
    list_plots <- list()
    list_plots[[1]] <- vals$plotgraph1
    imgs_url <- sapply( 1:length(list_plots), function(i) {
              file_n <- paste0(dirname(file),"/plot_" ,i,".jpeg")
              ggsave( 
                  filename =  file_n , 
                  plot = list_plots[[i]], 
                    width = input$w, height = input$h,units="in",dpi=input$ppi
                )
              return(file_n)
            }
        )
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    imgs <- image_read(imgs_url)
    im<-image_montage(imgs, tile = '2x3',geometry = "842x595")
    file_combined <- "plot_combined.jpeg"
    image_write(im, path = file_combined, format = "jpeg")
    zip(file,c(basename(imgs_url),file_combined))
  }
)

#TIFF
output$tiff <- downloadHandler(
  filename="summary_tiff.zip",
  content = function(file){
    list_plots <- list()
    list_plots[[1]] <- vals$plotgraph1
    imgs_url <- sapply( 1:length(list_plots), function(i) {
              file_n <- paste0(dirname(file),"/plot_" ,i,".tiff")
              ggsave( 
                  filename =  file_n , 
                  plot = list_plots[[i]], 
                    width = input$w, height = input$h,units="in",dpi=input$ppi
                )
              return(file_n)
            }
        )
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    imgs <- image_read(imgs_url)
    im<-image_montage(imgs, tile = '2x3',geometry = "842x595")
    file_combined <- "plot_combined.tiff"
    image_write(im, path = file_combined, format = "tiff")
    zip(file,c(basename(imgs_url),file_combined))
  }
)

    }
  )
}