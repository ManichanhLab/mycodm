##########################################################
## Project: FUNOMIC
## Script purpose: Shiny App analysis
## Date: 18/11/22
## Organization: Vall Hebron Institut de Recerca (VHIR)
## Author: Xavier Martinez
##########################################################

##########################################################
###################### REFERENCES ########################
##########################################################

# Zhao T and Wang Z (2022) GraphBio: A shiny web app to easily perform popular visualization analysis for omics data. Front. Genet. 13:957317. doi: 10.3389/fgene.2022.957317
# AbdulMajedRaja RS; https://medium.com/free-code-camp/build-your-first-web-app-dashboard-using-shiny-and-r-ec433c9f3f6c
# Mastering Shiny by Hadley Wickham; https://www.oreilly.com/library/view/mastering-shiny/9781492047377/ch01.html

##########################################################
####################### LIBRARIES ########################
##########################################################


## Environment:
# install.packages("renv")
# install.packages("shiny")
# install.packages("rsconnect")
# install.packages("usethis")
# renv::snapshot()

## Complements:
# install.packages("shinyjs")
# install.packages("shinyWidgets")
# install.packages("shinydashboard")
# install.packages("shinycssloaders")
# install.packages("r2symbols")
# install.packages("tools")
# install.packages("rlang")

## Complements for printing out plots:
# install.packages("magick")
# install.packages("patchwork") #not used
# install.packages("png") #not used
# install.packages("RColorBrewer")

## Analysis:
# install.packages("BiocManager")
# install.packages("DT")
# install.packages("flextable")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("phyloseq")
# install.packages("ggpubr")
# install.packages("stringr")
# install.packages("edgeR")
# install.packages("ape")
# install.packages("vegan")
# install.packages("microbiome")
# install.packages("ggprism")

## Debugging:
# install.packages("profvis") #debug
# > library(profvis)
# >library(profvis); profvis(shiny::runApp()) #debug
# Example:
# data <- as.data.frame( matrix(rnorm(4e5*150,mean=5),ncol=150) )
# profvis({
# normCols <- function(d) {
# means <- apply(d,2,mean)
# for(i in seq_along(means)) {
# d[,i] <- d[,i] - means[i]
# }
# d}
# normCols(data)
# })

# Shiny
library(shiny)
library(shinyjs)
library(shinyWidgets)
require(shinydashboard)
library(shinycssloaders)

# Bioconductor
library(BiocManager)
library(microbiome)
library(ggprism)

# Complements
library(r2symbols)   # symbol.setting(font.size=30,font.weight ="bold",font.color = "brown")
library(magick)
library(tools)
library(RColorBrewer) # Color Palettes



##########################################################
#################### GLOBAL OBJECTS ######################
##########################################################
source("global.R")


##########################################################
#################### LOAD MODULES ########################
##########################################################

source("module/alpha.diveristy.R",encoding = "utf-8")
source("module/abundance.R",encoding = "utf-8")
source("module/beta.diversity.R",encoding = "utf-8")


##########################################################
######################### WEB CLIENT #####################
##########################################################


# Define UI for application
ui <- dashboardPage(


              title = "AnalysisPlatform",   #Title shown in the internet browser web-tab.
              # titlePanel("AnalysisPlatform"), 

              skin='green',

              # header <- dashboardHeader(title = "Analysis Platform"),
              header <- dashboardHeader( 
                                            title = tags$a("Analysis platform",href='', style="color: white")  # ,titleWidth = 450 
                                        ),

              sidebar <- dashboardSidebar(    
                
                                            width = 130,    
                                            
                                            # Symbols: "alpha","beta","mu","chi" - "\u03B1","\u03B2","\u03BC","\u03C7"
                                            sidebarMenu(
                                              tags$style(".fa-chart-column {color:#00a65a}"),  # Change icon color
                                              menuItem("Home", tabName = "home" , icon = icon("home")),
                                              # menuItem(strong("\u03B1-diversity"), tabName = "alpha"),
                                              menuItem(HTML(paste(symbol("alpha",font.size=20,font.weight ="bold",font.color = "#00a65a"),"-diversity")), tabName = "alpha"),
                                              menuItem("Abundance", tabName = "abundance", icon = icon("chart-column")),
                                              menuItem(HTML(paste(symbol("beta",font.size=20,font.weight ="bold",font.color = "#00a65a"),"-diversity")), tabName = "beta"),
                                              # menuItem("Alpha Diversity", tabName = "alpha", icon = icon("circle")),
                                              menuItem("References", tabName = "references", icon = icon("circle-info")),
                                              menuItem("Contact us", tabName = "contact", icon = icon("location-dot"))
                                              #menuItem("Contact", icon = icon("send",lib='glyphicon'), 
                                              #menuItem("Contact", icon = ionicon(name="arrow-forward")
                                            )
                
                                      ),
                  

              # browser(), #debug

              body <- dashboardBody(

                          # JS support:
                          # tags$head(includeScript("./js/<FILE>.js")),
                          useShinyjs(),
                          
                          # Connects menuitems with modules
                          tabItems(
                             tabItem(tabName = "home",
                                        frow1 <- fluidRow( 

                                                          # htmltools::HTML("<script>$(function(){ $('.btn').mouseup(function() { this.blur(); console.log('hi'); }) });</script>"),
                                                          # The button class is "bttn-pill" (F12). Adds pop-back action:
                                                          htmltools::HTML("<script>$(function(){ $('.bttn-pill').mouseup(function() { this.blur(); }) });</script>"),

                                                          # Add datatable row selection style:
                                                          tags$style(HTML(css_row_select)),

                                                          ## Control box
                                                          box(
                                                              width=3,
                                                              
                                                              # Use built-in files (may last up to 3m to produce the 6 plots)
                                                              div(class = "buttonupload", 
                                                                HTML('&nbsp;'), HTML('&nbsp;'),
                                                                actionBttn(
                                                                  inputId = "loadbd",
                                                                  label = "Use built-in dataset",
                                                                  style = "pill", 
                                                                  color = "success",
                                                                  icon = icon("server"),
                                                                  size = "md",
                                                                )),
                                                              tags$style(".buttonupload .bttn-primary{ background-color: green; }"),
                                                              tags$h4(HTML('&nbsp;'),HTML('&nbsp;'),"Or upload your own 2 files:"),
                                                              tags$br(),
                                                              tags$ul(
                                                                        # Use uploaded files
                                                                        div(tags$li(strong("Metadata")," file:", style = "font-size:1.5em")),
                                                                        tags$h5("(format: csv,txt,xls,xlsx; <5MB)"),
                                                                        actionBttn(
                                                                          inputId = "show1",
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
                                                                        fileInput("file1",NULL,
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
                                                                          inputId = "show2",
                                                                          label = "Built-in fungi taxa",
                                                                          style = "pill", 
                                                                          color = "success",
                                                                          size = "sm",
                                                                          icon = icon("eye")

                                                                        ),
                                                                        tags$br(),
                                                                        tags$br(),
                                                                        fileInput("file2",NULL,
                                                                                  multiple = FALSE,
                                                                                  buttonLabel = tags$div(HTML('<i class="fa fa-cloud-upload" style = "color:#28b78d;"></i> Browse...')),
                                                                                  accept = c("text/csv",
                                                                                             "text/comma-separated-values,text/plain",
                                                                                             ".csv"))
                                                              ),
                                                              
                                                              
                                                              
                                                              # Download plot.
                                                              # HTML specs default screen resolution 96ppi and optimal for Win-PC 
                                                              # Internet and Mac optimal resolution is 72ppi (for printing: 300dpi)
                                                              # A4 (WxH): 8.3 x 11.7 in 
                                                              # A3 : 11.7 x 16.5 in
                                                              # A5 : 5.8 x 8.3 in
                                                              # Book size: 6 x 9 in (octavo) (9 x 6 in in landscape mode)
                                                              numericInput("w", label = "Figure Width (in) (def: A4 landscape)", value = 11.7),
                                                              numericInput("h", label = "Figure Height (in) (def: A4 landscape)", value = 8.3),
                                                              numericInput("ppi", label = "Figure Resolution (ppi)", value = 72), 
                                                              dropdownButton(
                                                                  downloadBttn(
                                                                    outputId = "pdf",
                                                                    label="PDF figure",
                                                                    style = "pill",
                                                                    color = "success",
                                                                    size='sm',
                                                                    block=TRUE
                                                                  ),
                                                                  downloadBttn(
                                                                    outputId = "png",
                                                                    label="PNG figure",
                                                                    style = "pill",
                                                                    color = "success",
                                                                    size='sm',
                                                                    block=TRUE
                                                                  ),
                                                                  downloadBttn(
                                                                    outputId = "jpeg",
                                                                    label="JPEG figure",
                                                                    style = "pill",
                                                                    color = "success",
                                                                    size='sm',
                                                                    block=TRUE
                                                                  ),
                                                                  downloadBttn(
                                                                    outputId = "tiff",
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
                                                          box(
                                                                id="intro_text",                                                                
                                                                width=9,  #width=10, #cant fit with control box
                                                                height="100%",                                                          
                                                                title = HTML("<h4 style='margin-top: 0px;margin-bottom: 0px;'>About</h4>"), 
                                                                status = "success", 
                                                                solidHeader = TRUE,
                                                                collapsible = TRUE,
                                                                HTML("<h4 style='text-align: justify'>This a user-friendly web platform for easy statistical analysis and visualization of microbiome data. You can choose from our built-in dataset or upload your own, and get detailed results like taxa bar plots, alpha diversity box plots, and beta diversity PCoA plots. Our goal is to empower researchers to explore and discover insights from their microbiome data in a simple and intuitive way.</h4>"),
                                                                HTML("<h4 style='text-align: justify'><br/><i>Quick start<br/><br/>Click <strong>'Use built-in dataset'</strong> to begin exploring the analysis with our curated dataset.</i></h4>"),
                                                                HTML("<h4 style='text-align: justify'><i>Use left-side menu to perform independent analysis. Use 'Browse' buttons to upload your own dataset following the format shown in our 'Metadata' and 'Taxon abundance' files.</i></h4>")
                                                          ),

                                                          ## Plots box
                                                          box (
                                                                # style = 'padding: 0px',
                                                                # background = "black",
                                                                id="sum_plot",
                                                                width=9,  #width=10, #cant fit with control box
                                                                height="100%",                                                          
                                                                title = HTML("Data Analysis Summary<br/><span style='font-size: 14px;'>(please,note that processing big tables may last several minutes)</span>"), 
                                                                status = "success", 
                                                                solidHeader = TRUE,
                                                                collapsible = TRUE,
                                                                collapsed = TRUE,
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
                                                                    splitLayout(cellArgs = list(style = "padding: 1px"), cellWidths = c("100%"),
                                                                          plotOutput("plotgraph5") %>% withSpinner(color="green",size=0.5), 
                                                                          ),
                                                                    splitLayout(cellArgs = list(style = "padding: 1px"),cellWidths = c("100%"),
                                                                          plotOutput("plotgraph1") %>% withSpinner(color="green",size=0.5), 
                                                                          ),
                                                                    splitLayout(cellArgs = list(style = "padding: 1px"),cellWidths = c("100%"),
                                                                          plotOutput("plotgraph3") %>% withSpinner(color="green",size=0.5), 
                                                                          )
                                                                                                                           
                                                          ),
                                                          # Hide "summary plot" at startup.
                                                          # NOTE: when hiding from client and not from server, then, when the "sum plot" 
                                                          # is later shown, it doesn't work properly and does not show the plots.
                                                          # htmltools::HTML("<script>$('#sum_plot').parent().hide()</script>")

                                                      )

                                ),
                               tabItem(tabName="alpha",alphaClient("alpha")),
                               tabItem(tabName="abundance",abundanceClient("abundance")),
                               tabItem(tabName="beta",betaClient("beta")),
                               tabItem(tabName = "contact",
                                      fluidRow(
                                        box(width=12,
                                            title="Contact us",solidHeader=TRUE,status='success',height=800,
                                            tags$p("Lab description..."),
                                      ))
                                    ),
                               tabItem(tabName = "references",
                                    fluidRow(
                                      box(width=12,
                                          title="References",solidHeader=TRUE,status='success',height=800,
                                          tags$h2("References:"),
                                          tags$hr(),
                                          tags$p("1. Zhao T and Wang Z (2022), GraphBio: A shiny web app to easily perform popular visualization analysis for omics data. Front. Genet. 13:957317.doi: 10.3389/fgene.2022.957317"),
                                       )
                                    )
                                )
                                     
                          ) #tabItems
                
              ) #dashboardBody
) #end ui




##########################################################
######################### SERVER #########################
##########################################################

# Define server logic required to draw a histogram
server <- function(input, output, session) {

## Shiny session timeout
#If no user interaction, though the app might be working - creating a plot f.i -, the timeout pops up and session gets disconnected.
#The default value for shiny.server.sessionTimeout is 30 minutes (1800 seconds)
# By default, the application timeout is set to 300 seconds (5 minutes) on shinyapps.io2
options(shiny.server.sessionTimeout = 1800) # If not, session stops at 5m
options(shiny.server.applicationTimeout = 1800) # If not, session stops at 5m


        # runjs("$('#sum_plot').parent().hide()")

        # Set initial blank plots to UI.
        output$plotgraph5 = renderPlot({ })
        output$plotgraph1 = renderPlot({ })
        output$plotgraph3 = renderPlot({ })
 

          #Analysis types (sidebar menus):
          alphaServer("alpha")
          abundanceServer("abundance")
          betaServer("beta")

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

            showBuiltinBacteriaData <- function(failed = FALSE) {
              modalDialog(
                easyClose=TRUE,
                size="l",
                footer = tagList(modalButton("Dismiss")),
                span(HTML(paste(tags$h1("Built-in abundance bacterial file"),tags$strong('(NOTE: File has to be tab delimited.)')))),
                tags$hr(),
                withSpinner( tagList( renderDT({DT_style(BDtaxa_micro)}, server = FALSE) ), color="green") 
              )
            }

          # React to show built-in dataset files.
          observeEvent(input$show1, {
            showModal(showBuiltinMetadata())
          })
          observeEvent(input$show2, {
            showModal(showBuiltinFungiData())
          })

        # vals=reactiveValues()


# new - old : sorting
# 3 - 1
# 4 - 2
# 5 - 3
# 6 - 4
# 1 - 5
# 2 - 6

# Make all "plotgraph" update when any other module update these plots.
# "vals$plotgraphX" are the plots "shared" across modules.
observeEvent(vals$plotgraph5, { output$plotgraph5 <- renderPlot({ vals$plotgraph5 }) })
observeEvent(vals$plotgraph1, { output$plotgraph1 <- renderPlot({ vals$plotgraph1 }) })
observeEvent(vals$plotgraph3, { output$plotgraph3 <- renderPlot({ vals$plotgraph3 }) })


 # ####### PLOT Alpha diversity #######
 # #### Alpha diversity in fungi partition mycobiome is significantly higher than that in bacteria and control partitions
summary_plot5 <- reactive({
          p=plot_alpha_diversity(physeq = vals$taxa_myco_phy, method = "Shannon",attribute = colnames(metadata_df)[1])  
          vals$plotgraph5<<-p[["plot"]]
          vals$wtable<<-p[["wtable"]]
          p[["plot"]]
})


###### PLOT Taxa Barplots #######
summary_plot1 <- reactive({
          p=plot_high_abundance(physeq = vals$taxa_myco_phy, level="Phylum", title = "Fungal composition", legend = "bottom", facet=colnames(metadata_df)[1])
          vals$plotgraph1<<-p
          p          
          })

# ####### PLOT Beta diversity #######
# ### Fungi   
summary_plot3 <- reactive({
          p=plot_beta_diversity(physeq =  vals$taxa_myco_phy, method = "bray", weighted = T, variable = colnames(metadata_df)[1]) 
           vals$plotgraph3 <<- p[["plot"]]
           vals$pvalue <<- p[["pvalue"]]
           p[["plot"]]
 })


# Make all plots using Built-in Dataset:
observeEvent(input$loadbd, {

        shinyjs::hide("intro_text",anim=T) # Collapses the box, don't hide
        shinyjs::showElement("sum_plot",anim=T)

        # toggle("intro_text") # Collapses the box, don't hide
        # toggle("sum_plot")

        # runjs("$('#intro_text').parent().hide()")
        # runjs("$('#sum_plot').parent().show()")

        metadata_df <<- BDmetadata_df
        taxa_myco <<- BDtaxa_myco

        vals$taxa_myco_phy <<- BDtaxa_myco_phy

        preloaded <<- T

        vals$plotgraph5 <<- BDplot5
        vals$wtable <<- BDwtable
        diversity_df <<- BDdiversity
        updateSelectInput(session, inputId="alpha-splitby", choices=colnames(metadata_df), selected="country")
        updateSelectInput(session, inputId="alpha-method", selected="Shannon")


        vals$plotgraph1 <<- BDplot1
        updateSelectInput(session, inputId="abundance-splitby", choices=c("none", colnames(metadata_df)), selected="disease")
        updateSelectInput(session, inputId="abundance-taxafill", selected="Phylum")
        


        vals$plotgraph3 <<- BDplot3
        vals$pvalue<<-BDpvalue
        updateSelectInput(session, inputId="beta-splitby", choices=c("none", colnames(metadata_df)), selected="country")

  })


# Make all plots if all 2 files are provided:
observeEvent(input$file2, {

        shinyjs::hide("intro_text") # Collapses the box, don't hide
        shinyjs::showElement("sum_plot")

        # toggle("intro_text") # Collapses the box, don't hide
        # toggle("sum_plot")

        # runjs("$('#intro_text').parent().hide()")
        # runjs("$('#sum_plot').parent().show()")


        # In order to make plots, file1 (metadata) and file2 (fungi abundance), must be provided first.
        
        req(input$file1) 

        # Updates metadata_df, taxa_myco and taxa_myco_phy
        # "taxa_myco_phy" should trigger "observeEvent(c(input$method,input$splitby,vals$taxa_myco_phy)"
        # to update the plot
        readData(input,session) # Update global vars with new input files
        
        output$plotgraph5 <- renderPlot({ summary_plot5() })
        updateSelectInput(session, inputId="alpha-splitby", choices=colnames(metadata_df), selected=colnames(metadata_df)[1])

        output$plotgraph1 <- renderPlot({ summary_plot1() })
        updateSelectInput(session, inputId="abundance-splitby", choices=c("none",colnames(metadata_df)), selected=colnames(metadata_df)[1])

        output$plotgraph3 <- renderPlot({ summary_plot3() })
        updateSelectInput(session, inputId="beta-splitby", choices=c("none", colnames(metadata_df)), selected=colnames(metadata_df)[1])

        preloaded <<- F

  })


# Download plots hanlders:

#PDF (up to 15s)
output$pdf <- downloadHandler(
    filename="summary.pdf",
    content = function(file){
      pdf(file,width=input$w,height=input$h,onefile = TRUE) #units: in
          print(vals$plotgraph5)
          print(vals$plotgraph1)
          print(vals$plotgraph3)
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
    list_plots[[1]] <- vals$plotgraph5
    list_plots[[2]] <- vals$plotgraph1
    list_plots[[3]] <- vals$plotgraph3
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
    list_plots[[1]] <- vals$plotgraph5
    list_plots[[2]] <- vals$plotgraph1
    list_plots[[3]] <- vals$plotgraph3
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
    list_plots[[1]] <- vals$plotgraph5
    list_plots[[2]] <- vals$plotgraph1
    list_plots[[3]] <- vals$plotgraph3
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

} #end server


##########################################################
######################### RUN APP ########################
##########################################################

# Run the application 
shinyApp(ui = ui, server = server)

