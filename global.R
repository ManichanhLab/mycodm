##########################################################
################## GLOBAL BUILT-IN DATASETS ##############
##########################################################

### DEBUGGING ###
# ## Shiny Options -  https://shiny.rstudio.com/reference/shiny/1.0.1/shiny-options.html
# #install.packages("reactlog")
# library(reactlog)
# options(shiny.reactlog = TRUE) #debug #Use Ctrl+F3 within the app to show reactive values. Ref: https://shiny.rstudio.com/reference/shiny/1.0.1/showreactlog
# options(shiny.trace=TRUE) #debug # display in the console the messages between server and ui
# options(shiny.fullstacktrace=TRUE) #debug
options(shiny.reactlog = FALSE)
options(shiny.trace=FALSE)
library(scriptName) #debug
DEBUGGING=T
# getCodeLine <- function ()  #Ref: https://stackoverflow.com/questions/59537482/how-to-get-line-number-of-a-function-call-in-r
# {
#     x <- .traceback(x = 1)

#     srcloc <- if (!is.null(srcref <- attr(x[[1]], "srcref"))) {
#         srcfile <- attr(srcref, "srcfile")
#         paste0("Called from ", x[[2]], ", at ", basename(srcfile$filename), "#", srcref[1])
#     }

#     cat(srcloc, "\n")
# }

### LIBRARIES REQUIRED FOR THE ANALYSIS ###
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(DT)        # format interactive tables
library(flextable) # format for stats tables 
library(rlang) # "parse_quo" function
library(vegan) # "adonis2" tests


# Shiny limitation:
# MAX UPLOAD FILE size: 5MB.
# (Ref: https://mastering-shiny.org/action-transfer.html)
# (Ref. File Upload: https://shiny.rstudio.com/gallery/file-upload.html)


##########################################################
#################### GLOBAL VARIABLES ####################
##########################################################

######### CURRENT DATA VALUES (non-reactive) ######### 

metadata_df <- NULL
taxa_myco   <- NULL

diversity_df <- NULL
preloaded <- F

######### LIVE DATA VALUES (reactive) ######### 
#"Reactive values"(input$ are reactive values): trigger reactive expressions ("reactive({})") or
# outputs to update when they change. "reactive expressions" can be called from whitn "renderXXX" functions.
vals <- reactiveValues(  # live/reactive variables that will trigger a render/plot update.
          taxa_myco_phy = NULL,

          plotgraph1 = NULL,
          plotgraph1Legend = NULL,

          plotgraph5 = NULL,
          wtable = NULL,

          plotgraph3 = NULL,
          pvalue = NULL,
        )

##########################################################
#################### GLOBAL FUNCTIONS ####################
##########################################################
source("module/helpers.R",encoding = "utf-8")


##########################################################
######################  PRELOADING #######################
##########################################################

######### PRELOAD DATASET ######### 
      BDplot1 <- NULL

      BDplot5 <- NULL
      BDwtable <- NULL
      BDdiversity <- NULL

      BDplot3 <- NULL
      BDpvalue <- NULL

      ######### METADATA INPUT FILE ######### 
      #It has to be comma "," delimited.
      BDmetadata_df <- read.csv("input/metadata_1156samples.csv", header = T,row.names = 1,sep = ",") #125 KiB
      BDmetadata_df %>% select_if(is.character) -> BDmetadata_df  # Filter out non-string cols (non-categorical)

      ######### TAXA INPUT FILES ######### 
      BDinput <- "input/DF_fungi_buglist_nonStratified.csv"

      #It has to be tab "\t" delimited.
      BDtaxa_myco  <- read.csv(BDinput, header = T,row.names = 1,sep = "\t") # 424 KiB

      #Import buglists (relies on BDinput and BDmetadata_df)
      BDtaxa_myco_phy <- import_buglist_phyloseq(path_buglist = BDinput,
                                               meta = BDmetadata_df, delimiter = '\t')

######### PRELOAD PLOTS ######### 

      printd("...PLOTTING BD_SUMMARY_PLOTS...")

      # ####### PLOT Alpha diversity #######
       # #### Alpha diversity in fungi partition mycobiome is significantly higher than that in bacteria and control partitions
      BD_summary_plot5 <- function() {                
                p=plot_alpha_diversity(physeq = BDtaxa_myco_phy, method = "Shannon",attribute = "country")
                BDplot5 <<- p[["plot"]]
                BDwtable <<- p[["wtable"]]
                BDdiversity <<- p[["diversity"]]
                assign("bdsp5", BDplot5, envir = globalenv()) #debug
      }

      ###### PLOT Taxa Barplots #######
      BD_summary_plot1 <- function() {
                p=plot_high_abundance(physeq = BDtaxa_myco_phy, level="Phylum", title = "Fungal composition", legend = "bottom", facet="disease")
                BDplot1 <<- p
                assign("bdsp1", p, envir = globalenv()) #debug

      }

      # ####### PLOT Beta diversity #######
      # ### Fungi   
      BD_summary_plot3 <- function() {
                printd("INITIAL BETA-DIV BD_PLOT ")
                p=plot_beta_diversity(physeq =  BDtaxa_myco_phy, method = "bray", weighted = T, variable = "country")
                BDplot3 <<- p["plot"]
                BDpvalue <<- p["pvalue"]
                assign("bdsp3", p, envir = globalenv()) #debug
       }


      ## PRELOAD PLOTS
      BD_summary_plot5()
      BD_summary_plot1()
      BD_summary_plot3()

