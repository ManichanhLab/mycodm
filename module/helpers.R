##########################################################
#################### GLOBAL FUNCTIONS ####################
##########################################################
# Non-reactive functions available to each user session (they are created only once):


# DEBUG function
printd <- function(text,DEB=T) {
  if(DEBUGGING & DEB) {
    print(paste("## DEBUG [",Sys.time(),"] ## ",text))
    # scriptName::current_filename()
    # scriptName::current_source_filename()
    # scriptName::current_cli_filename()
    # getCodeLine()
  }
}

# DEBUG function
# Display debugging messages in R (if local) and in the console log (if running in shiny)
# Ref: https://debruine.github.io/shinyintro/debugging.html
debug_msg <- function(...) {
  is_local <- Sys.getenv('SHINY_PORT') == ""
  in_shiny <- !is.null(shiny::getDefaultReactiveDomain())
  txt <- toString(list(...))
  if (is_local) message(txt)
  if (in_shiny) shinyjs::runjs(sprintf("console.debug(\"%s\")", txt))
}

# DEBUG. Code to get shiny app variables into the R console after ShinyApp finish
# Ex: assign("test_p_table", p_table, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.

# DEBUG.:
# observeEvent(input$browser,{
#   browser()
# })


import_buglist_phyloseq <- function(path_buglist, meta, delimiter) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO GENERATE A PHYLOSEQ OBJECT FOR TAXA ABUNDANCE TABLE 
  #
  # parameters:
  #   - path_buglist: path of the taxa ab table 
  #       rownames are full leanage of taxonomies (species level), 
  #         (e.g. k__Fungi|p__Ascomycota|c__Dothideomycetes|o__Botryosphaeriales|f__Botryosphaeriaceae|g__Neofusicoccum|s__Neofusicoccum parvum); 
  #       colnames are sampleID
  #   - meta: metadata table
  #   - delimiter: delimiter of taxa ab table
  #
  # returns: a phyloseq object that contains otu table (in relative abundance), taxa table, metadata, phy tree 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  library(phyloseq)
  library(ape)
  buglist <- read.csv(path_buglist,header = T, sep = delimiter, row.names = 1, check.names = F)
  buglist[is.na(buglist)] <- 0
  meta[is.na(meta)] <- "NA"
  buglist <- buglist %>%  sweep(2,colSums(.),"/")
  buglist <- buglist*100
  if(exists("taxa_Table")) rm(taxa_Table)
  taxa_Table <- extract_taxaTable(buglist)
  taxa_Table = tax_table(taxa_Table)
  colnames(taxa_Table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxa_names(taxa_Table) <- rownames(buglist)
  otable <- otu_table(buglist,taxa_are_rows = T)
  mdata <- sample_data(meta)
  phyObject <- phyloseq(otable, mdata, taxa_Table)
  random_tree = rtree(ntaxa(phyObject), rooted=TRUE, tip.label=taxa_names(phyObject))
  phyObject = merge_phyloseq(phyObject,random_tree, mdata)
  return(phyObject)
}

extract_taxaTable <- function(buglist){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO EXTRACT A TAXA TABLE FROM TAXA AB TABLE
  #
  # parameters:
  #   - buglist: taxa ab table, rownames are full leanage of taxonomies (species level), colnames are sampleID
  #
  # returns: a taxa table for creating phyloseq object 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(stringr)
  if(exists("taxa_Table")) rm(taxa_table)
  for (i in row.names(buglist)) {
    tax_list <- str_split(i,"\\|", simplify = T)
    kingdom <- str_split(tax_list[1], pattern = "k__", simplify = T)[2]
    phylum<- str_split(tax_list[2], pattern = "p__", simplify = T)[2]
    class<- str_split(tax_list[3], pattern = "c__", simplify = T)[2]
    order<- str_split(tax_list[4], pattern = "o__", simplify = T)[2]
    family<- str_split(tax_list[5], pattern = "f__", simplify = T)[2]
    genus <- str_split(tax_list[6], pattern = "g__", simplify = T)[2]
    specie <- str_split(tax_list[7], pattern = "s__", simplify = T)[2]
    
    if (!exists("taxa_table")){
      taxa_table <- data.frame(kingdom, phylum, class, order, family, genus, specie, stringsAsFactors = F)
    } 
    else if (exists("taxa_table")){
      taxa_table <- rbind(taxa_table, c(kingdom, phylum, class, order, family, genus, specie))
    }
  }
  rownames(taxa_table) <- rownames(buglist)
  return(taxa_table)
}


filter_prevalence_abundance_dataframe <- function(unfilter_dataframe, abundance = 0.00, prevalence = 0.1){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO FILTER AN ABUNDANCE TABLE BASED ON RELATIVE ABUNDANCE AND PREVALANCE
  #
  # parameters:
  #   - unfilter_dataframe: unfiltered ab table, rownames are items, colnames are sampleID
  #   - abundance: minimum relative abundance to keep a row
  #   - prevalence: minimum prevalence to keep a row
  #
  # returns: a filted abundance table
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  rows_to_keep=c()
  number_of_samples=dim(unfilter_dataframe)[2]
  colsums_vector = colSums(unfilter_dataframe)
  for (i in 1:dim(unfilter_dataframe)[1]) {
    row_vector = unfilter_dataframe[i,]
    relabun_row_vector = row_vector/colsums_vector
    num_over_abundance = sum((relabun_row_vector > abundance) == TRUE, na.rm = TRUE)
    if (num_over_abundance/number_of_samples > prevalence) {
      rows_to_keep = c(rows_to_keep, i)
    }
  }
  filtered_dataframe <- unfilter_dataframe[rows_to_keep,]
  filtered_dataframe <- filtered_dataframe[colSums(abs(filtered_dataframe), na.rm = TRUE) > 0]
  return(filtered_dataframe)
}

import_pwylist_phyloseq <- function(path_pwylist, meta, delimiter, ab, prev){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO GENERATE A PHYLOSEQ OBJECT FOR FUNCTION ABUNDANCE TABLE 
  #
  # parameters:
  #   - path_pwylist: path of the function ab table (raw counts) 
  #       rownames are pathways,
  #         ( e.g. Fatty acid biosynthesis); 
  #       colnames are sampleID
  #   - meta: metadata table
  #   - delimiter: delimiter of the function ab table
  #   - ab: minimum relative abundance
  #   - prev: minimum prevalance
  #
  # returns: a phyloseq object that contains TMM normalized function table, metadata, and a random tree
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  library(phyloseq)
  library(edgeR)
  library(ape)
  pwylist <- read.csv(path_pwylist,header = T, sep = delimiter, row.names = 1, check.names = F)
  pwylist[is.na(pwylist)] <- 0
  flt_pwy <- filter_prevalence_abundance_dataframe(pwylist,abundance = ab,prevalence = prev)
  pwy_dge <- DGEList(flt_pwy,lib.size = colSums(flt_pwy))
  pwy_dge <- calcNormFactors(pwy_dge, method = "TMM")
  pwy_tmm <- as.data.frame(cpm(pwy_dge, log = F))
  meta[is.na(meta)] <- "NA"
  mdata <- sample_data(meta)
  ptable <- otu_table(pwy_tmm, taxa_are_rows = T)
  phy_function <- phyloseq(ptable, mdata)
  random_tree = rtree(ntaxa(phy_function), rooted=TRUE, tip.label=taxa_names(phy_function))
  phy_function = merge_phyloseq(phy_function, mdata, random_tree)
  return(phy_function)
}

get_core <- function(table, abundance = 0.001, prevalence = 0.9){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO get core taxa/function
  #
  # parameters:
  #   - table: abundance table
  #   - abundance: minimum relative abundance, default is 0.001
  #   - prevalence: minimum prevalence, default is 0.9
  #
  # returns: three tables: shannon index, chao1 index, observed species numbers
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  rows_to_keep=c()
  number_of_samples=dim(table)[2]
  colsums_vector = colSums(table)
  for (i in 1:dim(table)[1]) {
    row_vector = table[i,]
    relabun_row_vector = row_vector/colsums_vector
    num_over_abundance = sum((relabun_row_vector > abundance) == TRUE, na.rm = TRUE)
    if (num_over_abundance/number_of_samples > prevalence) {
      rows_to_keep = c(rows_to_keep, i)
    }
    table[i,"prevalence"] <- num_over_abundance/number_of_samples*100
  }
  filtered_dataframe <- table[rows_to_keep,]
  return(list(filtered_dataframe, rows_to_keep,table))
}

plot_high_abundance <- function (physeq, x = "Sample", y = "Abundance", level = NULL, 
                                 title = NULL, lab = FALSE, n_taxa=20, legend = c("right", "bottom"), facet=NULL) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO CREATE A TAXONOMY BARPLOT FROM A PHYLOSEQ OBJECT EXCLUDING
  # LOW ABUNDANT GROUPS
  #
  # parameters:
  #   - physeq: phyloseq object
  #   - x: x axis variable. Default: "Sample"
  #   - y: y axis variable. Default: "Abundance"
  #   - fill: taxa level to display. For example: "Species"
  #   - title: plot title if needed
  #   - lab: display x axis labels. Default: FALSE
  #   - ab: relative abundance threshold. Default: 1.0
  #   - legend: legend position
  #
  # returns: ggplot barplot 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  palette_c25 <- c(
    "#bfd1d0","#f4d7de","#f4ceb8","#ebe1b9","#e2cbd9","#a4d4dc","#c2a2c2","#f3f3ab","#d7e5ec","#e9a7b8",
    "#fccfb8","#cabed8","#738b8a","#e0f0e3","#fcecc0","#da9076","#ffdada","#647c8b","#caa26b","#f4bbad",
    "#E3971F","#5A7C4D","#8B5357","#123B57","#3B6F7C", "#A27554",
    "#CAC3BD","#646F83","#A08D75","#CDA632","#4F7792","#3A530D",
    "#E2C48D","#88A0B5","#444349","#C4BBBE","#72939E","#C25E7B",
    "#9C918E","#887434","#D48B28","#636B83","#5B8BAF","#E5A93C",
    "#D7D6D7","#BB2649",
    "#CE9C9D", "#F6F3EE", "#DFD8AB","#EDECEB",
    "#ADB5BE","#EAE8EB",  "#FBEFE3", "#FCE4E2",
    "#FAE8E8", "#F8E6E6", "#FAE9E2","#D5A1A3",
    "#E9BDBE", "#FAF3EB", "#EAEEE0","#B5A28A", 
    "#D9E6EC", "#EBF6FA", "#A4C8D5", "#80ADBC",
    "#AED9EA", "#C0D3D8", "#E7E7DB",
    "#ECDCDC", "#EEE8E8", "#E9CCC4",
    "#EFDFD5", "#DEE9EB", "#F0EDE8", "#C2B7B1", "#E1D7CD",
    "#D7E0E5", "#D4C8B6", "#D9BD9C",
    "#CAB08B", "#DEC8BD", "#B4C6DC"
  ) #https://www.tadmint.com/blog/8-pastel-color-palettes-inspired-by-nature
  
  
  # palette_c25 <- c( "#E3971F","#5A7C4D","#8B5357","#123B57","#3B6F7C", "#A27554",
  # "#CAC3BD","#646F83","#A08D75","#CDA632","#4F7792","#3A530D",
  # "#E2C48D","#88A0B5","#444349","#C4BBBE","#72939E","#C25E7B",
  # "#9C918E","#887434","#D48B28","#636B83","#5B8BAF","#E5A93C",
  # "#D7D6D7","#BB2649")
  
  rel_ab <- function(x){
    # compute relative abundance
    if (sum(x) == 0){
      return(x)
    } else {
      return(100 * x/sum(x))
    }
  }
  
  # Transform to relative abundance
  AID_norm <- transform_sample_counts(physeq, rel_ab)
  # print(otu_table(AID_norm))
  top20_species <- top_taxa(AID_norm,n=n_taxa)
  
  #Compile taxa by rank (filtering out low abundance taxa)
  AID_Rank <- AID_norm  %>%
    tax_glom(taxrank = level) %>%   # agglomerate taxa at level
    psmelt() %>%     # Melt phyloseq object to long format for producing graphics with ggplot2
    mutate(tax = replace(get(level), ! OTU %in% top20_species, 'Others'))    # Filter out taxa below ab threshold in each sample
  
  
  ## order taxa levels by abundance 
  
  taxa_ordered <- AID_Rank %>%  
    select(Abundance, tax) %>% 
    group_by(tax) %>% 
    summarise(Total=sum(Abundance)) %>% 
    ungroup() %>% 
    arrange(-Total) %>% 
    filter(tax != 'Others') %>%
    pull('tax')
  
  AID_Rank$tax <- factor(AID_Rank$tax, levels=append(taxa_ordered, 'Others'))
  
  # order x axis by abundance
  
  # print(sum(is.na(AID_Rank[, x])))
  
  most_abun <- taxa_ordered[1]
  
  x_ordered <- AID_Rank %>%
    select(Abundance, tax, x) %>%
    filter(tax == most_abun) %>%
    arrange(Abundance) %>% 
    pull(get(x))
  
  x_missing <- as.vector(AID_Rank[!(AID_Rank[, x] %in% x_ordered), x])
  # print(x_missing)
  
  a <- append(x_ordered, x_missing)
  a <- a[!duplicated(a)]
  
  AID_Rank[, x] <- factor(AID_Rank[, x], levels=a)
  
  
  # rename level colname with the level name
  AID_Rank <- AID_Rank %>%
    select(-c(ncol(AID_Rank)-1)) 
  
  colnames(AID_Rank)[ncol(AID_Rank)] <- level
  
  # print(AID_Rank)
  p = ggplot(AID_Rank, aes_string(x = x, y = y, fill = level))
  p = p + geom_bar(stat = "identity", position = position_stack(reverse = TRUE))
  p = p + theme_prism(base_fontface = "bold", base_line_size = 1, base_family = "Arial",base_size = 14)+ theme(
    legend.position = legend,
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14),
    legend.spacing.x = unit(0, "pt"),
    legend.text = element_text(margin = margin(r = 20))
  ) + guides(fill = guide_legend(override.aes = list(size = 3)))+ scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    breaks = seq(0, 100, 25),
    guide = "prism_offset"
  ) + scale_fill_manual(values=rep(palette_c25, 10)) 
  
  
  if (lab == TRUE){
    p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
  } else {
    p = p + scale_x_discrete(labels = NULL)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (legend == "bottom"){
    p = p + theme(legend.position="bottom") + guides(fill = guide_legend(ncol=2))
  }

  p <- p + {if(!is.null(facet)) facet_wrap(as.formula(paste0("~ factor(",facet,")")), scales= "free_x", nrow=1)}

  return(p) 
}
  

plot_beta_diversity <- function(physeq, method="bray", weighted=T, variable=NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FUNCTION TO FILTER AN ABUNDANCE TABLE BASED ON RELATIVE ABUNDANCE AND PREVALANCE
  #
  # parameters:
  #   - physeq: phyloseq object
  #   - method: distance metric, "unifrac" or "bray"
  #   - weighted: when "unifrac" is chosen, you have to select weighted or no, default is TRUE
  #   - variable: the metadata variable that you want to use to color the data points
  #
  # returns: a filtered abundance table
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # if (method=="unifrac") {
  #   if (weighted == T) {
  #     wunifrac_dist <- phyloseq::distance(physeq, method="unifrac", weighted=T)
  #   } else if (weighted == F) {
  #     wunifrac_dist <- phyloseq::distance(physeq, method="unifrac", weighted=F)
  #   }
  #   ordination <- ordinate(physeq, method="PCoA", distance=wunifrac_dist)
  # } else {
  #   ordination <- ordinate(physeq, method="PCoA", distance=method)
  # }

  printd("--INSIDE plot_beta_diversity--")

  ordination <- ordinate(physeq, method="PCoA", distance=method)

  p <- plot_ordination(physeq, ordination, color=variable)
  assign("beta.vals_plotgraph3.plot_ordination", p, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.

  p <- p +  theme(aspect.ratio=1)+ 
            theme_classic()+
            scale_fill_brewer(palette = "Dark2")+
            scale_color_brewer(palette = "Dark2") +
            theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
            {if(!is.null(variable)) stat_ellipse(aes_string(fill = variable), geom="polygon",level=0.95,alpha=0.2)}
            # stat_ellipse(aes(fill = .data[[input$splitby]]), geom="polygon",level=0.95,alpha=0.2)

  ### Statistics
  # adonis2 example (rows are samples; dune.env: metadata (see colnames(dune.env) for vars to use in adonis) ):
  # library(vegan)
  # data(dune)
  # data(dune.env)
  # adonis2(dune ~ Management, data = dune.env)

  ares <- "null or 'none' group selected"
  printd("RUNNING plot_beta_diversity # BEFORE IF STATEMENT")
  print(variable) #debug
  printd("#--++#") 

  # "&&" evaluates the conditions from left to right and stops as soon as one of them is FALSE
  if(!is.null(variable) && variable != "none") {
      # otu.tb <- data.frame(otu_table(physeq))
      # printd(paste("meta dim:",dim(metadata_df)))
      otu.tb <- data.frame(otu_table(physeq))
        # metadata <- as(sample_data(phyloseq), "data.frame")

      # parse_var <- parse(text = variable)[[1]]
      # adonis_res <- adonis2(as.data.frame(t(otu.tb)) ~ eval(parse_var), data = data.frame(sample_data(physeq)))
      adonis_res <- adonis2(as.formula(paste0(" as.data.frame(t(otu.tb)) ~ ",variable)), data = data.frame(sample_data(physeq)))

      # adonis_res <- adonis2(as.data.frame(t(otu.tb)) ~ variable, data = meta)
      # head(metadata_df,4)
      # adonis_res <- adonis2(as.data.frame(t(otu.tb)) ~ country, data = metadata_df)
      # variable<-NULL
      # adonis_res <- adonis2(as.formula(paste0("as.data.frame(t(otu.tb)) ~ ",variable)), data=metadata_df)
      # adonis_res <- adonis2(as.data.frame(t(otu.tb)) ~ country, data=meta)
      printd(paste("p-value: ", adonis_res$`Pr(>F)`[1]))
      ares<-adonis_res$`Pr(>F)`[1]
  }

    assign("pbd_plot", p, envir = globalenv()) #debug
    assign("pbd_ares", ares, envir = globalenv()) #debug

  return( list("plot"=p, "pvalue"=ares) )
  # return(p)
}

plot_alpha_diversity <- function(physeq, method, attribute) {
            
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # FUNCTION TO CALCULATE ALPHA DIVERSITIES OF A PHYLOSEQ OBJECT
                #
                # parameters:
                #   - physeq: phyloseq object
                #   - method: alpha diversity index, "Shannon", "Chao1", "Obs"
                #   - attribute: the metadata variable that you want to use to color the data points
                #
                # returns: three tables: shannon index, chao1 index, observed species numbers
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                printd("--INSIDE plot_alpha_diversity--")

                round_meta <- data.frame(sample_data(physeq))
                round_otu <- data.frame(otu_table(physeq))
                round_otu <- round_otu %>% mutate_if(is.numeric,round)
                round_otu <- round_otu[colSums(round_otu) > 0]
                round_otu <- round_otu[rowSums(round_otu) > 0,]
                otable <- otu_table(round_otu,taxa_are_rows = T)
                mdata <- sample_data(round_meta)
                round_taxa <- extract_taxaTable(round_otu) %>% as.matrix() %>% tax_table()
                phy_round <- phyloseq(otable,mdata,round_taxa)
                random_tree = rtree(ntaxa(phy_round), rooted=TRUE, tip.label=taxa_names(phy_round))
                phy_round = merge_phyloseq(phy_round, mdata, random_tree)
                
                diversity <- plot_richness(phy_round, measures = method)$data
                
                printd("Diversity table:")
                printd(dim(diversity))
                print(head(diversity)) #debug

                # parse_var <- parse(text = attribute)[[1]]
                # Use "!!" and "sym()" to use strings as symbols/colnames
                
                # Get more palette colors in case grouping includes up to 20 categories:
                # Define the number of colors you want
                ncolors = 20
                # To get the total number of palette colors, check: ?ColorBrewer
                getPalette1 = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(ncolors)
                getPalette1 = union(RColorBrewer::brewer.pal(8, "Dark2"), getPalette1) # Keep Dark2 original colors first

                # p <- ggplot(diversity,aes(x=reorder(eval(parse_var),value,na.rm=T),y=value)) + 
                p <- ggplot(diversity,aes(x=reorder(!!sym(attribute),value,na.rm=T),y=value)) + 
                  geom_boxplot(width=0.5,lwd=1.5) +
                  # geom_jitter(width=0.3,aes(color=eval(parse_var))) +
                  geom_jitter(width=0.3,aes(color=!!sym(attribute))) +
                  #scale_color_brewer(palette="Dark2")+theme_classic()+
                  scale_colour_manual(values=getPalette1)+theme_classic()+
                  theme(text = element_text(size=25),
                        axis.text.x = element_text(angle=0, hjust=1)) +
                  xlab(attribute) +
                  ylab(method) +
                  theme(axis.text.x = element_text(angle = 30))

                ### statistics
                ls_grps <- unique(diversity[,attribute])
                p_table <- data.frame(matrix(NA, nrow = length(ls_grps), ncol = length(ls_grps)))
                colnames(p_table) <- rownames(p_table) <- ls_grps
                for (i in 1:length(ls_grps)) {
                  for (j in 1:length(ls_grps)) {
                    grp1 <- diversity[diversity[,attribute]==ls_grps[i],][,"value"]
                    grp2 <- diversity[diversity[,attribute]==ls_grps[j],][,"value"]
                    pvalues <- wilcox.test(grp1,grp2)$p.value
                    p_table[i,j] <- pvalues
                  }
                }
  
                assign("test_p_table", p_table, envir = globalenv()) #debug. Save Shiny var outside Shiny env, in the console.
                printd("raw p-table:")
                print(head(p_table)) #debug

                # p_table[!lower.tri(p_table, diag = FALSE)] <- "" # This operation converts all numeric columns to string
                p_table[!lower.tri(p_table, diag = FALSE)] <- NA

                printd("p-table after filtering diagonal:")
                print(head(p_table)) #debug

                return(list("plot"=p, "wtable"=p_table, "diversity"=diversity))
}


### Change datatable row selection color:
css_row_select <- "table.dataTable tr.selected td, table.dataTable td.selected {
                          box-shadow: inset 0 0 0 9999px #00BE9A !important;
                         }"

### Styling datatables:
# metadata_df %>%  datatable(., extensions = 'Buttons', ...)
# now:  dt_style(metadata_df)
DT_style <- function(df) {
              datatable(df, extensions = 'Buttons',
                              options = list(
                                scrollX = T,
                                #dom = 'Bfrtip', # without row number select
                                dom = 'Blfrtip',
                                buttons = list(

                                  # I('colvis'),
                                  # list(extend = 'copy', text = "Clipboard", title = NULL), 
                                  # list(extend = 'csv'), # title = "Title"), 
                                  # list(extend = 'excel'), #, title = "Title"), 
                                  # list(extend = 'pdf', size='A4',orientation='landscape') #, title = "Title")
                                  'colvis', list(extend = 'copy',text="Clipboard"),'csv','excel', 
                                  list(
                                    extend= 'collection',
                                    text= 'Print',
                                    className= "btnArrow",
                                     buttons= list(
                                      list(
                                        extend= "print" #,
                                        # exportOptions= list(columns= ':visible', rows= ':visible')
                                        ),
                                     list(
                                        extend= "pdf", text="PDF portrait A4", orientation="portrait", pageSize='A4' #,
                                          ),
                                     list(
                                        extend= "pdf", text="PDF landscape A4", orientation="landscape", pageSize='A4'#,
                                          ),
                                     list(
                                        extend= "pdf", text="PDF landscape A3", orientation="landscape", pageSize='A3' #,
                                          ),
                                     list(
                                        extend= "pdf", text="PDF landscape A0", orientation="landscape", pageSize='A0' #,
                                          )
                                     )
                                     )


                                ),
                                deferRender = TRUE,
                                # scroller = TRUE,
                                # 'background-color': '#00a65a', 
                                # 'color': '#fff'
                                initComplete = JS(
                                                    "function(settings, json) {",
                                                    "$(this.api().table().header()).css({
                                                          'background-color': '#009966',
                                                          'color': 'white',
                                                          'font-weight': 'bold'
                                                          });",
                                                    "}"
                                                  )
                              )
                    )
}


### Update input files
# Req: Make sure this function "readData()" is called after "vals <- reactiveValues()" is declared in "global.R"
readData <- function(input,session){

        ntasks<-4
        itask=-1

        conv <- function(n) {
            return(as.double(format(round(n, 2), nsmall = 2)))
        }

 ntasks<-ntasks-1; progress <- Progress$new(session)

## Update "metadata_df" ##
printd(conv(itask/ntasks))
itask<-itask+1; progress$set(value = conv(itask/ntasks), message = 'Loading file1...(1/3)'); Sys.sleep(0.5)
        #FILE1 ("metadata_df") DATA
        if(file_ext(input$file1$datapath) == "csv"){
          metadata_df <<- read.csv(input$file1$datapath, header = T,row.names = 1,sep = ",") 
        }else if(file_ext(input$file1$datapath) == "txt"){
          metadata_df <<- read.csv(input$file1$datapath, header = T,row.names = 1,sep = ",") 
        }else if(file_ext(input$file1$datapath) == "xls"){
            metadata_df <<- readxl::read_xls(input$file1$datapath)
            metadata_df <<- as.data.frame(metadata_df)
            rownames(metadata_df) <<- metadata_df[,1]
            metadata_df <<- metadata_df[,-1]
        }else if(file_ext(input$file1$datapath) == "xlsx"){
            metadata_df <<- readxl::read_xlsx(input$file1$datapath)
            metadata_df <<- as.data.frame(d)
            rownames(metadata_df) <<- metadata_df[,1]
            metadata_df <<- metadata_df[,-1]
        }
        metadata_df %>% select_if(is.character) ->> metadata_df

## Update "taxa_myco" ##
printd(conv(itask/ntasks))
itask<-itask+1; progress$set(value = conv(itask/ntasks), message = 'Loading file2...(2/3)'); Sys.sleep(0.5)
        #FILE2 ("taxa_myco") DATA
        if(file_ext(input$file2$datapath) == "csv"){
          taxa_myco <<- read.csv(input$file2$datapath, header = T,row.names = 1,sep = "\t")
        }else if(file_ext(input$file2$datapath) == "txt"){
          taxa_myco <<- read.csv(input$file2$datapath, header = T,row.names = 1,sep = "\t")
        }else if(file_ext(input$file2$datapath) == "xls"){
            taxa_myco <<- readxl::read_xls(input$file2$datapath)
            taxa_myco <<- as.data.frame(taxa_myco)
            rownames(taxa_myco) <<- taxa_myco[,1]
            taxa_myco <<- taxa_myco[,-1]
        }else if(file_ext(input$file2$datapath) == "xlsx"){
            taxa_myco <<- readxl::read_xlsx(input$file2$datapath)
            taxa_myco <<- as.data.frame(d)
            rownames(taxa_myco) <<- taxa_myco[,1]
            taxa_myco <<- taxa_myco[,-1]
        }

## Update "taxa_myco_phy" ##  
# Req: "vals$taxa_myco_phy" must be first declared in global variables
# Use of "metadata_df" and "taxon input file(input$file2)"
printd(conv(itask/ntasks))
itask<-itask+1; progress$set(value = conv(itask/ntasks), message = 'Importing buglists...(3/3)'); Sys.sleep(0.5)
        #Import buglists
        vals$taxa_myco_phy <<- import_buglist_phyloseq(path_buglist = input$file2$datapath,
                                                 meta = metadata_df, delimiter = '\t')

itask<-itask+1; progress$set(value = conv(itask/ntasks), message = 'Data uploaded successfully! ...plotting...'); Sys.sleep(1)  
        progress$close()
}