#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(gridExtra)
library(circlize)
library(DT)
library(ComplexHeatmap)
library(reshape2)
library(shinycssloaders)

table <- read.csv("all_panc_lookup_scores_2.tsv_uniprot_revel.tsv", sep = '\t', stringsAsFactors = FALSE)
table <- table[table$Chromosome != "error searching variant", ]
clinvar <- gsub(table$Clinvar2.accessions, pattern = "\\[|\\'|\\]", replacement = "")
clinvar <- strsplit(clinvar, ", ")
a <- lapply(clinvar, function(x){
  if(length(x) > 1){
    x[1]
  }
  else{
    x
  }
})
table$Clinvar2.accessions <- unlist(a)
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Benign", replacement = "B")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Likely Benign", replacement = "LB")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "drug response", replacement = "dr")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Likely pathogenic", replacement = "LP")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "not provided", replacement = "None")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "other", replacement = "None")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Pathogenic", replacement = "P")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Pathogenic", replacement = "LP")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Conflicting interpretations of pathogenicity", replacement = "cip")
# table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Uncertain significance", replacement = "VUS")

numerictests <- c("Phylop.Prediction", "Dann.score", "aggregatescore", "revel", "MutationTaster.Score","Sift.Score","MutationAssessor.Score","Fathmm.mkl.Score","Fathmm.Score","Metasvm.Score","Metalr.Score","Provean.Score", "lrt.Score")
cattests <- c("Clinvar2.accessions", "interpretation", "MutationTaster.Prediction","Sift.Prediction","MutationAssessor.Prediction","Fathmm.mkl.Prediction","Fathmm.Prediction","Metasvm.Prediction","Metalr.Prediction","Provean.Prediction", "lrt.Prediction")


circostable <- table[, c(6,9,9,2,3,4,5,7,8,10:ncol(table))]
circostable <- circostable[circostable[,1] != "error searching variant", ]
circostable[,c(2,3,7)] <- sapply(circostable[,c(2,3,7)], as.numeric)
circostable[circostable == "None"] <- NA

#initialize chr list for checkboxes
title <- paste("chr", c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "X", "Y"), sep = "")
chr <- as.list(c(1:22,"X","Y"))
names(chr) <- title

mt <- c("['A', 'D']" = 1, "['A']" = 1, "['D', 'N']" = 0, "['D']" = 1, "['N']" = -1, "['P']" = -1)
sift <- c("Damaging" = 1, "Tolerated" = -1)
ma <- c("['H']" = 1, "['L']" = .33, "['M']" = .66, "['N']" = 0)
fathmm <- c("['D']" = 1,"['N']" = -1)
metasvm <- c("['D']" = 1,"['T']" = -1)
metalr <- c("['D']" = 1,"['T']" = -1)
provean <- c("['D']" = 1,"['N']" = -1, "['D', 'N']" = 0)
lrt <- c("['D']" = 1,"['N']" = -1, "['U']" = 0)
tests <- c("MutationTaster.Prediction","Sift.Prediction","MutationAssessor.Prediction","Fathmm.Prediction","Metasvm.Prediction","Metalr.Prediction","Provean.Prediction", "lrt.Prediction")
           #, "Dann.score", "GERP.RS")

#extend function for zoom functionality
extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  if(length(zoom_bed[[1]]) == 0){
    return(bed)
  }
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  bed1 <- rbind(bed, zoom_bed)
  bed1[!bed1[[1]] %in% chromosome, ]}

# Define UI for application
ui <- navbarPage(
   
   # Application title
   "Variant summary",
  
   
   tabPanel("Upload data",
            sidebarLayout(
              sidebarPanel(
                
                fileInput("file", "Upload variant list",
                          multiple = FALSE,
                          accept = c("text/csv/tsv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv",
                                     ".tsv")),
                
                checkboxInput("lookup", label = "Lookup variants?", value = FALSE),
                helpText(h6("Optional: Sample variant list already loaded"))
                ),
                
              # Show a plot of the generated distribution
              mainPanel(
                textOutput("upload_message")
              )
            )
   ),
   
   
   # Circos Tab
   tabPanel("Circos",
            sidebarLayout(
              sidebarPanel(
                conditionalPanel(
                  condition = 'input.circosplots === "Chromosome"',
                  selectInput('circhr', "Chromosome to zoom:",
                              choices = chr, 
                              multiple=TRUE, 
                              selectize=TRUE),
                  actionButton("updatechr", "Plot")
                ),
                conditionalPanel(
                  condition = 'input.circosplots === "Gene"',
                  selectInput('cirgene', label = "Gene to zoom:",
                              choices = unique(table$marker), 
                              multiple=TRUE, 
                              selectize=TRUE
                              ),
                  actionButton("updategene", "Plot")
                ),

                
                checkboxInput("circheckbox", label = "Show points?", value = TRUE),
                br(),
                selectInput("circostest",
                                   "Select a test:",
                                   c("Interpretation"="interpretation", "Mutation Taster"="MutationTaster.Prediction","Sift"="Sift.Prediction", "Mutation Assessor"="MutationAssessor.Prediction","Fathmm"="Fathmm.Prediction","Metasvm"="Metasvm.Prediction","Metalr"="Metalr.Prediction","Provean"="Provean.Prediction","lrt"="lrt.Prediction")
                ),
                radioButtons("circosy",
                                  "y-value",
                                  c("Aggregate score" = "aggregatescore", "Patient Count" = "patient_count")),
                conditionalPanel(
                  condition = 'input.circosplots === "Chromosome"',
                  selectInput('heatcircoschr', "Gene:",
                              choices = c("All",unique(as.character(table$marker))), 
                              multiple=TRUE, 
                              selectize=TRUE),
                  actionButton("heatbuttonchr", "Plot heatmap")
                ),
                conditionalPanel(
                  condition = 'input.circosplots === "Gene"',
                  selectInput('heatcircosgene', "Gene:",
                              choices = c("All",unique(as.character(table$marker))), 
                              multiple=TRUE, 
                              selectize=TRUE),
                  actionButton("heatbuttongene", "Plot heatmap")
                )
              ),
              
              mainPanel(
                tabsetPanel(
                  id = "circosplots",
                  tabPanel("Chromosome",withSpinner(plotOutput("circoschrview"), 4)),
                  tabPanel("Gene",withSpinner(plotOutput("circosgeneview"), 4))
                )
              )
            )
          ),
   
   
   # Heatmap Tab
   tabPanel("Concordance",
            sidebarLayout(
              sidebarPanel(
                selectInput("heattests",
                                   "Select 2 tests to compare:",
                                    multiple=TRUE,
                                    selectize=TRUE,
                                    list("Score" = numerictests,
                                          "Prediction" = c("Interpretation"="interpretation","Mutation Taster"="MutationTaster.Prediction","Sift"="Sift.Prediction", "Mutation Assessor"="MutationAssessor.Prediction","Fathmm"="Fathmm.Prediction","Metasvm"="Metasvm.Prediction","Metalr"="Metalr.Prediction","Provean"="Provean.Prediction","lrt"="lrt.Prediction"))
                ),
                selectInput("heatgene",
                            "Gene:",
                            c("All",unique(as.character(table$marker)))
                ),
                downloadButton('downloadheat', 'Download heatmap')
              ),
              
              # Show a plot of the heatmap
              mainPanel(
                plotOutput("heat")
                )
              )
           ),
   
   
   #Barcharts
   tabPanel("Barcharts",
            sidebarLayout(
              sidebarPanel(
                selectInput("bartests",
                                   "Select a test:",
                                   c("Interpretation"="interpretation","Mutation Taster"="MutationTaster.Prediction","Sift"="Sift.Prediction", "Mutation Assessor"="MutationAssessor.Prediction","Fathmm"="Fathmm.Prediction","Metasvm"="Metasvm.Prediction","Metalr"="Metalr.Prediction","Provean"="Provean.Prediction","lrt"="lrt.Prediction")
                ),
                tags$hr(),
                selectInput('bargene', "Gene:",
                            choices = c(unique(as.character(table$marker))), 
                            multiple=TRUE, 
                            selectize=TRUE),
                sliderInput("barthreshold", label = "Aggregate Score Threshold", min = 0, 
                            max = 8, value = c(0, 8), step = .5),
                downloadButton('downloadbar', 'Download barchart')
              ),
              
              # Show a plot of the heatmap
              mainPanel(
                fluidRow(
                  plotOutput("zoomdist"))
              )
            )
   ),
   
   tabPanel("Prediction Heatmap",
            sidebarLayout(
              sidebarPanel(
                selectInput('predheatgene', "Gene:",
                            choices = c("All",unique(as.character(table$marker))), 
                            multiple=TRUE, 
                            selectize=TRUE),
                actionButton("updatepredheat", "Plot"),
                p("\n"),
                sliderInput("predheatthreshold", label = "Aggregate Score Threshold", min = 0, 
                            max = 8, value = c(0, 8), step = .5),
                checkboxInput("genefacet",label = "Split by gene?", value = FALSE),
                checkboxGroupInput("interp", label = "Clinvar Interpretation",
                               choices = c("All", unique(table$Clinvar2.accessions)), selected = c("All")),
                downloadButton('downloadpredheat', 'Download heatmap')
              ),
              
              # Show a plot of the heatmap
              mainPanel(
                withSpinner(uiOutput("predheatcontainer"), 4)
              )
            )
   ),
   
  
   tabPanel("Table view",
            sidebarLayout(
              sidebarPanel(
                sliderInput("tablethreshold", label = "Aggregate Score Threshold", min = 0, 
                            max = 8, value = c(0, 8), step = .5),
                sliderInput("revelthreshold", label = "Revel Threshold", min = 0, 
                            max = 1, value = c(0, 1), step = .01),
                checkboxGroupInput("table",
                            "Select columns to view:",
                            c("All","None","aggregatescore",colnames(table)),
                            selected = colnames(table)[1:40]
                ),
                downloadButton('downloadtable', 'Download table')
              ),
              
              # Show a plot of the heatmap
              mainPanel(
                dataTableOutput("table")
              )
            )
   )
   
)
# Define server logic
server <- function(input, output, session) {
  
  read_table <- reactive({
    if(!is.null(input$file)){
      if (input$lookup){
        system(paste0("python3 variant_workflow.py -f ",input$file$datapath, " -o ", input$file$datapath, "_lookup.tsv"))
        table <- read.csv(paste0(input$file$datapath, "_lookup.tsv"), sep = '\t', stringsAsFactors = FALSE)
      }
      else{
      table <- read.csv(input$file$datapath, sep = '\t', stringsAsFactors = FALSE)
      }
    }
    else{
      table <- read.csv("all_panc_lookup_scores_2.tsv_uniprot_revel.tsv", sep = '\t', stringsAsFactors = FALSE)
    }
    table <- table[table$Chromosome != "error searching variant", ]
    scores <- predtable(table)
    scores[,tests] <- data.frame(sapply(scores[,tests], function(x){
      x[is.na(x)] <- .5
      x
    }))
    table$aggregatescore <- rowSums(scores[,tests])
    clinvar <- gsub(table$Clinvar2.accessions, pattern = "\\[|\\'|\\]", replacement = "")
    clinvar <- strsplit(clinvar, ", ")
    a <- lapply(clinvar, function(x){
          if(length(x) > 1){
            x[1]
          }
          else{
            x
          }
    })
    table$Clinvar2.accessions <- unlist(a)
    
    for(col in numerictests[numerictests != "revel"]){
      cur <- gsub(table[[col]], pattern = "\\[|\\'|\\]", replacement = "")
      cur <- strsplit(cur, ", ")
      a <- lapply(cur, function(x){
        if(length(x) > 1){
          x[1]
        }
        else{
          x
        }
      })
      table[[col]] <- unlist(a)
      table[[col]] <- as.numeric(as.character(table[[col]]))
    }
    
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Benign", replacement = "B")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Likely benign", replacement = "LB")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "drug response", replacement = "dr")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Likely pathogenic", replacement = "LP")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "not provided", replacement = "None")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "other", replacement = "None")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Pathogenic", replacement = "P")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Pathogenic", replacement = "LP")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Conflicting interpretations of pathogenicity", replacement = "cip")
    # table$Clinvar2.accessions <- gsub(table$Clinvar2.accessions, pattern = "Uncertain significance", replacement = "VUS")
    table
  })
  circostable <- reactive({
    table <- read_table()
    circostable <- table[, c(6,9,9,2,3,4,5,7,8,10:ncol(table))]
    circostable <- circostable[circostable[,1] != "error searching variant", ]
    circostable[,c(2,3,7)] <- sapply(circostable[,c(2,3,7)], as.character)
    circostable[,c(2,3,7)] <- sapply(circostable[,c(2,3,7)], as.numeric)
    circostable[circostable == "None"] <- NA
    circostable
  })
  
  plots <- reactiveValues()
  
#-----------------------------------------------------------------------------------
  #render circos plot
  
  genelist <- reactive({
    genelist <- table
    genelist$aalength <- as.numeric(as.character(gsub("aa [0-9]+ of ", "", genelist$Coding.Location)))
    genelist$aalength[is.na(genelist$aalength)] <- 0
    genelist <- genelist[order(genelist$marker, -rank(genelist$aalength)),]
    genelist <- genelist[!duplicated(genelist$marker),]
    genelist$start <- 0
    genelist <- genelist[,c("marker","start","aalength")]
    genelist <- genelist[genelist$aalength>0,]
    colnames(genelist) <- c("name","start","end")
    if(!is.null(zoomopt$cirgene)){
      genelist[genelist$name %in% zoomopt$cirgene, ]
    }
    else{
      genelist
    }
  })
  circosgenelist <- reactive({
    variants <- circostable()
    variants$Coding.Location <- gsub("aa ", "", variants$Coding.Location)
    variants$position <- as.numeric(as.character(gsub(" of [0-9]+", "", variants$Coding.Location)))
    variantbed <- variants[, c("marker","position","position")]
    variants <- cbind(variantbed, variants[5:ncol(variants)])
    variants <- variants[!is.na(variants$position),]
    if(!is.null(zoomopt$cirgene)){
      variants[variants$marker %in% zoomopt$cirgene, ]
    }
    else{
      variants
    }
  })
  
  zoom_param_chr <- reactive({
    cytoband = read.cytoband()
    cytoband_df = cytoband$df
    chromosome = cytoband$chromosome
    xrange = c(cytoband$chr.len)
    zoom = zoomopt$circhr
    zoom <- gsub("X","23",zoom)
    zoom <- gsub("Y","24",zoom)
    zoom <- as.numeric(zoom)
    normal_chr_index = (1:24)[!(1:24) %in% zoom]
    zoomed_chr_index = zoom
    # normalize in normal chromsomes and zoomed chromosomes separately
    # sector_width
    c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                     xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index]))
  })
  buttons <- reactiveValues()
  observeEvent(input$heatbuttonchr,{
    buttons$heatcircoschr <- input$heatcircoschr
  })
  observeEvent(input$heatbuttongene,{
    buttons$heatcircosgene <- input$heatcircosgene
  })
  
  zoomopt <- reactiveValues()
  
  observeEvent(input$updatechr,{
    zoomopt$circhr <- input$circhr
  })
  observeEvent(input$updategene,{
    zoomopt$cirgene <- input$cirgene
  })

  #render plot
  output$circoschrview <- renderPlot({
    input$update
    circostable <- circostable()
    circos.clear()
    
    if (!is.null(zoomopt$circhr)){
      sector.width <- zoom_param_chr()
      circos.par(start.degree = 90)
      circos.initializeWithIdeogram(extend_chromosomes(read.cytoband()$df, paste("chr",zoomopt$circhr,sep="")),
                                    sector.width = sector.width)
    }
    else{
      circos.initializeWithIdeogram()
      #circos.genomicLabels(labels, labels.column = 4, side = "inside")
    }
    if(input$circheckbox){
      if (!is.null(zoomopt$circhr)){
        t <- extend_chromosomes(circostable, paste0("chr",zoomopt$circhr))
      }
      else{
        t <- circostable
      }
      t <- data.frame(lapply(t, function(x) gsub(pattern = "\\[|\\]|\\'", replacement = "", x = x)))
      test <- input$circostest
      lv <- levels(factor(t[,test]))
      t <- lapply((lv), function(x){
        t[t[[test]] == x, ] 
      })
      t <- lapply(t, function(x){
        x <- cbind(x[,1:3], x[[input$circosy]])
        x[complete.cases(x),]
      })
      t <- lapply(t, function(x){
        x[,2:4] <- sapply(x[,2:4],as.character)
        x[,2:4] <- sapply(x[,2:4],as.numeric)
        x <- data.frame(x)
      })
      circos.genomicTrack(t, track.index = 3,
                          panel.fun = function(region, value, ...) 
                          {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16, cex = .5, col = i, ...)
                          })
      if (!is.null(buttons$heatcircoschr)){
        heattable <- predtable(read_table())
        heattable <- heattable[, c(6,9,9,2,3,4,5,7,8,10:ncol(heattable))]
        if(!"All" %in% buttons$heatcircoschr){
        heattable <- heattable[heattable$marker %in% buttons$heatcircoschr,]
        }
        if (!is.null(zoomopt$circhr)){
          heattable <- extend_chromosomes(heattable, paste0("chr",zoomopt$circhr))
        }
        heattable[,c(2,3,7)] <- sapply(heattable[,c(2,3,7)], function(x) as.numeric(as.character(x)))
        heattable <- heattable[, c("Chromosome","Position","Position.1",tests)]
        heattable[is.na(heattable)] <- 0
        
        col_fun = colorRamp2(c(0, 0.5, 1), c("green", "black", "red"))
        circos.genomicHeatmap(heattable, col = col_fun, side = "inside")
        circos.clear()
      }
      
      lgd_points = Legend(at = lv, type = "points", 
                          legend_gp = gpar(col = 1:length(lv)), title_position = "topleft", 
                          title = test)
      pushViewport(viewport(x = unit(2, "mm"), y = unit(4, "mm"), 
                            width = grobWidth(lgd_points), 
                            height = grobHeight(lgd_points), 
                            just = c("left", "bottom")))
      grid.draw(lgd_points)
      upViewport()
    }
  output$circosgeneview <- renderPlot({
    genelist <- genelist()
      circos.clear()
      circos.par(gap.degree = 0)
      if(!is.null(zoomopt$cirgene)){
        circos.genomicInitialize(genelist, plotType = c("axis", "labels"))
      }
      else{
        circos.genomicInitialize(genelist, plotType = NULL)
      }
      
      circos.par(cell.padding = c(0.02, 0, 0.02, 0))
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        gene = genelist$marker
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        circos.rect(xlim[1], 0, xlim[2], 1, col = (rand_color(1, transparency = .8)))
      }, track.height = 0.15, bg.border = NA)
    if(input$circheckbox){
      genepoints <- circosgenelist()
      scores <- predtable(genepoints)
      scores[,tests] <- data.frame(sapply(scores[,tests], function(x){
        x <- (x+1)/2.0
        x[is.na(x)] <- .5
        x
      }))
      genepoints$aggregatescore <- rowSums(scores[,tests])
      
      test <- input$circostest
      lv <- levels(factor(genepoints[,test]))
      ylim <- genepoints[[input$circosy]]
      genepoints <- lapply((lv), function(x){
        genepoints[genepoints[[test]] == x, ]
      })
      genepoints <- lapply(genepoints, function(x){
        x <- x[, c("marker","position","position",input$circosy)]
        x[complete.cases(x),]
      })
      genepoints <- lapply(genepoints, function(x){
        x[,2:4] <- sapply(x[,2:4],as.character)
        x[,2:4] <- sapply(x[,2:4],as.numeric)
        x <- data.frame(x)
      })

      if (min(ylim, na.rm = TRUE)==max(ylim, na.rm = TRUE)){
        circos.genomicTrack(genepoints, ylim = c(min(ylim),max(ylim)+1),
                            panel.fun = function(region, value, ...) {
                              
                              i = getI(...)
                              circos.genomicPoints(region, value, pch = 16, cex = .5, col = i, ...)
                            })
      }
      else{
        circos.genomicTrack(genepoints,
                            panel.fun = function(region, value, ...) {
                              
                              i = getI(...)
                              circos.genomicPoints(region, value, pch = 16, cex = .5, col = i, ...)
                            })
      }
      lgd_points = Legend(at = lv, type = "points", 
                          legend_gp = gpar(col = 1:length(lv)), title_position = "topleft", 
                          title = test)
      pushViewport(viewport(x = unit(2, "mm"), y = unit(4, "mm"), 
                            width = grobWidth(lgd_points), 
                            height = grobHeight(lgd_points), 
                            just = c("left", "bottom")))
      grid.draw(lgd_points)
      upViewport()
      #circos.track(factor(genepoints[[1]]), genepoints[[2]], genepoints$patient_count, ylim = c(0, max(genepoints$patient_count)), panel.fun = function(x, y) {
      #circos.genomicPoints(x, y, pch = 16, cex = .5)
    }
      if(!is.null(buttons$heatcircosgene)){
        
        if (!identical(intersect(buttons$heatcircosgene, zoomopt$cirgene), character(0))){
        heatgenetable <- predtable(circosgenelist())
        if(!"All" %in% buttons$heatcircosgene){
          heatgenetable <- heatgenetable[heatgenetable$marker %in% buttons$heatcircosgene,]
        }
        heatgenetable <- heatgenetable[,c("marker","position","position.1",tests)]
        heatgenetable[,c(2:ncol(heatgenetable))] <- sapply(heatgenetable[,c(2:ncol(heatgenetable))], function(x) as.numeric(as.character(x)))
        
        heatgenetable[is.na(heatgenetable)] <- 0
        col_fun = colorRamp2(c(0, 0.5, 1), c("green", "black", "red"))
        circos.genomicHeatmap(heatgenetable, col = col_fun, side = "inside")
        circos.clear()
      }
    }
  })
  })
#-----------------------------------------------------------------------------------
  
  #render heatmap
  output$heat <- renderPlot({
    table <- read_table()
    if(length(input$heattests) == 2){
      #filter for selected gene
      
      if (input$heatgene == "All"){
        t <- table
      }
      else{
        t <- table[table$marker == input$heatgene, ]
      }
      if(sum(input$heattests %in% cattests) == 2){
      # select those columns
      t <- t[,c("result",input$heattests)]
      colnames(t)[2:3] <- c(gsub(".Prediction", "", colnames(t)[2:3]))
      #removes brackets and quotes from entries
      t <- data.frame(lapply(t, function(x) gsub(pattern = "\\[|\\]|\\'", replacement = "", x = x)))
      #replaces "-" with NAs
      t[t=="None"] <- NA
      #removes rows that have NA
      t <- t[complete.cases(t), ]
      #removes any row with an entry with multiple values
      #t <- t[grep(pattern = ",", x = t[ , 2], invert = TRUE ) , ]
      #t <- t[grep(pattern = ",", x = t[ , 3], invert = TRUE ) , ]
      t <- data.frame(lapply(t, function(x) gsub(pattern = ", ", replacement = "/", x = x)))
      #Count combinations of predictions
      t$comb <- paste(t[ , 2], t[ , 3])
      combs <- expand.grid(levels(factor(t[ , 2])), levels(factor(t[ , 3])),stringsAsFactors = FALSE)
      combs$Count <- table(t$comb)[paste(combs[ , 1], combs[ , 2])]
      combs$Count<- as.numeric(as.character(combs$Count))
      # draw the heatmap with the specified test results
      title <- paste(paste(colnames(t[ , 2:3])[[1]], colnames(t[ , 2:3])[[2]], sep = " vs. "), "\n", "")
      plots$heatplot <- ggplot(combs, aes(x=Var1, y=Var2, fill = Count)) + 
        geom_tile() + 
        ggtitle(title) + 
        xlab(colnames(t[ , 2:3])[[1]]) + 
        ylab(colnames(t[ , 2:3])[[2]]) +
        theme(axis.title = element_text(face="bold", size=20)) +
        theme(plot.title = element_text(face="bold", size=20)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text =  element_text(size=15))
      plots$heatplot
      }
    else if (sum(input$heattests %in% cattests) == 1){
      numtest <- input$heattests[!input$heattests %in% cattests][1]
      cattest <- input$heattests[input$heattests %in% cattests][1]
      #replaces "-" with NAs
      t[t=="None"] <- NA
      t <- data.frame(lapply(t, function(x) gsub(pattern = "\\[|\\]|\\'", replacement = "", x = x)))
      t <- data.frame(lapply(t, function(x) gsub(pattern = ", ", replacement = "/", x = x)))
      t <- t[ ,c(input$heattests)]
      t <- t[complete.cases(t), ]
      t[[numtest]] <- as.numeric(as.character(t[[numtest]]))
      plots$heatplot <- ggplot(t, aes(x = t[[cattest]], y = t[[numtest]])) + geom_violin() +
        xlab(cattest) + ylab(numtest)
      plots$heatplot
    }
      else if (sum(input$heattests %in% cattests) == 0){
        numtest1 <- input$heattests[1]
        numtest2 <- input$heattests[1]
        #replaces "-" with NAs
        t[t=="None"] <- NA
        t <- data.frame(lapply(t, function(x) gsub(pattern = "\\[|\\]|\\'", replacement = "", x = x)))
        t <- t[ ,c(input$heattests)]
        t <- t[complete.cases(t), ]
        t[[input$heattests[1]]] <- as.numeric(as.character(t[[input$heattests[1]]]))
        t[[input$heattests[2]]] <- as.numeric(as.character(t[[input$heattests[2]]]))
        if ("aggregatescore" %in% input$heattests){
          plots$heatplot <- ggplot(t, aes(x = t[[input$heattests[1]]], y = t[[input$heattests[2]]])) + geom_smooth(method = "lm") + geom_jitter()+
            xlab(input$heattests[1]) + ylab(input$heattests[2])
        }
        else {
        plots$heatplot <- ggplot(t, aes(x = t[[input$heattests[1]]], y = t[[input$heattests[2]]])) + geom_smooth(method = "lm") + geom_point()+
          xlab(input$heattests[1]) + ylab(input$heattests[2])
        }
        plots$heatplot
        }
    }
  })
  
  output$downloadheat <- downloadHandler(
    filename = function() { paste(gsub(pattern = " vs. ", replacement = "_vs_", gsub(pattern = " \n ", replacement = "", plots$heatplot$labels$title)), '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = plots$heatplot, device = "png")
    }
  )
  
#-----------------------------------------------------------------------------------
  #test classfication distribution

  output$zoomdist <- renderPlot({
    table <- read_table()
    table <- table[(table$aggregatescore <= input$barthreshold[2]), ]
    table <- table[(table$aggregatescore >= input$barthreshold[1]), ]
  
    t <- table
    
    if (!is.null(input$bargene)){
      t <- table[table$marker %in% input$bargene, ]
    }
    if(length(t[[1]])==0){
      return()
    }
    #replaces "None" with NAs
    t[t=="None"] <- NA
    t <- data.frame(lapply(t, function(x) gsub(pattern = "\\[|\\]|\\'", replacement = "", x = x)))
    
    test <- input$bartests
    t <- t[!is.na(t[test]),c("marker", test)]
    plots$barplot <- ggplot(t, aes(x = marker, fill = t[[test]])) +
      geom_bar() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(fill=guide_legend(title=input$bartests))
    if (is.null(input$bargene) || length(input$bargene) > 1){
      plots$barplot <- ggplot(t, aes(x = marker, fill = t[[test]])) +
        geom_bar() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        guides(fill=guide_legend(title=input$bartests))
      plots$barplot
    }
    else{
      plots$barplot <- ggplot(t, aes(x = marker, fill = t[[test]])) +
        geom_bar() +
        coord_polar("y", start=0) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        guides(fill=guide_legend(title=input$bartests))
      plots$barplot
    }
    })
  
  output$downloadbar <- downloadHandler(
    filename = function() { paste(paste(plots$barplot$guides$fill$title,"_barchart_dist", sep = ""), '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = plots$barplot, device = "png")
    }
  )

#-----------------------------------------------------------------------------------
  #table view
  
  output$table <- renderDataTable(
    datatable({
      table <- read_table()
      table <- table[((table$aggregatescore <= input$tablethreshold[2]) || table$aggregatescore == "None"), ]
      table <- table[((table$aggregatescore >= input$tablethreshold[1] || table$aggregatescore == "None")), ]
      table <- table[((table$revel <= input$revelthreshold[2] || table$aggregatescore == "None")), ]
      table <- table[((table$revel >= input$revelthreshold[1] || table$aggregatescore == "None")), ]
      if(is.null(input$table)){
        plots$table <- table
        plots$table
      }
      else if("None" %in% input$table) {
        plots$table <- NULL
        plots$table
      }
      else if("All" %in% input$table){
        plots$table <- table
        plots$table
      }
      else{
        plots$table <- table[, input$table]
        plots$table
      }
    }),
    options = list(scrollX = TRUE)
  )
  
  output$downloadtable <- downloadHandler(
    filename = function() { paste("variant_table", '.csv', sep='') },
    content = function(file) {
      write.csv(plots$table, file)
    }
  )



#-----------------------------------------------------------------------------------
#Prediction score heatmap
  predtable <- function(t){
    mt <- c("['A', 'D']" = 1, "['A']" = 1, "['D', 'N']" = 0.5, "['D']" = 1, "['N']" = 0, "['P']" = 0)
    sift <- c("Damaging" = 1, "Tolerated" = 0)
    ma <- c("['H']" = 1, "['L']" = .33, "['M']" = .66, "['N']" = 0)
    fathmm <- c("['D']" = 1,"['N']" = 0)
    metasvm <- c("['D']" = 1,"['T']" = 0)
    metalr <- c("['D']" = 1,"['T']" = 0)
    provean <- c("['D']" = 1,"['N']" = 0, "['D', 'N']" = 0.5)
    lrt <- c("['D']" = 1,"['N']" = 0, "['U']" = 0.5)
    
    t$MutationTaster.Prediction <- as.vector(sapply(as.character(t$MutationTaster.Prediction), function(x) mt[x]))
    t$Sift.Prediction <- as.vector(sapply(as.character(t$Sift.Prediction), function(x) sift[x]))
    t$MutationAssessor.Prediction <- as.vector(sapply(as.character(t$MutationAssessor.Prediction), function(x) ma[x]))
    t$Fathmm.Prediction <- as.vector(sapply(as.character(t$Fathmm.Prediction), function(x) fathmm[x]))
    t$Metasvm.Prediction <- as.vector(sapply(as.character(t$Metasvm.Prediction), function(x) metasvm[x]))
    t$Metalr.Prediction <- as.vector(sapply(as.character(t$Metalr.Prediction), function(x) metalr[x]))
    t$Provean.Prediction <- as.vector(sapply(as.character(t$Provean.Prediction), function(x) provean[x]))
    t$lrt.Prediction <- as.vector(sapply(as.character(t$lrt.Prediction), function(x) lrt[x]))
    
    t
  }
  
  predheatval <- reactiveValues()
  
  observeEvent(input$updatepredheat,{
    predheatval$predheatgene <- input$predheatgene
  })
  
  output$predheat <- renderPlot({
    table <- predtable(read_table())
    if (!is.null(predheatval$predheatgene)){
      if(!"All" %in% predheatval$predheatgene){
      t <- table[table$marker %in% predheatval$predheatgene, ]
      }
      else{
        t <- table
      }
      t <- t[(t$aggregatescore <= input$predheatthreshold[2]), ]
      t <- t[(t$aggregatescore >= input$predheatthreshold[1]), ]
      t <- t[order(as.numeric(as.character(t$Position))),]
      if (!"All" %in% input$interp){
        t <- t[t$Clinvar2.accessions %in% input$interp,]
      }
      
      if(length(t[[1]])==0){
        return()
      } 
      t <- t[,c("marker","result","interpretation","Clinvar2.accessions", tests)]
      t <- melt(t)
      # vus <- t[t$interpretation == "VUS",]
      # ben <- t[t$interpretation == "Benign",]
      # path <- t[t$interpretation == "Path",]
      # if (length(vus[[1]])!=0){
      # plots$predheatvus <- ggplot(vus, aes(x=result, y=variable, fill = value)) + 
      #   geom_tile() + 
      #   scale_fill_gradientn(colours = c("green", "black", "red"), breaks = c(0,1), labels = c("Neutral", "Damaging"), limits =c(0,1)) +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #   guides(fill=guide_legend(title="")) +
      #   xlab("Variants") +
      #   ylab("Tests") +
      #   ggtitle("VUS") +
      #   facet_grid(.~Clinvar2.accessions, scale="free", space="free")
      # }
      # else{
      #   plots$predheatvus <- NULL
      # }
      # if (length(ben[[1]])!=0){
      # plots$predheatbenign <- ggplot(ben, aes(x=result, y=variable, fill = value)) + 
      #   geom_tile() + 
      #   scale_fill_gradientn(colours = c("green", "black", "red"), breaks = c(0,1), labels = c("Neutral", "Damaging"), limits =c(0,1)) +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #   guides(fill=guide_legend(title="")) +
      #   xlab("Variants") +
      #   ylab("Tests") +
      #   ggtitle("Benign") +
      #   facet_grid(.~Clinvar2.accessions, scale="free", space="free")
      # }
      # else{
      #   plots$predheatben <- NULL
      # }
      # if (length(path[[1]])!=0){
      # plots$predheatpath <- ggplot(path, aes(x=result, y=variable, fill = value)) + 
      #   geom_tile() + 
      #   scale_fill_gradientn(colours = c("green", "black", "red"), breaks = c(0,1), labels = c("Neutral", "Damaging"), limits =c(0,1)) +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #   guides(fill=guide_legend(title="")) +
      #   xlab("Variants") +
      #   ylab("Tests") +
      #   ggtitle("Pathogenic") +
      #   facet_grid(.~Clinvar2.accessions, scale="free", space="free")
      # }
      # else{
      #   plots$predheatpath <- NULL
      # }
      
      lvl <- levels(factor(t$Clinvar2.accessions))
      ptlist <- lapply(lvl, function(x){
        cur <- t[t$Clinvar2.accession == x, ]
        curplot <- ggplot(cur, aes(x=result, y=variable, fill = value)) + 
          geom_tile() + 
          scale_fill_gradientn(colours = c("green", "black", "red"), breaks = c(0,1), labels = c("Neutral", "Damaging"), limits =c(0,1)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          theme(axis.text.y = element_blank()) +
          guides(fill=guide_legend(title="")) +
          xlab("Variants") +
          ylab("Tests") +
          ggtitle(x)
        if (input$genefacet){
          curplot <- curplot + facet_grid(marker~interpretation, scale="free", space="free")
        }
        else{
          curplot <- curplot + facet_grid(.~interpretation, scale="free", space="free")
        }
        curplot
      })
      
      #to_delete <- !sapply(ptlist,is.null)
      #ptlist <- ptlist[to_delete] 
      #if (length(ptlist)==0) return(NULL)
      plots$predheat <- grid.arrange(grobs=ptlist)
      plots$predheat
      
      
  
      
      
    }
  }, height = 1000, width = 1000)
  output$predheathelp <- renderText({"Select a gene and click Plot!"})
  output$predheatcontainer <- renderUI({
    if (!is.null(predheatval$predheatgene)){
      plotOutput("predheat")
    }
    else{
      h4(textOutput("predheathelp"))
    }
  })
  output$downloadpredheat <- downloadHandler(
    filename = function() { paste("predictionheatmap", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = plots$predheat, device = "png")
    }
  )
}


# Run the application 
shinyApp(ui = ui, server = server)

