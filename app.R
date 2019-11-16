library(shinydashboard)
library(data.table)
library(DT)
library(shiny)
library(viridis)
library(ComplexHeatmap) # normalization figures
library(DAtest) # DAtest
library(impute) # knn imputation
library(pcaMethods) # bpca imputation
library(imputeLCMD) # QRILC, MinDet imputation
library(sva) # combat
library(shinyjs) # hiding button
library(rhandsontable) # edibale table


source("normFunctions.R")
source("functions.R")

# Sets maximum upload size to 100MB
options(shiny.maxRequestSize = 100*1024^2)
options(stringsAsFactors=FALSE)

tweaks <- 
  list(tags$head(tags$style(HTML("
                                 .multicol { 
                                 height: 150px;
                                 -webkit-column-count: 3; /* Chrome, Safari, Opera */ 
                                 -moz-column-count: 3;    /* Firefox */ 
                                 column-count: 3; 
                                 -moz-column-fill: auto;
                                 -column-fill: auto;
                                 } 
                                 ")) 
  ))


DAtestTests = eval(formals(testDA)$tests)
allChecks <- seq_along(DAtestTests)
names(allChecks) <- DAtestTests
controls <-
  list(h3("Test to be EXCLUDED:"), 
       tags$div(align = 'left', 
                class = 'multicol', 
                checkboxGroupInput(inputId  = "checkboxDAtestTests", 
                                   label    = NULL, 
                                   choices  = allChecks,
                                   selected = which(DAtestTests %in% "per"),
                                   inline   = FALSE))) 

header <- dashboardHeader(title = "proteiNorm")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Data", tabName = "data", icon = icon("table")),
    menuItem("Filters", tabName = "filters", icon = icon("filter")),
    menuItem("Normalization", tabName = "norm", icon = icon("poll")),
    menuItem("DAtest", tabName = "DAtest", icon = icon("file-signature")),
    menuItem("DAtest Power", tabName = "DAtestPower", icon = icon("chart-line"))
  )
)

body <- dashboardBody(
  tabItems(
    # Loading data
    tabItem(
      tabName = "data",
      fluidRow(
        # Files
        box(
          title="Inputs", width=12,
          fileInput("peptFile", "Peptides"),
          
          useShinyjs(),
          textInput("saveProteinFilename", "Save filtered Proteins to", value="proteinGroups_filtered.txt"),
          actionButton("saveButton", "Save filtered Proteins"),
          shinyjs::hidden(p(id = "textFilteringPeptides", "Processing...")),
          
          fileInput("protFile", "Proteins")
          # fileInput("metaFile", "Metadata")
        )
      ),
      
      fluidRow(
        rHandsontableOutput("metaData")#,
        
        # tabBox(
        #   title="Raw Data (upload to view)", width=12, id="rawdata",
        #   tabPanel(
        #     "Peptides",
        #     div(DTOutput("rawPeptTable"),
        #         style = "overflow-y: scroll;overflow-x: scroll;")
        #   ),
        #   tabPanel(
        #     "Proteins",
        #     div(DTOutput("rawProtTable"),
        #         style = "overflow-y: scroll;overflow-x: scroll;")
        #   )#,
        #   # tabPanel(
        #   #   "Metadata",
        #   #   div(DTOutput("rawMetaTable"),
        #   #       style = "overflow-y: scroll;overflow-x: scroll;")
        #   # )
        # )
      ),
      textOutput("selected_var")
    ),
    
    # Filters
    tabItem(
      tabName = "filters",
      fluidRow(
        box(
      #     title="Proteins",
      #     uiOutput('protId'),
      #     checkboxGroupInput('protFilt', "Filter proteins by",
      #                        choices=c("Count", "Score"), inline=TRUE),
          radioButtons('protMDSColor', 'Color MDS by',
                       choices=c("Group", "Quantity missing"), inline=TRUE)
      #     uiOutput('protCountSlide'),
      #     uiOutput('protScoreSlide')
        ),
        box(
      #     title="Peptides",
      #     uiOutput('peptId'),
      #     checkboxGroupInput('peptFilt', "Filter peptides by",
      #                        choices=c("Count", "Score"), inline=TRUE),
          radioButtons('peptMDSColor', 'Color MDS by',
                       choices=c("Group", "Quantity missing"), inline=TRUE)
      #     uiOutput('peptCountSlide'),
      #     uiOutput('peptScoreSlide')
        )
      ),
      fluidRow(
        tabBox(
          title="Protein Plots", width=6,
          tabPanel(
            "Boxplot",
            plotOutput("filtProtBoxplot")
          ),
          tabPanel(
            "Histogram",
            plotOutput("filtProtHist")
          ),
          tabPanel(
            "MDS",
            plotOutput("filtProtMDS")
          )
        ),
        tabBox(
          title="Peptide Plots", width=6,
          tabPanel(
            "Boxplot",
            plotOutput("filtPeptBoxplot")
          ),
          tabPanel(
            "Histogram",
            plotOutput("filtPeptHist")
          ),
          tabPanel(
            "MDS",
            plotOutput("filtPeptMDS")
          )
        )
      )#,
      # fluidRow(
      #   tabBox(
      #     title="Cleaned Data", width=12, id="cleandata",
      #     tabPanel(
      #       "Proteins",
      #       div(DTOutput("filtProtTable"),
      #           style = "overflow-y: scroll;overflow-x: scroll;")
      #     ),
      #     tabPanel(
      #       "Peptides",
      #       div(DTOutput("filtPeptTable"),
      #           style = "overflow-y: scroll;overflow-x: scroll;")
      #     )
      #   )
      # )
    ),
    
    
    
    # Normalization
    tabItem(
      tabName = "norm",
      
      h2("Proteins"), # cleanProteins()
      fluidRow(
        tabBox(
          width=12,
          tabPanel(
            "Total Intensity",
            plotOutput("totalIntensity_barplot")
          ),
          tabPanel(
            "PCV",
            plotOutput("PCV_boxplot")
          ),
          tabPanel(
            "PMAD",
            plotOutput("PMAD_boxplot")
          ),
          tabPanel(
            "PEV",
            plotOutput("PEV_boxplot")
          ),
          tabPanel(
            "Intragroup Correlation",
            plotOutput("cor_boxplolt")
          ),
          tabPanel(
            "Missing Values",
            plotOutput("NA_heatmap")
          ),
          tabPanel(
            "Correlation heatmap",
            plotOutput("cor_heatmap")
          ),
          tabPanel(
            "LogFC density",
            plotOutput("logFC_density")
          )
        )
      )
    ),
    
    # DAtest
    tabItem(
      tabName = "DAtest",
      
      fluidRow(inputPanel(
        selectInput("normMethod", "Normalization Method:", 
                    c("loggedInt", "medianNorm", "meanNorm", "vsnNorm",
                      "quantNorm", "cycLoessNorm", "rlrNorm", "giNorm")),
        selectInput("imputationMethod", "Imputation Method:", 
                    c("No Imputation", "KNN", "QRILC", "MinDet", "MinProb", "min", "zero"))
      )),
      
      tweaks,
      fluidRow(column(width = 4, controls)),
      
      shiny::actionButton(inputId = "goButtonDAtest", label = "Go!", width = '100%', style='font-size:150%'),
      
      div(DTOutput("DAtestResults"), style = "overflow-y: scroll;overflow-x: scroll;"),
      
      plotOutput(outputId = "DAtestFigure")
    ),
    
    # DAtest Power
    tabItem(
      tabName = "DAtestPower",
      fluidRow(inputPanel(
        selectInput(inputId = "DAtest4power", label = "Estimating statistical power for:", choices = eval(formals(testDA)$tests))
      )),
      
      shiny::actionButton(inputId = "goButtonDAPower", label = "Go!", width = '100%', style='font-size:150%'),
      
      div(DTOutput("DAtestPowerResults"), style = "overflow-y: scroll;overflow-x: scroll;"),
      
      plotOutput(outputId = "DAtestPowerFigure")
      
    )
  )
)



ui <- dashboardPage(header, sidebar, body)

server <- function(input, output, session) {
  peptides <- reactive({
    peptFile <- input$peptFile
    if(is.null(peptFile)) return(NULL)
    # fread(peptFile$datapath)
    
    tmpRawData = data.frame(fread(peptFile$datapath))
    tmpData = extractPeptideData(rawData = tmpRawData)
    ifelse(tmpData[["isTMT"]], print("Peptides: TMT detected"), print("Peptides: Label-Free detected"))
    tmpData[["data"]]
  })
  
  observeEvent(input$saveButton, {
    validate(
      need(input$saveProteinFilename != "", message="Please enter a valid filename")
    )
    
    shinyjs::disable("saveButton")
    shinyjs::show("textFilteringPeptides")
    
    peptides = peptides()
    cat("Filtering peptides (can take a few minutes) ...")
    filteredProteins = filterPeptides(peptides)
    write.table(filteredProteins, file = input$saveProteinFilename, sep = "\t", row.names = F) #, col.names=NA
    cat("done\n")
    shinyjs::enable("saveButton")
    shinyjs::hide("textFilteringPeptides")
  })
  
  proteins <- reactive({
    protFile <- input$protFile
    if(is.null(protFile)) return(NULL)
    # fread(protFile$datapath)
    
    tmpRawData = data.frame(fread(protFile$datapath))
    tmpData = extractProteinData(rawData = tmpRawData)
    ifelse(tmpData[["isTMT"]], print("Proteins: TMT detected"), print("Proteins: Label-Free detected"))
    tmpData[["data"]]
  })
  
  # meta <- reactive({
  #   metaFile <- input$metaFile
  #   if(is.null(metaFile)) return(NULL)
  #   fread(metaFile$datapath)
  # })
  
  ### Creates a Dataframe that will be used to edit the metaData ####
  ### Creates handsontable where metaData1() will be edited ####
  metaData1 <- reactive({
    req(peptides())
    req(proteins())
    tempTable = data.frame(Peptide.Sample.Names = colnames(peptides())[-1],
                           Protein.Sample.Names = colnames(proteins())[-1],
                           Custom.Sample.Names = "",
                           Group = "",
                           Batch = "")
    
    rhandsontable(tempTable) %>% 
      hot_col("Peptide.Sample.Names", readOnly = T) %>% 
      hot_col("Protein.Sample.Names", readOnly = T)
  })
  #### Outputs the Rhandsontable ####
  output$metaData<- renderRHandsontable({metaData1()})
  
  #### Changes the handsontable back into a dataframe ####
  metaData2 <- reactive({
    # req(peptides())
    # req(proteins())
    hot_to_r(input$metaData) 
  })
  
  
  # output$rawProtTable <- renderDT({ proteins() })
  # output$rawPeptTable <- renderDT({ peptides() })
  # output$rawMetaTable <- renderDT({ meta() }, editable='cell')
  
  # output$protId <- renderUI({
  #   selectInput('protIdCol', 'Numerical ID column in proteins',
  #               c("Select none to make one"="", "None", names(proteins())),
  #               selected="None")
  # })
  
  # output$peptId <- renderUI({
  #   selectInput('peptIdCol', 'Numerical ID column in peptides',
  #               c("Select none to make one"="", "None", names(peptides())),
  #               selected="None")
  # })
  
  # scoreFilter <- function(dat, scoreCutoff = 5) {
  #   dat[dat$score >= scoreCutoff, ]
  # }
  
  # countFilter <- function(dat, meta, countCutoff = 8, numReplicates = 2) {
  #   dat <- as.data.frame(dat)
  #   for (group in levels(as.factor(meta$group))) {
  #     tmp <- meta[meta$group==group, ]
  #     countCols <- paste0(tmp$sample, ".counts")
  #     sampleCols <- tmp$sample
  #     dat[rowSums(dat[countCols]>countCutoff) < numReplicates, sampleCols] = NA
  #   }
  #   as.data.table(dat)
  # }
  
  # cleanProteins <- reactive({
  #   proteins <- proteins()
  #   meta <- meta()
  #   
  #   if(is.null(proteins) | is.null(meta)) return(NULL)
  #   
  #   idCol <- input$protIdCol
  #   if(idCol=="None") {
  #     proteins <- cbind(proteins, idCol=seq(1, nrow(proteins)))
  #     idCol <- "idCol"
  #   }
  #   
  #   colNames <- c(paste0("Reporter intensity corrected ", meta$ionNum, " TMT", meta$batch),
  #                 paste0("Reporter intensity count ", meta$ionNum, " TMT", meta$batch)) 
  #   tmpNames <- c(idCol, colNames, "Score")
  #   cleanProt <- proteins[, ..tmpNames]
  #   colnames(cleanProt) <- c("id", meta$sample, paste0(meta$sample, ".counts"), "score")
  #   
  #   # Filters
  #   if("Score" %in% input$protFilt & !is.null(input$protScoreSlider)) {
  #     cleanProt <- scoreFilter(cleanProt, scoreCutoff = input$protScoreSlider)
  #   }
  #   
  #   if("Count" %in% input$protFilt & !is.null(input$protCountSlider)) {
  #     cleanProt <- countFilter(cleanProt, meta, countCutoff = input$protCountSlider)
  #   }
  #   
  #   cleanProt
  # })
  
  # # cleanPeptides <- reactive({
  #   peptides <- peptides()
  #   meta <- meta()
  #   
  #   if(is.null(peptides) | is.null(meta)) return(NULL)
  #   
  #   idCol <- input$peptIdCol
  #   if(idCol=="None") {
  #     peptides <- cbind(peptides, idCol=seq(1, nrow(peptides)))
  #     idCol <- "idCol"
  #   }
  #   
  #   colNames <- c(paste0("Reporter intensity corrected ", meta$ionNum, " TMT", meta$batch),
  #                 paste0("Reporter intensity count ", meta$ionNum, " TMT", meta$batch)) 
  #   tmpNames <- c(idCol, colNames, "Score")
  #   cleanPept <- peptides[, ..tmpNames]
  #   colnames(cleanPept) <- c("id", meta$sample, paste0(meta$sample, ".counts"), "score")
  #   
  #   # Filters
  #   if("Score" %in% input$peptFilt & !is.null(input$peptScoreSlider)) {
  #     cleanPept <- scoreFilter(cleanPept, scoreCutoff = input$peptScoreSlider)
  #   }
  #   
  #   if("Count" %in% input$peptFilt & !is.null(input$peptCountSlider)) {
  #     cleanPept <- countFilter(cleanPept, meta, countCutoff = input$peptCountSlider)
  #   }
  #   
  #   cleanPept
  # })
  
  # wideProteins <- reactive({
  #   proteins <- cleanProteins()
  #   meta <- meta()
  #   
  #   if(is.null(proteins) | is.null(meta)) return(NULL)
  #   
  #   colNames <- meta$sample
  #   wideProt <- proteins[, ..colNames]
  #   row.names(wideProt) <- proteins$id
  #   
  #   wideProt
  # })
  # 
  # widePeptides <- reactive({
  #   peptides <- cleanPeptides()
  #   meta <- meta()
  #   
  #   if(is.null(peptides) | is.null(meta)) return(NULL)
  #   
  #   colNames <- meta$sample
  #   widePept <- peptides[, ..colNames]
  #   row.names(widePept) <- peptides$id
  #   
  #   widePept
  # })
  
  normProteins <- reactive({
    sampleCols = metaData2()$Peptide.Sample.Names
    proteins <- proteins()[,sampleCols]
    if(is.null(proteins)) return(NULL)
    
    normList <- vector("list", 8)
    names(normList) <- c("loggedInt", "medianNorm", "meanNorm", "vsnNorm",
                         "quantNorm", "cycLoessNorm", "rlrNorm", "giNorm")
    normList[["loggedInt"]] <- logNorm(proteins)
    normList[["medianNorm"]] <- medianNorm(normList[["loggedInt"]])
    normList[["meanNorm"]] <- meanNorm(normList[["loggedInt"]])
    normList[["vsnNorm"]] <- vsnNorm(proteins)
    normList[["quantNorm"]] <- quantNorm(normList[["loggedInt"]])
    normList[["cycLoessNorm"]] <- cycLoessNorm(normList[["loggedInt"]])
    normList[["rlrNorm"]] <- rlrNorm(normList[["loggedInt"]])
    normList[["giNorm"]] <- giNorm(normList[["loggedInt"]])
    
    normList
  })
  
  logProteins <- reactive({
    normList <- normProteins()
    
    as.data.table(normList[["loggedInt"]])
  })
  
  output$filtProtTable <- renderDT({
    logProteins()
  })
  
  output$filtPeptTable <- renderDT({
    cleanPeptides()
  })
  
  # output$protCountSlide <- renderUI({
  #   if("Count" %in% input$protFilt) {
  #     sliderInput("protCountSlider", "Count Cutoff", 1, 50, value=8)
  #   }
  # })
  # 
  # output$protScoreSlide <- renderUI({
  #   if("Score" %in% input$protFilt) {
  #     sliderInput("protScoreSlider", "Score Cutoff", 1, 100, value=5)
  #   }
  # })
  # 
  # output$peptCountSlide <- renderUI({
  #   if("Count" %in% input$peptFilt) {
  #     sliderInput("peptCountSlider", "Count Cutoff", 1, 50, value=8)
  #   }
  # })
  # 
  # output$peptScoreSlide <- renderUI({
  #   if("Score" %in% input$peptFilt) {
  #     sliderInput("peptScoreSlider", "Score Cutoff", 1, 200, value=5)
  #   }
  # })
  
  # ---
  # findDensityCutoff
  # ---
  # Finds upper limit of the histogram so each box contains at least cutoffPercent of the data
  #
  findDensityCutoff <- function(longdat, cutoffPercent = .001) {
    densityCutoff <- 0
    newDensityCutoff <- max(longdat, na.rm=TRUE)
    totalNum <- length(longdat)
    while(newDensityCutoff != densityCutoff) {
      densityCutoff <- newDensityCutoff
      breakSeq <- seq(max(min(longdat, na.rm=TRUE), 0), densityCutoff,
                      by=(densityCutoff - max(min(longdat, na.rm=TRUE), 0)) / 30)
      freqs <- hist(longdat[longdat < densityCutoff], breaks = breakSeq, plot=FALSE)
      newDensityCutoff <- freqs$breaks[which((freqs$counts / totalNum) < cutoffPercent)[1]+1]
    }
    return(densityCutoff)
  }
  
  # Boxplot Proteins
  output$filtProtBoxplot <- renderPlot({
    proteins <- proteins()
    meta <- metaData2()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    boxplot(proteins[, sampleCols], outline=FALSE, col=colorGroup(meta$Group)[meta$Group], main="Boxplot of Proteins",
            xlab="Sample", ylab="Corrected Intensity", names=sampleLabels)
    legend("top", inset=c(0, -.085), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
           col=colorGroup(meta$Group), pch=15, horiz=TRUE)
  })
  
  # Boxplot Peptides
  output$filtPeptBoxplot <- renderPlot({
    peptides <- peptides()
    meta <- metaData2()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    boxplot(peptides[, sampleCols], outline=FALSE, col=colorGroup(meta$Group)[meta$Group], main="Boxplot of Peptides",
            xlab="Sample", ylab="Corrected Intensity", names=sampleLabels)
    legend("top", inset=c(0, -.085), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
           col=colorGroup(meta$Group), pch=15, horiz=TRUE)
  })
  
  # Histogram Proteins
  output$filtProtHist <- renderPlot({
    proteins <- proteins()
    meta <- metaData2()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    longdat <- unlist(proteins[, sampleCols])
    densityCutoff <- findDensityCutoff(longdat)
    longdat <- longdat[longdat < densityCutoff]
    breakSeq <- seq(max(min(longdat, na.rm=TRUE), 0), max(c(longdat, densityCutoff), na.rm = TRUE),
                    length.out=31)
    hist(longdat, main="Histogram of Proteins", breaks=breakSeq, xlab="Corrected Intensity", ylab="Counts")
  })
  
  # Histogram Peptides
  output$filtPeptHist <- renderPlot({
    peptides <- peptides()
    meta <- metaData2()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    longdat <- unlist(peptides[, sampleCols])
    densityCutoff <- findDensityCutoff(longdat)
    longdat <- longdat[longdat < densityCutoff]
    breakSeq <- seq(max(min(longdat, na.rm=TRUE), 0), max(c(longdat, densityCutoff), na.rm = TRUE),
                    length.out=31)
    hist(longdat, main="Histogram of Peptides", breaks=breakSeq, xlab="Corrected Intensity", ylab="Counts")
  })
  
  # MDS Proteins
  output$filtProtMDS <- renderPlot({
    proteins <- proteins()
    meta <- metaData2()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    d <- cmdscale(dist(t(proteins[, sampleCols])))
    if(input$protMDSColor=="Group") {
      plot(d, type="n", xlab="", ylab="", main="MDS of Proteins, Colored by Group")
      text(d, labels = sampleLabels, col=colorGroup(meta$Group)[meta$Group])
      legend("top", inset=c(0, -.085), levels(as.factor(meta$Group)),
             bty="n", xpd=TRUE,
             col=colorGroup(meta$Group), pch=15, horiz=TRUE)
    } else {
      plot(d, type="n", xlab="", ylab="", main="MDS of Peptides, Colored by Missing")
      missing <- colSums(is.na(proteins[, sampleCols]))
      colPalette <- colorRampPalette(c('blue','red'))
      missingCol <- colPalette(5)[as.numeric(cut(missing,breaks = 5))]
      text(d, labels = sampleLabels, col=missingCol)
    }
  })
  
  # MDS Peptides
  output$filtPeptMDS <- renderPlot({
    peptides <- peptides()
    meta <- metaData2()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    d <- cmdscale(dist(t(peptides[, sampleCols])))
    if(input$peptMDSColor=="Group") {
      plot(d, type="n", xlab="", ylab="", main="MDS of Peptides, Colored by Group")
      text(d, labels = sampleLabels, col=colorGroup(meta$Group)[meta$Group])
      legend("top", inset=c(0, -.085), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
             col=colorGroup(meta$Group), pch=15, horiz=TRUE)
    } else {
      plot(d, type="n", xlab="", ylab="", main="MDS of Peptides, Colored by Missing")
      missing <- colSums(is.na(peptides[, sampleCols]))
      colPalette <- colorRampPalette(c('blue','red'))
      missingCol <- colPalette(5)[as.numeric(cut(missing,breaks = 5))]
      text(d, labels = sampleLabels, col=missingCol)
    }
  })
  
  output$totalIntensity_barplot <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()

    if(is.null(normList) | is.null(meta)) return(NULL)
    
    # groups <- meta$group
    # batch <- meta$batch
    # peptides = cleanPeptides()
    # proteins = cleanProteins()
    # save(normList, groups, batch, meta, peptides, proteins, file = "../../proteiNorm_data.Rdata")
    # print("data saved")
    
    par(mfrow = c(3,3), mar=c(4,8,3,1))
    for(i in names(normList)){
      barplot(colSums(normList[[i]], na.rm = T), main = i, las = 2, yaxt="n", cex.main = 1.5, col = plasma(ncol(normList[[i]])))
      axis(side = 2, cex.axis=1.5, las = 2)
      # axis(side = 1, at = seq_along(colnames(normList[[i]])), labels = colnames(normList[[i]]), cex.axis=1.5, las = 2)
      if(i == "vsnNorm") mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
      abline(h = max(colSums(normList[[i]], na.rm = T)), lty = 2)
    }
  })
  
  output$PCV_boxplot <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    par(mar=c(10,6,3,1))
    plotData = lapply(normList, function(x) PCV(x, groups))
    boxplot(plotData, main = "PCV", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
    axis(2,cex.axis=1.5, las = 2)
    axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
    mtext(side = 2, text = "Pooled Coefficient of Variation", line = 4.5, cex = 1.5)
    points(rep(seq_along(normList), each = 2), unlist(plotData), pch = "*", cex = 1.3)
  })
  
  output$PMAD_boxplot <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    
    par(mar=c(10,6,3,1))
    plotData = lapply(normList, function(x) PMAD(x, groups))
    boxplot(plotData, main = "PMAD", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
    axis(2,cex.axis=1.5, las = 2)
    axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
    mtext(side = 2, text = "Median Absolute Deviation", line = 4.5, cex = 1.5)
    points(rep(seq_along(normList), each = 2), unlist(plotData), pch = "*", cex = 1.3)
  })
  
  output$PEV_boxplot <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    par(mar=c(10,6,3,1))
    plotData = lapply(normList, function(x) PEV(x, groups))
    boxplot(plotData, main = "PEV", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
    axis(2,cex.axis=1.5, las = 2)
    axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
    mtext(side = 2, text = "Pooled Estimate of Variance", line = 4.5, cex = 1.5)
    points(rep(seq_along(normList), each = 2), unlist(plotData), pch = "*", cex = 1.3)
  })
  
  output$cor_boxplolt <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    par(mar=c(10,6,3,1))
    plotData = lapply(normList, function(x) unlist(COR(x, groups)))
    boxplot(plotData, main = "Cor", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
    axis(side = 2,cex.axis=1.5, las = 2)
    axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
    mtext(side = 2, text = "Intragroup Correlation", line = 4.5, cex = 1.5)
    points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
  })
  
  output$NA_heatmap <- renderPlot({
    normList <- normProteins()
    
    # peptides <- peptides()
    # save(peptides, file = "../../proteiNorm_data_peptide.Rdata")
    # print("Pep saved")
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    heatmapMissing(normList[["loggedInt"]], groups, batch)
  })
  
  output$cor_heatmap <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    heatmapCorr(normList[["loggedInt"]], groups, batch)
  })
  
  output$logFC_density <- renderPlot({
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    densityLogFC(normList, groups)
  })
  
  shiny::observeEvent(input$goButtonDAtest, {
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    groups <- meta$Group
    batch <- meta$Batch
    normData <- normList[[input$normMethod]]
    
    cat("Excluding: ", DAtestTests[as.numeric(input$checkboxDAtestTests)], "\n")
    DAtestResults = if(input$imputationMethod == "No Imputation"){
      DAtest(normData, groups, batch, imputed = FALSE, exludedTests = DAtestTests[as.numeric(input$checkboxDAtestTests)]) # specify effect size
    } else {
      tempImpute = impute(normData, input$imputationMethod)
      DAtest(tempImpute, groups, batch, imputed = TRUE, exludedTests = DAtestTests[as.numeric(input$checkboxDAtestTests)])
    }
    
    output$DAtestResults <- renderDT({ summary(DAtestResults) })
    output$DAtestFigure <- renderPlot({ plot(DAtestResults) })
    
    updateSelectInput(session, "DAtest4power", choices = rownames(DAtestResults$run.times))
  })
  
  
  shiny::observeEvent(input$goButtonDAPower, {
    normList <- normProteins()
    meta <- metaData2()
    if(is.null(normList) | is.null(meta)) return(NULL)
    groups <- meta$Group
    batch <- meta$Batch
    normData <- normList[[input$normMethod]]
    if(input$imputationMethod != "No Imputation") normData = impute(normData, input$imputationMethod)
    
    DAtestPowerResults = powerDA(2^normData, predictor = as.character(groups), test = input$DAtest4power, cores = 10, relative = FALSE)
    output$DAtestPowerResults <- renderDT({ summary(DAtestPowerResults) })
    output$DAtestPowerFigure <- renderPlot({ plot(DAtestPowerResults) })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
