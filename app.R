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
library(shinyFiles) # save to button

source("normFunctions.R")
source("functions.R")

# Sets maximum upload size to 100MB
options(shiny.maxRequestSize = 100*1024^2)
options(stringsAsFactors=FALSE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


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
          shinySaveButton("savePep2Pro", "Save file", "Save file as ...", filetype=list(txt="txt")),
          shinyjs::hidden(p(id = "textFilteringPeptides", "Processing...")),
          
          # useShinyjs(),
          # textInput("saveProteinFilename", "Save filtered Proteins to", value="proteinGroups_filtered.txt"),
          # actionButton("saveButton", "Save filtered Proteins"),
          # shinyjs::hidden(p(id = "textFilteringPeptides", "Processing...")),
          
          fileInput("protFile", "Proteins")
        )
      ),
      
      fluidRow(rHandsontableOutput("metaData"))
    ),
    
    # Filters
    tabItem(
      tabName = "filters",
      fluidRow(
        box(
          radioButtons('protMDSColor', 'Color MDS by', choices=c("Group", "Quantity missing"), inline=TRUE)
        ),
        box(
          radioButtons('peptMDSColor', 'Color MDS by', choices=c("Group", "Quantity missing"), inline=TRUE)
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
      ),
      
      checkboxGroupInput(inputId = "sampleCheckbox", label = "Samples to be removed")
      # actionButton("updateSampleFilterButton", "Updates Samples"),
      # textOutput("selected_var")
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
        ),
        selectInput("normMethodCorrelationHeatmap", "Normalization Method:", 
                    c("loggedInt", "medianNorm", "meanNorm", "vsnNorm",
                      "quantNorm", "cycLoessNorm", "rlrNorm", "giNorm"))
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
      ),
      shinySaveButton("saveNormProtein", "Save file", "Save file as ...", filetype=list(txt="txt"))
      
      ),
      
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
  peptideAnnotationColums = c("id", "Protein.group.IDs")
  proteinAnnotationColums = c("id")
  # shinyjs::hide("normMethodCorrelationHeatmap")
  
  peptides <- reactive({
    peptFile <- input$peptFile
    if(is.null(peptFile)) return(NULL)
    
    tmpRawData = data.frame(fread(peptFile$datapath))
    tmpData = extractPeptideData(rawData = tmpRawData, peptideAnnotationColums = peptideAnnotationColums)
    cat(ifelse(tmpData[["isTMT"]], "Peptides: TMT detected\n", "Peptides: Label-Free detected\n"))
    tmpData[["data"]]
  })
  
  
  observe({
    volumes <- c("UserFolder" = getwd())
    shinyFileSave(input, "savePep2Pro", roots=volumes, session=session)
    fileinfo <- parseSavePath(volumes, input$savePep2Pro)
    if (nrow(fileinfo) > 0) {
      shinyjs::disable("savePep2Pro")
      shinyjs::show("textFilteringPeptides")
      peptides = peptides()
      cat("Filtering peptides (can take a few minutes) ...")
      # save(peptides, file = "peptides.Rdata")
      filteredProteins = filterPeptides(peptides)
      write.table(filteredProteins, file = as.character(fileinfo$datapath), sep = "\t", row.names = F) #, col.names=NA
      cat(" done\n")
      shinyjs::enable("savePep2Pro")
      shinyjs::hide("textFilteringPeptides")
    }
  })
  
  proteins <- reactive({
    protFile <- input$protFile
    if(is.null(protFile)) return(NULL)
    
    tmpRawData = data.frame(fread(protFile$datapath))
    tmpData = extractProteinData(rawData = tmpRawData, proteinAnnotationColums = proteinAnnotationColums)
    cat(ifelse(tmpData[["isTMT"]], "Proteins: TMT detected\n", "Proteins: Label-Free detected\n"))
    tmpData[["data"]]
  })
  
  
  ### Creates a Dataframe that will be used to edit the metaData ####
  ### Creates handsontable where metaData1() will be edited ####
  metaData1 <- reactive({
    req(peptides())
    req(proteins())
    tempTable = data.frame(Peptide.Sample.Names = colnames(peptides())[-seq_along(peptideAnnotationColums)],
                           Protein.Sample.Names = colnames(proteins())[-seq_along(proteinAnnotationColums)],
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
    hot_to_r(input$metaData)
  }) 
  
  observe({
    req(peptides())
    req(proteins())
    req(metaData1())
    req(input$metaData)
    meta = metaData2()
    
    updateCheckboxGroupInput(session, inputId = "sampleCheckbox",
                             choices = #c("1","2","3")
                               if(meta$Custom.Sample.Name[1] == ""){
                                 meta$Peptide.Sample.Names
                               } else {
                                 meta$Custom.Sample.Names
                               }
    )
  })
  
  peptidesSampleFiltered <- reactive({
    req(peptides())
    req(metaData1())
    req(input$metaData)
    peptides = peptides()
    meta = metaData2()
    if(meta$Custom.Sample.Name[1] == ""){
      filteredPeptides = peptides[,c(rep(TRUE, length(peptideAnnotationColums)), !meta$Peptide.Sample.Names %in% input$sampleCheckbox)] # TRUE to keep annotation
    } else {
      filteredPeptides = peptides[,c(rep(TRUE, length(peptideAnnotationColums)), !meta$Custom.Sample.Names %in% input$sampleCheckbox)]
    }
    filteredPeptides
  })
  
  proteinsSampleFiltered <- reactive({
    req(proteins())
    req(metaData1())
    req(input$metaData)
    proteins = proteins()
    meta = metaData2()
    if(meta$Custom.Sample.Name[1] == ""){
      filteredProteins = proteins[,c(rep(TRUE, length(proteinAnnotationColums)), !meta$Peptide.Sample.Names %in% input$sampleCheckbox)] # TRUE to keep annotation
    } else {
      filteredProteins = proteins[,c(rep(TRUE, length(proteinAnnotationColums)), !meta$Custom.Sample.Names %in% input$sampleCheckbox)]
    }
    filteredProteins
  })
  
  metaDataFiltered <- reactive({
    req(proteins())
    req(metaData1())
    req(input$metaData)
    meta = metaData2()
    if(meta$Custom.Sample.Name[1] == ""){
      filteredMetadata = meta[!meta$Peptide.Sample.Names %in% input$sampleCheckbox,]
    } else {
      filteredMetadata = meta[!meta$Custom.Sample.Names %in% input$sampleCheckbox,]
    }
    filteredMetadata
  })
  
  
  
  
  # Boxplot Proteins
  output$filtProtBoxplot <- renderPlot({
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    boxplot(proteins[, sampleCols], outline=FALSE, col=colorGroup(meta$Group)[meta$Group], main="Boxplot of Proteins",
            xlab="Sample", ylab="Corrected Intensity", names=sampleLabels, las = 2)
    legend("top", inset=c(0, -.085), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
           col=colorGroup(meta$Group), pch=15, horiz=TRUE)
  })
  
  # Boxplot Peptides
  output$filtPeptBoxplot <- renderPlot({
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    
    sampleCols = meta$Peptide.Sample.Names
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    boxplot(peptides[, sampleCols], outline=FALSE, col=colorGroup(meta$Group)[meta$Group], main="Boxplot of Peptides",
            xlab="Sample", ylab="Corrected Intensity", names=sampleLabels, las = 2)
    legend("top", inset=c(0, -.085), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
           col=colorGroup(meta$Group), pch=15, horiz=TRUE)
  })
  
  # Histogram Proteins
  output$filtProtHist <- renderPlot({
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
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
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
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
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
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
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
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
  
  
  
  
  
  
  
  
  
  normProteins <- reactive({
    sampleCols = metaDataFiltered()$Peptide.Sample.Names
    proteins <- proteinsSampleFiltered()[,sampleCols]
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
  
  
  
  
  
  
  
  
  output$totalIntensity_barplot <- renderPlot({
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    meta <- metaDataFiltered()
    
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    
    par(mfrow = c(3,3), mar=c(4,8,3,1))
    for(i in names(normList)){
      barplot(colSums(normList[[i]], na.rm = T), main = i, las = 2, yaxt="n", cex.main = 1.5, col = plasma(ncol(normList[[i]])), names.arg = sampleLabels)
      axis(side = 2, cex.axis=1.5, las = 2)
      # axis(side = 1, at = seq_along(colnames(normList[[i]])), labels = colnames(normList[[i]]), cex.axis=1.5, las = 2)
      if(i == "vsnNorm") mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
      abline(h = max(colSums(normList[[i]], na.rm = T)), lty = 2)
    }
  })
  
  output$PCV_boxplot <- renderPlot({
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    meta <- metaDataFiltered()
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
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    meta <- metaDataFiltered()
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
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    meta <- metaDataFiltered()
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
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    meta <- metaDataFiltered()
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
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    
    heatmapMissing(normList[["loggedInt"]], groups, batch, sampleLabels)
  })
  
  output$cor_heatmap <- renderPlot({
    # shinyjs::show("normMethodCorrelationHeatmap")

    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    batch <- meta$Batch
    
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    
    heatmapCorr(normList[[input$normMethodCorrelationHeatmap]], groups, batch, sampleLabels)
  })
  
  output$logFC_density <- renderPlot({
    # shinyjs::hide("normMethodCorrelationHeatmap")
    
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    groups <- meta$Group
    densityLogFC(normList, groups)
  })
  
  observe({
    volumes <- c("UserFolder" = getwd())
    shinyFileSave(input, "saveNormProtein", roots=volumes, session=session)
    fileinfo <- parseSavePath(volumes, input$saveNormProtein)
    if (nrow(fileinfo) > 0) {
      normList <- normProteins()
      meta <- metaDataFiltered()
      if(is.null(normList) | is.null(meta)) return(NULL)
      groups <- meta$Group
      batch <- meta$Batch
      normData <- normList[[input$normMethod]]
      if(input$imputationMethod == "No Imputation"){
        write.table(normData, file = as.character(fileinfo$datapath), sep = "\t", row.names = F) 
      } else {
        write.table(impute(normData, input$imputationMethod), file = as.character(fileinfo$datapath), sep = "\t", row.names = F)
      }
    }
  })
  
  shiny::observeEvent(input$goButtonDAtest, {
    normList <- normProteins()
    meta <- metaDataFiltered()
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
    meta <- metaDataFiltered()
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
