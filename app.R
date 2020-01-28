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
library(parallel) # detect cores
library(rJava) # finding screen resolution
library(ggfortify) # PCA autoplot
library(NormalyzerDE) # PCA autoplot



source("normFunctions.R")
source("functions.R")

# Sets maximum upload size to 100MB
options(shiny.maxRequestSize = 5000*1024^2)
options(stringsAsFactors=FALSE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# finding screen resolution
.jinit()
toolkit <- J("java.awt.Toolkit")
default_toolkit <- .jrcall(toolkit, "getDefaultToolkit")
screenDim <- .jrcall(default_toolkit, "getScreenSize")
screenHeight <- .jcall(screenDim, "D", "getHeight")
screenWidth <- .jcall(screenDim, "D", "getWidth")


tweaks <- 
  list(tags$head(tags$style(HTML("
                                 .multicol { 
                                 height: 300px;
                                 -webkit-column-count: 2; /* Chrome, Safari, Opera */ 
                                 -moz-column-count: 2;    /* Firefox */ 
                                 column-count: 2; 
                                 -moz-column-fill: auto;
                                 -column-fill: auto;
                                 } 
                                 ")) 
  ))


DAtestTests = eval(formals(testDA)$tests)
allChecks <- seq_along(DAtestTests)
names(allChecks) <- DAtestTestNames[DAtestTests]
controlsDAtest <-
  list(h3("Test to be excluded:"), 
       tags$div(align = 'left', 
                class = 'multicol', 
                checkboxGroupInput(inputId  = "checkboxDAtestTests", 
                                   label    = NULL, 
                                   choices  = allChecks,
                                   selected = which(DAtestTests %in% "per"),
                                   inline   = FALSE))) 

controlsRmSamples <-
  list(h3("Samples to be removed:"), 
       tags$div(align = 'left', 
                class = 'multicol', 
                checkboxGroupInput(inputId  = "sampleCheckbox", 
                                   label    = NULL, 
                                   inline   = FALSE))) 

# checkboxGroupInput(inputId = "sampleCheckbox", label = "Samples to be removed")


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
          radioButtons('protPCAColor', 'Color by', choices=c("Group", "Batch"), inline=TRUE)
        ),
        box(
          radioButtons('peptPCAColor', 'Color by', choices=c("Group", "Batch"), inline=TRUE)
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
            "PCA",
            plotOutput("filtProtPCA")
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
            "PCA",
            plotOutput("filtPeptPCA")
          )
        ) 
      ),
      
      
      tweaks,
      fluidRow(column(width = 4, controlsRmSamples))
      
      
      # checkboxGroupInput(inputId = "sampleCheckbox", label = "Samples to be removed")
      # actionButton("updateSampleFilterButton", "Updates Samples"),
      # textOutput("selected_var")
    ),
    
    
    
    # Normalization
    tabItem(
      tabName = "norm",
      
      h2("Proteins"), # cleanProteins()
      fluidRow(
        tabBox(
          id = "normalizationTab",
          width = 12,
          height = round(0.8 * screenHeight),
          tabPanel(
            "Total Intensity",
            plotOutput("totalIntensity_barplot")
          ),
          tabPanel(
            "PCA",
            selectInput("normMethodPCA", "Normalization Method:", 
                        c("Log2", "Median", "Mean", "VSN",
                          "Quantile", "Cyclic Loess", "RLR", "Global Intensity")),
            radioButtons('groupBatchPCA', 'Color by', choices=c("Group", "Batch"), inline=TRUE),
            # selectInput("groupBatchPCA", "Color code Group or Batch:", 
            #             c("Group", "Batch")),
            plotOutput("pcaPlot")
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
            "Correlation heatmap",
            selectInput("normMethodCorrelationHeatmap", "Normalization Method:", 
                        c("Log2", "Median", "Mean", "VSN",
                          "Quantile", "Cyclic Loess", "RLR", "Global Intensity")),
            plotOutput("cor_heatmap")
          ),
          tabPanel(
            "LogRatio density",
            plotOutput("logRatio_density")
          ),
          tabPanel(
            "Missing Values",
            checkboxInput("showAllProtein_NA_heatmap", "Show all protein", value = FALSE, width = NULL),
            plotOutput("NA_heatmap")
          )
        )
      )
    ),
    
    # DAtest
    tabItem(
      tabName = "DAtest",
      
      fluidRow(inputPanel(
        selectInput("normMethod", "Normalization Method:", 
                    c("Log2", "Median", "Mean", "VSN",
                      "Quantile", "Cyclic Loess", "RLR", "Global Intensity")),
        selectInput("imputationMethod", "Imputation Method:", 
                    c("No Imputation", "KNN", "QRILC", "MinDet", "MinProb", "Min", "Zero"))
      ),
      shinySaveButton("saveNormProtein", "Save file", "Save file as ...", filetype=list(txt="txt"))
      
      ),
      
      numericInput(inputId = "DAtestNumberTests", 
                   label = "Number of times to run the tests", 
                   value = 20, min = 1, step = 1),
      numericInput(inputId = "DAtestEffectsize", 
                   label = "The effect size for the spike-ins", 
                   value = 5, min = 0, step = 0.1),
      numericInput(inputId = "DAtestCores", 
                   label = "Number of cores to use for parallel computing", 
                   value = min(detectCores()-1, 10), min = 1, step = 1),
      selectInput(inputId = "groupSelection1", 
                  label = "Group 1",
                  choices = ""),
      selectInput(inputId = "groupSelection2", 
                  label = "Group 2",
                  choices = ""),
      
      # tweaks,
      fluidRow(column(width = 7, controlsDAtest)),
      
      shiny::actionButton(inputId = "goButtonDAtest", label = "Go!", width = '100%', style='font-size:150%'),
      
      div(DTOutput("DAtestResults"), style = "overflow-y: scroll;overflow-x: scroll;"),
      
      plotOutput(outputId = "DAtestFigure")
    ),
    
    # DAtest Power
    tabItem(
      tabName = "DAtestPower",
      fluidRow(inputPanel(
        selectInput(inputId = "DAtest4power", label = "Estimating statistical power for:", choices = unname(DAtestTestNames[eval(formals(testDA)$tests)]))
      )),
      
      shiny::actionButton(inputId = "goButtonDAPower", label = "Go!", width = '100%', style='font-size:150%'),
      
      div(DTOutput("DAtestPowerResults"), style = "overflow-y: scroll;overflow-x: scroll;"),
      
      plotOutput(outputId = "DAtestPowerFigure")
      
    )
  )
)



ui <- dashboardPage(header, sidebar, body)

server <- function(input, output, session) {
  peptideAnnotationColums = c("id", "Protein.group.IDs", "Leading.razor.protein", "Gene.names")
  proteinAnnotationColums = c("id")
  
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
                               })
  })
  
  observe({
    req(peptides())
    req(proteins())
    req(metaData1())
    req(input$metaData)
    meta = metaData2()
    
    if(meta$Custom.Sample.Name[1] == ""){
      groups = unique(meta$Group[!meta$Peptide.Sample.Names %in% input$sampleCheckbox])
    } else {
      groups = unique(meta$Group[!meta$Custom.Sample.Names %in% input$sampleCheckbox])
    }
    updateSelectInput(session, inputId = "groupSelection1", choices = groups,
                      selected = groups[1])
    updateSelectInput(session, inputId = "groupSelection2", choices = groups,
                      selected = ifelse(length(groups)>1, groups[2], groups[1]))
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
    filteredPeptides[filteredPeptides == 0] = NA
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
    selectData = !colnames(filteredProteins) %in% proteinAnnotationColums
    filteredProteins[,selectData][filteredProteins[,selectData] == 0] = NA
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
  
  # PCA Proteins
  output$filtProtPCA <- renderPlot({
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    
    data = proteins[, meta$Peptide.Sample.Names]
    data = data[!apply(is.na(data), 1, any),]
    pca = stats::prcomp(t(data))
    
    if(meta$Custom.Sample.Name[1] == ""){
      rownames(meta) = meta$Peptide.Sample.Names
    } else {
      rownames(meta) = meta$Custom.Sample.Names
    }
    
    ggplot2::autoplot(pca, data = meta, colour = input$protPCAColor, x = 1, y = 2, label  = TRUE)
  })
  
  # PCA Peptides
  output$filtPeptPCA <- renderPlot({
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    
    data = peptides[, meta$Peptide.Sample.Names]
    data = data[!apply(is.na(data), 1, any),]
    pca = stats::prcomp(t(data))
    
    if(meta$Custom.Sample.Name[1] == ""){
      rownames(meta) = meta$Peptide.Sample.Names
    } else {
      rownames(meta) = meta$Custom.Sample.Names
    }
    ggplot2::autoplot(pca, data = meta, colour = input$peptPCAColor, x = 1, y = 2, label  = TRUE)    
  })
  
  
  
  
  
  
  
  
  
  normProteins <- reactive({
    sampleCols = metaDataFiltered()$Peptide.Sample.Names
    proteins <- proteinsSampleFiltered()[,sampleCols]
    if(is.null(proteins)) return(NULL)
    
    normList <- vector("list", 8)
    names(normList) <- c("Log2", "Median", "Mean", "VSN",
                         "Quantile", "Cyclic Loess", "RLR", "Global Intensity")
    normList[["Log2"]] <- logNorm(proteins)
    normList[["Median"]] <- medianNorm(normList[["Log2"]])
    normList[["Mean"]] <- meanNorm(normList[["Log2"]])
    normList[["VSN"]] <- vsnNorm(proteins)
    normList[["Quantile"]] <- quantNorm(normList[["Log2"]])
    normList[["Cyclic Loess"]] <- cycLoessNorm(normList[["Log2"]])
    normList[["RLR"]] <- rlrNorm(normList[["Log2"]])
    normList[["Global Intensity"]] <- giNorm(normList[["Log2"]])
    
    normList
  })
  
  
  observeEvent(input$normalizationTab, {
    choice = input$normalizationTab
    # print(choice)
    if(choice == "PCA"){
      shinyjs::show("normMethodPCA")
      shinyjs::show("groupBatchPCA")
    } else {
      shinyjs::hide("normMethodPCA")
      shinyjs::hide("groupBatchPCA")
    }
    if(choice == "Correlation heatmap"){
      shinyjs::show("normMethodCorrelationHeatmap")
    } else {
      shinyjs::hide("normMethodCorrelationHeatmap")
    }
    if(choice == "Missing Values"){
      shinyjs::show("showAllProtein_NA_heatmap")
    } else {
      shinyjs::hide("showAllProtein_NA_heatmap")
    }
  })
  
  
  
  
  
  output$totalIntensity_barplot <- renderPlot({
    
    normList <- normProteins()
    meta <- metaDataFiltered()
    
    if(is.null(normList) | is.null(meta)) return(NULL)
    # save(normList, meta, file = "savedData.Rdata")
    
    if(meta$Custom.Sample.Name[1] == ""){
      sampleLabels = meta$Peptide.Sample.Names
    } else {
      sampleLabels = meta$Custom.Sample.Names
    }
    
    par(mfrow = c(3,3), mar=c(8,8,3,1))
    for(i in names(normList)){
      barplot(colSums(normList[[i]], na.rm = T), main = i, las = 2, yaxt="n", cex.main = 1.5, col = plasma(ncol(normList[[i]])), names.arg = sampleLabels)
      axis(side = 2, cex.axis=1.5, las = 2)
      # axis(side = 1, at = seq_along(colnames(normList[[i]])), labels = colnames(normList[[i]]), cex.axis=1.5, las = 2)
      if(i == "VSN") mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
      abline(h = max(colSums(normList[[i]], na.rm = T)), lty = 2)
    }
  }, height = round(0.6 * screenHeight))
  
  output$pcaPlot <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)

    data = normList[[input$normMethodPCA]]
    data = data[!apply(is.na(data), 1, any),]
    pca = stats::prcomp(t(data))
    if(meta$Custom.Sample.Name[1] == ""){
      rownames(meta) = meta$Peptide.Sample.Names
    } else {
      rownames(meta) = meta$Custom.Sample.Names
    }
    ggplot2::autoplot(pca, data = meta, colour = input$groupBatchPCA, x = 1, y = 2, label  = TRUE)
  }, height = round(0.5 * screenHeight))
  
  output$PCV_boxplot <- renderPlot({
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
    points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
  }, height = round(0.6 * screenHeight))
  
  output$PMAD_boxplot <- renderPlot({
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
    points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
  }, height = round(0.6 * screenHeight))
  
  output$PEV_boxplot <- renderPlot({
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
    points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
  }, height = round(0.6 * screenHeight))
  
  output$cor_boxplolt <- renderPlot({
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
  }, height = round(0.6 * screenHeight))
  
  output$NA_heatmap <- renderPlot({
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
    
    heatmapMissing(normList[["Log2"]], groups, batch, sampleLabels, input$showAllProtein_NA_heatmap)
  }, height = round(0.7 * screenHeight))
  
  output$cor_heatmap <- renderPlot({
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
  }, height = round(0.6 * screenHeight))
  
  output$logRatio_density <- renderPlot({
    
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    groups <- meta$Group
    
    densityLog2Ratio(normList, groups)
  }, height = round(0.6 * screenHeight))
  
  shiny::observeEvent(input$saveNormProtein,{
    volumes <- c("UserFolder" = getwd())
    shinyFileSave(input, "saveNormProtein", roots=volumes, session=session)
    fileinfo <- parseSavePath(volumes, input$saveNormProtein)
    if (nrow(fileinfo) > 0) {
      normList <- normProteins()
      meta <- metaDataFiltered()
      anno <- proteinsSampleFiltered()[,proteinAnnotationColums]
      if(is.null(normList) | is.null(meta)) return(NULL)
      groups <- meta$Group
      batch <- meta$Batch
      normData <- normList[[input$normMethod]]
      if(input$imputationMethod == "No Imputation"){
        normDataAnno = cbind(anno, normData)
        colnames(normDataAnno)[seq_along(proteinAnnotationColums)] = proteinAnnotationColums
        write.table(normDataAnno, file = as.character(fileinfo$datapath), sep = "\t", row.names = F) 
      } else {
        normDataImp = impute(normData, input$imputationMethod)
        normDataImpAnno = cbind(anno, normDataImp)
        colnames(normDataImpAnno)[seq_along(proteinAnnotationColums)] = proteinAnnotationColums
        write.table(normDataImpAnno, file = as.character(fileinfo$datapath), sep = "\t", row.names = F)
      }
    }
  })
  
  shiny::observeEvent(input$goButtonDAtest, {
    normList <- normProteins()
    meta <- metaDataFiltered()
    
    if(is.null(normList) | is.null(meta)) return(NULL)
    normData <- normList[[input$normMethod]]
    idx = meta$Group %in% c(input$groupSelection1, input$groupSelection2)
    normData = normData[,idx]
    meta = meta[idx,]

    groups <- meta$Group
    batch <- meta$Batch
    
    cat("Excluding: ", DAtestTestNames[DAtestTests[as.numeric(input$checkboxDAtestTests)]], "\n")
    DAtestResults = if(input$imputationMethod == "No Imputation"){
      DAtest(normData, groups, batch, imputed = FALSE, 
             exludedTests = DAtestTests[as.numeric(input$checkboxDAtestTests)],
             R = input$DAtestNumberTests, cores = input$DAtestCores, effectSize = input$DAtestEffectsize)
    } else {
      tempImpute = impute(normData, input$imputationMethod)
      DAtest(tempImpute, groups, batch, imputed = TRUE, 
             exludedTests = DAtestTests[as.numeric(input$checkboxDAtestTests)],
             R = input$DAtestNumberTests, cores = input$DAtestCores, effectSize = input$DAtestEffectsize)
    }
    
    output$DAtestResults <- renderDT({ summary(DAtestResults) })
    output$DAtestFigure <- renderPlot({ plot(DAtestResults) })
    
    updatedTests = DAtestTestNames[rownames(DAtestResults$run.times)]; names(updatedTests) = NULL
    updateSelectInput(session, "DAtest4power", choices = updatedTests)
  })
  
  
  shiny::observeEvent(input$goButtonDAPower, {
    normList <- normProteins()
    meta <- metaDataFiltered()
    
    if(is.null(normList) | is.null(meta)) return(NULL)
    normData <- normList[[input$normMethod]]
    idx = meta$Group %in% c(input$groupSelection1, input$groupSelection2)
    normData = normData[,idx]
    meta = meta[idx,]
    
    groups <- meta$Group
    batch <- meta$Batch
    if(input$imputationMethod != "No Imputation") normData = impute(normData, input$imputationMethod)
    
    DAtestEffectsize = input$DAtestEffectsize
    DAtestPowerResults = powerDA(2^normData, 
                                 predictor = as.character(groups), 
                                 test = names(DAtestTestNames)[DAtestTestNames == input$DAtest4power], 
                                 cores = input$DAtestCores,
                                 relative = FALSE, 
                                 effectSizes = c(DAtestEffectsize/5, DAtestEffectsize/2, DAtestEffectsize, 
                                                 1.5*DAtestEffectsize, 2*DAtestEffectsize, 5*DAtestEffectsize, 
                                                 10*DAtestEffectsize),
                                 R = 5 # more than 5 don't work (bug with error)
                                 # R = input$DAtestNumberTests
    )
    output$DAtestPowerResults <- renderDT({ summary(DAtestPowerResults) })
    output$DAtestPowerFigure <- renderPlot({ plot(DAtestPowerResults) })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
