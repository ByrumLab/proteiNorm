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
library(bit64) # handle large intensities



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("normFunctions.R")
source("functions.R")

# Sets maximum upload size to 5GB
options(shiny.maxRequestSize = 5000*1024^2)
options(stringsAsFactors=FALSE)

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
  list(h3("Samples:"), 
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
          selectInput(inputId = "PeptideFilterMethod", label = "Filtering method Peptide to Protein", choices = "Top3"),
          useShinyjs(),
          shinySaveButton("savePep2Pro", "Filter peptides", "Save file as ...", filetype=list(txt="txt")),
          shinyjs::hidden(p(id = "textFilteringPeptides", "Processing...")),
          
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
      
      numericInput(inputId = "minSamplesPerXGroup", label = "Minimum number of samples with measurement in Y groups", value = 0, min = 0, step = 1),
      numericInput(inputId = "YGroup", label = "Number of Y groups", value = 0, min = 0, step = 1),
      
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
      
      actionButton(inputId = "saveFigures", label = "Save Figures"),
      
      
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
      checkboxGroupInput(inputId = "groupSelection", 
                         label = "Group Selection"),
      
      # tweaks,
      fluidRow(column(width = 7, controlsDAtest)),
      
      "Run DAtest (requires no NA values):",
      shiny::actionButton(inputId = "goButtonDAtest", label = "Go!", width = '100%', style='font-size:150%'),
      
      div(DTOutput("DAtestResults"), style = "overflow-y: scroll;overflow-x: scroll;"),
      
      plotOutput(outputId = "DAtestFigure"),
      actionButton(inputId = "saveDAtestFigure", label = "Save DAtest figure")
    ),
    
    # DAtest Power
    tabItem(
      tabName = "DAtestPower",
      fluidRow(inputPanel(
        selectInput(inputId = "DAtest4power", label = "Estimating statistical power for:", choices = unname(DAtestTestNames[eval(formals(testDA)$tests)]))
      )),
      
      "Run powerDA (requires no NA values):",
      shiny::actionButton(inputId = "goButtonDAPower", label = "Go!", width = '100%', style='font-size:150%'),
      
      div(DTOutput("DAtestPowerResults"), style = "overflow-y: scroll;overflow-x: scroll;"),
      
      plotOutput(outputId = "DAtestPowerFigure"),
      actionButton(inputId = "saveDAtestPower", label = "Save DAtest power figure")
    )
  )
)



ui <- dashboardPage(header, sidebar, body)

server <- function(input, output, session) {
  peptideAnnotationColums = c("id", "Protein.group.IDs", "Leading.razor.protein", "Gene.names")
  proteinAnnotationColums = c("id")
  DAtestResultsGlobal = reactiveValues(test = NULL, power = FALSE)

  shinyjs::disable("goButtonDAtest")
  shinyjs::disable("goButtonDAPower")
  shinyjs::hide("saveDAtestFigure")
  shinyjs::hide("saveDAtestPower")
  
  
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
      start = Sys.time()
      filteredProteins = filterPeptides(peptides = peptides, method = input$PeptideFilterMethod)
      write.table(filteredProteins, file = as.character(fileinfo$datapath), sep = "\t", row.names = F) #, col.names=NA
      cat(" done (")
      cat(difftime(Sys.time(), start)) 
      cat(") \n")
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
    tmp = tmpData[["data"]]
    tmpData[["data"]]
  })
  
  
  ### Creates a Dataframe that will be used to edit the metaData ####
  ### Creates handsontable where metaData1() will be edited ####
  metaData1 <- reactive({
    # req(peptides())
    req(proteins())
    tempTable = data.frame(Protein.Sample.Names = colnames(proteins())[-seq_along(proteinAnnotationColums)],
                           Custom.Sample.Names = "",
                           Group = "",
                           Batch = "")
    
    rhandsontable(tempTable) %>% 
      hot_col("Protein.Sample.Names", readOnly = T)
  })
  #### Outputs the Rhandsontable ####
  output$metaData<- renderRHandsontable({metaData1()})
  
  #### Changes the handsontable back into a dataframe ####
  metaData2 <- reactive({
    hot_to_r(input$metaData)
  }) 
  
  observe({
    # req(peptides())
    req(proteins())
    req(metaData1())
    req(input$metaData)
    meta = metaData2()
    
    updateCheckboxGroupInput(session, inputId = "sampleCheckbox",
                             choices = 
                               if(meta$Custom.Sample.Name[1] == ""){
                                 meta$Protein.Sample.Names
                               } else {
                                 meta$Custom.Sample.Names
                               },
                             selected = 
                               if(meta$Custom.Sample.Name[1] == ""){
                                 meta$Protein.Sample.Names
                               } else {
                                 meta$Custom.Sample.Names
                               })
  })
  
  observe({
    # req(peptides())
    req(proteins())
    req(metaData1())
    req(input$metaData)
    meta = metaData2()
    
    if(meta$Custom.Sample.Name[1] == ""){
      groups = unique(meta$Group[meta$Protein.Sample.Names %in% input$sampleCheckbox])
    } else {
      groups = unique(meta$Group[meta$Custom.Sample.Names %in% input$sampleCheckbox])
    }
    updateCheckboxGroupInput(session, inputId = "groupSelection", 
                             choices = groups,
                             selected = c(groups[1],groups[2]))
    
    if(meta$Custom.Sample.Name[1] == ""){
      selected = meta$Protein.Sample.Names %in% input$sampleCheckbox
    } else {
      selected = meta$Custom.Sample.Names %in% input$sampleCheckbox
    }
    updateNumericInput(session, inputId = "minSamplesPerXGroup", max = min(table(meta[selected,]$Group)))
    updateNumericInput(session, inputId = "YGroup", value = length(unique(meta[selected,]$Group)), max = length(unique(meta[selected,]$Group)))
  })
  
  
  peptidesSampleFiltered <- reactive({
    req(peptides())
    req(metaData1())
    req(input$metaData)
    peptides = peptides()
    meta = metaData2()
    if(meta$Custom.Sample.Name[1] == ""){
      filteredPeptides = peptides[,c(rep(TRUE, length(peptideAnnotationColums)), meta$Protein.Sample.Names %in% input$sampleCheckbox)] # TRUE to keep annotation
    } else {
      filteredPeptides = peptides[,c(rep(TRUE, length(peptideAnnotationColums)), meta$Custom.Sample.Names %in% input$sampleCheckbox)]
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
      selected = meta$Protein.Sample.Names %in% input$sampleCheckbox
      filteredProteins = proteins[,c(rep(TRUE, length(proteinAnnotationColums)), selected)] # TRUE to keep annotation
    } else {
      selected = meta$Custom.Sample.Names %in% input$sampleCheckbox
      filteredProteins = proteins[,c(rep(TRUE, length(proteinAnnotationColums)), selected)]
    }
    selectData = !colnames(filteredProteins) %in% proteinAnnotationColums
    filteredProteins[,selectData][filteredProteins[,selectData] == 0] = NA
    
    # filtering proteins with at least X protein measurements in Y groups 
    myGroups =  unique(meta[selected,]$Group)
    tmp = filteredProteins[,-seq_along(proteinAnnotationColums)]
    numberOfSamplesPerGroup = NULL
    for(group in myGroups){
      dat = tmp[, meta[selected,]$Group %in% group]
      if(is.numeric(dat)) dat = as.data.frame(dat)
      numberOfSamplesPerGroup = cbind(numberOfSamplesPerGroup, apply(dat, 1, FUN = function(x) sum(!is.na(x))))
    }
    # colnames(numberOfSamplesPerGroup) = unique(meta$Group[meta$Protein.Sample.Names %in% input$sampleCheckbox])
    filteredProteins = filteredProteins[apply(numberOfSamplesPerGroup, 1, FUN = function(x) sum(x >= input$minSamplesPerXGroup)) >= input$YGroup, ] # at least X samples in at least Y groups

    if(is.integer64(filteredProteins[1, length(proteinAnnotationColums)+1])){
      print("64")
      tmp = data.frame(filteredProteins[,seq_along(proteinAnnotationColums)], apply(filteredProteins[,-seq_along(proteinAnnotationColums)], 2, bit64::as.double.integer64))
      colnames(tmp)[seq_along(proteinAnnotationColums)] = proteinAnnotationColums
      return(tmp)
    } else if(is.numeric(filteredProteins[1, length(proteinAnnotationColums)+1])){
      return(filteredProteins)
    } else {
      stop("Unknown data class")
    }
  })
  
  metaDataFiltered <- reactive({
    req(proteins())
    req(metaData1())
    req(input$metaData)
    meta = metaData2()
    if(meta$Custom.Sample.Name[1] == ""){
      filteredMetadata = meta[meta$Protein.Sample.Names %in% input$sampleCheckbox,]
    } else {
      filteredMetadata = meta[meta$Custom.Sample.Names %in% input$sampleCheckbox,]
    }
    filteredMetadata
  })
  
  
  
  
  # Boxplot Proteins
  output$filtProtBoxplot <- renderPlot({
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    plotBoxplotProtein(proteins, meta)
    })
  
  # Boxplot Peptides
  output$filtPeptBoxplot <- renderPlot({
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    plotBoxplotPeptide(peptides, meta)

  })
  
  # Histogram Proteins
  output$filtProtHist <- renderPlot({
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    plotHistogramProtein(proteins, meta)
  })
  
  # Histogram Peptides
  output$filtPeptHist <- renderPlot({
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    plotHistogramPeptide(peptides, meta)
  })
  
  # PCA Proteins
  output$filtProtPCA <- renderPlot({
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(proteins) | is.null(meta)) return(NULL)
    plotPCAProtein(proteins, meta, col = input$protPCAColor)
  })
  
  # PCA Peptides
  output$filtPeptPCA <- renderPlot({
    peptides <- peptidesSampleFiltered()
    meta <- metaDataFiltered()
    
    if(is.null(peptides) | is.null(meta)) return(NULL)
    plotPCAPeptide(peptides, meta, col = input$peptPCAColor)
  })
  
  
  
  
  
  
  
  
  
  normProteins <- reactive({
    sampleCols = metaDataFiltered()$Protein.Sample.Names
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
    
    # enable DAtest and powerDA buttons if no NA's in data
    if(!any(is.na(normList[["Log2"]]))){
      shinyjs::enable("goButtonDAtest")
      shinyjs::enable("goButtonDAPower")
    }
    
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
  
  observeEvent(input$imputationMethod, {
    choice = input$imputationMethod
    req(normProteins())
    normList <- normProteins()
    if(any(is.na(normList[[input$normMethod]])) & choice == "No Imputation"){
      shinyjs::disable("goButtonDAtest")
      shinyjs::disable("goButtonDAPower")
    } else {
      shinyjs::enable("goButtonDAtest")
      shinyjs::enable("goButtonDAPower")
    }
  })
  
  
  output$totalIntensity_barplot <- renderPlot({
    
    normList <- normProteins()
    meta <- metaDataFiltered()
    
    if(is.null(normList) | is.null(meta)) return(NULL)
    
    # peptides <- peptidesSampleFiltered()
    # proteins <- proteinsSampleFiltered()
    # save(peptides, proteins, normList, meta, file = "savedData.Rdata")
    
    plotTotInten(normList, meta)
  }, height = round(0.6 * screenHeight))
  
  
  output$pcaPlot <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotPCA(normList, meta, method = input$normMethodPCA, col = input$groupBatchPCA)
  }, height = round(0.5 * screenHeight))
  
  
  output$PCV_boxplot <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotPCV(normList, meta)
  }, height = round(0.6 * screenHeight))
  
  output$PMAD_boxplot <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotPMAD(normList, meta)
  }, height = round(0.6 * screenHeight))
  
  
  output$PEV_boxplot <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotPEV(normList, meta)
  }, height = round(0.6 * screenHeight))
  
  output$cor_boxplolt <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotCOR(normList, meta)
  }, height = round(0.6 * screenHeight))
  
  output$NA_heatmap <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotNaHM(normList, meta, show = input$showAllProtein_NA_heatmap)
  }, height = round(0.7 * screenHeight))
  
  output$cor_heatmap <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotCorHM(normList, meta, method = input$normMethodCorrelationHeatmap)
  }, height = round(0.6 * screenHeight))
  
  
  output$logRatio_density <- renderPlot({
    normList <- normProteins()
    meta <- metaDataFiltered()
    if(is.null(normList) | is.null(meta)) return(NULL)
    plotLogRatio(normList, meta)
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
      if(meta$Custom.Sample.Name[1] != "") colnames(normData) = meta$Custom.Sample.Names
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
    idx = meta$Group %in% input$groupSelection
    normData = normData[,idx]
    meta = meta[idx,]
    
    groups <- meta$Group
    batch <- meta$Batch
    
    # save(norm=List, meta, normData, file = "DAtest.Rdata")
    
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
    
    DAtestResultsGlobal$test <- DAtestResults
    shinyjs::show("saveDAtestFigure")
    
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
    idx = meta$Group %in% input$groupSelection
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
    
    DAtestResultsGlobal$power <- DAtestPowerResults
    shinyjs::show("saveDAtestPower")
    
    output$DAtestPowerResults <- renderDT({ summary(DAtestPowerResults) })
    output$DAtestPowerFigure <- renderPlot({ plot(DAtestPowerResults) })
  })
  
  observeEvent(input$saveFigures, {
    peptides <- peptides()
    proteins <- proteinsSampleFiltered()
    meta <- metaDataFiltered()
    normList <- normProteins()

    if(!is.null(peptides) & !is.null(meta)){
      peptides <- peptidesSampleFiltered()
      png("BoxplotPeptide.png", width = 960, height = 960)
      plotBoxplotPeptide(peptides, meta)
      dev.off()

      png("HistogramPeptide.png", width = 960, height = 960)
      plotHistogramPeptide(peptides, meta)
      dev.off()
      
      plotPCAPeptide = plotPCAPeptide(peptides, meta, col = input$peptPCAColor)
      png("PCAPeptide.png", width = 960, height = 960)
      plot(plotPCAPeptide)
      dev.off()
    }
    
    if(!is.null(proteins) & !is.null(meta)){
      png("BoxplotProtein.png", width = 960, height = 960)
      plotBoxplotProtein(proteins, meta)
      dev.off()
      
      png("HistogramProtein.png", width = 960, height = 960)
      plotHistogramProtein(proteins, meta)
      dev.off()
      
      plotPCAProtein = plotPCAProtein(proteins, meta, col = input$protPCAColor)
      png("PCAProtein.png", width = 960, height = 960)
      plot(plotPCAProtein)
      dev.off()
    }


    
    if(!is.null(normList) & !is.null(meta)){
      png("TotalIntensity.png", width = 960, height = 960)
      plotTotInten(normList, meta)
      dev.off()
      
      plotPCA = plotPCA(normList, meta, method = input$normMethodPCA, col = input$groupBatchPCA)
      png("PCA.png", width = 960, height = 960)
      plot(plotPCA)
      dev.off()
      
      png("PCV.png", width = 960, height = 960)
      plotPCV(normList, meta)
      dev.off()
      
      png("PMAD.png", width = 960, height = 960)
      plotPMAD(normList, meta)
      dev.off()
      
      png("PEV.png", width = 960, height = 960)
      plotPEV(normList, meta)
      dev.off()
      
      png("COR.png", width = 960, height = 960)
      plotCOR(normList, meta)
      dev.off()
      
      png("NA-HM.png", width = 960, height = 960)
      plotNaHM(normList, meta, show = input$showAllProtein_NA_heatmap)
      dev.off()
      
      png("Cor-HM.png", width = 960, height = 960)
      plotCorHM(normList, meta, method = input$normMethodCorrelationHeatmap)
      dev.off()
      
      png("LogRatio.png", width = 960, height = 960)
      plotLogRatio(normList, meta)
      dev.off()
    }
  })
  
  observeEvent(input$saveDAtestFigure, {
    DAplot = DAtestResultsGlobal$test
    test = plot(DAplot)
    png("DAtest.png", width = 960, height = 960)
    plot(test)
    dev.off()
  })
  
  observeEvent(input$saveDAtestPower, {
    DAplot = DAtestResultsGlobal$power
    test = plot(DAplot)
    png("DAtestPower.png", width = 960, height = 960)
    plot(test)
    dev.off()
  })
  
   
}

# Run the application 
shinyApp(ui = ui, server = server)
