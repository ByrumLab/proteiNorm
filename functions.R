extractColumnNames = function(columnNames){
  # check data type
  isTMT = any(grepl("Reporter\\.intensity\\.corrected\\.[[:digit:]]+", columnNames))
  
  if(isTMT){
    ## TMT
    pattern = "Reporter\\.intensity\\.corrected\\.[[:digit:]]+\\.."
    rawDataCols = grep(pattern, columnNames, value = TRUE)
    if(identical(rawDataCols, character(0))){
      pattern = "Reporter\\.intensity\\.corrected\\.[[:digit:]]+"
      rawDataCols = grep(pattern, columnNames, value = TRUE)
    }
  } else {
    ## Label Free
    pattern = "Intensity\\.[[:digit:]]+"
    rawDataCols = grep(pattern, columnNames, value = TRUE)
  }
  
  return(list(isTMT = isTMT, rawDataCols = rawDataCols))
}

extractPeptideData = function(rawData, peptideAnnotationColums){
  # quality filter 
  remove = rawData$Reverse == "+" | rawData$Potential.contaminant == "+"
  rawData = rawData[!remove,]
  
  # extract relevant data
  tmpCol = extractColumnNames(colnames(rawData))
  isTMT = tmpCol[["isTMT"]]
  rawDataCols = tmpCol[["rawDataCols"]]
  extractedData = rawData[,c(peptideAnnotationColums, rawDataCols)]
  
  return(list(data = extractedData, isTMT = isTMT))
}

extractProteinData = function(rawData, proteinAnnotationColums){
  # quality filter 
  if(all(c("Reverse", "Potential.contaminant", "Only.identified.by.site") %in% colnames(rawData))){
    remove = rawData$Reverse == "+" | rawData$Potential.contaminant == "+" | rawData$Only.identified.by.site == "+"
    rawData = rawData[!remove,]
  }
  # extract relevant data
  tmpCol = extractColumnNames(colnames(rawData))
  isTMT = tmpCol[["isTMT"]]
  rawDataCols = tmpCol[["rawDataCols"]]
  extractedData = rawData[,c(proteinAnnotationColums, rawDataCols)]
  return(list(data = extractedData, isTMT = isTMT))
}


outlier = function(vector){
  vector < quantile(vector)[2] - 1.5 * IQR(vector) | quantile(vector)[4] + 1.5 * IQR(vector) < vector
}


filterPeptides = function(peptides){
  peptides = as.data.frame(peptides)
  # save(peptides, file = "peptidesPreFilter.Rdata")
  # print("saved")
  peptides$proteinIDs = unlist(lapply(strsplit(peptides$Protein.group.IDs, ";"), function(x) x[1]))
  
  colmnNames = extractColumnNames(colnames(peptides))[["rawDataCols"]]
  FilteredProteins = NULL
  for(i in unique(peptides$proteinIDs)){
    tempProtein = peptides[peptides$proteinIDs == i, colmnNames]
    if(nrow(tempProtein) == 1){
      FilteredProteins = rbind(FilteredProteins, tempProtein)
    } else {
      tempProtein_log = log2(tempProtein)
      tempOutliers = apply(tempProtein_log, 2, outlier)
      keepPeptide = !apply(tempOutliers, 1, any)
      tempFilteredProtein = colSums(tempProtein[keepPeptide,])
      FilteredProteins = rbind(FilteredProteins, tempFilteredProtein)
    }
  }
  FilteredProteins = cbind(id = unique(peptides$proteinIDs), FilteredProteins)
  return(FilteredProteins)
}




# ---
# findDensityCutoff
# ---
# Finds upper limit of the histogram so each box contains at least cutoffPercent of the data

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


PCV = function(data, groups){
  PCV = NULL
  for(group in unique(groups)){
    tempData = data[,groups %in% group]
    CVs = genefilter::rowSds(tempData, na.rm = FALSE) /
      rowMeans(tempData, na.rm = FALSE)
    PCV[group] = mean(CVs, na.rm = T)
  }
  return(PCV)
}

PMAD = function(data, groups){
  PMAD = NULL
  for(group in unique(groups)){
    tempData = data[,groups %in% group]
    MAD = matrixStats::rowMads(tempData, na.rm = FALSE)
    PMAD[group] = mean(MAD, na.rm = T)
  }
  return(PMAD)
}



PEV = function(data, groups){
  PEV = NULL
  for(group in unique(groups)){
    tempData = data[,groups %in% group]
    
    rowNonNACnt = rowSums(!is.na(tempData)) - 1
    EV = rowNonNACnt * matrixStats::rowVars(tempData, na.rm = FALSE)
    PEV[group] = sum(EV, na.rm = TRUE)/sum(rowNonNACnt, na.rm = TRUE)
  }
  return(PEV)
}




COR = function(data, groups){
  COR = NULL
  for(group in unique(groups)){
    corGroup = NULL
    tempData = data[,groups %in% group]
    
    corVals = stats::cor(tempData, use = "pairwise.complete.obs", method = "pearson")
    for (index in seq_len(ncol(corVals) - 1)) {
      corGroup <- c(corGroup, corVals[index, -(seq_len(index)), drop = "FALSE"])
    }
    COR[[group]] = corGroup
  }
  return(COR)
}



colorBatch = function(batch){
  batchCol =  unlist(ifelse(length(unique(batch)) == 1, rainbow(1, start = 0.5), list(rainbow(length(unique(batch))))))
  names(batchCol) = sort(unique(batch))
  return(batchCol)
}

colorGroup = function(groups){
  groupCol = viridis::viridis(length(unique(groups)))
  names(groupCol) = sort(unique(groups))
  return(groupCol)
}

heatmapMissing = function(data, groups, batch, sampleLabels){
  
  missing = !is.na(data)
  complete = apply(missing, 1, all)
  completeNA = apply(!missing, 1, all)
  missing = missing[!complete & !completeNA,]
  
  
  
  batchCol =  colorBatch(batch)
  groupCol = colorGroup(groups)
  
  hm_clust = heatmapClustered(groups, batch, missing, groupCol, batchCol, sampleLabels)
  hm_group = heatmapGroup(groups, batch, missing, groupCol, batchCol, sampleLabels)
  hm_batch = heatmapBatch(groups, batch, missing, groupCol, batchCol, sampleLabels)
  
  draw(hm_clust+hm_group+hm_batch, heatmap_legend_side = "right", annotation_legend_side = "right", ht_gap = unit(2, "cm"), column_title = "Missing values")
}


heatmapClustered = function(groups, batch, missing, groupCol, batchCol, sampleLabels){
  ColAnn <- HeatmapAnnotation(Sample = groups, Batch = batch, col = list(Sample = groupCol, Batch = batchCol), 
                              annotation_legend_param = list(
                                Sample = list(
                                  title = "Sample",
                                  at = sort(unique(groups)),
                                  labels = paste("Group", sort(unique(groups)))
                                ),
                                Batch = list(
                                  title = "Batch",
                                  at = sort(unique(batch)),
                                  labels = paste("Batch", sort(unique(batch)))
                                )
                              ))
  
  hm_clust = Heatmap(missing+0, col = c("white", "black"), column_names_side = "top",  column_title = "Clustered",
                     show_row_names = FALSE, show_column_names = TRUE, name = "Missing values pattern", 
                     column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(at = c(0, 1), labels = c("Missing value", "Valid value")), 
                     top_annotation = ColAnn, column_labels = sampleLabels)
  # draw(hm_clust, heatmap_legend_side = "right", annotation_legend_side = "right")
  return(hm_clust)
}




heatmapGroup = function(groups, batch, missing, groupCol, batchCol, sampleLabels){
  orderGroups = order(groups, batch)
  groups_groupSorted = groups[orderGroups]
  batch_groupSorted = batch[orderGroups]
  missing_groupSorted = missing[,orderGroups]
  ColAnn <- HeatmapAnnotation(Sample = groups_groupSorted, Batch = batch_groupSorted, col = list(Sample = groupCol, Batch = batchCol), 
                              annotation_legend_param = list(
                                Sample = list(
                                  title = "Sample",
                                  at = sort(unique(groups_groupSorted)),
                                  labels = paste("Group", sort(unique(groups_groupSorted)))
                                ),
                                Batch = list(
                                  title = "Batch",
                                  at = sort(unique(batch_groupSorted)),
                                  labels = paste("Batch", sort(unique(batch_groupSorted)))
                                )
                              ))
  hm_group = Heatmap(missing_groupSorted+0, col = c("white", "black"), column_names_side = "top",  column_title = "Sorted by Groups",
                     show_row_names = FALSE, show_column_names = TRUE, name = "Missing values pattern", 
                     column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(at = c(0, 1), labels = c("Missing value", "Valid value")), 
                     top_annotation = ColAnn, cluster_columns = FALSE, column_labels = sampleLabels)
  # draw(hm_group, heatmap_legend_side = "right", annotation_legend_side = "right")
  return(hm_group)
}

heatmapBatch = function(groups, batch, missing, groupCol, batchCol, sampleLabels){
  orderBatch = order(batch, groups)
  groups_batchSorted = groups[orderBatch]
  batch_batchSorted = batch[orderBatch]
  missing_batchSorted = missing[,orderBatch]
  ColAnn <- HeatmapAnnotation(Sample = groups_batchSorted, Batch = batch_batchSorted, col = list(Sample = groupCol, Batch = batchCol),
                              annotation_legend_param = list(
                                Sample = list(
                                  title = "Sample",
                                  at = sort(unique(groups_batchSorted)),
                                  labels = paste("Group", sort(unique(groups_batchSorted)))
                                ),
                                Batch = list(
                                  title = "Batch",
                                  at = sort(unique(batch_batchSorted)),
                                  labels = paste("Batch", sort(unique(batch_batchSorted)))
                                )
                              ))
  hm_batch = Heatmap(missing_batchSorted+0, col = c("white", "black"), column_names_side = "top",  column_title = "Sorted by Batch",
                     show_row_names = FALSE, show_column_names = TRUE, name = "Missing values pattern", 
                     column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(at = c(0, 1), labels = c("Missing value", "Valid value")), 
                     top_annotation = ColAnn, cluster_columns = FALSE, column_labels = sampleLabels)
  # draw(hm_group, heatmap_legend_side = "right", annotation_legend_side = "right")
  return(hm_batch)
}


heatmapCorr = function(data, groups, batch, sampleLabels){
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  
  ColAnn <- HeatmapAnnotation(Sample = groups, 
                              col = list(Sample = colorGroup(groups)), 
                              annotation_legend_param = list(
                                Sample = list(
                                  title = "Sample",
                                  at = sort(unique(groups)),
                                  labels = paste("Group", sort(unique(groups)))
                                )
                              ))
  RowAnn <- rowAnnotation(Batch = batch, 
                          col = list(Batch = colorBatch(batch)), 
                          annotation_legend_param = list(
                            Batch = list(
                              title = "Batch",
                              at = sort(unique(batch)),
                              labels = paste("Batch", sort(unique(batch)))
                            )
                          ))
  
  hm_corr = Heatmap(cor_mat, 
                    col = circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7)),RColorBrewer::brewer.pal(8, "Reds")), 
                    heatmap_legend_param = list(color_bar = "continuous", 
                                                legend_direction = "horizontal", 
                                                legend_width = unit(5, "cm"), 
                                                title_position = "topcenter"), 
                    name = "Pearson correlation", 
                    column_names_gp = gpar(fontsize = 12), 
                    row_names_gp = gpar(fontsize = 12), 
                    top_annotation = ColAnn, 
                    left_annotation = RowAnn,
                    column_labels = sampleLabels, row_labels = sampleLabels
  )
  
  draw(hm_corr, heatmap_legend_side = "top")
}



densityLogFC = function(normList, groups){
  logFC = NULL
  groupList = sort(unique(groups))
  for(method in names(normList)){
    for(g1 in 1:(length(groupList)-1)){
      for(g2 in (g1+1):length(groupList)){
        
        logFC[[method]] = rowMeans(normList[[method]][,groups == groupList[g1]], na.rm = T) - 
          rowMeans(normList[[method]][,groups == groupList[g2]], na.rm = T)
      }
    }
  }
  
  # minX = min(unlist(lapply(logFC, FUN = function(x) min(density(x, na.rm = T)$x))))
  # maxX = max(unlist(lapply(logFC, FUN = function(x) max(density(x, na.rm = T)$x))))
  maxY = max(unlist(lapply(logFC, FUN = function(x) max(density(x, na.rm = T)$y))))
  
  minX = 0.5 * min(unlist(min(density(logFC[["vsnNorm"]], na.rm = T)$x)))
  maxX = 0.5 * max(unlist(max(density(logFC[["vsnNorm"]], na.rm = T)$x)))
  
  plot(NA, xlim = c(minX, maxX), ylim = c(0,maxY), xlab = "logFC", ylab = "Density")
  abline(v = 0, lty = 2, col = "grey")
  
  for(method in names(normList)){
    lines(density(logFC[[method]], na.rm = T), 
          col = rainbow(length(logFC))[which(names(logFC) %in% method)], 
          lty = ifelse(which(names(logFC) %in% method) %% 2 == 0, 2, 1),
          lwd = 2)
  }
  legend("topright", names(logFC), lty=1:2, col = rainbow(length(logFC)), lwd = 2)
}





impute = function(normData, imputMethod){ 
  
  # seNormData = SummarizedExperiment(normData)
  # mssNormData = as(seNormData, "MSnSet")
  switch(imputMethod,
         # "bpca" = pcaMethods::bpca(Matrix = normData),
         "KNN" = impute.knn(normData, rowmax = 0.9)$data,
         "QRILC" = imputeLCMD::impute.QRILC(dataSet.mvs = normData)[[1]],
         # "MLE" = ,
         "MinDet" = imputeLCMD::impute.MinDet(dataSet.mvs = normData),
         "MinProb" = imputeLCMD::impute.MinProb(dataSet.mvs = normData),
         # "man" = ,
         "min" = {
           normData[is.na(normData)] = min(normData, na.rm = T)
           normData
         },
         "zero" = {
           normData[is.na(normData)] = 0
           normData
         },
         # "mixed" = ,
         # "nbavg" = ,
         stop("Imputation method not found")
  )
}


DAtest = function(normData, groups, batch, imputed, exludedTests, R, cores, effectSize){ # will need to include other covariates later
  predictor = as.character(groups)
  covar = list(batch = as.character(batch))
  
  availableTests  = eval(formals(testDA)$tests)
  performTests = availableTests[!availableTests %in% exludedTests]
  
  # Prevent error if only 1 batch
  if(length(unique(batch)) == 1){ 
    res <- testDA(data = 2^normData, 
                  predictor = predictor, 
                  cores = cores, 
                  R = R, 
                  relative = FALSE, 
                  tests = performTests,
                  effectSize = effectSize)
    return(res)
  }
  
  if(!imputed){
    ## DAtest
    res <- testDA(data = 2^normData, 
                  predictor = predictor, 
                  cores = cores, 
                  R = R, 
                  covars = covar, 
                  relative = FALSE, 
                  tests = performTests,
                  effectSize = effectSize)
    return(res)
  }
  
  if(imputed){
    ## combat
    tempCombat = ComBat(dat=normData, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
    
    ## DAtest
    res <- testDA(data = 2^tempCombat, 
                  predictor = predictor, 
                  cores = cores, 
                  R = R, 
                  relative = FALSE, 
                  tests = performTests,
                  effectSize = effectSize)
    return(res)
  } 
}





