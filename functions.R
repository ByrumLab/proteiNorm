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
  vector < quantile(vector, na.rm = T)[2] - 1.5 * IQR(vector, na.rm = T) | 
    quantile(vector, na.rm = T)[4] + 1.5 * IQR(vector, na.rm = T) < vector
}


filterPeptides = function(peptides, method){
  FilteredProteins = switch(method,
                            Top3 = Top3Filter(peptides),
                            stop("Method not found")
  )
  return(FilteredProteins)
}


Top3Filter = function(peptides){
  peptides = as.data.frame(peptides)
  peptides$proteinIDs = peptides$Leading.razor.protein
  
  colmnNames = extractColumnNames(colnames(peptides))[["rawDataCols"]]
  FilteredProteins = NULL
  for(i in unique(peptides$proteinIDs)){
    tempProtein = peptides[peptides$proteinIDs == i, colmnNames]
    if(nrow(tempProtein) == 1){
      FilteredProteins = rbind(FilteredProteins, tempProtein)
    } else {
      tempProtein[tempProtein == 0] = NA
      top3mean = apply(tempProtein, 2, function(x) mean(x[order(-x)][1:min(3, length(x))], na.rm = TRUE))
      FilteredProteins = rbind(FilteredProteins, top3mean)
    }
  }
  rownames(FilteredProteins) = unique(peptides$proteinIDs)
  FilteredProteins[FilteredProteins == 0] = NA
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

heatmapMissing = function(data, groups, batch, sampleLabels, showAllProtein){
  missing = !is.na(data)
  if(!showAllProtein){
    complete = apply(missing, 1, all)
    completeNA = apply(!missing, 1, all)
    missing = missing[!complete & !completeNA,]
  }
  
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
                     top_annotation = ColAnn, cluster_columns = FALSE, column_labels = sampleLabels[orderGroups])
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
                     top_annotation = ColAnn, cluster_columns = FALSE, column_labels = sampleLabels[orderBatch])
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



densityLog2Ratio = function(normList, groups){
  log2Ratio = NULL
  groupList = sort(unique(groups))
  for(method in names(normList)){
    for(g1 in 1:(length(groupList)-1)){
      for(g2 in (g1+1):length(groupList)){
        
        log2Ratio[[method]] = c(log2Ratio[[method]], rowMeans(normList[[method]][,groups == groupList[g1]], na.rm = T) - 
                                  rowMeans(normList[[method]][,groups == groupList[g2]], na.rm = T))
      }
    }
  }
  
  # minX = min(unlist(lapply(log2Ratio, FUN = function(x) min(density(x, na.rm = T)$x))))
  # maxX = max(unlist(lapply(log2Ratio, FUN = function(x) max(density(x, na.rm = T)$x))))
  maxY = max(unlist(lapply(log2Ratio, FUN = function(x) max(density(x, na.rm = T)$y))))
  
  minX = 0.5 * min(unlist(min(density(log2Ratio[["VSN"]], na.rm = T)$x)))
  maxX = 0.5 * max(unlist(max(density(log2Ratio[["VSN"]], na.rm = T)$x)))
  
  plot(NA, xlim = c(minX, maxX), ylim = c(0,maxY), xlab = "Log2 ratio", ylab = "Density", main = "Log2-ratio")
  abline(v = 0, lty = 2, col = "grey")
  
  for(method in names(normList)){
    lines(density(log2Ratio[[method]], na.rm = T), 
          col = rainbow(length(log2Ratio))[which(names(log2Ratio) %in% method)], 
          lty = ifelse(which(names(log2Ratio) %in% method) %% 2 == 0, 2, 1),
          lwd = 2)
  }
  legend("topright", names(log2Ratio), lty=1:2, col = rainbow(length(log2Ratio)), lwd = 2)
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
         "Min" = {
           normData[is.na(normData)] = min(normData, na.rm = T)
           normData
         },
         "Zero" = {
           normData[is.na(normData)] = 0
           normData
         },
         # "mixed" = ,
         # "nbavg" = ,
         stop("Imputation method not found")
  )
}


DAtestTestNames = c(
  "adx" = "ALDEx2 (adx)",
  "bay" = "baySeq (bay)",
  "ds2" = "DESeq2 man. geoMeans (ds2)",
  "ds2x" = "DESeq2 (ds2x)",
  "per" = "Permutation (per)",
  "znb" = "ZI-NegBin GLM (znb)",
  "zpo" = "ZI-Poisson GLM (zpo)",
  "msf" = "MgSeq Feature (msf)",
  "zig" = "MgSeq ZIG (zig)",
  "erq" = "EdgeR qll - TMM (erq)",
  "erq2" = "EdgeR qll - RLE (erq2)",
  "neb" = "Negbinom GLM (neb)",
  "qpo" = "Quasi-Poisson GLM (qpo)",
  "poi" = "Poisson GLM (poi)",
  "sam" = "SAMseq (sam)",
  "lrm" = "Linear regression (lrm)",
  "llm" = "Log Linear reg. (llm)",
  "llm2" = "Log Linear reg. 2 (llm2)",
  "lma" = "Linear model - ALR (lma)",
  "lmc" = "Linear model - CLR (lmc)",
  "ere" = "EdgeR exact - TMM (ere)",
  "ere2" = "EdgeR exact - RLE (ere2)",
  "pea" = "Pearson (pea)",
  "spe" = "Spearman (spe)",
  "wil" = "Wilcox (wil)",
  "kru" = "Kruskal-Wallis (kru)",
  "qua" = "Quade (qua)",
  "fri" = "Friedman (fri)",
  "ttt" = "t-test (ttt)",
  "ltt" = "Log t-test (ltt)",
  "ltt2" = "Log t-test2 (ltt2)",
  "tta" = "t-test - ALR (tta)",
  "ttc" = "t-test - CLR (ttc)",
  "aov" = "ANOVA (aov)",
  "lao" = "Log ANOVA (lao)",
  "lao2" = "Log ANOVA 2 (lao2)",
  "aoa" = "ANOVA - ALR (aoa)",
  "aoc" = "ANOVA - CLR (aoc)",
  "vli" = "LIMMA voom (vli)",
  "lim" = "LIMMA (lim)",
  "lli" = "Log LIMMA (lli)",
  "lli2" = "Log LIMMA 2 (lli2)",
  "lia" = "LIMMA - ALR (lia)",
  "lic" = "LIMMA - CLR (lic)"
)



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






################# figures
plotBoxplotProtein = function(proteins, meta){
  sampleCols = meta$Protein.Sample.Names
  if(meta$Custom.Sample.Name[1] == ""){
    sampleLabels = meta$Protein.Sample.Names
  } else {
    sampleLabels = meta$Custom.Sample.Names
  }
  par(mar=c(11,5,4,2))
  boxplot(proteins[, sampleCols], outline=FALSE, col=colorGroup(meta$Group)[meta$Group], main="Boxplot of Proteins",
          xlab="", ylab="", names=sampleLabels, las = 2)
  mtext(side = 2, text = "Intensity", line = 4)
  legend("top", inset=c(0, -.14), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
         col=colorGroup(meta$Group), pch=15, horiz=TRUE)
}

plotBoxplotPeptide = function(peptides, meta){
  sampleCols = meta$Protein.Sample.Names
  if(meta$Custom.Sample.Name[1] == ""){
    sampleLabels = meta$Protein.Sample.Names
  } else {
    sampleLabels = meta$Custom.Sample.Names
  }
  par(mar=c(11,5,4,2))
  boxplot(peptides[, sampleCols], outline=FALSE, col=colorGroup(meta$Group)[meta$Group], main="Boxplot of Peptides",
          xlab="", ylab="", names=sampleLabels, las = 2)
  mtext(side = 2, text = "Corrected Intensity", line = 4)
  legend("top", inset=c(0, -.14), levels(as.factor(meta$Group)), bty="n", xpd=TRUE,
         col=colorGroup(meta$Group), pch=15, horiz=TRUE)
}


plotHistogramProtein = function(proteins, meta){
  sampleCols = meta$Protein.Sample.Names
  longdat <- unlist(proteins[, sampleCols])
  densityCutoff <- findDensityCutoff(longdat)
  longdat <- longdat[longdat < densityCutoff]
  breakSeq <- seq(max(min(longdat, na.rm=TRUE), 0), max(c(longdat, densityCutoff), na.rm = TRUE),
                  length.out=31)
  hist(longdat, main="Histogram of Proteins", breaks=breakSeq, xlab="Corrected Intensity", ylab="Counts")
}

plotHistogramPeptide = function(proteins, meta){
  sampleCols = meta$Protein.Sample.Names
  longdat <- unlist(peptides[, sampleCols])
  densityCutoff <- findDensityCutoff(longdat)
  longdat <- longdat[longdat < densityCutoff]
  breakSeq <- seq(max(min(longdat, na.rm=TRUE), 0), max(c(longdat, densityCutoff), na.rm = TRUE),
                  length.out=31)
  hist(longdat, main="Histogram of Peptides", breaks=breakSeq, xlab="Corrected Intensity", ylab="Counts")
}

plotPCAProtein = function(proteins, meta, col){
  data = proteins[, meta$Protein.Sample.Names]
  data = data[!apply(is.na(data), 1, any),]
  pca = stats::prcomp(t(data))
  
  if(meta$Custom.Sample.Name[1] == ""){
    rownames(meta) = meta$Protein.Sample.Names
  } else {
    rownames(meta) = meta$Custom.Sample.Names
  }
  
  ggplot2::autoplot(pca, data = meta, colour = col, x = 1, y = 2, label  = TRUE)
}


plotPCAPeptide = function(proteins, meta, col){
  data = peptides[, meta$Protein.Sample.Names]
  data = data[!apply(is.na(data), 1, any),]
  pca = stats::prcomp(t(data))
  
  if(meta$Custom.Sample.Name[1] == ""){
    rownames(meta) = meta$Protein.Sample.Names
  } else {
    rownames(meta) = meta$Custom.Sample.Names
  }
  ggplot2::autoplot(pca, data = meta, colour = col, x = 1, y = 2, label  = TRUE)  
}


plotTotInten = function(normList, meta){
  if(meta$Custom.Sample.Name[1] == ""){
    sampleLabels = meta$Protein.Sample.Names
  } else {
    sampleLabels = meta$Custom.Sample.Names
  }
  
  par(mfrow = c(3,3), mar=c(10,8,3,1))
  for(i in names(normList)){
    barplot(colSums(normList[[i]], na.rm = T), main = i, las = 2, yaxt="n", cex.main = 1.5, col = plasma(ncol(normList[[i]])), names.arg = sampleLabels)
    axis(side = 2, cex.axis=1.5, las = 2)
    # axis(side = 1, at = seq_along(colnames(normList[[i]])), labels = colnames(normList[[i]]), cex.axis=1.5, las = 2)
    if(i == "VSN") mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
    abline(h = max(colSums(normList[[i]], na.rm = T)), lty = 2)
  }
}

plotPCA = function(normList, meta, method, col){
  data = normList[[method]]
  data = data[!apply(is.na(data), 1, any),]
  pca = stats::prcomp(t(data))
  if(meta$Custom.Sample.Name[1] == ""){
    rownames(meta) = meta$Protein.Sample.Names
  } else {
    rownames(meta) = meta$Custom.Sample.Names
  }
  ggplot2::autoplot(pca, data = meta, colour = col, x = 1, y = 2, label  = TRUE)
}


plotPCV = function(normList, meta){
  groups <- meta$Group
  batch <- meta$Batch
  
  par(mar=c(10,6,3,1))
  plotData = lapply(normList, function(x) PCV(x, groups))
  boxplot(plotData, main = "PCV", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
  axis(2,cex.axis=1.5, las = 2)
  axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
  mtext(side = 2, text = "Pooled Coefficient of Variation", line = 4.5, cex = 1.5)
  points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
}

plotPMAD = function(normList, meta){
  groups <- meta$Group
  batch <- meta$Batch
  
  par(mar=c(10,6,3,1))
  plotData = lapply(normList, function(x) PMAD(x, groups))
  boxplot(plotData, main = "PMAD", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
  axis(2,cex.axis=1.5, las = 2)
  axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
  mtext(side = 2, text = "Median Absolute Deviation", line = 4.5, cex = 1.5)
  points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
}


plotPEV = function(normList, meta){
  groups <- meta$Group
  batch <- meta$Batch
  
  par(mar=c(10,6,3,1))
  plotData = lapply(normList, function(x) PEV(x, groups))
  boxplot(plotData, main = "PEV", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
  axis(2,cex.axis=1.5, las = 2)
  axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
  mtext(side = 2, text = "Pooled Estimate of Variance", line = 4.5, cex = 1.5)
  points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
}


plotCOR = function(normList, meta){
  groups <- meta$Group
  batch <- meta$Batch
  
  par(mar=c(10,6,3,1))
  plotData = lapply(normList, function(x) unlist(COR(x, groups)))
  boxplot(plotData, main = "Cor", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
  axis(side = 2,cex.axis=1.5, las = 2)
  axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
  mtext(side = 2, text = "Intragroup Correlation", line = 4.5, cex = 1.5)
  points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)
}


plotNaHM = function(normList, meta, show){
  groups <- meta$Group
  batch <- meta$Batch
  
  if(meta$Custom.Sample.Name[1] == ""){
    sampleLabels = meta$Protein.Sample.Names
  } else {
    sampleLabels = meta$Custom.Sample.Names
  }
  heatmapMissing(normList[["Log2"]], groups, batch, sampleLabels, show)
}


plotCorHM = function(normList, meta, method){
  groups <- meta$Group
  batch <- meta$Batch
  
  if(meta$Custom.Sample.Name[1] == ""){
    sampleLabels = meta$Protein.Sample.Names
  } else {
    sampleLabels = meta$Custom.Sample.Names
  }
  
  heatmapCorr(normList[[method]], groups, batch, sampleLabels)
}


plotLogRatio = function(normList, meta){
  groups <- meta$Group
  densityLog2Ratio(normList, groups)
}









