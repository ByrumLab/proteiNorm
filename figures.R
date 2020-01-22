load("savedData.Rdata")

par(mfrow = c(3,3))
i = "Cyclic Loess"
groups <- meta$Group
batch <- meta$Batch
if(meta$Custom.Sample.Name[1] == ""){
  sampleLabels = meta$Peptide.Sample.Names
} else {
  sampleLabels = meta$Custom.Sample.Names
}


par(mar=c(12,8,3,0))
barplot(colSums(normList[[i]], na.rm = T), main = i, las = 2, yaxt="n", cex.main = 1.5, col = plasma(ncol(normList[[i]])), names.arg = sampleLabels)
axis(side = 2, cex.axis=1.5, las = 2)
# axis(side = 1, at = seq_along(colnames(normList[[i]])), labels = colnames(normList[[i]]), cex.axis=1.5, las = 2)
mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
abline(h = max(colSums(normList[[i]], na.rm = T)), lty = 2)


data = normList[[i]]
data = data[!apply(is.na(data), 1, any),]
pca = stats::prcomp(t(data))
if(meta$Custom.Sample.Name[1] == ""){
  rownames(meta) = meta$Peptide.Sample.Names
} else {
  rownames(meta) = meta$Custom.Sample.Names
}
ggplot2::autoplot(pca, data = meta, colour = "Group", x = 1, y = 2, label  = TRUE)



par(mar=c(10,6,3,1))
plotData = lapply(normList, function(x) PCV(x, groups))
boxplot(plotData, main = "PCV", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
axis(2,cex.axis=1.5, las = 2)
axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
mtext(side = 2, text = "Pooled Coefficient of Variation", line = 4.5, cex = 1.5)
points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)


par(mar=c(10,6,3,1))
plotData = lapply(normList, function(x) PMAD(x, groups))
boxplot(plotData, main = "PMAD", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
axis(2,cex.axis=1.5, las = 2)
axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
mtext(side = 2, text = "Median Absolute Deviation", line = 4.5, cex = 1.5)
points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)

par(mar=c(10,6,3,1))
plotData = lapply(normList, function(x) PEV(x, groups))
boxplot(plotData, main = "PEV", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
axis(2,cex.axis=1.5, las = 2)
axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
mtext(side = 2, text = "Pooled Estimate of Variance", line = 4.5, cex = 1.5)
points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)


par(mar=c(10,6,3,1))
plotData = lapply(normList, function(x) unlist(COR(x, groups)))
boxplot(plotData, main = "Cor", las = 2, border = rainbow(length(normList)), boxlwd = 2, yaxt="n", xaxt="n", cex.main = 1.5)
axis(side = 2,cex.axis=1.5, las = 2)
axis(side = 1, at = seq_along(names(normList)), labels = names(normList), cex.axis=1.5, las = 2)
mtext(side = 2, text = "Intragroup Correlation", line = 4.5, cex = 1.5)
points(rep(seq_along(normList), each = length(plotData[[1]])), unlist(plotData), pch = "*", cex = 1.3)


data = normList[["Log2"]]
showAllProtein = FALSE
missing = !is.na(data)
if(!showAllProtein){
  complete = apply(missing, 1, all)
  completeNA = apply(!missing, 1, all)
  missing = missing[!complete & !completeNA,]
}
batchCol =  colorBatch(batch)
groupCol = colorGroup(groups)
hm_clust = heatmapClustered(groups, batch, missing, groupCol, batchCol, sampleLabels)
draw(hm_clust, heatmap_legend_side = "right", annotation_legend_side = "right", ht_gap = unit(2, "cm"), column_title = "Missing values")


heatmapCorr(normList[[i]], groups, batch, sampleLabels)

densityLog2Ratio(normList, groups)


