logNorm <- function(dat) {
  loggedInt <- log2(dat)
  #loggedInt <- replace(is.infinite(loggedInt), NA)
  loggedInt[is.infinite(as.matrix(loggedInt))] <- NA
  return(as.matrix(loggedInt))
}

medianNorm <- function(loggedDat) {
  # Find medians of each sample
  # Divide by median
  # Multiply by mean of medians
  sampleMed <- apply(loggedDat, 2, median, na.rm=TRUE)
  meanMed <- mean(sampleMed, na.rm=TRUE)
  out <- t(t(loggedDat) / sampleMed)
  out <- out * meanMed
  return(as.matrix(out))
}

meanNorm <- function(loggedDat) {
  # Find means of each sample
  # Divide by mean
  # Multiply by mean of means
  sampleMean <- apply(loggedDat, 2, mean, na.rm=TRUE)
  meanMean <- mean(sampleMean, na.rm=TRUE)
  out <- t(t(loggedDat) / sampleMean)
  out <- out * meanMean
  return(as.matrix(out))
}

vsnNorm <- function(dat) {
  vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsnNormed) <- colnames(dat)
  row.names(vsnNormed) <- rownames(dat)
  return(as.matrix(vsnNormed))
}

quantNorm <- function(loggedDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(loggedDat), copy=FALSE)
  colnames(quantNormed) <- colnames(loggedDat)
  row.names(quantNormed) <- rownames(loggedDat)
  return(as.matrix(quantNormed))
}

cycLoessNorm <- function(loggedDat) {
  cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(loggedDat), method="fast")
  colnames(cycLoessNormed) <- colnames(loggedDat)
  row.names(cycLoessNormed) <- rownames(loggedDat)
  return(as.matrix(cycLoessNormed))
}

rlrNorm <- function(loggedDat) {
  rlrNormed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(loggedDat), noLogTransform=TRUE)
  colnames(rlrNormed) <- colnames(loggedDat)
  row.names(rlrNormed) <- rownames(loggedDat)
  return(as.matrix(rlrNormed))
}

giNorm <- function(loggedDat) {
  giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(loggedDat), noLogTransform=TRUE)
  colnames(giNormed) <- colnames(loggedDat)
  row.names(giNormed) <- rownames(loggedDat)
  return(as.matrix(giNormed))
}
