library(ggplot2)
library(reshape2)
library(ggpubr)
source("https://raw.githubusercontent.com/VanAndelInstitute/TenXTools/master/TenX_util.R")


plot_row_height <- 250

setwd("/mnt/data/scratch/christy/AGG/outs")

filter <- function(data, symbol, threshold) {
  ix <- which(fData(data)$symbol == symbol)[1]
  ix <- which(exprs(data)[ix,] > threshold)
  return(data[,ix])
}

dat <- readRDS("cycling_10x.rds")
#gbm <- load_cellranger_matrix("../")

ix <- which(fData(dat)$symbol == "TMEM88")[1]
ix2 <- which(fData(dat)$symbol == "TNNT2")[1]
dat <- dat[, which(exprs(dat)[ix,] > 0 | exprs(dat)[ix2,] > 0)]

cell_cycle <- c("CCNA2", "CCNB1", "CCNB2", "CDK1", "CDK4", "WEE1", "CDC20", "CDC5L", "CDK14", "CDC14B", "CDC14A", "CCND1", "CDH1")
cm <- c("TNNT2", "ISL1", "TMEM88", "MYH6", "NKX2-5", "ACTN2", "SIRPA")
mitosis <- c("ANAPC10", "ANLN", "KIF11", "BUB3", "SPG20", "SPAST", "AURKB")
apoptosis <- c("TP53", "CDKN1A", "BAX")
proliferation <- c("MKI67", "E2F1", "E2F3")

plot_genes <- function(gg, data, class, title = "Gene counts by class") {
  ix <- which(fData(data)$symbol %in% gg)
  slice <- exprs(data)[ix,]
  rownames(slice) <- fData(data)$symbol[ix]
  slice <- as.data.frame(t(as.matrix(slice)))
  slice$class <- class
  slice$class <- factor(slice$class, labels = c("G1-CDC20-", "G2/M-CDC20-", "G2/M-CDC20+"))
  slice_m <- melt(slice, id.vars = c("class"))
  dummy <- slice_m
  dummy$value <- 2 * dummy$value
  comparisons <- list( c("G2/M-CDC20+", "G2/M-CDC20-"),  c("G1-CDC20-", "G2/M-CDC20-"), c("G1-CDC20-", "G2/M-CDC20+"))
  p <- ggplot(slice_m, aes(class, value)) +
    geom_jitter(aes(color=class), height=0.01, width=0.1, alpha=0.5, shape=20, size=1) +
    geom_violin(fill="orange", alpha=0.5) +
    theme_bw() + 
    stat_compare_means(method = "wilcox.test", comparisons = comparisons, textsize=3, label="p.signif", hide.ns=TRUE) +
    scale_y_continuous(expand = c(0.2, 0)) +
    ggtitle(title) +
    facet_wrap(~variable, ncol=3, scales="free")
  print(p)
}

violins <- function() {
  png(type = "cairo", 
      filename = "cell_cycle.png", 
      width = 1200, 
      height = plot_row_height * ceiling(length(cell_cycle)/3), 
      units = "px",
      res = 100)
  plot_genes(cell_cycle, dat, pData(dat)$type_cdc20, "Cell cycle associated genes")
  dev.off()
  
  png(type = "cairo", 
      filename = "cardiomyocyte.png", 
      width = 1200, 
      height = plot_row_height * ceiling(length(cm)/3), 
      units = "px",
      res = 100)
  plot_genes(cm, dat, pData(dat)$type_cdc20, "Cardiomyocyte associated genes")
  dev.off()
  
  png(type = "cairo", 
      filename = "mitosis.png", 
      width = 1200, 
      height = plot_row_height * ceiling(length(mitosis)/3), 
      units = "px",
      res = 100)
  plot_genes(mitosis, dat, pData(dat)$type_cdc20, "Mitosis associated genes")
  dev.off()
  
  png(type = "cairo", 
      filename = "apoptosis.png", 
      width = 1200, 
      height = plot_row_height * ceiling(length(apoptosis)/3), 
      units = "px",
      res = 100)
  plot_genes(apoptosis, dat, pData(dat)$type_cdc20, "Apoptosis associated genes")
  dev.off()
  
  png(type = "cairo", 
      filename = "proliferation.png", 
      width = 1200, 
      height = plot_row_height * ceiling(length(proliferation)/3), 
      units = "px",
      res = 100)
  plot_genes(proliferation, dat, pData(dat)$type_cdc20, "Proliferation associated genes")
  dev.off()
  
}

tsne_plot <- function() {
  counts <- tenXsym(gbm, c("TMEM88", "CDC20"))
  analysis_results <- load_cellranger_analysis_results("./")
  tsne_proj <- analysis_results$tsne
  tenXoverlay(tsne_proj[c("TSNE.1","TSNE.2")],
              size = counts$CDC20, 
              color = pd$type_cdc20, 
              title = "TSNE Plot with CDC20 expression",
              marker_range = c(2,6))
}

plot_data <- function(data, class, title = "Data by class") {
  df <- data.frame(value=data, class=class)
  df$class <- factor(df$class, labels = c("G1-CDC20-", "G2/M-CDC20-", "G2/M-CDC20+"))
  slice_m <- melt(df, id.vars = c("class"))
  comparisons <- list( c("G2/M-CDC20+", "G2/M-CDC20-"),  c("G1-CDC20-", "G2/M-CDC20-"), c("G1-CDC20-", "G2/M-CDC20+"))
  p <- ggplot(slice_m, aes(class, value)) +
    geom_jitter(aes(color=class), height=0.01, width=0.1, alpha=0.5, shape=20, size=1) +
    geom_violin(fill="orange", alpha=0.5) +
    theme_bw() + 
    stat_compare_means(method = "t.test", comparisons = comparisons, textsize=3, label="p.signif", hide.ns=TRUE) +
    scale_y_continuous(expand = c(0.2, 0)) +
    ggtitle(title) +
    facet_wrap(~variable, ncol=3, scales="free")
  print(p)
}


auc_plots <- function(sig, file) {
  dat <- readRDS("cycling_10x.rds")
  ix <- which(fData(dat)$symbol == "TMEM88")[1]
  ix2 <- which(fData(dat)$symbol == "TNNT2")[1]
  ix <- which(exprs(dat)[ix,] > 0 | exprs(dat)[ix2,] > 0)
  dat <- dat[, ix]
  cells_AUC <- AUCell_calcAUC(txt2sig(sig), cells_rankings)
  png(type = "cairo", 
      filename = file, 
      width = 500, 
      height = 400, 
      units = "px",
      res = 100)
  plot_data(getAUC(cells_AUC)[1,], pData(dat)$type_cdc20, names(cells_AUC)[1])
  dev.off()
  
}

txt2sig <- function(file) {
  data <- readLines(file)
  res <- list()
  res[[data[1]]] <- c(data[-c(1:2)])
  ix <- which(res[[data[1]]] %in% ls(revmap(org.Hs.egSYMBOL)))
  res[[data[1]]] <- unlist(mget(res[[data[1]]][ix], revmap(org.Hs.egSYMBOL)))
  res[[data[1]]] <- unlist(mget(res[[data[1]]], org.Hs.egENSEMBL))
  res
}

