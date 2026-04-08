library(randomForest)
library(clue)
library(ggplot2)
library(reshape2)
library(ggalluvial)
library(ggdendro)
library(networkD3)
library(htmlwidgets)
library(keras)
library(nnet)
library(reticulate)
library(cluster)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(gridExtra)
library(ggrepel)

options(scipen=0)

setwd("./project/Analysis/Hi-C_differential_contacts/contacts_file/200bp/max_4MB/")

merge_chr <- function(search_path=search_path,search_pattern=search_pattern){
  
  # list the files 
  files <- list.files(path = search_path,
                      pattern = search_pattern,
                      full.names = TRUE)
  # adapt the code if needed for finding pattern of particular tyoe in case multiple replicates within same folder 
  
  data <- lapply(files, function(x) read.table(x,header = FALSE,stringsAsFactors = FALSE))
  dd <- do.call(rbind, data)
  df <- aggregate(V3 ~ V2 + V1,
                  data=dd,
                  FUN=sum)
  df <- as.data.frame(df)
  colnames(df) <- c("Source","Target","Counts")
  return(df)
}

nbc1 <- merge_chr(search_path = "NBC_rep1/",search_pattern = "^NBC_rep1.*\\.txt$")
mbc1 <- merge_chr(search_path = "MBC_rep1/",search_pattern = "^MBC_rep1.*\\.txt$")
gcbc2 <- merge_chr(search_path = "GCBC_rep2/",search_pattern = "^GCBC_rep2.*\\.txt$")
pc1 <- merge_chr(search_path = "PC_rep1/",search_pattern = "^PC_rep1.*\\.txt$")

#mat should be nbc1, ...
plot_circos <- function(mat, title_text, order = c("Tx_Wk", "Tx_Str", "Pr", "Enh", "Enh_Wk","Mix", "ReprPC", "Het_F", "Het_C", "NS")) {
  circos.clear()
  par(xpd = NA, mar = c(1, 1, 5, 1)) 
  
  circos.par(
    start.degree = 90,
    gap.degree = 4,
    track.margin = c(0.01, 0.01),
    points.overflow.warning = FALSE
  )
  
  chordDiagram(
    mat,
    grid.col = all_colors,
    transparency = 0.2,
    annotationTrack = "grid", 
    preAllocateTracks = list(track.height = uh(2, "mm")),  # ADJUST THE 2 TO 8 IF LABELS SHOULD GO IN!
    order = order
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      
      sector_width <- diff(xlim)
      
      ticks <- seq(1000000, 40000000, by = 2500000)
      #showing only ticks in field
      ticks <- ticks[ticks >= xlim[1] & ticks <= xlim[2]]
      print(ticks)
      for (t in ticks) {
        circos.lines(c(t, t), c(ylim[1], ylim[1] - mm_y(6)), lwd = 1.2)
        circos.text(
          x = t,
          y = ylim[1] + mm_y(5),
          labels = paste0(t / 1000000, "MB"),
          facing = "inside",
          cex = 1.3,
          adj = c(0.5, 1)
        )
      }
    },
    bg.border = NA
  )
  title(title_text, line = .5, cex.main = 3)
}