# list of packages to be loaded
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "kohonen","keras","nnet","reticulate","cluster",
                      "ComplexHeatmap","circlize","dendsort")

# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

# set the working directory as well as the environment
setwd("D:/Work/Epigenetics/IHEC/EpiSegMix/")
options(scipen = 999)

#  Mode selection 
# Set mode to one of: "CHMM" (ChromHMM), "ESM" (EpiSegMix), "ESMM" (EpiSegMixMeth)
mode <- "ESM"

# Mode-specific parameters 
if (mode == "CHMM") {
  feature_file     <- "Classifier/IHEC_ChromHMM_10N_all_states_samples_classifier.csv"
  annotation_file  <- "IHEC_Companion/Random_forest_classifier/2024_09_10_Annotations_improved/ChromHMM_annotation_without_features.csv"
  output_file      <- "IHEC_Companion/Figure2/ChromHMM_300dpi.png"
  plot_title       <- "ChromHMM"
  annotation_levels <- c("Tx_Str","Tx_Wk","Tx","Pr","Pr_Wk","Enh",
                          "Enh_Wk","ReprPC","Het_F","Het_C","Rep","NS")
  dna_methyl_label <- "Annotation"   # not a training feature in CHMM

} else if (mode == "ESM") {
  feature_file     <- "Classifier/IHEC_EpiSegMix_10N_all_states_samples_classifier.csv"
  annotation_file  <- "IHEC_Companion/Random_forest_classifier/2024_09_10_Annotations_improved/EpiSegMix_annotation_without_features.csv"
  output_file      <- "IHEC_Companion/Figure2/EpiSegMix_300dpi.png"
  plot_title       <- "EpiSegMix"
  annotation_levels <- c("Tx_Str","Tx_Wk","Pr","Enh","Enh_Wk",
                          "ReprPC","Mix","Het_F","Het_C","NS")
  dna_methyl_label <- "Annotation"   # not a training feature in ESM

} else if (mode == "ESMM") {
  feature_file     <- "Classifier/IHEC_EpiSegMix_meth_10N_all_states_samples_classifier.csv"
  annotation_file  <- "IHEC_Companion/Random_forest_classifier/2024_09_10_Annotations_improved/EpiSegMixMeth_annotation_without_features.csv"
  output_file      <- "IHEC_Companion/Figure2/EpiSegMixMeth_300dpi.png"
  plot_title       <- "EpiSegMixMeth"
  annotation_levels <- c("Tx_Str","Tx_Wk","Pr","Enh","Enh_Wk",
                          "ReprPC","Mix","Het_F","Het_C","NS")
  dna_methyl_label <- "TrainFeature" # DNA methylation is a training feature in ESMM

} else {
  stop(paste("Unknown mode:", mode, "-- must be one of: CHMM, ESM, ESMM"))
}

# Read the feature file and do min-max normalization 
df <- read.csv(feature_file, header = TRUE)

cols <- c(3:ncol(df))
epirrids <- unique(df$EpiRR)
for (epirrid in epirrids) {
  for (x in cols) {
    min_val <- min(df[df$EpiRR == epirrid, x], na.rm = TRUE)
    max_val <- max(df[df$EpiRR == epirrid, x], na.rm = TRUE)
    df[df$EpiRR == epirrid, x] <- (df[df$EpiRR == epirrid, x] - min_val) / (max_val - min_val)
  }
}

df$id <- paste0(df$EpiRR, "_", df$states)

# Read annotation file and merge 
annotation_df <- read.csv(annotation_file, header = TRUE)
annotation_df$id <- paste0(annotation_df$EpiRR, "_", annotation_df$states)
annotation_df <- subset(annotation_df, select = -c(EpiRR, states, annotated_probs, corrected_probs, annotated))

df <- merge(df, annotation_df, by = "id")
rownames(df) <- paste0(df$EpiRR, "_", df$states)
df <- df[order(row.names(df)), ]
df$id <- NULL

# Rename certain columns 
names(df)[3] <- "Coverage"
names(df)[8] <- "GeneExpression"

# ChromHMM-only: remap "Mix" → "Enh_Wk" 
if (mode == "CHMM") {
  df$corrected[df$corrected == "Mix"] <- "Enh_Wk"
}

# Column order and feature-group labels 
custom_order <- c("H3K4me3","H3K27ac","H3K4me1",
                  "H3K36me3","H3K27me3","H3K9me3",
                  "DNA_methyl","Coverage","CpG",
                  "GeneExpression","RefSeqGene",
                  "RefSeqTSS","RefSeqTSS2kb",
                  "RefSeqExon","RefSeqTES",
                  "PMDs_Common","SINE","RNA",
                  "Retroposon",
                  "LINE","LTR",
                  "RC","TandemRepeats")

custom_order_labels <- c(
  "TrainFeature","TrainFeature","TrainFeature",
  "TrainFeature","TrainFeature","TrainFeature",
  dna_methyl_label,                              # mode-specific (index 7)
  "Annotation","Annotation",
  "Annotation","Annotation",
  "Annotation","Annotation",
  "Annotation","Annotation",
  "Annotation","RepeatAnnotation","RepeatAnnotation",
  "RepeatAnnotation",
  "RepeatAnnotation","RepeatAnnotation",
  "RepeatAnnotation","RepeatAnnotation")

colors_custom_order_labels <- c("Yellow","Green","Red")
names(colors_custom_order_labels) <- unique(custom_order_labels)

# Annotation factor levels (mode-specific) 
df$annotation <- factor(df$corrected, levels = annotation_levels)

#  Color palette 
col_fun <- colorRamp2(c(0, 1), c("white","#3F7F93"))

set.seed(123)
all_states <- c("NS","Mix","Het_C","Het_F","ReprPC","Enh_Wk","Enh",
                "Pr","Tx_Wk","Tx_Str","Tx","Pr_Wk","Rep")
all_colors <- c("#EAEAEA","#8A91D0","#92C5DE","#11C1FF","#808080",
                "#FFFF00","#FFC34D","#FF0000","#006400","#008000",
                "#006400","#FF4500","#66CDAA")
names(all_colors) <- all_states

# Heatmap annotations 
text_fontsize <- 11

col_ha <- HeatmapAnnotation(
  Features = custom_order_labels,
  col = list(Features = colors_custom_order_labels),
  annotation_name_gp = gpar(fontsize = text_fontsize),
  annotation_legend_param = list(
    Features = list(grid_height = unit(1, "mm"), grid_width = unit(5, "mm"),
                    title_gp = gpar(fontsize = text_fontsize, fontface = "bold"),
                    legend_height = 2,
                    labels_gp = gpar(fontsize = text_fontsize))
  ))

row_ha <- rowAnnotation(
  States = df$annotation,
  col = list(States = all_colors),
  annotation_name_gp = gpar(fontsize = text_fontsize),
  annotation_legend_param = list(
    States = list(grid_height = unit(0.35, "mm"), grid_width = unit(5, "mm"),
                  title_gp = gpar(fontsize = text_fontsize, fontface = "bold"),
                  labels_gp = gpar(fontsize = text_fontsize))
  ))

#  Main heatmap (annotated, ordered) 
png(output_file, width = 150, height = 150, res = 300, units = "mm")
features2 <- as.matrix(df[, custom_order])

ht_opt(legend_gap = unit(1, "mm"))
a <- Heatmap(features2,
             show_row_names = FALSE,
             show_row_dend = FALSE,
             cluster_columns = FALSE, cluster_rows = FALSE,
             row_title_rot = 0,
             left_annotation = row_ha, bottom_annotation = col_ha,
             row_split = as.factor(df$annotation),
             col = col_fun,
             cluster_row_slices = FALSE, column_title_rot = 30,
             column_title_gp = gpar(fontsize = text_fontsize),
             column_names_gp = gpar(fontsize = text_fontsize),
             row_title_gp = gpar(fontsize = text_fontsize),
             row_names_gp = gpar(fontsize = text_fontsize),
             row_gap = unit(1, "mm"),
             heatmap_legend_param = list(
               grid_height = unit(5, "mm"), grid_width = unit(2, "mm"),
               title = "Min-Max", direction = "horizontal",
               title_gp = gpar(fontsize = text_fontsize, fontface = "bold"),
               labels_gp = gpar(fontsize = text_fontsize)))

draw(a, merge_legend = TRUE,
     column_title = plot_title,
     column_title_gp = gpar(fontsize = text_fontsize))
dev.off()


#  RFC input heatmap (clustered, no annotation) 
rfc_output_file <- sub("_300dpi\\.png$", "_RFC_input_300dpi.png", output_file)

png(rfc_output_file, width = 150, height = 150, res = 300, units = "mm")
features2 <- as.matrix(df[, custom_order])
text_fontsize <- 14

ht_opt(legend_gap = unit(1, "mm"))
a <- Heatmap(features2,
             show_row_names = FALSE,
             show_row_dend = FALSE, show_column_dend = FALSE,
             cluster_columns = TRUE, cluster_rows = TRUE,
             row_title_rot = 0,
             col = col_fun,
             cluster_row_slices = FALSE, column_title_rot = 45,
             column_title_gp = gpar(fontsize = text_fontsize),
             column_names_gp = gpar(fontsize = text_fontsize, rot = 45),
             row_title_gp = gpar(fontsize = text_fontsize),
             row_names_gp = gpar(fontsize = text_fontsize),
             row_gap = unit(1, "mm"),
             heatmap_legend_param = list(
               grid_height = unit(5, "mm"), grid_width = unit(2, "mm"),
               title = "Min-Max", direction = "horizontal",
               title_gp = gpar(fontsize = text_fontsize, fontface = "bold"),
               labels_gp = gpar(fontsize = text_fontsize)))

draw(a, merge_legend = TRUE,
     column_title = paste0(plot_title, " - RFC Input"),
     column_title_gp = gpar(fontsize = text_fontsize))
dev.off()
