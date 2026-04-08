# set the environment
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "keras","nnet","reticulate","cluster",
                      "ComplexHeatmap","circlize","dplyr","readr","gridExtra","ggrepel","jsonlite")
# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

set.seed(123)
# set the working directory
setwd("./project/")

epirr_order <- read.table("Analysis/pygenome/159_samples_order.txt",header=FALSE)
epirr_order$id <- c(1:nrow(epirr_order))

annotation <- read.csv("metadata/IHEC_metadata_harmonization.v1.1.csv",header = TRUE)

json_data <- fromJSON("metadata/IHEC_EpiATLAS_IA_colors_Mar18_2024.json",flatten=TRUE)

df <- merge(epirr_order,annotation,by.x = "V1",by.y="epirr_id_without_version",all.x = TRUE)
df.resorted <- df[order(df$id), ]



colors.cell_type <- c(json_data$harmonized_sample_ontology_intermediate[[4]])

colors.hex.cell_type <- sapply(colors.cell_type,function(x) {
  rgb_vals <- as.numeric(unlist(strsplit(x, ",")))
  rgb_str <- sprintf("#%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  return(rgb_str)
}) 



row_ha <- rowAnnotation(cell_types= df.resorted$harmonized_sample_ontology_intermediate,
                       col = list(cell_types = colors.hex.cell_type),
                       annotation_name_gp = gpar(fontsize=20),
                       annotation_legend_param = list(
                         cell_types=list(
                           title_gp=gpar(fontsize=20,fontface="bold"),
                           labels_gp = gpar(fontsize = 20))
                       ))

png(file = "Analysis/Figures/Figure0/154_samples_order.png",width = 0.5,height = 40,res = 200,units = "cm")
draw(row_ha)
dev.off()


png(file = "Analysis/Figures/Figure0/154_samples_order_legends.png",width = 25,height = 40,res = 200,units = "cm")
rowAnnotation(cell = 1:154) + rowAnnotation(cell_types= df.resorted$harmonized_sample_ontology_intermediate,
                                           col = list(cell_types = colors.hex.cell_type),
                                           annotation_name_gp = gpar(fontsize=20),
                                           annotation_legend_param = list(
                                             cell_types=list(
                                               title_gp=gpar(fontsize=20,fontface="bold"),
                                               labels_gp = gpar(fontsize = 20))
                                           ))


dev.off()

#ggsave(filename = "Analysis/sankey_plot_Bcell_full_reduced/NBC_I_ChromHMM_Sankey.png",plot =gplot,
#       dpi = 200,
#       width = 10,height = 10)
