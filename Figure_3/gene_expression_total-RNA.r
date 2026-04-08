# Load necessary library
# set the environment
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "keras","nnet","reticulate","cluster","readr","caret","uwot","gprofiler2",
                      "ComplexHeatmap","circlize","dplyr","readr","gridExtra","ggrepel","scales","ggfortify","ggpubr","ggsignif")
# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)
options(scipen=999)
set.seed(123)


# set the working directory 
setwd("./project/")
metadata_ihec ="./project/metadata/IHEC_metadata_harmonization.v1.2_with_coverage_avg_meth.extended.csv"
metadata <- read.csv(metadata_ihec,header=T)

json_color ="./project/metadata/IHEC_EpiATLAS_IA_colors_Mar18_2024.json"

json_data <- fromJSON(json_color,flatten=TRUE)
# convert the rgb color codes to hex codes 
colors.harmonized_sample_ontology_intermediate <- c(json_data$harmonized_sample_ontology_intermediate[[4]])
colors.hex.harmonized_sample_ontology_intermediate <- sapply(colors.harmonized_sample_ontology_intermediate,function(x) {
  rgb_vals <- as.numeric(unlist(strsplit(x, ",")))
  rgb_str <- sprintf("#%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  return(rgb_str)
}) 
# Set the path to the directory containing the CSV files
folder_path <- "./scratch/gene_expression_predictions/total-RNA_200bp/"

# List all CSV files in the directory
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
list_of_data_frames <- lapply(csv_files, read.csv)


combined_data_frame <- bind_rows(list_of_data_frames)
combined_data_frame$Model <- factor(combined_data_frame$Model ,levels = c("EpiSegMixMeth","EpiSegMix","ChromHMM")) 
metadata2 <- read.csv("./project/metadata/epiatlas_chipseq_qc_summary.csv",header=TRUE)

df <- combined_data_frame
df <- merge(df,metadata,by.x="EpiRR",by.y="epirr_id_without_version")
metadata2$samples <- sapply(metadata2$epirr_id, function(x) strsplit(x,".",fixed = TRUE)[[1]][1])
h3k27me3 <- metadata2[which(metadata2$antibody=="H3K27me3"),c(91:94)]
h3k9me3 <- metadata2[which(metadata2$antibody=="H3K9me3"),c(91:94)]
df <- merge(df,h3k9me3,by.y="samples",by.x="EpiRR")

dd <- aggregate(lm_r2~Model + harmonized_sample_ontology_intermediate,data=df,median,na.rm=TRUE)

p <- ggplot(df, aes(x=Model, y=lm_r2))  + geom_point(aes(color=harmonized_sample_ontology_intermediate)) + 
  geom_boxplot() + xlab("Models") + ylab( bquote(R^2)) + 
  geom_path(aes(group = EpiRR,color=harmonized_sample_ontology_intermediate),alpha=0.7) +  
  theme_minimal() +  scale_color_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + 
  scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate)   + 
  facet_wrap(~harmonized_sample_ontology_intermediate) 
  
ggplot(dd, aes(x=Model, y=lm_r2))  + geom_point(aes(color=harmonized_sample_ontology_intermediate)) + 
  xlab("Models") + ylab( bquote(R^2)) + 
  geom_path(aes(group = harmonized_sample_ontology_intermediate,color=harmonized_sample_ontology_intermediate),alpha=0.7) +  
  theme_minimal() +  scale_color_manual(values=colors.hex.harmonized_sample_ontology_intermediate) 


p <- p + theme(axis.text = element_text(size = 20)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 20)) # change axis titles

p <- p + theme(text = element_text(size = 20))    
p  
options(scipen=0)
p + geom_boxplot() 
ggplot(df, aes(x=harmonized_sample_ontology_intermediate, y=lm_r2,fill=Model)) + geom_boxplot() + theme_minimal()


get_box_stats <- function(y, upper_limit = max(df$lm_r2) * 1.15) {
  return(data.frame(
    y = 0.93 * upper_limit,
    label = paste(
      "Median", "\n",round(median(y), 3), "\n"
    )
  ))
}

# check the comparison for p-values
compare_means(lm_r2 ~ Model,  data = df)

# my_comparisons <- list(c("EpiSegMixMeth ", "ChromHMM"), c("EpiSegMix", "ChromHMM"), c("EpiSegMixMeth", "EpiSegMix"))

p <- ggplot(df, aes(x=Model, y=lm_r2,fill=Model))+ 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1,binwidth = 1/159) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9,size=6) + 
  theme_minimal() +   xlab("Models") + ylab( bquote(R^2)) #+ scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate)
p <- p + geom_signif(y_position = c(.66, .64,.62),
                     comparisons = list(c("EpiSegMixMeth", "ChromHMM"), c("EpiSegMix", "ChromHMM"),c("EpiSegMixMeth","EpiSegMix")),
                     map_signif_level = function(x) paste("p =", scales::pvalue(x)))


p <- p + theme(axis.text = element_text(size = 16,color = "black")) # changes axis labels

p <- p + theme(axis.title = element_text(size = 16,color = "black")) # change axis titles

p <- p + theme(text = element_text(size = 16,color = "black"))    
p <- p + guides(fill="none")

p

ggsave(filename = "./scratch/gene_expression_predictions/total-RNA_200bp/total-RNA_predictions_200bp.png",plot=p,height = 6, width=7, dpi=300)
ggsave(filename = "./scratch/gene_expression_predictions/total-RNA_200bp/total-RNA_predictions_200bp.pdf",plot=p,height = 8, width=9)
# supplementory slide

dd <- df 
dd$harmonized_sample_ontology_intermediate[which(df$harmonized_sample_ontology_intermediate == "melanocyte" | 
                                                   df$harmonized_sample_ontology_intermediate == "mucosa" |
                                                   df$harmonized_sample_ontology_intermediate == "eosinophil" |
                                                   df$harmonized_sample_ontology_intermediate == "myeloid cell" )] <- "Others"
unique(df$harmonized_sample_ontology_intermediate)
ggplot(dd, aes(x = Model, y = lm_r2, fill = Model)) +
  theme_minimal() +
  labs(title = "Gene expression prediction", x = "Model", y = "Value") + 
  geom_path(aes(group = EpiRR,color=harmonized_sample_ontology_intermediate),alpha=0.7) +
  scale_color_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + 
  facet_wrap(~harmonized_sample_ontology_intermediate)

