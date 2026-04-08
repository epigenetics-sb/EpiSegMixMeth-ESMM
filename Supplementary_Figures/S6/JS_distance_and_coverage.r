
# set the environment
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "keras","nnet","reticulate","cluster",
                      "ComplexHeatmap","circlize","dplyr","readr")

# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

set.seed(123)
# set the working directory
setwd("./project/")


chmm <- read.table("./scratch/gene_expression_predictions/ChromHMM_total-RNA-predictions.v2.txt",header = F)
epi <- read.table("./scratch/gene_expression_predictions/EpiSegMix_total-RNA-predictions.v2.txt",header=F)
epim <- read.table("./scratch/gene_expression_predictions/EpiSegMixMeth_total-RNA-predictions.v2.txt",header=F)
test <- read.table("Analysis/genomic_coverage/EpiSegMixMeth_EpiSegMix_Het-F_qc.tab",header = FALSE)
metadata <- read.csv("./project/metadata/IHEC_metadata_harmonization.v1.2_with_coverage_avg_meth.extended.10x_filtered.sorted.csv")
metadata2 <- read.csv("./project/metadata/epiatlas_chipseq_qc_summary.csv",header=TRUE)

unique(metadata$harmonized_sample_ontology_intermediate)
metadata$epirr_id_without_version <- factor(metadata$epirr_id_without_version,levels=c(metadata$epirr_id_without_version))

metadata2$epirr_id_without_version <- sapply(metadata2$epirr_id, function(x) { strsplit(x,'.',fixed=TRUE)[[1]][1]})



pruned_metdata <- metadata2[which(metadata2$antibody =="H3K27me3"),c(3,4,5,7,92,93,94)]



colnames(test) <- c("epirr","EpiSegMixMeth","EpiSegMix","jsd_H3K9me3","jsd_H3K27me3")

test$jsd_H3K9me3 <- test$jsd_H3K9me3*100
test$jsd_H3K27me3 <- test$jsd_H3K27me3*100

temp <- merge(test,metadata,by.x = "epirr",by.y = "epirr_id_without_version")

# ordered_temp <- temp %>% 
#   arrange(jsd_H3K27me3,jsd_H3K9me3,flagstat_qc.total,flagstat_qc.mapped_pct) %>% 
#   mutate(sample_ordered = factor(epirr, levels = unique(epirr)))

unique(melt_df$harmonized_sample_ontology_intermediate)

temp$epirr = with(temp, reorder(epirr, jsd_H3K27me3))

df <- melt(data = temp[,c(1,2,3,4,5)],id.vars = c("epirr", "jsd_H3K27me3","jsd_H3K9me3"),measure.vars =  c("EpiSegMixMeth","EpiSegMix"))


epirr_colors <- list(
  IHECRE00004617 = "pink",
  IHECRE00004616 = "black",
  IHECRE00004614 = "red",
  IHECRE00004607 = "blue",
  IHECRE00004606 = "green",
  IHECRE00004603 = "orange",
  IHECRE00000231 = "darkgreen",
  IHECRE00000223 = "yellow"
)

d_plot <- ggplot(df, aes(x = jsd_H3K27me3, y = value)) +
  geom_point(size=4,color="grey", aes(color=variable)) +  # Assuming 'color' is intended to show 'variable'
  labs(color = "Model") +
  facet_wrap(~variable) +
  geom_smooth(method = "lm", aes(x = jsd_H3K27me3, y = value), color = "blue") +  # Fit line for 'jsd_H3K27me3'
  theme_minimal()  + 
  geom_text_repel(
    data= subset(df, epirr %in% c("IHECRE00004617", "IHECRE00004616", "IHECRE00004614", "IHECRE00004607",
                                  "IHECRE00004606", "IHECRE00004603", "IHECRE00000231", "IHECRE00000223")),
    aes(label=sub("IHECRE0000", "", paste("(",round(jsd_H3K27me3,0),",",round(value,0),")",sep = ""))),
    max.overlaps = 8,  # Adjusted for potential more overlaps
    # point.padding = 0.5,  # Increased padding
    nudge_x = -0.5,  # Optional: nudge labels to the right
    nudge_y = -0.5   # Optional: nudge labels upward
  ) + 
  geom_point(data=subset(df, epirr %in% c("IHECRE00004617", "IHECRE00004616", "IHECRE00004614", "IHECRE00004607",
                                     "IHECRE00004606", "IHECRE00004603", "IHECRE00000231", "IHECRE00000223")),
             
             aes(x = jsd_H3K27me3, y = value,color=epirr),size=4) + scale_color_manual(values=epirr_colors) 

p <- d_plot 
#p <- ggplot(test, aes(x=V2, y=V3,fill=V2)) +  geom_boxplot() + theme_minimal() 

p <- p +   labs(y="Het_F % ",x="H3K27me3 Quality \n JS Distance", fill="Model",color="EpiRR") 

p <- p + theme(axis.text = element_text(size = 20)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 20)) # change axis titles

p <- p + theme(text = element_text(size = 20)) #+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

p
# p +  theme(axis.text.x=element_blank(),
#            axis.ticks.x=element_blank())
ggsave(filename = "Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/luminal_epithelial_cell_of_mammary_gland/js_coverage.v4.png",plot=p,height = 9, width=14, dpi=200)
ggsave(filename = "Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/luminal_epithelial_cell_of_mammary_gland/js_coverage.v4.pdf",plot=p,height = 9, width=14)
