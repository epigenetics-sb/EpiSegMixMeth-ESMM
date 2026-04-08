# set the environment
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "keras","nnet","reticulate","cluster",
                      "ComplexHeatmap","circlize","dplyr","readr","gridExtra")
# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

set.seed(123)
# set the working directory
setwd("./project/")



# reading annotation file to relabel the states
annotation_list <- c("Analysis/classifier/Random_forest_classifier/ChromHMM_annotation_without_features.csv",
                     "Analysis/classifier/Random_forest_classifier/EpiSegMix_annotation_without_features.csv",
                     "Analysis/classifier/Random_forest_classifier/EpiSegMixMeth_annotation_without_features.csv")


# file contain the genomic coverage they were obtained from the overlap enrichment file generated for each epirr 
# during the classifier training data 
files_list <- c("Analysis/genomic_coverage/ChromHMM_state_coverage.txt",
                "Analysis/genomic_coverage/EpiSegMix_state_coverage.txt",
                "Analysis/genomic_coverage/EpiSegMixMeth_state_coverage.txt")




coverage <- function(annotation_path,df_path,identifier){
  # reading annotation file
  annotation <- read.csv(annotation_path,header = TRUE)
  # temporary id
  annotation$id <- paste(annotation$EpiRR,annotation$states,sep = "_") 
  
  # read data file
  df <- read.table(df_path,header = FALSE)
  colnames(df) <- c("epirr","states","coverage")
  # temporary id
  df$id <- paste(df$epirr,df$states,sep = "_")
  
  # merging the annotation and the genomic coverage based on the temporary id created
  df.merged <- merge(df,annotation,by = "id",all.x = TRUE)
  df <- df.merged[,c("coverage","corrected")]
  df$Model <- rep(identifier,nrow(df))
  return(df)
}

chmm_df <- coverage(annotation_path = annotation_list[1],df_path = files_list[1],identifier = "ChromHMM")
chmm_df <- rbind(chmm_df, data.frame(coverage = as.numeric(0), corrected = "Enh_Wk", Model = "ChromHMM"))
episegmix_df <- coverage(annotation_path = annotation_list[2], df_path = files_list[2], identifier = "EpiSegMix")
#episegmix_df <- rbind(episegmix_df, data.frame(coverage = as.numeric(0), corrected = "Rep", Model = "EpiSegMix"))

episegmixmeth_df <- coverage(annotation_path = annotation_list[3],df_path = files_list[3],identifier = "EpiSegMixMeth")

# color code for the states
all_states <- c("NS",
                "Mix",
                "Het_C",
                "Het_F",
                "ReprPC",
                "Enh_Wk",
                "Enh",
                "Pr",
                "Tx_Wk",
                "Tx_Str",
                "Tx","Pr_Wk","Rep") # replace "Other" <- "Mix"
all_colors <- c("#EAEAEA",
                "#8A91D0",
                "#92C5DE",
                "#11C1FF",
                "#808080",
                "#FFFF00",
                "#FFC34D",
                "#FF0000",
                "#006400",
                "#008000",
                "#006400","#FF4500","#66CDAA") # "#0000CC" <- Het_F

names(all_colors) <- all_states

# merge all dataframes from three Models 
df <- as.data.frame(rbind(chmm_df,episegmix_df,episegmixmeth_df))
# reorder the levels
#
df$Model <- factor(df$Model, levels = c("EpiSegMixMeth", "EpiSegMix", "ChromHMM"))

text_size=16
# boxplot
p <- ggplot(df, aes(x=corrected, y=coverage,fill = Model)) + 
  geom_boxplot() + xlim(c("NS","Het_C","Het_F","ReprPC","Mix","Enh_Wk","Enh","Pr","Tx_Str","Tx_Wk")) + 
  theme_minimal() +   xlab("States") + ylab("Coverage (%)")

p <- p + theme(axis.text = element_text(size = text_size,color = "black")) # changes axis labels

p <- p + theme(axis.title = element_text(size = text_size,color = "black")) # change axis titles

p <- p + theme(text = element_text(size = text_size,color = "black")) + guides(fill="none")

p 

jaccard <- function(annotation_path,jaccard_path,identifier){
  # reading annotation file
  annotation <- read.csv(annotation_path,header = TRUE)
  # temporary id
  annotation$id <- paste(annotation$EpiRR,annotation$states,sep = "_") 
  
  # read data file
  jaccard_max_df <- read.table(jaccard_path,header = FALSE)
  colnames(jaccard_max_df) <- c("epirr","states","final_states","ji")
  # temporary id
  jaccard_max_df$id <- paste(jaccard_max_df$epirr,jaccard_max_df$states,sep = "_")
  
  # merging the annotation and the genomic coverage based on the temporary id created
  jaccard_max_df.merged <- merge(jaccard_max_df,annotation,by = "id",all.x = TRUE)
  df <- jaccard_max_df.merged[,c("ji","corrected")]
  df$Model <- rep(identifier,nrow(df))
  return(df)
}


# Jaccard Index files 
jaccard_dfs <- c("Analysis/jaccard_index_full_reduced/ChromHMM_full_reduced_max_si.txt",
                 "Analysis/jaccard_index_full_reduced/EpiSegMix_full_reduced_max_si.txt",
                 "Analysis/jaccard_index_full_reduced/EpiSegMixMeth_full_reduced_max_si.txt")




jaccard_chmm <- jaccard(annotation_path = annotation_list[1],jaccard_path = jaccard_dfs[1],identifier = "ChromHMM")
jaccard_chmm <- rbind(jaccard_chmm, data.frame(ji = as.numeric(0), corrected = "Enh_Wk", Model = "ChromHMM"))

jaccard_episegmix <- jaccard(annotation_path = annotation_list[2], jaccard_path = jaccard_dfs[2], identifier = "EpiSegMix")
jaccard_episegmixmeth <- jaccard(annotation_path = annotation_list[2], jaccard_path= jaccard_dfs[3], identifier = "EpiSegMixMeth")






  
ji <- rbind(jaccard_chmm,jaccard_episegmix,jaccard_episegmixmeth)
colnames(ji) <- c("SI","full","Model")
head(ji)
ji$Model <- factor(ji$Model, levels = c("EpiSegMixMeth", "EpiSegMix", "ChromHMM"))

j <- ggplot(ji, aes(x=full, y=SI,fill = Model)) + 
  geom_boxplot() + xlim(c("NS","Het_C","Het_F","ReprPC","Mix","Enh_Wk","Enh","Pr","Tx_Str","Tx_Wk")) + 
  theme_minimal() +   xlab("States") + ylab("Similarity index")


j <- j + theme(axis.text = element_text(size = text_size,color = "black")) # changes axis labels

j <- j + theme(axis.title = element_text(size = text_size,color = "black")) # change axis titles

j <- j + theme(text = element_text(size = text_size,color = "black")) + guides(fill="none")


#pdf(file = "Analysis/Figures/Figure_II/jaccard_index_coverage.pdf",width = 10,height = 10)
png(file = "Analysis/Figures/Figure2/jaccard_index_coverage.v3.png",width = 9,height = 5,res = 300,units = "in")
grid.arrange(j,                             # First row with one plot spaning over 2 columns
             p, ncol = 1, # Second row with 2 plots in 2 different columns
             nrow = 2) 

dev.off()
ggsave(filename = "Analysis/genomic_coverage/coverage_Models.png",plot =p,
       dpi = 200,
       width = 15,height = 5)


# finding max jaccard index for each states and epirr : 
ji_max <- read.table("Analysis/jaccard_index_full_reduced/EpiSegMixMeth_full_reduced_max_si.txt",
                      header = FALSE)

epirrids=unique(ji_max$V1)
deep_annotation <- read.csv("Analysis/classifier/Random_forest_classifier/EpiSegMixMeth_annotation_without_features.csv",header = TRUE)


for (epirrid in epirrids) {
  temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid),c(2,6)]
  replacement_vector <- setNames(temp$corrected, temp$states)
  ji_max$V2[which(ji_max$V1==epirrid)] <- replacement_vector[as.character(ji_max$V2[which(ji_max$V1==epirrid)])]  
}

ji_max$V5 <- rep("EpiSegMixMeth",nrow(ji_max))
colnames(ji_max) <- c("epirr","full","reduced","SI","Model")
#ji_max$Model <- factor(ji$Model, levels = c("EpiSegMixMeth", "EpiSegMix", "ChromHMM"))

j <- ggplot(ji_max, aes(x=full, y=SI,fill = Model)) + 
  geom_boxplot() + xlim(c("NS","Het_C","Het_F","ReprPC","Mix","Enh_Wk","Enh","Pr","Tx_Str","Tx_Wk")) + 
  theme_minimal() +   xlab("States") + ylab("Similarity index")

j <- j + theme(axis.text = element_text(size = 20)) # changes axis labels

j <- j + theme(axis.title = element_text(size = 20)) # change axis titles

j <- j + theme(text = element_text(size = 20)) 

j



