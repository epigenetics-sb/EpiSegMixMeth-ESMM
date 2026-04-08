# set the environment
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "keras","nnet","reticulate","cluster",
                      "ComplexHeatmap","circlize","dplyr","readr","gridExtra","ggrepel")
# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

# environment 
options(scipen=0) 

# set the directory 
setwd("./data/shared/IHEC/Analysis/heatmap_matrix/emission_matrix/")

# final_df is read from ./project/Analysis/jaccard_index_full_reduced/all_compared/All_models_max_jaccard_coverage_full_model_coverage.tab
df <- read.table("./project/Analysis/jaccard_index_full_reduced/all_compared/All_models_max_jaccard_coverage_full_model_coverage.tab",
                   header = TRUE)



# colors schemes for heatmap 
col_fun = colorRamp2(c(0,1), c("white","#3F7F93"))


# color according to facets
all_states <- c("NS","Mix","Het_C","Het_F","ReprPC",
                "Enh_Wk","Enh","Pr","Tx_Wk","Tx_Str",
                "Tx","Pr_Wk","Rep") # replace "Other" <- "Mix"
all_colors <- c("#EAEAEA","#8A91D0","#92C5DE","#11C1FF","#808080",
                "#FFFF00","#FFC34D","#FF0000","#006400","#008000",
                "#006400","#FF4500","#66CDAA") # "#0000CC" <- Het_F

names(all_colors) <- all_states

# colors for train features
train_cols <- c("Absent"="red","Present"="green")

# colors for histone marks 
epigenetics_rgb <- c("236,232,56","238,30,37","13,127,66","184,185,189","144,206,219","225,148,37","30,30,30","81,89,168","30,176,75")
epigenetics_marks <-c("H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K9me3","H3K27ac","Input","WGBS","RNA-seq")
names(epigenetics_rgb) <- epigenetics_marks
epigenetics_hex <- sapply(strsplit(epigenetics_rgb, ","), function(x) {
  rgb(red = as.integer(x[1]), 
      green = as.integer(x[2]), 
      blue = as.integer(x[3]),
      maxColorValue = 255)
})
names(epigenetics_hex) <- epigenetics_marks


merge_epirr <- function(search_path=search_path,search_pattern=search_pattern,jaccard_coverage=jaccard_coverage,mark=mark,model=model){
  
  # list the files 
  files <- list.files(path = search_path,
                      pattern = search_pattern,
                      full.names = TRUE)
  
  # adapt the code if needed for finding pattern of particular tyoe in case multiple replicates within same folder 
  data <- lapply(files, function(x) read.table(x,header = TRUE,stringsAsFactors = FALSE))
  dd <- do.call(rbind, data)
  
  #mim-max-normalization, adapt column numbers for min-max
  cols <- c(3,4,5,6,7,8,10)
  epirrids <- unique(dd$epirr)
  for (epirrid in epirrids){
    for (x in cols ){
      min_val <- min(dd[dd$epirr == epirrid, x], na.rm = TRUE)
      max_val <- max(dd[dd$epirr == epirrid, x], na.rm = TRUE)
      dd[dd$epirr==epirrid,x]=(dd[dd$epirr==epirrid,x] - min_val)/(max_val-min_val)
    }
  }
  
  dd$id <- paste(dd$epirr,dd$states,sep="_")
  jaccard_coverage <- jaccard_coverage[which(jaccard_coverage$mark==mark & jaccard_coverage$model == model),]
  jaccard_coverage$id <- paste(jaccard_coverage$EpiRR,jaccard_coverage$target,sep="_")
  dd_merge <- merge(dd,jaccard_coverage,by="id",all.x=TRUE)
  
  # convert the target covertage between 0 and 1 
  dd_merge$target_coverage <- round(dd_merge$target_coverage/100,4)
  
  # set the levels for the heatmap 
  ## !!!! Change here between chromHMM and EpiSegMix and EpiSegMixMeth
  dd_merge$source <- factor(dd_merge$source,levels = c("Tx_Str","Tx_Wk","Tx","Pr","Pr_Wk",
                                                           "Enh","Enh_Wk","ReprPC","Mix","Het_F","Het_C",
                                                           "Rep","NS"))
  
  
  return(dd_merge)
}


## !!!! Change here 
data <- merge_epirr(search_path = "H3K27me3-H3K9me3/",
                    search_pattern = "^EpiSegMix_",
                    jaccard_coverage = df,
                    mark = "-H3K27-9me3",
                    model = "EpiSegMix")

## !!!! Change here for ChromHMM only uncomment
# data$source[which(data$source == "Mix")] <- "Enh_Wk"

data <- data[complete.cases(data), ]
data <- data[order(data$jaccard,data$DNA_methyl),]


# add DNA methyl for EpiSegMixMeth
## !!!! Change here
features_to_plot <- c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K27ac","H3K4me3","DNA_methyl")
features_to_plot <- c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K27ac","H3K4me3")

features <- as.matrix(data[,features_to_plot])
train_features <- c(rep("Present",7)) # ## !!!! Change here 7 to 6 for ChromHMM and EpiSegMix
train_features <- c(rep("Present",6))
## !!!! Change here 1- for -H3K9me3 2 - -H3K27me3 
train_features[c(1,2)] <- "Absent"

# metadata$jsd_qc.syn_jsd for the plotting quality of the data
# chipseq quality data 
metadata <- read.csv("./project/metadata/epiatlas_chipseq_qc_summary.csv",header=TRUE)
metadata$epirr_id <- sapply(metadata$epirr_id, function(x) strsplit(x,".",fixed = TRUE)[[1]][1])
epirrs <- unique(data$epirr)
metadata <- metadata[metadata$epirr_id %in% epirrs, c("epirr_id","antibody","jsd_qc.syn_jsd")]
metadata_wide <- reshape(metadata, idvar = "epirr_id", timevar = "antibody", direction = "wide")
colnames(metadata_wide) <- c(sub("jsd_qc.syn_jsd.", "", colnames(metadata_wide)))
metadata_wide <- metadata_wide[,features_to_plot[1:6]] # use same column order as in the main dataframe where the heatmap is generated
# metadata_wide$DNA_methyl <- 0 # !!! Chnage here Dummy for DNA methyl model 

# column annotation
col_ha = HeatmapAnnotation(ChIP_Quality = anno_boxplot(metadata_wide,height = unit(3, "cm"),
                                                   gp=gpar(fontsize=20,fill=epigenetics_hex[colnames(metadata_wide)]),
                                                   axis_param=list(gp = gpar(fontsize = 20))),
                       Train_features=train_features,col=list(Train_features=train_cols),
                       annotation_legend_param = list(
                         Train_features=list(
                           title_gp=gpar(fontsize=20,fontface="bold"),
                           labels_gp = gpar(fontsize = 20))))






# row annotation # Ideally for WGBS color is  > "#5159A8"
row_ha = rowAnnotation(States=data$source,
                       `Meth. %` = anno_barplot(bar_width = 1,
                                                  width = unit(2, "cm"),
                                                  data$DNA_methyl,
                                                  gp=gpar(fill="#5159A8")),
                       `States %` = anno_barplot(bar_width = 1,
                                                 width = unit(2, "cm"),
                                                 data$target_coverage,
                                                 gp=gpar(fill="#3F7F93")),
                       Jaccard = anno_barplot(bar_width = 1,
                                              width = unit(2, "cm"),
                                              data$jaccard,
                                              gp=gpar(fill=all_colors[factor(data$source,levels = names(all_colors))])),
                       annotation_name_gp = gpar(fontsize=20),
                       annotation_name_rot = 90,
                       gap = unit(c(3, 3, 3), "mm"),
                       col = list(States = all_colors),
                       annotation_legend_param = list(
                         States=list(
                           title_gp=gpar(fontsize=20,fontface="bold"),
                           labels_gp = gpar(fontsize = 20))))


## !!!! Change here
a <- Heatmap(features,show_row_names = FALSE,show_column_names=TRUE,
             show_row_dend = FALSE,
             cluster_columns = FALSE ,cluster_rows = FALSE,
             right_annotation = row_ha,
             col = col_fun,
             column_gap = unit(3, "mm"),
             border = TRUE,
             cluster_row_slices = FALSE,
             column_split = train_features,
             top_annotation = col_ha,
             row_split = data$source,
             column_title = "EpiSegMixMeth: -H3K27me3-H3K9me3",
             column_title_gp = gpar(fontsize = 20),  
             column_names_gp = gpar(fontsize = 20),   
             row_title_gp = gpar(fontsize=20),
             row_names_gp = gpar(fontsize = 20),row_gap = unit(5, "mm"),
             heatmap_legend_param = list(
               title = "Min-Max",direction="horizontal",                     
               title_gp = gpar(fontsize = 20,fontface="bold"),       
               labels_gp = gpar(fontsize = 20)        
             ))
a
# ChromHMM - H3K27me3, H3K9me3 and removing both 
# EpiSegMix - H3K27me3, H3K9me3 and removing both 
# EpiSegMixMeth - H3K27me3, H3K9me3 and removing both 


## !!!! Change here
png(file = "EpiSegMixMeth_H3K27me3-H3K9me3_sorted_jaccard.png",
    width = 15,
    height = 17,
    res = 200,
    units = "in")
a
dev.off()


#### Code for the comparison of the full model emission matrix
annotation_list <- c("./project/Analysis/classifier/classifier_input/IHEC_EpiSegMix_meth_10N_all_states_samples_classifier.tab",
                     "./project/Analysis/classifier/classifier_input/IHEC_EpiSegMix_10N_all_states_samples_classifier.tab",
                     "./project/Analysis/classifier/classifier_input/IHEC_ChromHMM_10N_all_states_samples_classifier.tab")
data <- read.table(annotation_list[3],
                   header = T)
data <- data[,c(1:8,10,12)]
data$EpiRR <- sapply(data$EpiRR, function(x) strsplit(x,".",fixed = TRUE)[[1]][1])
# split the dataframe for each epirr
data_list <- lapply(epirrs, function(x) {temp <- data[data$EpiRR==x,] ; 
                                          colnames(temp) <- NULL; 
                                          column_names <- temp[1,2:ncol(temp)];
                                          colnames(temp) <- c("EpiRR",column_names); 
                                          temp <- temp[2:11,];
                                          return(temp)})
dd <- do.call(rbind, data_list)

# Add the metdata information annotated state and more information CHIP quality infromation 
# add the metadata for the the states annotation 
annotation_list_states <- c("./project/Analysis/classifier/Random_forest_classifier/EpiSegMixMeth_annotation_without_features.csv",
                            "./project/Analysis/classifier/Random_forest_classifier/EpiSegMix_annotation_without_features.csv",
                            "./project/Analysis/classifier/Random_forest_classifier/ChromHMM_annotation_without_features.csv")
state_annotation_metadata <- read.csv(annotation_list_states[3],
                                      header=TRUE)
state_annotation_metadata$EpiRR <- sapply(state_annotation_metadata$EpiRR, function(x) strsplit(x,".",fixed = TRUE)[[1]][1])
state_annotation_metadata$corrected[which(state_annotation_metadata$corrected=="Mix")] <- "Enh_Wk"
for (epirrid in epirrs) {
  temp <- state_annotation_metadata[which(state_annotation_metadata$EpiRR==epirrid),c(2,6)]
  replacement_vector <- setNames(temp$corrected, temp$states)
  dd$states[which(dd$EpiRR==epirrid)] <- replacement_vector[as.character(dd$states[which(dd$EpiRR==epirrid)])]
}
dd[3:10] <- lapply(dd[3:10], as.numeric)
# min-max normalization: 
cols <- c(3:8)
for (epirrid in epirrs){
  for (x in cols ){
    min_val <- min(dd[dd$EpiRR == epirrid, x], na.rm = TRUE)
    max_val <- max(dd[dd$EpiRR == epirrid, x], na.rm = TRUE)
    dd[dd$EpiRR==epirrid,x]=(dd[dd$EpiRR==epirrid,x] - min_val)/(max_val-min_val)
  }
}

dd$states <- factor(dd$states,levels = c("Tx_Str","Tx_Wk","Tx","Pr","Pr_Wk",
                                                     "Enh","Enh_Wk","ReprPC","Mix","Het_F","Het_C",
                                                     "Rep","NS"))
dd <- dd[order(dd$DNA_methyl),]
features <- as.matrix(dd[,features_to_plot])
train_features <- c(rep("Present",6))

col_ha = HeatmapAnnotation(ChIP_Quality = anno_boxplot(metadata_wide,height = unit(3, "cm"),
                                                       gp=gpar(fontsize=20,fill=epigenetics_hex[colnames(metadata_wide)]),
                                                       axis_param=list(gp = gpar(fontsize = 20))),
                           Train_features=train_features,col=list(Train_features=train_cols),
                           annotation_legend_param = list(
                             Train_features=list(
                               title_gp=gpar(fontsize=20,fontface="bold"),
                               labels_gp = gpar(fontsize = 20))))






# row annotation # Ideally for WGBS color is  > "#5159A8"
row_ha = rowAnnotation(States=dd$states,
                       `Meth. %` = anno_barplot(bar_width = 1,
                                                width = unit(2, "cm"),
                                                dd$DNA_methyl,
                                                gp=gpar(fill="#5159A8")),
                       `States %` = anno_barplot(bar_width = 1,
                                                 width = unit(2, "cm"),
                                                 round(dd$coverage/100,4),
                                                 gp=gpar(fill="#3F7F93")),
                       annotation_name_gp = gpar(fontsize=20),
                       annotation_name_rot = 90,
                       gap = unit(c(3, 3, 3), "mm"),
                       col = list(States = all_colors),
                       annotation_legend_param = list(
                         States=list(
                           title_gp=gpar(fontsize=20,fontface="bold"),
                           labels_gp = gpar(fontsize = 20))))


## !!!! Change here
b <- Heatmap(features,show_row_names = FALSE,show_column_names=TRUE,
             show_row_dend = FALSE,
             cluster_columns = FALSE ,cluster_rows = FALSE,
             right_annotation = row_ha,
             col = col_fun,
             column_gap = unit(3, "mm"),
             border = TRUE,
             cluster_row_slices = FALSE,
             column_split = train_features,
             top_annotation = col_ha,
             row_split = dd$states,
             column_title = "ChromHMM",
             column_title_gp = gpar(fontsize = 20),  
             column_names_gp = gpar(fontsize = 20),   
             row_title_gp = gpar(fontsize=20),
             row_names_gp = gpar(fontsize = 20),row_gap = unit(5, "mm"),
             heatmap_legend_param = list(
               title = "Min-Max",direction="horizontal",                     
               title_gp = gpar(fontsize = 20,fontface="bold"),       
               labels_gp = gpar(fontsize = 20)        
             ))

b
png(file = "ChromHMM.png",
    width = 15,
    height = 17,
    res = 200,
    units = "in")
b
dev.off()
