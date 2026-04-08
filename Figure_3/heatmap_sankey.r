# set the environment
packages_to_load <- c("randomForest","clue","ggplot2","reshape2",
                      "ggalluvial","ggdendro","networkD3","htmlwidgets",
                      "keras","nnet","reticulate","cluster",
                      "ComplexHeatmap","circlize","dplyr","readr","gridExtra","ggrepel","scales")
# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

set.seed(123)
# set the working directory
setwd("./project/")

file_to_read <- list(repimeth = "results/EpiSegMix_meth_4_marks/10N/IHECRE00000354.2/IHECRE00000354.2_EpiSegMix_10N_classifier.txt",
                     fepimeth ="results/EpiSegMix_meth_10N/IHECRE00000354.2/IHECRE00000354.2_EpiSegMix_10N_classifier.txt",
                     repi ="results/EpiSegMix_4_marks/10N/IHECRE00000354.2/IHECRE00000354.2_EpiSegMix_10N_classifier.txt",
                  fepi ="results/EpiSegMix_10N/IHECRE00000354.2/IHECRE00000354.2_EpiSegMix_10N_classifier.txt",
                  rchmm = "results/ChromHMM_4marks/10N/IHECRE00000354.2/IHECRE00000354.2_ChromHMM_10N_classifier.txt",
                  fchmm ="results/ChromHMM_10N/IHECRE00000354.2/IHECRE00000354.2_ChromHMM_10N_classifier.txt")  

# 0298 IHEC epirr for nbc II ReprPC
nbc_II <- c("ReprPC","NS","Pr","Enh","Mix","Tx_Str","Het_F","Het_C","Enh_Wk","Tx_Wk")
nbc_II_file <- "results/EpiSegMix_10N/IHECRE00000298.2/IHECRE00000298.2_EpiSegMix_10N_classifier.txt"
rnbc_II_file <- "results/EpiSegMix_4_marks/10N/IHECRE00000298.2/IHECRE00000298.2_EpiSegMix_10N_classifier.txt"
hm_full <- Heatmap_episegmix(model_file = nbc_II_file, 
                             features = c(custom_order[1:6]),
                             state_annotation = nbc_II,
                             splitting_columns = FALSE)
hm_reduced <- Heatmap_episegmix(model_file = rnbc_II_file, 
                                features = c(custom_order[1:6]),
                                state_annotation = states$rstates,
                                splitting_columns = TRUE)

hm_reduced
hm_full
psankey$normal
# State sequences of the models
# After observation in EpiSegMix the states were switched between Het_F and Enh_Wk as this made more sense
states <- list(rstates = c(1:10),
            episegmixmeth= c("Het_F","NS","Pr","ReprPC","Tx_Str","Enh","Enh_Wk","Het_C","Mix","Tx_Wk"),
            episegmix = c("Enh_Wk","NS","Pr","ReprPC","Tx_Str","Enh","Mix","Het_F","Het_C","Tx_Wk"),
            chmm = c("Het_C","NS","Het_F","ReprPC","Pr_Wk","Enh","Mix","Tx_Wk","Tx_Str","Pr"))



# Heatmap features to be used for EpiSegMixMeth and first 6 for ChromHMM and EpiSegMix
custom_order <- c("H3K4me3","H3K27ac","H3K4me1","H3K36me3","H3K27me3","H3K9me3","DNA_methyl")



#function to generate the heatmap
Heatmap_episegmix <- function(model_file,features,state_annotation,splitting_columns){
  # read the file
  df <-  read.table(model_file,header = TRUE)
  # select the columns to be normalized
  cols <- c(2:ncol(df)) 
  
  
  #mim-max-normalization
  for (x in cols ){
    min_val <- min(df[, x], na.rm = TRUE)
    max_val <- max(df[, x], na.rm = TRUE)
    df[,x]=(df[,x] - min_val)/(max_val-min_val)
  }
  
  
  # adding rownames 
  rownames(df) <- state_annotation
  
  # adding color annotation 
  col_fun = colorRamp2(c(0,1), c("white","#3F7F93"))
  
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
                  "Tx","Pr_Wk","Rep",
                  c(1:10)) # replace "Other" <- "Mix"
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
                  "#006400","#FF4500","#66CDAA",
                  "#CFEAE4","#B01C49","#A5B6C5","#A6F936","#77213D","#969B2D","#71AF38","#CB2AAB","#1D29F8","#79C69E")
  
  names(all_colors) <- all_states
  
  # column annotation
  col_ha = HeatmapAnnotation(Marks=features,
                             annotation_name_gp = gpar(fontsize=20),
                             annotation_legend_param = list(
                               Marks=list(
                                 title_gp=gpar(fontsize=20,fontface="bold"),
                                 labels_gp = gpar(fontsize = 20))
                             ))
  
  # row annotation
  row_ha = rowAnnotation(States=c(levels(factor(state_annotation))),
                         col = list(States = all_colors),
                         annotation_name_gp = gpar(fontsize=20),
                         annotation_legend_param = list(
                           States=list(
                             title_gp=gpar(fontsize=20,fontface="bold"),
                             labels_gp = gpar(fontsize = 20))
                         ))
  
  # selecting features for the matrix
  features2 <- as.matrix(df[c(levels(factor(rownames(df)))),c(features)])
  
  
  #png(file = image_file,width = 15,height = 10,res = 200,units = "in")
  
  
  a <- Heatmap(features2,show_row_names = FALSE,show_column_names=TRUE,
               show_row_dend = FALSE,
               cluster_columns = FALSE ,cluster_rows = FALSE,
               right_annotation = row_ha,
               col = col_fun,
               cluster_row_slices = FALSE,
               border = TRUE,
               column_title_gp = gpar(fontsize = 20),  
               column_names_gp = gpar(fontsize = 20),   
               row_title_gp = gpar(fontsize=20),
               row_names_gp = gpar(fontsize = 20),row_gap = unit(5, "mm"),
               heatmap_legend_param = list(
                 title = "Min-Max",direction="horizontal",                     
                 title_gp = gpar(fontsize = 20,fontface="bold"),       
                 labels_gp = gpar(fontsize = 20)        
               ))
  
  if(splitting_columns==TRUE){
    state_annotation <- factor(state_annotation, levels = c("1", "2","3","4","5","6","7","8","9","10"))
    row_ha = rowAnnotation(States=c(levels(factor(state_annotation))),
                           col = list(States = all_colors),
                           annotation_name_gp = gpar(fontsize=20),
                           annotation_legend_param = list(
                             States=list(at = c("1", "2","3","4","5","6","7","8","9","10"),
                               title_gp=gpar(fontsize=20,fontface="bold"),
                               labels_gp = gpar(fontsize = 20))
                           ))
    train_features <- c(rep("True",ncol(features2)))
    train_features[which(colnames(features2)=="H3K27me3")] <- "False"
    train_features[which(colnames(features2)=="H3K9me3")] <- "False"
    train_features <- factor(train_features, levels = c("True", "False"))
    a <- Heatmap(features2,show_row_names = FALSE,show_column_names=TRUE,
                 show_row_dend = FALSE,
                 cluster_columns = FALSE ,cluster_rows = FALSE,
                 right_annotation = row_ha,
                 col = col_fun,
                 column_gap = unit(3, "mm"),
                 border = TRUE,
                 cluster_row_slices = FALSE,
                 column_split = train_features,
                 column_title_gp = gpar(fontsize = 20),  
                 column_names_gp = gpar(fontsize = 20),   
                 row_title_gp = gpar(fontsize=20),
                 row_names_gp = gpar(fontsize = 20),row_gap = unit(5, "mm"),
                 heatmap_legend_param = list(
                   title = "Min-Max",direction="horizontal",                     
                   title_gp = gpar(fontsize = 20,fontface="bold"),       
                   labels_gp = gpar(fontsize = 20)        
                 ))
    }
  
  # plot the heatmap 
  a
  
  # save the plot and close the plot
  return(a)
  
}




hm_reduced <- Heatmap_episegmix(model_file = file_to_read$repimeth, 
                  features = c(custom_order[1:7]),
                  state_annotation = states$rstates,
                  splitting_columns = TRUE)


hm_full <- Heatmap_episegmix(model_file = file_to_read$fepi, 
                                        features = c(custom_order[1:6]),
                                        state_annotation = states$episegmix,
                                        splitting_columns = FALSE)

hm_full

png(file = "Analysis/sankey_plot_Bcell_full_reduced/NBC_II_EpiSegMix_Heatmap_reduced.png",
    width = 7.5,
    height = 5,
    res = 200,
    units = "in")
hm_reduced
dev.off()

pdf(file = "Analysis/sankey_plot_Bcell_full_reduced/NBC_I_EpiSegMix_Heatmap_reduced.pdf",
    width = 7.5,
    height = 5)
hm_reduced
dev.off()

png(file = "Analysis/sankey_plot_Bcell_full_reduced/NBC_II_EpiSegMix_Heatmap_full.png",
    width = 7.5,
    height = 5,
    res = 200,
    units = "in")
hm_full
dev.off()

pdf(file = "Analysis/sankey_plot_Bcell_full_reduced/NBC_I_EpiSegMix_Heatmap_full.pdf",
    width = 7.5,
    height = 5)
hm_full
dev.off()

# sankey part 


sankey_files <- list(sepimeth = "Analysis/sankey_plot_Bcell_full_reduced/EpiSegMixMeth_reduced_full_counts.txt",
                  sepi = "Analysis/sankey_plot_Bcell_full_reduced/EpiSegMix_reduced_full_counts.txt",
                  schmm = "Analysis/sankey_plot_Bcell_full_reduced/ChromHMM_reduced_full_counts.txt")

sankey_episegmix <- function(model_file,state_annotation){
  
  sdf <- read.table(model_file,header = FALSE)
  replacement_vector <- setNames(state_annotation, c(1:10))
  sdf$V1 <- replacement_vector[as.character(sdf$V1)]
  colnames(sdf) <- c("Full","Reduced","Counts")
  sorder <- c(levels(factor(sdf$Full)),levels(factor(sdf$Reduced)))
  sdf$Counts <- sdf$Counts * 200 #  convert bins to base pairs
  
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
                  "Tx","Pr_Wk","Rep",
                  c(1:10)) # replace "Other" <- "Mix"
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
                  "#006400","#FF4500","#66CDAA",
                  "#CFEAE4","#B01C49","#A5B6C5","#A6F936","#77213D","#969B2D","#71AF38","#CB2AAB","#1D29F8","#79C69E")
  
  names(all_colors) <- all_states
  
  
  pplot <-  ggplot(sdf,aes(axis1=Full,axis2=Reduced,y=Counts)) +
    scale_x_discrete(limits=c("Full","Reduced"),expand=c(0.0,0.25)) +
    xlab("Models") + 
    geom_alluvium(aes(fill=Full),width=1/6) +  
    geom_stratum(fill= c(rev(all_colors[sorder[1:10]]),rev(all_colors[sorder[11:20]])))  +
    theme_minimal() + 
    scale_fill_manual(values = all_colors) + 
    guides(fill="none")  + 
    labs(y = "Base pairs") +
    scale_y_continuous(labels = label_comma(scale = 1e-6, suffix = "Mb", accuracy = 1)) + 
    geom_text_repel(aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),stat = "stratum", size = 4, direction = "y", nudge_x = -.4) + 
    geom_text_repel(aes(label = ifelse(after_stat(x)  == 2, as.character(after_stat(stratum)), "")), stat = "stratum", size = 4, direction = "y", nudge_x = .4)
  
  pplot <- pplot + theme(axis.text = element_text(size = 20)) # changes axis labels
  pplot <- pplot + theme(axis.title = element_text(size = 20)) # change axis titles
  pplot <- pplot + theme(text = element_text(size = 20)) 
  
  gplot <- ggplot(sdf,aes(axis1=Full,axis2=Reduced,y=Counts)) +
    scale_x_discrete(limits=c("Full","Reduced"),expand=c(0.0,0.25)) +
    xlab("Models") + 
    geom_alluvium(aes(fill=Full),width=1/3) +  
    geom_stratum(fill= c(rev(all_colors[sorder[1:10]]),rev(all_colors[sorder[11:20]])))  +
    theme_minimal() + 
    scale_fill_manual(values = all_colors) + 
    guides(fill="none")  + 
    labs(y = "Base pairs") +
    scale_y_continuous(labels = label_comma(scale = 1e-9, suffix = "Gb", accuracy = 1))
  
    #scale_y_continuous(labels = label_number(scale_cut = cut_long_scale(),suffix = "b"))
  
  gplot <- gplot + theme(axis.text = element_text(size = 20)) # changes axis labels
  gplot <- gplot + theme(axis.title = element_text(size = 20)) # change axis titles
  gplot <- gplot + theme(text = element_text(size = 20)) 
  
  sankey_plots <-list(pplot,gplot)
  names(sankey_plots) <- c("text_repel","normal")
  return(sankey_plots)
  
}

nbc_II_sankey_file <- "Analysis/sankey_plot_Bcell_full_reduced/EpiSegMix_reduced_full_counts_IHECRE00000298.txt"
psankey <- sankey_episegmix(model_file = nbc_II_sankey_file,nbc_II)
psankey$normal


psankey <- sankey_episegmix(model_file = sankey_files$sepi,states$episegmix)
psankey$normal
psankey$text_repel


ggsave(filename = "Analysis/sankey_plot_Bcell_full_reduced/NBC_II_EpiSegMix_Sankey.png",plot =psankey$normal,
       dpi = 200,
       width = 5,height = 5)


ggsave(filename = "Analysis/sankey_plot_Bcell_full_reduced/NBC_I_EpiSegMix_Sankey_ggrepel.png",plot = psankey$text_repel,
       dpi = 200,
       width = 5,height = 5)


