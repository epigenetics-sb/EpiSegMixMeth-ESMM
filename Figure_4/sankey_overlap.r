# loading necessary libraries
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggalluvial)


# set the working directory
setwd("./resources/")

# mapping file for the states to biological funtional states
deep_annotation <- read.csv("IHEC/Analysis/classifier/Random_forest_classifier/EpiSegMixMeth_annotation_without_features.csv",header = TRUE)

# NBC-I <- IHECRE00000354.2
# GCBC-I <- IHECRE00000332.2
# MBC-I <- IHECRE00000335.2 
# PC-I <- IHECRE00000341.2

epirrid <- c("IHECRE00000354.2","IHECRE00000332.2","IHECRE00000335.2","IHECRE00000341.2")
# replace columns for selected epirrs with biolgical states

# color mapping
states <- c("NS","Mix","Het_C",
            "Het_F","ReprPC","Enh_Wk",
            "Enh","Pr","Tx_Wk","Tx_Str",
            "Tx","Rep","Pr_Wk")

colors <- c("234,234,234","138,145,208","146,197,222",
            "17,193,255","128,128,128","255,255,0",
            "255,195,77","255,0,0","0,100,0","0,128,0",
            "0,100,0","102,205,170","255,69,0")

names(colors) <- states
# convert colors to hex codes
colors.hex <- sapply(colors, function(x) {
  rgb_vals <- as.numeric(unlist(strsplit(x, ",")))
  rgb_str <- sprintf("#%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  return(rgb_str)
})

#colors for cell types 
# celltypes_colors <- c(pc="#ff00e2ee", mbc = "#0031ffee",gcbc ="#ff6442ee" , nbc="#3ee709ee")
# celltypes_colors

######################### Memory b cell ###################################################

# read the counts label
df <- read.table("IHEC/Analysis/overlap/NBC_GCBC_MBC_I_overlap_counts.txt",header = TRUE)

# calculate basepairs
df$Counts <- df$Counts*200

# NBC
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[1]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df$NBC <- replacement_vector[as.character(df$NBC)]

# GCBC
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[2]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df$GCBC <- replacement_vector[as.character(df$GCBC)]

# MBC 
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[3]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df$MBC <- replacement_vector[as.character(df$MBC)]

states <- c(levels(factor(df$NBC)),levels(factor(df$GCBC)),levels(factor(df$MBC)))


# plotting the overlap
mbc_plot <- ggplot(df,aes(axis1=NBC,axis2=GCBC,axis3=MBC,y=Counts)) +
  scale_x_discrete(limits=c("NBC","GCBC","MBC"),expand=c(0.1,0.1)) + 
  xlab("B cells") + ylab("Base pairs") +
  geom_alluvium(aes(fill=NBC)) + 
  geom_stratum(fill=c(rev(colors.hex[states])))  + 
  theme_minimal() + 
  scale_fill_manual("States",values = colors.hex)  #
  # geom_text_repel(aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),stat = "stratum", size = 4, direction = "y", nudge_x = -.4) + 
  # geom_text_repel(aes(label = ifelse(after_stat(x)  == 3, as.character(after_stat(stratum)), "")), stat = "stratum", size = 4, direction = "y", nudge_x = .4)

mbc_plot <- mbc_plot + scale_y_continuous(labels = label_comma(scale = 1e-9, suffix = "Gb", accuracy = 1)) 

mbc_plot <- mbc_plot +  theme(axis.text.x = element_text(face="bold",colour = c("#3ee709ee","#ff6442ee","#0031ffee")))

mbc_plot <- mbc_plot + theme(axis.text = element_text(size = 20)) # changes axis labels

mbc_plot <- mbc_plot + theme(axis.title = element_text(size = 20)) # change axis titles

mbc_plot <- mbc_plot + theme(text = element_text(size = 20)) + guides(fill = "none")   # remove guide false to see legends


mbc_plot
ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_MBC_I_overlap_counts.txt.pdf",plot =  mbc_plot,
       width = 8,height = 8)
ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_MBC_I_overlap_counts.txt.png",plot =  mbc_plot,
       dpi = 300,
       width = 5,height = 5)

ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_MBC_I_overlap_counts.txt_repel.pdf",plot =  mbc_plot,
       width = 12,height = 10)
ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_MBC_I_overlap_counts.txt_repel.png",plot =  mbc_plot,
       dpi = 200,
       width = 12,height = 10)

############################################### plasma cell ################################################
df <- read.table("IHEC/Analysis/overlap/NBC_GCBC_PC_I_overlap_counts.txt",header = TRUE)

# calculate basepairs
df$Counts <- df$Counts*200

# NBC
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[1]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df$NBC <- replacement_vector[as.character(df$NBC)]

# GCBC
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[2]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df$GCBC <- replacement_vector[as.character(df$GCBC)]

# PC 
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[4]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df$PC <- replacement_vector[as.character(df$PC)]

states <- c(levels(factor(df$NBC)),levels(factor(df$GCBC)),levels(factor(df$PC)))

# plot for PC
pc_plot <- ggplot(df,aes(axis1=NBC,axis2=GCBC,axis3=PC,y=Counts)) +
  scale_x_discrete(limits=c("NBC","GCBC","PC"),expand=c(0.1,0.1)) + 
  xlab("B cells") + ylab("Base pairs") + 
  geom_alluvium(aes(fill=NBC)) + 
  geom_stratum(fill=c(rev(colors.hex[states])))  + 
  theme_minimal() + 
  scale_fill_manual("States",values = colors.hex) #+
  #geom_text_repel(aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),stat = "stratum", size = 4, direction = "y", nudge_x = -.4) + 
  #geom_text_repel(aes(label = ifelse(after_stat(x)  == 3, as.character(after_stat(stratum)), "")), stat = "stratum", size = 4, direction = "y", nudge_x = .4) 

pc_plot <-  pc_plot + scale_y_continuous(labels = label_comma(scale = 1e-9, suffix = "Gb", accuracy = 1))

pc_plot <- pc_plot +  theme(axis.text.x = element_text(face="bold",colour = c("#3ee709ee","#ff6442ee","#ff00e2ee")))

pc_plot <- pc_plot + theme(axis.text = element_text(size = 20)) # changes axis labels

pc_plot <- pc_plot + theme(axis.title = element_text(size = 20))  + guides(fill = "none")  # change axis titles



ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_PC_I_overlap_counts.txt.pdf",plot =pc_plot,
       width = 5,height = 5)
ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_PC_I_overlap_counts.txt.png",plot =pc_plot,
       dpi = 300,
       width = 5,height = 5)


ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_PC_I_overlap_counts.txt_repel.png",plot =pc_plot,
       dpi = 300,
       width = 5,height = 5)
ggsave(filename = "IHEC/Analysis/overlap/NBC_GCBC_PC_I_overlap_counts.txt_repel.pdf",plot =pc_plot,
       dpi = 300,
       width = 12,height = 10) 



