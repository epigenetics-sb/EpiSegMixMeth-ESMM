library(dplyr)

library(tidyverse) 



setwd("./resources/")


# possible 
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


# read the dataframes
df1 <- read.table("IHEC/Analysis/hic_compartments_overlap/NBC_I_EpiSegMixMeth_10N_Hi-C.ovelap_counts.txt",header = FALSE)
df2 <- read.table("IHEC/Analysis/hic_compartments_overlap/GCBC_I_EpiSegMixMeth_10N_Hi-C.ovelap_counts.txt",header = FALSE)
df3 <- read.table("IHEC/Analysis/hic_compartments_overlap/MBC_I_EpiSegMixMeth_10N_Hi-C.ovelap_counts.txt",header = FALSE)
df4 <- read.table("IHEC/Analysis/hic_compartments_overlap/PC_I_EpiSegMixMeth_10N_Hi-C.ovelap_counts.txt",header = FALSE)
deep_annotation <- read.csv("IHEC/Analysis/classifier/Random_forest_classifier/EpiSegMixMeth_annotation_without_features.csv",header = TRUE)

df1$Type <-  rep("NBC",nrow(df1))
df2$Type <-  rep("GCBC",nrow(df2))
df3$Type <-  rep("MBC",nrow(df3))
df4$Type <-  rep("PC",nrow(df4))

epirrid <- c("IHECRE00000354.2","IHECRE00000332.2","IHECRE00000335.2","IHECRE00000341.2")
# replace columns for selected epirrs with biolgical states

# NBC
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[1]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df1$V3 <- replacement_vector[as.character(df1$V3)]

# GCBC
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[2]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df2$V3 <- replacement_vector[as.character(df2$V3)]

# MBC 
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[3]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df3$V3 <- replacement_vector[as.character(df3$V3)]

# PC 
temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid[4]),c(2,6)]
replacement_vector <- setNames(temp$corrected, temp$states)
df4$V3 <- replacement_vector[as.character(df4$V3)]


df <- rbind(df1,df2,df3,df4)
levels(df$Hi_C)

names(df) <- c("Counts","Hi_C","States","Type")
df$Hi_C <- factor(df$Hi_C, levels = c("A","B","I"))
df$Type <- factor(df$Type , levels=c("NBC","GCBC","MBC","PC"))

# When using group_by and summarize functon ensure to use dplyr package by speicifically calling via dplyr:: 
# Otherwise the plyr package if loaded overwrites the functions called and the following piece of code does not work
# as expected
df_summary <- df %>% 
  dplyr::group_by(Type,Hi_C) %>% 
  dplyr::summarize(Count=sum(Counts)) %>% 
  dplyr::mutate(Percentage = (Count)/sum(Count))

colors.hex.title <- c("#3ee709ee","#ff6442ee","#0031ffee","#ff00e2ee")
names(colors.hex.title) <- c("NBC","GCBC","MBC","PC")
colors.hex.title
p <- ggplot(df) + 
     geom_bar(mapping=aes(x=Hi_C,y=Counts,fill=States),position="fill",stat="identity") + 
     scale_fill_manual(values=colors.hex) + 
     facet_wrap(vars(Type),strip.position = "top",scales='free',ncol=4)  + 
     geom_text(data=df_summary,aes(x=Hi_C,y=1,label=round(Percentage,2)),color="black",vjust = -0.25, size = 6) +
     coord_cartesian(ylim = c(0, 1.1)) + 
     theme_minimal() + ylab("Counts %")
p <- p + theme_bw() + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), text = element_text(size = 20), plot.title = element_text(size = 20),legend.title = element_text(size = 20), legend.text = element_text(size = 20)) 

p
g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
codes <- levels(df$Type)
replacement_vector <- setNames(codes, stript)
replacement_vector
k <- 1
for (i in stript) {
  mark <- as.character(replacement_vector[as.character(i)])
  print(mark)
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <-colors.hex.title[mark]
}

g
pdf("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.pdf",height = 6, width=12)
grid.draw(g)
dev.off()

png("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.png",height = 6, width=12,units="in",res=200)
grid.draw(g)
dev.off()

ggsave("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.pdf",device="pdf",plot=p,dpi = 200,width = 12,height = 10)
ggsave("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.png",device="png",plot=p,dpi = 200,width = 12,height = 10)


df_summary <- df %>% 
  dplyr::group_by(Type,Hi_C) %>% 
  dplyr::summarize(Count=sum(Counts)) %>% 
  dplyr::mutate(Percentage = (Count)/sum(Count))

abi_colors.hex.title <- c("red","blue","green")
names(abi_colors.hex.title) <- c("A","B","I")
abi_colors.hex.title
p <- ggplot(df) + 
  geom_bar(mapping=aes(x=Type,y=Counts,fill=States),position="fill",stat="identity") + 
  scale_fill_manual(values=colors.hex) + 
  facet_wrap(vars(Hi_C),strip.position = "top",scales='free',ncol=4)  + 
  geom_text(data=df_summary,aes(x=Type,y=1,label=round(Percentage,2)),color="black",vjust = -0.25, size = 6) +
  coord_cartesian(ylim = c(0, 1.1)) + scale_y_continuous(breaks = seq(0, 1, 0.2)) + 
  theme_minimal() + ylab("Counts %") + xlab("Cell type")
plot_text_size=20
plot_text_color="black"
p <- p + theme_bw() + 
  theme(axis.text = element_text(size = plot_text_size,colour = plot_text_color), 
        axis.title = element_text(size = plot_text_size,colour = plot_text_color), 
        text = element_text(size = plot_text_size,colour = plot_text_color),
        plot.title = element_text(size = plot_text_size,colour = plot_text_color),
        legend.title = element_text(size = plot_text_size,colour = plot_text_color), 
        legend.text = element_text(size = plot_text_size,colour = plot_text_color)) 

p
g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
codes <- levels(df$Hi_C)
replacement_vector <- setNames(codes, stript)
replacement_vector
k <- 1
for (i in stript) {
  mark <- as.character(replacement_vector[as.character(i)])
  print(mark)
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <-abi_colors.hex.title[mark]
}

png("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.v2.png",height = 6, width=14,units="in",res=300)
grid.draw(g)
dev.off()


pdf("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.v2.pdf",height = 8, width=16)
grid.draw(g)
dev.off()


ggsave("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.v2.pdf",device="pdf",plot=p,dpi = 200,width = 12,height = 10)
ggsave("IHEC/Analysis/hic_compartments_overlap/HiC_compartments_relabelled.v2.png",device="png",plot=p,dpi = 200,width = 12,height = 10)