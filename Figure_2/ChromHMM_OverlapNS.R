library(ggrepel)
library(ggalluvial)
library(scales)
library(cowplot)
setwd("./project/Analysis/state_overlap/")
overlap <- read.csv("windows_laptop/EpiSegMixMeth_ChromHMM_NS-overlap_counts.csv",header = T)
names(overlap) <- c("ChromHMM","EpiSegMixMeth","Counts")
overlap$ChromHMM <- c(rep("ChromHMM-NS",10))
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
                "Tx","Pr_Wk","Rep","ChromHMM-NS") # replace "Other" <- "Mix"
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
                "#006400","#FF4500","#66CDAA","#EAEAEA") # "#0000CC" <- Het_F

names(all_colors) <- all_states
#overlap$EpiSegMixMeth <- factor(overlap$EpiSegMixMeth , 
#                        levels=c("Tx_Str","Tx_Wk","Pr","Enh","Enh_Wk",
#                                 "ReprPC","Mix","Het_F","Het_C","NS"))
overlap$EpiSegMixMeth <- factor(overlap$EpiSegMixMeth,levels=c(overlap$EpiSegMixMeth[order(overlap$Counts)]))

states <- c(levels(factor(overlap$ChromHMM)),levels(factor(overlap$EpiSegMixMeth)))

overlap$Counts <- as.numeric(overlap$Counts)

overlap$Percent <- round((overlap$Counts/sum(overlap$Counts))*100,2)



d <- ggplot(overlap[order(overlap$Counts),],aes(axis1=ChromHMM,axis2=EpiSegMixMeth,y=Percent)) +
  scale_x_discrete(limits=c("ChromHMM","EpiSegMixMeth"),expand=c(0.0,0.0)) + 
  xlab("Model") + ylab("Overlap %") + labs(fill="States") + 
  geom_alluvium(aes(fill=EpiSegMixMeth),width = 0.2,alpha=0.8,reverse=FALSE,curve_type = "quintic")  + 
  geom_stratum(alpha=0.9,fill=c(all_colors[states]),width = 0.2,reverse=FALSE) + 
  theme_minimal() + scale_fill_manual(values = all_colors) #+ ggtitle("ChromHMM - NS distribution") + 
  # geom_text_repel(aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),stat = "stratum", size = 5, direction = "y", nudge_x = 0) 
  #+  geom_text_repel(aes(label = ifelse(after_stat(x)  == 2, as.character(after_stat(stratum)), "")), stat = "stratum", size = 5, direction = "y", nudge_x = .2,reverse) 
text_size=34

d <- d + theme_minimal_grid() + 
  theme(axis.text = element_text(size = text_size,color = "black")) +
  theme(axis.title = element_text(size = text_size,color="black")) +
  theme(text = element_text(size = text_size,color = "black")) + theme(axis.title.x = element_blank())
d
# d <- d + coord_flip()

ggsave("windows_laptop/ChromHMM_NS.png",width = 8,height = 8,plot=d,dpi = 300) 
#scale_y_continuous(labels = unit_format(unit = "Mbp", scale = 1e-6)) + 

reorder(overlap,decreasing = overlap$Counts)
overlap[order(overlap$Counts),]

ggplot(overlap[order(overlap$Counts),],aes(axis1=ChromHMM,axis2=EpiSegMixMeth,y=percent)) + 
  scale_x_discrete(limits=c("ChromHMM","EpiSegMixMeth"),expand=c(0.2,0.5))  + 
  xlab("Model") + ylab("Overlap percentage") + 
  geom_alluvium(aes(fill=EpiSegMixMeth),width = 0.2,alpha=0.8)  
  geom_stratum(alpha=0.9,fill=c(rev(all_colors[states])),width = 0.2) +
  theme_minimal()  + 
  ggtitle("ChromHMM - NS distribution") + 
  scale_fill_manual(values = all_colors) + 
  geom_text_repel(aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),stat = "stratum", size = 5, direction = "y", nudge_x = 0) + 
  geom_text_repel(aes(label = ifelse(after_stat(x)  == 2, as.character(after_stat(stratum)), "")), stat = "stratum", size = 5, direction = "y", nudge_x = .2) # + 
