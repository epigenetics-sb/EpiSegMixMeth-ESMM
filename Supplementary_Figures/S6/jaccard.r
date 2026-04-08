#loading necessary library

library("ggExtra")
library("patchwork")

# Comparing the Het-F as well as other states jaccard index of bad quality chip to states from similar cell type of better quality 
epi_09_13 <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/EpiSegMix_1909-1913.csv",header=FALSE)
epi_09_16 <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/EpiSegMix_1909-1916.csv",header=FALSE)
epi_13_16 <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/EpiSegMix_1913-1916.csv",header=FALSE)
epim_09_13 <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/EpiSegMixMeth_1909-1913.csv",header=FALSE)
epim_09_16 <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/EpiSegMixMeth_1909-1916.csv",header=FALSE)
epim_13_16 <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/EpiSegMixMeth_1913-1916.csv",header=FALSE)

# reading the summary file
epim <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/luminal_epithelial_cell_of_mammary_gland/jaccard_index/EpiSegMixMeth_ji.csv",header=FALSE)
epi <- read.csv("./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/luminal_epithelial_cell_of_mammary_gland/jaccard_index/EpiSegMix_ji.csv",header=FALSE)

# adding model information
epim$V5 <- rep("EpiSegMixMeth",nrow(epim))
epi$V5 <- rep("EpiSegMix",nrow(epi))


# changing the levels order to fit the upper traingular matrix
df$V2 <- factor(df$V2,levels = c("IHECRE00004616.1","IHECRE00004617.2","IHECRE00004614.1","IHECRE00004607.2",
                                   "IHECRE00004606.2","IHECRE00004603.3",
                                   "IHECRE00000223.3","IHECRE00000231.3"))

df$V1 <- factor(df$V1,levels = c("IHECRE00004616.1","IHECRE00004617.2","IHECRE00004614.1","IHECRE00004607.2",
                                   "IHECRE00004606.2","IHECRE00004603.3",
                                   "IHECRE00000231.3","IHECRE00000223.3"))

# combine the dataframe
df <- rbind(epim,epi)

# change the factor level order for episegmixmeth to come first
df$V5 <- factor(df$V5,levels=c("EpiSegMixMeth","EpiSegMix"))

# calulating difference between the EpiSegMixMeth and EpiSegMix / state
dd <- data.frame()
for (x in 1:280) {
    a=as.character(df[x,1])
    b=as.character(df[x,2])
    c=as.character(df[x,3])
    d=round(df$V4[which(df$V5=="EpiSegMix" & df$V1==a & df$V2==b & df$V3==c)],2)
    e=round(df[x,4],2)
    rchange=e-d
    dd[x,1] <- a
    dd[x,2] <- b
    dd[x,3] <- c
    dd[x,4] <- as.numeric(rchange)
    print(c(a,b,c,rchange))
  }

# change the factor order of staes to match figure 3 of main figures
dd$V3 <- factor(dd$V3,levels = c("NS","Het_C","Het_F","ReprPC","Mix",
                                 "Enh_Wk","Enh","Pr","Tx_Wk","Tx_Str"))
dd$V2 <- factor(dd$V2,levels = c("IHECRE00004616.1","IHECRE00004617.2","IHECRE00004614.1","IHECRE00004607.2",
                                 "IHECRE00004606.2","IHECRE00004603.3",
                                 "IHECRE00000223.3","IHECRE00000231.3"))

dd$V1 <- factor(dd$V1,levels = c("IHECRE00004616.1","IHECRE00004617.2","IHECRE00004614.1","IHECRE00004607.2",
                                 "IHECRE00004606.2","IHECRE00004603.3",
                                 "IHECRE00000231.3","IHECRE00000223.3"))

# plotting the differences
p <- ggplot(dd, aes(x=V1, y=V2,fill=V4)) + 
     geom_tile() + facet_grid(~V3) + theme_minimal() + coord_fixed(ratio=1) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        limits = c(-.7,.7)) # geom_point(shape=21,size=4)

p <- p +   labs(y="Source EpiRR",x="Target EpiRR", fill=expression(paste(Delta,JI[EpiSegMixMeth-EpiSegMix]))) 

p <- p + theme(axis.text = element_text(size = 20)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 20)) # change axis titles

p <- p + theme(text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p <- p + theme(strip.text = element_text(size = 20)) 

p

  
j <- ggplot(df, aes(x=V3, y=V4,fill = V5)) + 
  geom_boxplot(width=0.5) + xlim(c("NS","Het_C","Het_F","ReprPC","Mix","Enh_Wk","Enh","Pr","Tx_Wk","Tx_Str")) + 
  theme_minimal() +   xlab("States") + ylab("Jaccard Index \n Similarity") + labs(fill="Model")

j <- j + theme(axis.text = element_text(size = 20)) # changes axis labels

j <- j + theme(axis.title = element_text(size = 20)) # change axis titles

j <- j + theme(text = element_text(size = 20)) + theme(axis.title.x=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank())



pdf(file = "./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/luminal_epithelial_cell_of_mammary_gland/jaccard_index/JI_compared.pdf",width = 25,height = 10)
(j | p) +
  plot_layout(ncol = 1, nrow = 2, heights = c(0.75, 1))

dev.off()
png(file = "./project/Analysis/jaccard_index_full_reduced_bad_qualityChIP_cell_types/luminal_epithelial_cell_of_mammary_gland/jaccard_index/JI_compared.png",width = 25,height = 10,res = 200,units = "in")
(j | p) +
  plot_layout(ncol = 1, nrow = 2, heights = c(0.75, 1))

dev.off()


find_max_JI <- function(state,df){
  colnames(df) <- c("source","target","JI")
  # subset the data for particular state
  subset_data <- df[which(df$source==state),]
  
  # maximum jaccard index state
  max_ji <- max(subset_data$JI)
  
  # Step 3: Subset to get the row where JI is maximum
  result_row <- subset_data[subset_data$JI == max_ji,]
  
  #print(which(df$source==state & df$JI == max_ji))
}
