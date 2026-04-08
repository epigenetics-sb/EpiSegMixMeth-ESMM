# list of packages to be loaded
packages_to_load <- c("dplyr","stringr","ggplot2","reshape2",
                      "scales","jsonlite","grid","patchwork")

# load the packages ! important character.only
lapply(packages_to_load, require, character.only=TRUE)

# reading metadatas
metadata <- read.csv("./data/metadata/IHEC_metadata_harmonization.v1.2_with_coverage_avg_meth.extended.10x_filtered.sorted.csv")
metadata2 <- read.csv("./data/metadata/epiatlas_chipseq_qc_summary.csv",header=TRUE)
metadata3 <- read.table("./data/metadata/ihec_companion_erirrs.txt",header=FALSE)

# color scheme for cell type 
json_color ="./data/metadata/IHEC_EpiATLAS_IA_colors_Mar18_2024.json"

json_data <- fromJSON(json_color,flatten=TRUE)
# convert the rgb color codes to hex codes 
colors.harmonized_sample_ontology_intermediate <- c(json_data$harmonized_sample_ontology_intermediate[[4]])
colors.hex.harmonized_sample_ontology_intermediate <- sapply(colors.harmonized_sample_ontology_intermediate,function(x) {
  rgb_vals <- as.numeric(unlist(strsplit(x, ",")))
  rgb_str <- sprintf("#%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  return(rgb_str)
}) 

colors.experiment <- c(json_data$experiment[[1]])
colors.hex.experiment <- sapply(colors.experiment,function(x) {
  rgb_vals <- as.numeric(unlist(strsplit(x, ",")))
  rgb_str <- sprintf("#%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  return(rgb_str)
}) 


# set the levels order as in IHEC flagship and figure 1 in main figures
metadata$epirr_id_without_version <- factor(metadata$epirr_id_without_version,levels=c(metadata$epirr_id_without_version))

# epirr without version and sorting levels
metadata2$epirr_id_without_version <- sapply(metadata2$epirr_id, function(x) { strsplit(x,'.',fixed=TRUE)[[1]][1]})
metadata2$epirr_id_without_version <- factor(metadata2$epirr_id_without_version , levels=c(metadata$epirr_id_without_version))

# same for metadata3
metadata3$V1 <- factor(metadata3$V1 , levels=c(metadata$epirr_id_without_version))
colnames(metadata3) <- "epirr_id_without_version"

# plotting cell types 

# merge to get cell types of metadata3
df <- merge(metadata3,metadata,by="epirr_id_without_version")


type="cell_type"
if(type=="cell_type"){
  dd <- data.frame(table(df$harmonized_cell_type))
  
  dd$V2 <- sapply(dd$Var1, function(x) unlist(df$harmonized_sample_ontology_intermediate[df$harmonized_cell_type==x])[1])
  dd$Var1 = with(dd, reorder(Var1, -Freq))
  
  dd <- dd %>%
    mutate(V3 = sapply(str_split(Var1, " "), function(words) str_remove(paste(words[1:min(length(words), 3)], collapse=" "), ",\\s*$")))
  
  # levels(dd) <- NULL
  # levels(dd) <- NULL
  dd <- dd %>%
    arrange(V2)
  dd$V5 <- seq(1:nrow(dd))  
  dd$V3 = with(dd, reorder(V3, V5))
  p <- ggplot(data=dd,aes(x=V3,y=Freq,fill=V2)) +
    geom_col() + #scale_y_continuous(breaks = seq(0, 14, by = 2)) + 
    scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + 
    labs(y="Counts",x= "Cell types",fill="Cell types \n (Intermediate) ") + coord_flip() + 
    theme_minimal() + scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2))
  
  p <- p + theme(axis.text = element_text(size = 12)) # changes axis labels
  
  p <- p + theme(axis.title = element_text(size = 12)) # change axis titles
  
  p <- p + theme(text = element_text(size = 12)) 
  
  pplot1 <- p
  
  #ggsave(filename = "./output/analysis/cell_type.v2.png",plot=p,height = 6, width=7, dpi=300)
}
dd <- data.frame(table(df$harmonized_sample_ontology_intermediate))
dd$V2 <- sapply(dd$Var1, function(x) unlist(df$harmonized_sample_ontology_intermediate[df$harmonized_cell_type==x])[1])
dd$Var1 = with(dd, reorder(Var1, -Freq))

dd <- dd %>%
  mutate(V3 = sapply(str_split(Var1, " "), function(words) paste(words[1:min(length(words), 3)], collapse=" ")))

# levels(dd) <- NULL
dd <- dd %>%
  arrange(V2)
dd$V5 <- seq(1:nrow(dd))  
dd$V3 = with(dd, reorder(V3, -V5))

p <- ggplot(data=dd,aes(x=Var1,y=Freq,fill=Var1)) +
  geom_col() + #scale_y_continuous(breaks = seq(0, 14, by = 2)) + 
  scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + 
  labs(y="Counts",x= "Cell types \n (Intermediate)",fill="Cell types \n (Intermediate) ") + coord_flip() + 
  theme_minimal() 

p <- p + theme(axis.text = element_text(size = 12)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 12)) # change axis titles

p <- p + theme(text = element_text(size = 12)) 

pplot2 <- p


combined <- (pplot1 + (pplot2 + guides(fill = "none"))) + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.justification = "center")

combined
#ggsave(filename = paste(out, "cell_type_combined_plot.png", sep="/"),
#       plot = combined, height = 9, width = 16, units = "in", dpi = 300)

#ggsave(filename = "./output/analysis/cell_type.png",plot=p,height = 6, width=7, dpi=300)
#ggsave(filename = "./output/analysis/cell_types.pdf",plot=p,height = 15, width=15)


df <- merge(metadata3,metadata2,by="epirr_id_without_version")
df <- df[,c(1,5,6,8,92,93,94)]
colnames(df) <- c("epirr_id_without_version","antibody","Total_Reads","Mapped_Pct",
                  "Syn_AUC","Syn_Elbow_Pt","Syn_JSD")  
df <- merge(df,metadata,by="epirr_id_without_version")
df <- df[,c(1:8)]
df$epirr_id_without_version <- factor(df$epirr_id_without_version , levels=c(metadata$epirr_id_without_version))
df$antibody <- factor(df$antibody,levels = c("H3K4me3","H3K27ac","H3K4me1","H3K36me3","H3K27me3","H3K9me3"))

melt_df <- reshape2::melt(df,id.vars = c("epirr_id_without_version","harmonized_sample_ontology_intermediate","antibody"))


label_fun <- function(x){
  if (max(x, na.rm = TRUE) >= 1e6) {
    paste0(x / 1e6, "M")
  } else {
    label_number_auto()(x)
  }
}

p <- ggplot(data = melt_df, aes(x = epirr_id_without_version, y = value, fill = harmonized_sample_ontology_intermediate)) +
  geom_col(width = 1, position = "dodge") +
  facet_grid(variable ~ antibody, scales = "free_y") +
  theme_light() +
  labs(x = "EpiRR", fill = "Cell types") +
  scale_fill_manual(values = colors.hex.harmonized_sample_ontology_intermediate) +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  scale_y_continuous(labels = label_fun)
p



p <- p + theme(axis.text = element_text(size = 12)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 12)) # change axis titles

p <- p + theme(text = element_text(size = 12)) 

colors.hex.rand <- c("#367588","#F4C430","#9966CC","#F88379","#708238")

g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
codes <- levels(melt_df$antibody)
replacement_vector <- setNames(levels(melt_df$antibody), stript)
for (i in stript) {
  mark <- as.character(replacement_vector[as.character(i)])
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <-colors.hex.experiment[mark]
}

stripr <- which(grepl('strip-r', g$layout$name))
k <- 1
for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colors.hex.rand[(k)]
  k <- k+1
}

grid.draw(g)

png("./output/analysis/ChIP_QC.v3.png",height = 7, width=10, units = "in",res = 300)
grid.draw(g)
dev.off()
#pdf("./output/analysis/ChIP_QC.v2.pdf",height = 20, width=30)
grid.draw(g)
#dev.off()



# Average methylation and coverage
df <- metadata[,c(1:6)]
df <- merge(metadata3,df,by="epirr_id_without_version")
df$epirr_id_without_version <- factor(df$epirr_id_without_version , levels=c(metadata$epirr_id_without_version))
melt_df <- reshape2::melt(df,id.vars = c("epirr_id_without_version","harmonized_sample_ontology_intermediate"))


for(x in levels(df$epirr_id_without_version)){
  a=df$harmonized_sample_ontology_intermediate[df$epirr_id_without_version==x]
  print(a)
  
}



p <- ggplot(df,aes(x=factor(epirr_id_without_version),y=mean_methylation,fill=harmonized_sample_ontology_intermediate)) +
  geom_col(width = 1) + theme_minimal() + labs(x="EpiRR",y="Average methylation +/- \n Standard deviation",fill="Cell types") + 
  scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + scale_y_continuous(breaks=c(seq(0,100,10))) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_errorbar(aes(ymin=mean_methylation-standard_devation, 
                    ymax=mean_methylation+standard_devation), 
                width=0.15)

p <- p + theme(axis.text = element_text(size =12)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 12)) # change axis titles

p <- p + theme(text = element_text(size = 12)) 

p1 <- p 


out="./output/analysis/"
#ggsave(filename = paste(out, "WGBS_average_methylation.v2.pdf", sep="/"),
#       plot = p, height = 4, width = 6.3, units = "in", dpi = 300)
#ggsave(filename = "./output/analysis/WGBS_average_methylation.pdf",plot=p,height = 10, width=15)


p <- ggplot(df,aes(x=factor(epirr_id_without_version),y=total_CpGs,fill=harmonized_sample_ontology_intermediate)) +
  geom_col(width=1) + theme_minimal() + labs(x="EpiRR",y="Total CpGs",fill="Cell types") + 
  scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6),breaks=c(0,seq(26000000,28000000,500000))) + coord_cartesian(ylim = c(26000000,28000000)) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

p <- p + theme(axis.text = element_text(size = 12)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 12)) # change axis titles

p <- p + theme(text = element_text(size = 12)) 

p2 <- p

#ggsave(filename = "./output/analysis/WGBS_total_CpGs.png",plot=p,height = 10, width=15, dpi=200)
#ggsave(filename = "./output/analysis/WGBS_total_CpGs.pdf",plot=p,height = 10, width=15)


# Coverage 
p <- ggplot(df,aes(x=factor(epirr_id_without_version),y=CG_coverage,fill=harmonized_sample_ontology_intermediate)) +
  geom_col(width=1) + theme_minimal() + labs(x="EpiRR",y="Average coverage \n (Reads/CpG)",fill="Cell types") + 
  scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + 
  scale_y_continuous(breaks=c(0,seq(0,100,10))) + 
  coord_cartesian(ylim = c(0,100)) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

p <- p + theme(axis.text = element_text(size = 12)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 12)) # change axis titles

p <- p + theme(text = element_text(size = 12)) 

p3 <- p
#ggsave(filename = "./output/analysis/WGBS_CpG_coverage.png",plot=p,height = 10, width=15, dpi=200)
#ggsave(filename = "./output/analysis/WGBS_CpG_coverage.pdf",plot=p,height = 10, width=15)

library(patchwork)

combined <- (p1 / (p2 + p3)) +
  plot_layout(guides = "collect", heights = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

combined

ggsave(filename = paste(out, "combined_plot.png", sep="/"),
       plot = combined, height = 8, width = 12.5, units = "in", dpi = 300)
