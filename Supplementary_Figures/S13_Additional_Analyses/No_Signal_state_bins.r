# load all the necessary packages
library(randomForest)
library(clue)
library(ggplot2)
library(reshape2)
library(ggalluvial)
library(ggdendro)
library(networkD3)
library(htmlwidgets)

# set the working directory and environment
options(scipen = 999)
setwd("./resources/src/windows/EpiSegMix/")
getwd()
rm(list = ls())
# read the input file 

df <- read.csv("Classifier/IHEC_ChromHMM_10N_all_states_samples_classifier.csv",header = T)
df <- read.csv("Classifier/IHEC_EpiSegMix_meth_10N_all_states_samples_classifier.csv",header = T)
df <- read.csv("Classifier/IHEC_EpiSegMix_10N_all_states_samples_classifier.csv",header = T)

# getting  samples  epirr
all_samples <- sort(unique(df$EpiRR))


# coloring scheme
json_color ="./project/metadata/IHEC_EpiATLAS_IA_colors_Mar18_2024.json"

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
                "#006400","#FF4500","#66CDAA") 


names(all_colors) <- all_states


# select columns
cols <- c(3:ncol(df))
epirrids <- unique(df$EpiRR)


# extract the training labels
deep_annotation <- read.csv("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_chmm.csv",header = T)
deep_annotation <- read.csv("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix_meth.csv",header = T)
deep_annotation <- read.csv("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.csv",header = T)
ihec.bak <- deep_annotation

df$id <- paste(df$EpiRR,df$states,sep = "_")
deep_annotation$id  <-  paste(deep_annotation$EpiRR,deep_annotation$states,sep = "_")
ihec <- merge(df,deep_annotation,by = "id")
# create id for anotation for merging later
#deep_annotation_chmm$id <- paste(deep_annotation_chmm$EpiRR,deep_annotation_chmm$states,sep = "_")
#deep_annotation$id <- paste(deep_annotation$EpiRR,deep_annotation$states,sep = "_")

epirrs <- ihec$EpiRR.y[which(ihec$corrected=="NS")]
epirrs <- sapply(ihec$EpiRR.y[which(ihec$corrected=="NS")], 
                 function(x) strsplit(x,".",fixed = TRUE)[[1]][1])

# extract the total coverage i.e. bins from counts for ChromHMM
# counts <- read.table("./project/counts/counts.bed") 
# genomic_coverage=sum(counts$V3-counts$V2) 
# genomic_coverage -- 3088262600
# Based on the chunk code to avoid reading files again and again 
genomic_coverage <- 3088262600


ihec$coverage_bp <- (ihec$coverage * genomic_coverage)/100 

ihec_ns <- subset(ihec, corrected == "NS")
ihec_ns$epirr_id_without_version <-  sapply(ihec_ns$EpiRR.y, 
                                         function(x) strsplit(x,".",fixed = TRUE)[[1]][1])

# metadata 
metadata <- read.csv("./project/metadata/IHEC_metadata_harmonization.v1.2_with_coverage_avg_meth.extended.10x_filtered.sorted.csv")
metadata$epirr_id_without_version <- factor(metadata$epirr_id_without_version,levels=c(metadata$epirr_id_without_version))

# set the levels for ihec ns to same as metadata 

df <- merge(ihec_ns,metadata,by="epirr_id_without_version" ,all.x = TRUE)
df$epirr_id_without_version <-  factor(df$epirr_id_without_version,levels=c(metadata$epirr_id_without_version))

df <- df[complete.cases(df), ]

# bins 


# Coverage 
p <- ggplot(data=subset(df,!is.na(harmonized_sample_ontology_intermediate)),
            aes(x=factor(epirr_id_without_version),
                y=coverage_bp,
                fill=harmonized_sample_ontology_intermediate)) +
  geom_col(width=1) + theme_minimal() + labs(x="EpiRR",y="coverage (bp)",fill="Cell types",title = "ChromHMM: No Signal (NS) state genomic proportion") + 
  scale_fill_manual(values=colors.hex.harmonized_sample_ontology_intermediate) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

p <- p + theme(axis.text = element_text(size = 20)) # changes axis labels

p <- p + scale_y_continuous(labels = label_comma(scale = 1e-9, suffix = "Gb", accuracy = 1))  +
  theme(axis.text = element_text(size = 20)) # changes axis labels

p <- p + theme(axis.title = element_text(size = 20)) # change axis titles

p <- p + theme(text = element_text(size = 20)) 

p 

print(paste("Median of NS state in ChromHMM: ",median(df$coverage), "%",sep = ""))
print(paste("Median of NS state in ChromHMM: ",round(median(df$coverage_bp),0), " bps",sep = ""))

# overlap with ESMM states 
dd <- read.table("../../../IHEC/Analysis/state_overlap/EpiSegMixMeth_overlap_ChromHMM-NS.states.txt")
table(dd$V2)

state_labels <- c("NS", "Mix", "Het_C", "Het_F", "ReprPC", "Enh_Wk", "Enh", "Pr", "Tx_Wk", "Tx_Str", "Tx")

dd$V2 <- state_labels[dd$V2]


dd$V2 <- reorder(dd$V2, dd$V1, FUN = mean, decreasing = TRUE)

g <- ggplot(dd, aes(x = V2, y = V1*200,fill = V2 )) +
  geom_boxplot() + 
  scale_fill_manual(values = all_colors) + 
  labs(x="ESMM states",
       y="Overlap with ChromHMM NS (bps)  ",
       fill="States",
       title = "Reclassification of ChromHMM No Signal (NS) state by ESMM ") +
  theme_minimal()

g <- g + scale_y_continuous(labels = label_comma(scale = 1e-9, suffix = "Gbp"))  +
  theme(axis.text = element_text(size = 12)) 

g <- g + theme(axis.title = element_text(size = 12)) # change axis titles

g <- g + theme(text = element_text(size = 12)) 

g 

ggsave(filename = "./project/IHEC/IHEC_companion/NS_overlap/Supplementary_Stacked_Boxplot.png",
       plot=g,height = 4, width=8.5)




# Diagnostic steps 
states <- c("Enh", "Enh_Wk", "Het_C", "Het_F", "Mix", "NS", "Pr", "Pr_Wk",
            "Rep", "ReprPC", "Tx", "Tx_Str", "Tx_Wk")
# Generate custom heatmaps
nstates=length(unique(ihec.bak$corrected))
colors_s <- c("#4a94c2","#cab0d4","#eb4f35","#fa9999","#eb8f45","#739c4f","#a6cce3","#fd9423",
              "#fdbd6e","#eb4245","#b0de8a","#5cb34a","#4c99a6")




test <- read.table("IHEC_Companion/Random_forest_classifier/ChromHMM_NS_overlap.txt",header = F)
colnames(test) <- c("Counts","States")
test <- as.data.frame(cbind(rep("NS",nrow(test)),test$States,test$Counts))
colnames(test) <- c("ChromHMM","EpiSegMix","Counts")
test[2,2] <- "Tx_Str"
ggplot(as.data.frame(test),aes(y=Counts,axis1 = ChromHMM, axis2 = EpiSegMix )) +
  geom_alluvium(aes(fill=as.factor(EpiSegMix)), width= 1/20) +
  geom_stratum(width=1/20) 

scale_x_continuous(breaks = 1:2, labels = c("ChromHMM", "EpiSegMix")) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold")) +
  labs(fill="States")


# Load necessary libraries
library(ggplot2)
library(ggalluvial)

# Create the dataframe
df <- data.frame(
  V1 = c("NS", "NS", "NS", "NS", "NS", "NS", "NS", "NS", "NS", "NS", "NS"),
  V2 = c("Pr", "Str", "NS", "Enh", "Tx_Wk", "ReprPC", "Het_C", "Rep", "Het_F", "Mix"),
  V3 = c(419206800, 1110102400, 1459845800, 1554500000, 2416219600, 2520243400, 2620735600, 2644782400, 2690852600, 2724604200)
)

# Create the Sankey plot

