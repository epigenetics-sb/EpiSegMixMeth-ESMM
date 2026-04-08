library(ggplot2)
library(reshape2)
args = commandArgs(trailingOnly=TRUE)
b <- read.table(args[1],header=TRUE) 
file_name=args[2]
output=args[3]
nstates = as.integer(args[4])
histone_marks_number = as.integer(args[5])
b$DNA_methyl <- round(b$Meth/b$Cov,2) 
b$Cov <- NULL
b$Meth <- NULL
data_frame <- b 

data_frame$length <- data_frame$end - data_frame$start
data_frame <- data_frame[!is.na(as.numeric(as.character(data_frame$state))), ]
temp <- data_frame[, c(4:dim(data_frame)[2])]
  
# Below code only used for normalization when we are adding across marks and then penalizing 
# for state length
# Uncomment for the case sum used instead of mean/median when generating input matrix  
#for (x in 2:(2 + histone_marks_number - 1)) {
    #temp[, x] <- log(1 + temp[, x] / temp$length)
#}

# Normalize CpG number for the length of states 
temp$CpG <- log(1 + temp$CpG/ temp$length)
temp$length <- NULL
heatmap <- aggregate(. ~ state, data = temp, mean, na.rm = TRUE)
colnames(heatmap) <- c("states", colnames(heatmap)[2:ncol(heatmap)])

write.table(heatmap,file=paste(output,"_classifier.txt",sep=""),sep="\t",quote=FALSE,row.names=F)
for (x in 2:dim(heatmap)[2]) {
  heatmap[, x] <- (heatmap[, x] - min(heatmap[, x])) / (max(heatmap[, x] - min(heatmap[, x])))
}


heatmap_melt <- melt(data = heatmap, id.vars = "states")
heatmap <- ggplot(heatmap_melt, aes(x = variable, y = as.factor(states), fill = value)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "#3F7F93") +
   # geom_text(aes(label = round(value, 2))) +
   ggtitle(file_name) +
   theme(axis.text = element_text(size = 16), axis.title = element_text(size = 14, face = "bold")) +
   ylab("States") +
   xlab("Marks") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
   scale_y_discrete(limits = rev(levels(as.factor(heatmap_melt$states)))) + scale_x_discrete(limits = c("H3K27me3","H3K9me3","H3K36me3","H3K4me1","H3K27ac","H3K4me3","DNA_methyl","CpG"))
    
ggsave(filename = paste(output,"_classifier.png",sep=""), plot = heatmap, device = "png")
ggsave(filename = paste(output,"_classifier.pdf",sep=""), plot = heatmap, device = "pdf")
  
# Average methylation
data_frame <- subset(b,select=c("state","DNA_methyl"))
colnames(data_frame)=c("state","raw_meth_perc")
data_frame = data_frame[!is.na(as.numeric(as.character(data_frame$state))),] # converting states to numeric so states other then numeric (LMR.UMR,outliers etc) would be converted to NAs and then finally removed 
data_frame = na.omit(data_frame)
methyl_plot = ggplot(data_frame, aes(x=as.factor(state), y=raw_meth_perc, fill="#3F7F93")) + geom_boxplot(fill="#3F7F93") + scale_x_discrete(limits = rev(levels(as.factor(heatmap_melt$states)))) + coord_flip() +
    ggtitle(file_name) + theme(axis.text=element_text(size=18), axis.title=element_text(size=16)) + labs(x= "state", y = "average methylation") +  stat_summary(geom = "point",fun = "mean",col = "black",size = 2,shape = 25,fill = "red")
ggsave(filename=paste(output,"_average_methylation.pdf",sep=""),plot=methyl_plot,device="pdf")
ggsave(filename=paste(output,"_average_methylation.png",sep=""),plot=methyl_plot,device="png")