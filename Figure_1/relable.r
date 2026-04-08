args = commandArgs(trailingOnly=TRUE)
names(args) <- c("file_to_read","output","annotation","epirrid","tag")
print(args)

if (args["tag"]=="ChromHMM") {
   file=read.table(args["file_to_read"],header=F,skip=1)
} else {
   file=read.table(args["file_to_read"],header=F)
}

deep_annotation = read.csv(args["annotation"],header=T)


#labels <- deep_annotation$corrected[which(deep_annotation$ihec_epirrid==epirrid)]
#temp <- deep_annotation[which(deep_annotation$EpiRR==epirrid),c(2,20)]
temp <- deep_annotation[which(deep_annotation$EpiRR==args["epirrid"]),c(2,6)]
#replacement_vector <- setNames(temp$clusters, temp$states)
replacement_vector <- setNames(temp$corrected, temp$states)
file$V4 <- replacement_vector[as.character(file$V4)]
file$V4[which(file$V4=="Mix")] <- "Enh_Wk"
#file$V4 <- factor(file$V4, levels = 1:10, labels = labels)
states <- c("NS","Mix","Het_C",
            "Het_F","ReprPC","Enh_Wk",
            "Enh","Pr","Tx_Wk","Tx_Str",
            "Tx","Rep","Pr_Wk")

colors <- c("255,255,255","138,145,208","146,197,222",
            "17,193,255","128,128,128","255,255,0",
            "255,195,77","255,0,0","0,100,0","0,128,0",
            "0,100,0","102,205,170","255,69,0")

names(colors) <- states
file$V9 <- colors[as.character(file$V4)]
# add colors coding as well later to have uniform colors
write.table(file,paste(args["output"],"/",args["epirrid"],"_",args["tag"],"_10N_relabelled.bed",sep=""),quote=FALSE,row.names=FALSE,sep="\t",col.names = FALSE)