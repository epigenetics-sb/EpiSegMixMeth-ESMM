library(readr)

# Read the dataframes ChromHMM-NS and EpiSegMixMeth all relabelled 200bp bins
setwd("./project/Analysis/state_overlap/")

# ESMM files that are relabelled and merged for all samples with state information
df <- read_table("EpiSegMixMeth/EpiSegMixMeth_merged_10N_relabelled.tab")

# All ChromHMM-NS states/bins across samples 
dd <- read_table("../relabelled_bed_files/random_forest_classifier/ChromHMM_NS.bed",col_names = FALSE)

df <- read.table("EpiSegMixMeth_overlap_ChromHMM-NS.states.txt",header=FALSE)
df$V1 <- df$V1*200

df
ggplot(df,aes(x=factor(V2),y=V1)) + geom_boxplot() + coord_flip() + theme_minimal() +
  scale_y_continuous(labels = function(x) paste0(x / 1000000000, "GB"))

dd <- aggregate(V1~V2,data=df,mean,na.rm=TRUE)
# Assuming dd is your dataframe and V2 is the column you're working with
dd$V2 <- factor(dd$V2,
                levels = 1:10,
                labels = c("NS", "Mix", "Het_C", "Het_F", "ReprPC", 
                           "Enh_Wk", "Enh", "Pr", "Tx_Wk", "Tx_Str"))
# Assuming dd is your dataframe and you want to sort it by columns 'Column1' and 'Column2'
# in decreasing order
dd <- dd[order(-dd$V1), ]
dd$p <- (dd$V1 - min(dd$V1)) / (max(dd$V1)- min(dd$V1))

# create id of bin coordinates 
dd$id <- paste(dd$X1,dd$X2,dd$X3,sep = "_")
df$id <- paste(df$Chr,df$Start,df$End,sep = "_")
head(dd$id)

# merge the dataframe to filter bins only for ChromHMM NS state
df_ns <- merge(dd,df,by="id")

# filter for column containing only states
ds <- df_ns[,c(8:ncol(df_ns))]

# remove all other dataframes  
rm(dd,df,df_ns)

# write the dataframe to a file and refer to the file in the future for any downstream analysis
write.table(x=ds,file="EpiSegMixMeth_overlap_ChromHMM-NS.txt",quote=FALSE,row.names=FALSE,sep="\t",col.names=TRUE)

# create an empty dataframe to store overlap counts
df <- data.frame()
df[c(1:10),1] <- table(ds$IHECRE00000004)

# for each sample calulate the counts of each state 
for(i in 1:2){
  counts <- table(ds[,i])
  print(counts)
}




# read files with matching pattern after moving the directory
files <- list.files(pattern = "\\.txt$")

# get all the epirr without version for all the files in the list (in the folder)
epirrs <- sapply(files, function(x){a=strsplit(x,'.',fixed=TRUE)[[1]][1];b=strsplit(a,'_',fixed=TRUE)[[1]][2];return(b)})
#  epirrs <- sapply(files, function(x){a=strsplit(x,'_',fixed=TRUE)[[1]][3];b=strsplit(a,".",fixed=TRUE)[[1]][1];return(b)})

# remove the names of the list since it is redudant information
names(epirrs) <- NULL

# add name to the list, perhaps not needed
names(files) <- epirrs

# read all the files which are in the list and save it as list
# works well when working with only one column otherwise columnames have to be given before merging
files_read = lapply(files, function(x){read.table(x,as.is=FALSE)})

# merge all the files columnwise from the list using do.call function
# do.call function is quite fast and powerful 
merged_df <- as.data.frame(do.call(cbind,files_read))

# set the column names since the order of the list and the merging by do.call is same as in the list
names(merged_df) <- epirrs

# read the counts bed files to add the coordinate for the states 
df <- read.table("../../counts/counts.bed",header=F)

# add the column names to the recently read dataframe
names(df) <- c("Chr","Start","End")

# merged the dataframe 
df_mod <- cbind(df,merged_df)

# write the dataframe to a file and refer to the file in the future for any downstream analysis
write.table(x=df_mod,file="EpiSegMixMeth_merged_10N_ov.tab",quote=FALSE,row.names=FALSE,sep="\t",col.names=TRUE)
savehistory("EpiSegMixMeth_merged_10N_ov_history.txt")