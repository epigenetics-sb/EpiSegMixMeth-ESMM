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

deep_annotation=read.table("IHEC_Companion/DEEP_samples_annotation_summary.tab",header = T)
ihec <- df

# getting  samples  epirr
DEEP_samples <- sort(unique(deep_annotation$ihec_epirrid[which(deep_annotation$cohort=="DEEP")]))
DEEP_samples <- c("IHECRE00001913.1",DEEP_samples)
all_samples <- sort(unique(df$EpiRR))


# perfrom min-max normalization 

# select columns
cols <- c(3:ncol(ihec))
epirrids <- unique(ihec$EpiRR)

#min-max-normalization
for (epirrid in epirrids){
  for (x in cols ){
    min_val <- min(ihec[ihec$EpiRR == epirrid, x], na.rm = TRUE)
    max_val <- max(ihec[ihec$EpiRR == epirrid, x], na.rm = TRUE)
    ihec[ihec$EpiRR==epirrid,x]=(ihec[ihec$EpiRR==epirrid,x] - min_val)/(max_val-min_val)
  }
}


# select the training epirrs
train_epirr <- DEEP_samples
train_epirr <- c(DEEP_samples[1:20],"IHECRE00000283.2",
                 "IHECRE00000298.2","IHECRE00000332.2",
                 "IHECRE00000335.2","IHECRE00000341.2",
                 "IHECRE00000347.2","IHECRE00000354.2",
                 "IHECRE00001375.1","IHECRE00001486.1")
train_epirr <- c(DEEP_samples,"IHECRE00000283.2",
                 "IHECRE00000298.2","IHECRE00000332.2",
                 "IHECRE00000335.2","IHECRE00000341.2",
                 "IHECRE00000347.2","IHECRE00000354.2",
                 "IHECRE00001375.1","IHECRE00001486.1")
# extract the training data
train_data <- ihec[ihec$EpiRR %in% train_epirr,]




# create an id for merging with the labels 
train_data$id <- paste(train_data$EpiRR,train_data$states,sep = "_")

# extract the training labels
deep_annotation <- read.csv("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_chmm.csv",header = T)
deep_annotation <- read.csv("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix_meth.csv",header = T)
deep_annotation <- read.csv("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.csv",header = T)
ihec.bak <- deep_annotation

ihec$id <- paste(ihec$EpiRR,ihec$states,sep = "_")
ihec.bak$id  <-  paste(ihec.bak$EpiRR,ihec.bak$states,sep = "_")
test <- merge(ihec,ihec.bak,by = "id")
# create id for anotation for merging later
deep_annotation_chmm$id <- paste(deep_annotation_chmm$EpiRR,deep_annotation_chmm$states,sep = "_")
deep_annotation$id <- paste(deep_annotation$EpiRR,deep_annotation$states,sep = "_")

#states_labels <- deep_annotation_chmm$labels[ihec$EpiRR %in% train_epirr]
#states_labels <- ihec.bak$corrected[ihec$EpiRR %in% train_epirr]

# merge the label and train dataframe 
df.merge <- merge(train_data,deep_annotation_chmm[,c("id","labels")],by = "id",all.x = TRUE )
df.merge <- merge(train_data,deep_annotation[,c("id","corrected")],by = "id",all.x = TRUE )

train_data <- df.merge


#deep_annotation_chmm <- train_data
#write.csv(x = deep_annotation_chmm,file = "DEEP_Bcells_EpiRR_annotation_summary_chmm.csv",quote = FALSE,row.names = FALSE)

#set the seed and train the model 

# extract the training features 
train_features <- train_data[train_data$EpiRR %in% train_epirr,c(4:29)]
names(train_features)

train_features <- train_features[,-c(3,5,26)]
names(train_features)

states_labels <-train_data$labels
states_labels <-train_data$corrected

set.seed(123)
# model <- randomForest(train_features, as.factor(states_labels), ntree=100 )
model <- randomForest(train_features, as.factor(states_labels), ntree=100, importance = TRUE )

model


# temp code to test gene expression contribution to the model. 
gene_exp_model <- model

## ------------------------------------------------------ RFC feature contribution -----------------------------------------------

# Extract importance
imp <- importance(model)

# --- Plot 1: Built-in plot (shows both metrics side by side) ---
varImpPlot(model, sort = TRUE, n.var = 23, main = "Feature Importance")

# --- Plot 2: MeanDecreaseAccuracy (ggplot2) ---
imp_acc <- data.frame(
  Feature    = rownames(imp),
  Accuracy   = imp[, "MeanDecreaseAccuracy"],
  Gini       = imp[, "MeanDecreaseGini"]
)
imp_acc <- imp_acc[order(imp_acc$Accuracy, decreasing = TRUE), ]



# Top 20 by MeanDecreaseAccuracy
ggplot(head(imp_acc, 23), aes(x = reorder(Feature, Accuracy), y = Accuracy)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  labs(x = "Feature", y = "Mean Decrease Accuracy", 
       title = "Top 20 Features by Accuracy Importance") +
  theme_minimal()

# --- Plot 3: MeanDecreaseGini (ggplot2) ---
imp_gini <- imp_acc[order(imp_acc$Gini, decreasing = TRUE), ]

d <- ggplot(head(imp_gini, 23), aes(x = reorder(Feature, Gini,decreasing = TRUE), y = Gini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Feature", y = "Mean Decrease Gini", 
       title = "Features by \"Gini\" Importance")
  



text_size <- 8
d <- d + theme_minimal() + 
  theme(axis.text = element_text(size = text_size,color = "black")) +
  theme(axis.title = element_text(size = text_size,color="black")) +
  theme(text = element_text(size = text_size,color = "black")) + theme(axis.title.x = element_blank()) + 
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
d

ggsave("./project/manuscript/NAR/Scripts/Supplementary_Figures/Reviewers_comments/Random_forest_classifier_feature_importance.pdf",
       width = 4, height = 4, plot = d,
       device = cairo_pdf)
## !!!! Make sure to use only feature which are used in the paper 23 features and then save the plot


## -----------------------------------------------------------------------------------------------------
# claude code for shap values 
library(fastshap)
library(shapviz)

# ── 1. Prediction wrapper — matches your binary label setup ──────────────────
pfun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, 2]
}

# ── 2. Compute SHAP values using your existing train_features ────────────────
set.seed(123)
shap_values <- explain(
  object       = model,
  X            = train_features,
  pred_wrapper = pfun,
  nsim         = 100,
  adjust       = TRUE
)

# ── 3. Build shapviz object ──────────────────────────────────────────────────
sv <- shapviz(shap_values, X = train_features)

# ── 4. Beeswarm — equivalent to Python SHAP summary plot ────────────────────
p_bee <- sv_importance(sv, kind = "beeswarm") +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"))

# ── 5. Bar plot — mean |SHAP| per feature (mirrors your imp_acc bar plot) ───
p_bar <- sv_importance(sv, kind = "bar") +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"))

# ── 6. Save ──────────────────────────────────────────────────────────────────
ggsave("path/to/output/SHAP_beeswarm.pdf",
       plot = p_bee, width = 4, height = 4, device = cairo_pdf)

ggsave("path/to/output/SHAP_bar.pdf",
       plot = p_bar, width = 4, height = 4, device = cairo_pdf)


# save the model for prediction on all the test features 
deep_trained <- model 
deep_trained_chmm_meth <- model
deep_trained_episegmix_meth <- gene_exp_model
# predict and save the predictions on other samples/test features 
ihec.features <- ihec[,c(3:ncol(ihec))]
ihec.bak <- ihec[,c(1:2)]
samples_list <- all_samples
ihec.bak$annotated_probs =rep(0,nrow(ihec.bak))
ihec.bak$corrected_probs =rep(0,nrow(ihec.bak))
predicted_corrected_states  <- as.data.frame(matrix(NA, nrow = 10, ncol = length(samples_list)))
colnames(predicted_corrected_states) = c(as.character(samples_list))
for (x in 1:length(samples_list)){
  epirr <- samples_list[x]
  print(epirr)
  test_features <- ihec.features[which(ihec$EpiRR==epirr),]
  predicted_labels <- predict(deep_trained_episegmix_meth, test_features)
  predicted_labels
  temp <- c(as.character(predicted_labels))
  predicted_probs <- predict(deep_trained_episegmix_meth, test_features, type="prob")
  cost_matrix <- 1 - predicted_probs
  assignment <- solve_LSAP(cost_matrix)
  label_names <- colnames(predicted_probs)
  final_labels <- label_names[assignment]
  final_labels
  ihec.bak$annotated[which(ihec.bak$EpiRR==epirr)] <- c(as.character(predicted_labels))
  ihec.bak$corrected[which(ihec.bak$EpiRR==epirr)] <- final_labels
  
  ann_indices <- match(predicted_labels, colnames(predicted_probs))
  corr_indices <- match(final_labels, colnames(predicted_probs))  
  ann_matrix <- cbind(1:length(predicted_labels), ann_indices)
  corr_matrix <- cbind(1:length(final_labels), corr_indices)
  ihec.bak$annotated_probs[which(ihec$EpiRR==epirr)] <- predicted_probs[ann_matrix]
  ihec.bak$corrected_probs[which(ihec$EpiRR==epirr)] <- predicted_probs[corr_matrix]
  
  ihec$corrected_probs 
  for (z in 1:length(temp)){
    if (temp[z] != final_labels[z]){
      predicted_corrected_states[z,x] = paste(temp[z],unlist(final_labels)[z],sep = "-->")
    } 
    else{
      predicted_corrected_states[z,x] = temp[z]
    }
    
  }
  
}


#####  temp code started #############
ihec.bak$id <- paste(ihec.bak$EpiRR,ihec.bak$states,sep="_")


df <- merge(deep_annotation[,c(4,6,7)],ihec.bak[,c(4,6,7)],by="id")

colnames(df) <- c("Id","pGE_probs","pGE","mGE_probs","mGE")

df_melt <- melt(df,ids=c("Id","pGE","mGE"))

table(df[which(df$pGE != df$mGE),c(3,5)])

ggplot(df_melt[which(df_melt$pGE == df_melt$mGE),], aes(x=pGE, y=value,fill=variable)) + 
  geom_boxplot(outlier.colour="red",
               outlier.size=2)

same_prediction <- length(df$Id[which(df$pGE == df$mGE)])
diff_prediction <- length(df$Id[which(df$pGE != df$mGE)])


print(paste("Number of states predicted same as with or without gene expresssion: ", same_prediction))

print(paste("Number of states predicted different as with and without gene expresssion: ",diff_prediction))


#####  temp code finished #############
train <- ihec.bak[ihec.bak$EpiRR %in% train_epirr,]
#save the predictions in  a file
#write.csv(x =train,file = "IHEC_Companion/Random_forest_classifier/DEEP_Bcells_EpiRR_annotation_summary_chmm.csv",quote = FALSE,row.names = FALSE)
#write.csv(x =ihec.bak,file = "IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_chmm.csv",quote = FALSE,row.names = FALSE)




# Diagnostic steps 
states <- c("Enh", "Enh_Wk", "Het_C", "Het_F", "Mix", "NS", "Pr", "Pr_Wk",
            "Rep", "ReprPC", "Tx", "Tx_Str", "Tx_Wk")
# Generate custom heatmaps
nstates=length(unique(ihec.bak$corrected))
colors_s <- c("#4a94c2","#cab0d4","#eb4f35","#fa9999","#eb8f45","#739c4f","#a6cce3","#fd9423",
              "#fdbd6e","#eb4245","#b0de8a","#5cb34a","#4c99a6")
              
names(colors_s) <- states

# entire dataframe after prediction
df_matrix <- as.matrix(test[,c(4,6,7,9:15,17,18,20,21,23,26,27,28)])
rownames(df_matrix) = paste0(test$states.x,test$EpiRR)
row_ha = rowAnnotation(clusters=test$corrected,col=list(clusters=colors_s))
row_ha2 = rowAnnotation(states=as.factor(test$states))
col_fun = colorRamp2(c(0, 1), c("white", "#3F7F93"))
Heatmap(df_matrix,name="enrichment",show_row_names = FALSE,col = col_fun,cluster_columns = FALSE ,cluster_rows = TRUE,right_annotation = row_ha,row_split = c(as.factor(test$corrected)))
dev.off()
# save the plot
png("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.png")
pdf("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.pdf",width = 10,height =10)
svg("IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.svg",width = 10,height =10)
Heatmap(df_matrix,name="enrichment",show_row_names = FALSE,col = col_fun,cluster_columns = FALSE ,cluster_rows = TRUE,right_annotation = row_ha,row_split = c(as.factor(ihec.bak$corrected)))
dev.off()


# train features only 
df_matrix <- as.matrix(train_data[,4:29])
rownames(df_matrix) = paste0(train_data$id)
row_ha = rowAnnotation(clusters=train_data$corrected)
row_ha2 = rowAnnotation(states=as.factor(train_data$states))
col_fun = colorRamp2(c(0, 1), c("white", "#3F7F93"))
Heatmap(df_matrix,name="enrichment",show_row_names = FALSE,col = col_fun,cluster_columns = FALSE ,cluster_rows = TRUE,right_annotation = row_ha,row_split = c(as.factor(train_data$corrected)))


# save the plot
pdf("DEEP_Bcells_EpiRR_annotation_summary_chmm.pdf",width = 10,height = 15)
Heatmap(df_matrix,show_row_names = FALSE,col = col_fun,cluster_columns = TRUE ,cluster_rows = TRUE,right_annotation = row_ha,row_split = c(as.factor(train_data$labels)))
dev.off()


# check the probabilites for the annotated states: 
# plot for probabilites 
complete_df <- as.data.frame(cbind(ihec.bak$EpiRR,ihec.bak$annotated,ihec.bak$corrected,ihec.bak$annotated_probs,ihec.bak$corrected_probs))
colnames(complete_df) <- c("epirrid","annotated","corrected","annotated_probs","corrected_probs")
complete_df_melt <- melt(data = complete_df,id.vars = "annotated",measure.vars = "corrected_probs")
t_plot <- ggplot(complete_df_melt,aes(x=reorder(annotated,as.numeric(value)),y=as.numeric(value))) + geom_boxplot(width=0.7,fill="#3F7F93") + theme_minimal() + coord_flip() + labs(x= "States", y = "Average probability")
t_plot <- t_plot + theme(axis.text = element_text(size=28))
t_plot <- t_plot + theme(axis.title = element_text(size=28))
t_plot <- t_plot + theme(text = element_text(size=24))
t_plot

ggsave(filename = "IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_chmm_probs.png",plot = t_plot,device = "png",width = 10, height = 10) 
ggsave(filename = "IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_chmm_probs.pdf",plot = t_plot,device = "pdf",width = 10, height = 10) 
ggsave(filename = "IHEC_Companion/Random_forest_classifier/DEEP_Bcells_trained_EpiRR_annotation_summary_chmm_probs.svg",plot = t_plot,device = "svg",width = 10, height = 10) 


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

