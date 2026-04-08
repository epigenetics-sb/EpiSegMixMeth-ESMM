
# set the environment
library(gprofiler2)
library(tidyr)
library(ggplot2)
library(patchwork)
set.seed(123)
setwd("./data/projects/")

# load the environment if needed
# load("IHEC/state_transitions/output/epiregios_overlap/4_4_6/gprofiler.RData")

# read the genes file 
dd <- read.table("IHEC/state_transitions/output/epiregios_overlap/4_4_6/transition_interaction_4_4_6_NBC_GCBC_MBC_sorted_significant.tab",header=F)
length(unique(dd$V4))

gostres_all <- gost(query=unique(dd$V4),
                organism="hsapiens",
                significant = TRUE, 
                correction_method = "fdr" ,
                measure_underrepresentation=F,
                evcodes=TRUE) #,
                #sources=c("GO:BP"))

gostres <- gost(query=unique(dd$V4),
                    organism="hsapiens",
                    significant = TRUE, 
                    correction_method = "fdr" ,
                    measure_underrepresentation=F,
                    evcodes=TRUE,sources=c("GO:BP"))

names(gostres$result)

genes <- read.table("IHEC/state_transitions/output/epiregios_overlap/4_4_6/transition_interaction_4_4_6_NBC_GCBC_MBC_sorted_significant.tab",header=F)

# read the header for the file: 
genes_header <- read.table(gzfile("./resources/references/epiregio/REMAnnotationModelScore_1.csv.gz"),nrows=1,sep=",",header=T)
colnames(genes) <- colnames(genes_header)

# add a column for predicted functions 
genes$predicted_function <- rep("Activating",nrow(genes))
genes$predicted_function[which(genes$regressionCoefficient < 0) ]  <- "Repressing"

# create a named list/dictionary to map the geneID to predicted_function
gene_functions <- c(genes$predicted_function)
names(gene_functions) <- genes$geneID


gene_regression <- data.frame(geneID=unique(genes$geneID))
gene_regression$regression_coef <- lapply(gene_regression$geneID, function(x) {
  a <- c(genes[which(genes$geneID==x),c("regressionCoefficient")])
  return(a)}
  )


# filter top 20 or 30 based on order
temp <- gostres_all$result
temp
temp <- temp[order(temp$p_value,-temp$precision),]
temp
temp <- temp[c(1:20),]
temp$regression_coeff <- unlist(lapply(temp$intersection, function(x) {
  a=strsplit(x,",",fixed=TRUE); 
  b=unlist(a[1]);
  c=unlist(lapply(b, function(y) {
    a=gene_regression$regression_coef[which(gene_regression$geneID==y)]
  }))
  return(a)})) 
lapply(temp$intersection[20], function(x) {
  a=strsplit(x,",",fixed=TRUE); 
  b=unlist(a[1]);
  c=lapply(b, function(y) {
    d=gene_regression$regression_coef[which(gene_regression$geneID==y)]
    return(d)}
    )
  return(c)})

temp$reg_coeff <- lapply(temp$intersection, function(x) {
  a=unlist(strsplit(x,",",fixed=TRUE)[1])
    b=unlist(lapply(a, function(y){
     gene_regression$regression_coef[which(gene_regression$geneID==y)]
    }))
  return(b)
  }
)
temp$length <- unlist(lapply(temp$reg_coeff, function(x) length(unlist(x))))

  df <- data.frame(id=c(1:temp$length[4]),
                   reg_coef =unlist(temp$reg_coeff[4]))
  ggplot(data=df,aes(y=reg_coef,fill = reg_coef > 0)) + 
    geom_density() + 
    theme_minimal() +  
    scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "red"))

df <- map_df(1:nrow(temp), function(i) {
    data.frame(id = temp$term_name[i],
               reg_coef = unlist(temp$reg_coeff[i]))
  })
  
df$type <- "Activating"
df$type[which(df$reg_coef<0)] <- "Repressing"
# cutoff for the regression coefficent
cutoff=0.050
df.filt <- df[which(abs(df$reg_coef) >= cutoff),]

 
  # Plotting
g_plot <- ggplot(data = df.filt, aes(x = abs(reg_coef), group = type,fill=type)) +
    geom_histogram(alpha = 0.75) + xlim(0,1.5) + 
    facet_wrap(~id, scales = "free_x", ncol = 4) +
    theme_minimal() +  scale_fill_manual(values = c("Activating" = "green", "Repressing" = "red")) + 
    labs(x = "Regression Coefficient",
         y = "Frequency",
         fill="Function") 

g_plot <- g_plot + theme(axis.text = element_text(size = 28),
                          axis.title = element_text(size = 28),
                          text = element_text(size = 28)) 
 
g_plot
ggsave(plot=g_plot,file = "./project/Analysis/gene_ontology/gene_ontology_BP_goprofiler_4_4_6.epiregio.png",
        width = 30,
        height = 15,dpi = 200)

# map the geneIDS of intersection size from the gprofiler
# unlist(lapply(gostres_bp$result$intersection, function(x) {
#                                               a=strsplit(x,",",fixed=TRUE); 
#                                               b=unlist(a[1]); 
#                                               gene_list <- c(gene_functions[b]);
#                                               total_a=length(grep("Activating",gene_list));
#                                               total_b=length(grep("Repressing",gene_list));
#                                               total_c = length(gene_list) - (total_a +total_b);
#                                               return(paste(total_a,total_b,total_c))}))
#save.image("./project/Analysis/gene_ontology/Rsession.RData")
gostres_all$result$activating <- unlist(lapply(gostres_all$result$intersection, function(x) {
  a=strsplit(x,",",fixed=TRUE); 
  b=unlist(a[1]);
  c=unlist(lapply(b, function(y) {
    a=length(grep("Activating",genes$predicted_function[which(genes$geneID == y)]))
  }))
  return(sum(c))})) 
gostres_all$result$repressing <- unlist(lapply(gostres_all$result$intersection, function(x) {
  a=strsplit(x,",",fixed=TRUE); 
  b=unlist(a[1]);
  c=unlist(lapply(b, function(y) {
    a=length(grep("Repressing",genes$predicted_function[which(genes$geneID == y)]))
  }))
  return(sum(c))})) 


# create a new dataframe from the selected columns
columns_to_have <- c("term_name","significant","p_value","term_size","activating","repressing","precision","query_size","intersection_size","source")
df <- gostres_all$result[which(gostres_all$result$source == "GO:BP"),]
names(df)

df$minus_log10_adjp <- -log10(df$p_value)
df <- df[order(-df$minus_log10_adjp,-df$precision),]
df <- df[c(1:20),]

df_long <- df %>%
  select(term_name, activating, repressing) %>%
  pivot_longer(cols = c(activating, repressing), names_to = "type", values_to = "count")



r_plot <- ggplot(df,aes(y = term_name, x = intersection_size,color = minus_log10_adjp, size = precision)) +   
                geom_point(alpha = 1.0) +
                theme_classic() + scale_size_continuous(range = c(7, 14)) + 
                labs(x="Genes", y="GO: BP Term \n <- Rank ", size="Gene ratio", 
                        col="-log10(adj. p-value)") + 
                theme(legend.title = element_text(margin=margin(b=15)))  + 
                geom_hline(aes(yintercept = term_name ),linewidth = 0.25,color="gray") 

r_plot <- r_plot + theme(axis.text = element_text(size = 20,colour = "black"),
                         axis.title = element_text(size = 20,colour = "black"),
                         text = element_text(size = 20,colour = "black"))
r_plot
ggsave(plot=r_plot,file = "./project/Analysis/gene_ontology/gene_ontology_BP_goprofiler_4_4_6.only.png",
       width = 15,
       height = 7,dpi = 300)


g_plot <- ggplot(data = df_long, aes(y = term_name, x = count, fill = type)) + 
                  geom_bar( stat = "identity") + 
                  scale_fill_manual(values=c("activating"="darkgreen","repressing"="red")) + 
                  theme_classic() + 
                  labs(x="# Reg. Elements (EpiRegio)", y="GO: BP Term \n Rank <-",fill="EpiRegio predicted funtions")
g_plot <-  g_plot + theme(axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.line.y=element_blank(),
               axis.ticks.y=element_blank()) + 
               geom_hline(aes(yintercept = term_name ),linewidth = 0.25,color="gray")
g_plot <- g_plot + theme(axis.text = element_text(size = 28),
                         axis.title = element_text(size = 28),
                         text = element_text(size = 28)) # changes axis labels text size


png(file = "./project/Analysis/gene_ontology/gene_ontology_BP_goprofiler_4_4_6.v2.png",width = 35,height = 15,res = 200,units = "in")
(r_plot + theme(plot.margin = unit(c(0,30,0,0), "pt"))) + (g_plot) +  plot_layout(guides = 'collect')  
# new_plot <- grid.arrange(r_plot,                             # First row with one plot spaning over 2 columns
#              g_plot, ncol = 2, # Second row with 2 plots in 2 different columns
#              nrow = 1) 
dev.off()

