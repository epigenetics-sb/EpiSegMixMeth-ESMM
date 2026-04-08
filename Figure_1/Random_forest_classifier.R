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
setwd("./project/manuscript/NAR/Scripts/Figure_1/")
rm(list = ls())

# mode configurations 
modes <- list(
  chmm = list(
    data_file       = "IHEC_ChromHMM_10N_all_states_samples_classifier.csv",
    annotation_file = "DEEP_Bcells_trained_EpiRR_annotation_summary_chmm.csv",
    label_col       = "labels",
    feature_cols    = 4:29,
    output_suffix   = "chmm"
  ),
  episegmix_meth = list(
    data_file       = "IHEC_EpiSegMix_meth_10N_all_states_samples_classifier.csv",
    annotation_file = "DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix_meth.csv",
    label_col       = "corrected",
    feature_cols    = 4:29,
    output_suffix   = "episegmix_meth"
  ),
  episegmix = list(
    data_file       = "IHEC_EpiSegMix_10N_all_states_samples_classifier.csv",
    annotation_file = "DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.csv",
    label_col       = "corrected",
    feature_cols    = 4:29,
    output_suffix   = "episegmix"
  )
)

# ── annotation & sample setup (shared across modes) ───────────────────────────
deep_annotation_meta <- read.table("DEEP_samples_annotation_summary.tab", header = TRUE)
DEEP_samples <- sort(unique(deep_annotation_meta$ihec_epirrid[deep_annotation_meta$cohort == "DEEP"]))
DEEP_samples <- c("IHECRE00001913.1", DEEP_samples)

# DEEP and B cells 
train_epirr <- c(DEEP_samples,
                 "IHECRE00000283.2", "IHECRE00000298.2", "IHECRE00000332.2",
                 "IHECRE00000335.2", "IHECRE00000341.2", "IHECRE00000347.2",
                 "IHECRE00000354.2", "IHECRE00001375.1", "IHECRE00001486.1")

states <- c("Enh", "Enh_Wk", "Het_C", "Het_F", "Mix", "NS", "Pr", "Pr_Wk",
            "Rep", "ReprPC", "Tx", "Tx_Str", "Tx_Wk")
colors_s <- c("#4a94c2","#cab0d4","#eb4f35","#fa9999","#eb8f45","#739c4f","#a6cce3","#fd9423",
              "#fdbd6e","#eb4245","#b0de8a","#5cb34a","#4c99a6")
names(colors_s) <- states

out_dir <- "IHEC_Companion/Random_forest_classifier"

# classifier function 
run_classifier <- function(cfg, train_epirr, colors_s, out_dir) {
  suffix <- cfg$output_suffix
  cat("\n=== Running classifier:", suffix, "===\n")

  # normalize 
  ihec <- read.csv(cfg$data_file, header = TRUE)
  all_samples <- sort(unique(ihec$EpiRR))
  cols <- seq(3, ncol(ihec))
  epirrids <- unique(ihec$EpiRR)

  for (epirrid in epirrids) {
    for (x in cols) {
      min_val <- min(ihec[ihec$EpiRR == epirrid, x], na.rm = TRUE)
      max_val <- max(ihec[ihec$EpiRR == epirrid, x], na.rm = TRUE)
      ihec[ihec$EpiRR == epirrid, x] <- (ihec[ihec$EpiRR == epirrid, x] - min_val) / (max_val - min_val)
    }
  }

  # training data 
  train_data <- ihec[ihec$EpiRR %in% train_epirr, ]
  train_data$id <- paste(train_data$EpiRR, train_data$states, sep = "_")
  ihec$id <- paste(ihec$EpiRR, ihec$states, sep = "_")

  deep_annotation <- read.csv(cfg$annotation_file, header = TRUE)
  deep_annotation$id <- paste(deep_annotation$EpiRR, deep_annotation$states, sep = "_")

  train_data <- merge(train_data, deep_annotation[, c("id", cfg$label_col)], by = "id", all.x = TRUE)
  states_labels <- train_data[[cfg$label_col]]

  # train model 
  train_features <- train_data[train_data$EpiRR %in% train_epirr, cfg$feature_cols]
  set.seed(123)
  model <- randomForest(train_features, as.factor(states_labels), ntree = 100, importance = TRUE)
  print(model)

  # --- predict on all samples -------------------------------------------------
  ihec_features <- ihec[, seq(3, ncol(ihec))]
  ihec_bak <- ihec[, c(1, 2)]
  ihec_bak$annotated_probs <- 0
  ihec_bak$corrected_probs <- 0

  predicted_corrected_states <- as.data.frame(
    matrix(NA, nrow = 10, ncol = length(all_samples))
  )
  colnames(predicted_corrected_states) <- as.character(all_samples)

  for (x in seq_along(all_samples)) {
    epirr <- all_samples[x]
    cat(epirr, "\n")
    test_features    <- ihec_features[ihec$EpiRR == epirr, ]
    predicted_labels <- predict(model, test_features)
    predicted_probs  <- randomForest::predict(model, test_features, type = "prob")

    cost_matrix  <- 1 - predicted_probs
    assignment   <- solve_LSAP(cost_matrix)
    label_names  <- colnames(predicted_probs)
    final_labels <- label_names[assignment]

    ihec_bak$annotated[ihec_bak$EpiRR == epirr] <- as.character(predicted_labels)
    ihec_bak$corrected[ihec_bak$EpiRR == epirr] <- final_labels

    ann_indices  <- match(predicted_labels, label_names)
    corr_indices <- match(final_labels, label_names)
    ihec_bak$annotated_probs[ihec$EpiRR == epirr] <- predicted_probs[cbind(seq_along(predicted_labels), ann_indices)]
    ihec_bak$corrected_probs[ihec$EpiRR == epirr] <- predicted_probs[cbind(seq_along(final_labels),    corr_indices)]

    temp <- as.character(predicted_labels)
    for (z in seq_along(temp)) {
      predicted_corrected_states[z, x] <- if (temp[z] != final_labels[z]) {
        paste(temp[z], final_labels[z], sep = "-->")
      } else {
        temp[z]
      }
    }
  }

  # save predictions 
  train_out <- ihec_bak[ihec_bak$EpiRR %in% train_epirr, ]
  # uncommant to write 
  
  # write.csv(train_out, file = file.path(out_dir, paste0("DEEP_Bcells_EpiRR_annotation_summary_", suffix, ".csv")),
  #          quote = FALSE, row.names = FALSE)
  # write.csv(ihec_bak,  file = file.path(out_dir, paste0("DEEP_Bcells_trained_EpiRR_annotation_summary_", suffix, ".csv")),
  #          quote = FALSE, row.names = FALSE)

  # heatmap: all predictions merged with labels
  ihec_bak_id <- ihec_bak
  ihec_bak_id$id <- paste(ihec_bak_id$EpiRR, ihec_bak_id$states, sep = "_")  # states col may not exist; kept for compatibility
  test_merged <- merge(ihec, deep_annotation, by = "id")

  df_matrix  <- as.matrix(test_merged[, c(4, 6, 7, 9:15, 17, 18, 20, 21, 23, 26, 27, 28)])
  rownames(df_matrix) <- paste0(test_merged$states.x, test_merged$EpiRR)
  col_fun  <- colorRamp2(c(0, 1), c("white", "#3F7F93"))
  row_ha   <- rowAnnotation(clusters = test_merged$corrected, col = list(clusters = colors_s))

  # uncomment to save heatmap (all predictions)
  # for (dev_fn in list(
  #   list(fn = png, args = list(filename = file.path(out_dir, paste0("DEEP_Bcells_trained_EpiRR_annotation_summary_", suffix, ".png")))),
  #   list(fn = pdf, args = list(file = file.path(out_dir, paste0("DEEP_Bcells_trained_EpiRR_annotation_summary_", suffix, ".pdf")), width = 10, height = 10)),
  #   list(fn = svg, args = list(filename = file.path(out_dir, paste0("DEEP_Bcells_trained_EpiRR_annotation_summary_", suffix, ".svg")), width = 10, height = 10))
  # )) {
  #   do.call(dev_fn$fn, dev_fn$args)
  #   draw(Heatmap(df_matrix, name = "enrichment", show_row_names = FALSE, col = col_fun,
  #                cluster_columns = FALSE, cluster_rows = TRUE,
  #                right_annotation = row_ha, row_split = as.factor(test_merged$corrected)))
  #   dev.off()
  # }

  # heatmap: training features only 
  df_matrix_train  <- as.matrix(train_data[, cfg$feature_cols])
  rownames(df_matrix_train) <- train_data$id
  col_fun_train    <- colorRamp2(c(0, 1), c("white", "#3F7F93"))
  row_ha_train     <- rowAnnotation(clusters = train_data[[cfg$label_col]])

  # uncomment to save heatmap (training features only)
  # pdf(paste0("DEEP_Bcells_EpiRR_annotation_summary_", suffix, ".pdf"), width = 10, height = 15)
  # draw(Heatmap(df_matrix_train, show_row_names = FALSE, col = col_fun_train,
  #              cluster_columns = TRUE, cluster_rows = TRUE,
  #              right_annotation = row_ha_train,
  #              row_split = as.factor(train_data[[cfg$label_col]])))
  # dev.off()

  # probability boxplot 
  complete_df <- data.frame(
    epirrid         = ihec_bak$EpiRR,
    annotated       = ihec_bak$annotated,
    corrected       = ihec_bak$corrected,
    annotated_probs = as.numeric(ihec_bak$annotated_probs),
    corrected_probs = as.numeric(ihec_bak$corrected_probs)
  )
  complete_df_melt <- melt(complete_df, id.vars = "annotated", measure.vars = "corrected_probs")

  t_plot <- ggplot(complete_df_melt, aes(x = reorder(annotated, as.numeric(value)), y = as.numeric(value))) +
    geom_boxplot(width = 0.7, fill = "#3F7F93") +
    theme_minimal() + coord_flip() +
    labs(x = "States", y = "Average probability") +
    theme(axis.text = element_text(size = 28),
          axis.title = element_text(size = 28),
          text = element_text(size = 24))

  # uncomment to save probability boxplot
  # for (fmt in c("png", "pdf", "svg")) {
  #   ggsave(
  #     filename = file.path(out_dir, paste0("DEEP_Bcells_trained_EpiRR_annotation_summary_", suffix, "_probs.", fmt)),
  #     plot = t_plot, device = fmt, width = 10, height = 10
  #   )
  # }

  invisible(list(model = model, predictions = ihec_bak))
}

# run all three modes 
results <- lapply(modes, run_classifier,
                  train_epirr = train_epirr,
                  colors_s    = colors_s,
                  out_dir     = out_dir)
