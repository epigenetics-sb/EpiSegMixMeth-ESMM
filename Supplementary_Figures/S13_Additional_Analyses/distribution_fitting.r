library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Required order
dist_levels <- c("ZINBI", "ZISICHEL", "ZIBNB", "DPO", "DEL",
                 "ZAPIG", "SICHEL", "WARING", "ZIPIG")

tissue_levels <- c("BlTA", "SFTR", "BlTN", "SFMa", "SFTM4", "BlEM", "BlCM", "CoMu",
                   "LiHe", "BlMa", "BlMo", "LiHR", "LiHG", "midbrain", "heart", "lung",
                   "kidney", "stomach", "liver", "forebrain", "hindbrain")

bin_levels <- c(200, 500, 1000, 1500, 2000)

# ---- Legend labels/order: add "m-" for mouse tissues and place them last ----
mouse_tissues <- c("midbrain", "heart", "lung", "kidney",
                   "stomach", "liver", "forebrain", "hindbrain")

legend_labels <- setNames(tissue_levels, tissue_levels)
legend_labels[mouse_tissues] <- paste0("m-", mouse_tissues)

legend_order <- c(setdiff(tissue_levels, mouse_tissues), mouse_tissues)

# ---- Distinct 21 colors from RColorBrewer ----
tissue_cols <- c(
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set1"),
  brewer.pal(5, "Paired")
)
names(tissue_cols) <- tissue_levels

# ---- Prepare data (include histonemark for faceting) ----
plot_df <- dist_fit %>%
  filter(bins %in% bin_levels, org %in% c("human", "mouse")) %>%
  mutate(
    bins = as.integer(bins),
    org = factor(org, levels = c("human", "mouse")),
    histonemark = factor(trimws(histonemark)),
    distribution = factor(trimws(distribution), levels = dist_levels),
    tissue = factor(trimws(tissue), levels = tissue_levels)
  ) %>%
  count(bins, org, histonemark, distribution, tissue, .drop = FALSE)

# ---- Function: one plot per bin; facets = histonemark x org ----
make_bin_plot <- function(bin_value, proportional = TRUE) {
  dat <- plot_df %>%
    filter(bins == bin_value) %>%
    group_by(org, histonemark, distribution) %>%
    filter(sum(n) > 0) %>%  # avoids warning from all-zero groups in position="fill"
    ungroup()
  
  ggplot(dat, aes(x = distribution, y = n, fill = tissue)) +
    geom_col(position = if (proportional) "fill" else "stack", width = 0.9) +
    facet_grid(histonemark ~ org) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(
      values = tissue_cols,
      breaks = legend_order,
      labels = legend_labels[legend_order],
      drop = FALSE
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    labs(
      title = paste("Bin size:", bin_value),
      x = "Distribution",
      y = if (proportional) "Proportion of Cell types" else "Count of tissue",
      fill = "Tissue"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
}

# ---- Named list of 5 plots ----
plots_by_bin <- setNames(
  lapply(bin_levels, make_bin_plot, proportional = TRUE),
  paste0("bin_", bin_levels)
)

# Example
plots_by_bin$bin_1500

# Print all
for (p in plots_by_bin) print(p)
# for (p in plots_by_bin) print(p)

# Create output folder
#dir.create("./project/EpiSegMix/distribution_fitting/", showWarnings = FALSE, recursive = TRUE)

# Save each plot
for (nm in names(plots_by_bin)) {
  ggsave(
    filename = file.path("./project/EpiSegMix/distribution_fitting/", paste0(nm, ".png")),
    plot = plots_by_bin[[nm]],
    width = 8,
    height = 8,
    dpi = 300
  )
}

