library(readr)
library(patchwork)
library(jsonlite)

# ── Paths ─────────────────────────────────────────────────────────────────────
path_overlap  <- "./project/IHEC/IHEC_companion/NS_overlap/ChromHMM_NS_overlap_withESMM_states.csv"
path_metadata <- "./project/metadata/IHEC_metadata_harmonization.v1.2_with_coverage_avg_meth.extended.10x_filtered.sorted.csv"
path_colors   <- "./project/metadata/IHEC_EpiATLAS_IA_colors_Mar18_2024.json"

setwd("./project/Analysis/state_overlap/")

# ── Helper: RGB string → hex ───────────────────────────────────────────────
rgb_to_hex <- function(rgb_str) {
  vals <- as.numeric(unlist(strsplit(rgb_str, ",")))
  sprintf("#%02X%02X%02X", vals[1], vals[2], vals[3])
}

# ── Color schemes ──────────────────────────────────────────────────────────
all_states <- c("NS", "Mix", "Het_C", "Het_F", "ReprPC",
                "Enh_Wk", "Enh", "Pr", "Tx_Wk", "Tx_Str",
                "Tx", "Pr_Wk", "Rep")
all_colors <- setNames(
  c("#EAEAEA", "#8A91D0", "#92C5DE", "#11C1FF", "#808080",
    "#FFFF00", "#FFC34D", "#FF0000", "#006400", "#008000",
    "#006400", "#FF4500", "#66CDAA"),
  all_states
)

json_data <- fromJSON(path_colors, flatten = TRUE)
colors_ontology <- sapply(json_data$harmonized_sample_ontology_intermediate[[4]], rgb_to_hex)
colors_experiment <- sapply(json_data$experiment[[1]], rgb_to_hex)

# ── Load data ──────────────────────────────────────────────────────────────
metadata <- read.csv(path_metadata)
metadata$epirr_id_without_version <- factor(
  metadata$epirr_id_without_version,
  levels = metadata$epirr_id_without_version
)

df <- read.table(path_overlap, sep = ",", header = FALSE)
df$epirr_id_without_version <- factor(
  sapply(df$V1, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1]),
  levels = levels(metadata$epirr_id_without_version)
)
df$V4 <- factor(df$V4, levels = all_states)
df <- df[complete.cases(df), ]

head(df)


# ── Annotation data (NS rows only, merged with metadata) ───────────────────
ihec <- merge(
  df[df$V4 == "NS", ],
  metadata,
  by  = "epirr_id_without_version",
  all.x = TRUE
)
ihec$epirr_id_without_version <- factor(
  ihec$epirr_id_without_version,
  levels = levels(metadata$epirr_id_without_version)
)

# Legend order = left-to-right appearance in plot
ihec$harmonized_sample_ontology_intermediate <- factor(
  ihec$harmonized_sample_ontology_intermediate,
  levels = unique(ihec$harmonized_sample_ontology_intermediate[
    order(ihec$epirr_id_without_version)
  ])
)

# ── Plots ──────────────────────────────────────────────────────────────────
p_anno <- ggplot(ihec, aes(x = epirr_id_without_version, y = 1,
                           fill = harmonized_sample_ontology_intermediate)) +
  geom_tile() +
  scale_fill_manual(values = colors_ontology, name = "Cell type") +
  theme_void() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.text     = element_text(size = 8)
  )

p_main <- ggplot(df, aes(x = V1, y = V5, fill = V4)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(labels = function(x) paste0(x / 1e9, " Gb")) +
  labs(x = "Sample", y = "Base pairs", fill = "State") +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

combined <- p_anno / p_main + plot_layout(heights = c(0.25, 12))
 
outfile = "./project/IHEC/IHEC_companion/NS_overlap/Supplementary_Stacked_Barplot.png"
ggsave(outfile,
       plot   = combined,
       width  = 250,
       height = 100,
       units  = "mm",
       dpi    = 300)


# sankey plot 
temp <- subset(df,select = c("V2","V3","V4","V5"))
head(temp)

sum(tapply(temp$V5, temp$V4, median),na.rm = TRUE)
overlap <- data.frame("Counts"=tapply(temp$V5, temp$V4, median))
overlap$EpiSegMixMeth <- rownames(overlap)
overlap$ChromHMM <- "NS"
overlap <- overlap[complete.cases(overlap), ]

overlap$EpiSegMixMeth <- factor(overlap$EpiSegMixMeth,levels=c(overlap$EpiSegMixMeth[order(overlap$Counts)]))
states <- c(levels(factor(overlap$ChromHMM)),levels(factor(overlap$EpiSegMixMeth)))


# ── ChromHMM collapsed to single "ChromHMM-NS" stratum ────────────────────
overlap$ChromHMM <- factor("ChromHMM-NS", levels = "ChromHMM-NS")

# ── EpiSegMixMeth factor: ascending so largest renders at top ──────────────
overlap$EpiSegMixMeth <- factor(
  overlap$EpiSegMixMeth,
  levels = overlap$EpiSegMixMeth[order(overlap$Counts)]
)

# ── Extend all_colors to include ChromHMM-NS mapped to NS color ───────────
all_colors["ChromHMM-NS"] <- all_colors["NS"]

# ── Legend: EpiSegMixMeth only, descending ─────────────────────────────────
legend_levels <- rev(levels(overlap$EpiSegMixMeth))

d <- ggplot(overlap[order(overlap$Counts), ],
            aes(axis1 = ChromHMM, axis2 = EpiSegMixMeth, y = Counts)) +
  scale_x_discrete(limits = c("ChromHMM", "EpiSegMixMeth"), expand = c(0.0, 0.0)) +
  ylab("Counts (median)") + labs(fill = "States") +
  geom_alluvium(aes(fill = EpiSegMixMeth), width = 0.2, alpha = 0.8,
                reverse = FALSE, curve_type = "quintic") +
  geom_stratum(aes(fill = after_stat(stratum)),
               alpha = 0.9, width = 0.2, reverse = FALSE) +
  scale_fill_manual(
    values  = all_colors,
    breaks  = legend_levels    # ChromHMM-NS excluded from legend
  ) +
  scale_y_continuous(labels = label_comma(scale = 1e-9, suffix = "Gb"), limits = c(0, 2e9))

text_size <- 34

d <- d +
  theme_minimal_grid() +
  theme(
    axis.text    = element_text(size = text_size, color = "black"),
    axis.title   = element_text(size = text_size, color = "black"),
    text         = element_text(size = text_size, color = "black"),
    axis.title.x = element_blank()
  )
d

outfile = "./project/IHEC/IHEC_companion/NS_overlap/ChromHMM_NS.png"

ggsave(outfile,width = 8.5,height = 8,plot=d,dpi = 300) 
ggsave(filename = "./project/Analysis/metadata_plots/cell_type.v2.png",plot=p,height = 6, width=7, dpi=300)
# fix the y axis for GB 