library(circlize)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(grid)

#read in the fisher file
fisher_summary <-  read.csv("G7_MG_hap_level_fisher_results_withenrichment.csv")
head(fisher_summary)
#hap mat file
hap_mat = readRDS("hap_mat_G7.rds")
# 1) Make sure hap_mat is numeric matrix
hap_bin <- as.matrix(hap_mat)
storage.mode(hap_bin) <- "numeric"
# remove haplotypes with zero variance, if any
keep <- apply(hap_bin, 2, function(x) length(unique(x)) > 1)
hap_bin <- hap_bin[, keep, drop = FALSE]
# 2) Phi coefficient function for binary data
phi_sim <- function(x, y) {
  a <- sum(x == 1 & y == 1)
  b <- sum(x == 1 & y == 0)
  c <- sum(x == 0 & y == 1)
  d <- sum(x == 0 & y == 0)
  
  den <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
  if (den == 0) return(NA_real_)
  (a * d - b * c) / den
}

haps <- colnames(hap_bin)

phi_mat <- matrix(
  NA_real_,
  nrow = length(haps),
  ncol = length(haps),
  dimnames = list(haps, haps)
)

for (i in seq_along(haps)) {
  for (j in seq_along(haps)) {
    phi_mat[i, j] <- phi_sim(hap_bin[, i], hap_bin[, j])
  }
}
diag(phi_mat) <- 1

head(phi_mat)

saveRDS(phi_mat, "Phi_g7.rds")

phi_mat <- readRDS("Phi_g7.rds")

hap_annot <- fisher_summary %>%
  filter(hap %in% haps) %>%
  distinct(hap, .keep_all = TRUE) %>%
  mutate(
    chr = str_extract(hap, "^Gm[0-9]+"),
    chr = gsub("^Gm", "", chr),
    logp = -log10(p_value),
    direction = case_when(
      enriched_phenotype == "R" & p_value < 0.05 ~ "SMV-G7 Resistant",
      enriched_phenotype == "N" & p_value < 0.05 ~ "SMV-G7 Necrotic",
      enriched_phenotype == "S" & p_value < 0.05 ~ "SMV-G7 Susceptible",
      TRUE ~ "NS"
    )
  )

hap_annot <- hap_annot[match(haps, hap_annot$hap), ]

chr_levels <- sort(unique(hap_annot$chr))

muted_chr_cols <- c(
  "#5A8FB8",  # soft blue
  "#D6C98A",  # muted khaki
  "#9C8762",  # soft brown
  "#B07AA1",  # dusty mauve
  "#707353",  # pale pink
  "#4E454F",  # muted rose
  "#D8A24A",  # soft orange
  "#7FB3D5",  # light blue
  "#B9A7D9",  # lavender
  "#6F5A96",  # muted purple
  "#A9744A",  # warm brown
  "#6D869C"   # slate blue
)

if (length(chr_levels) > length(muted_chr_cols)) {
  # fallback if you ever have more chromosomes than the fixed palette covers
  chr_cols <- setNames(
    colorRampPalette(muted_chr_cols)(length(chr_levels)),
    chr_levels
  )
} else {
  chr_cols <- setNames(muted_chr_cols[seq_along(chr_levels)], chr_levels)
}

phi_col <- colorRamp2(
  c(-1, 0, 1),
  c("#2166ac", "#f7f7f7", "#b2182b")  
)

top_ha <- HeatmapAnnotation(
  Chr = hap_annot$chr,
  col = list(
    Chr = chr_cols
  ),
  annotation_name_gp = gpar(fontsize = 9)
)

left_ha <- rowAnnotation(
  Chr = hap_annot$chr,
  col = list(
    Chr = chr_cols
  ),
  annotation_name_gp = gpar(fontsize = 9),
  show_annotation_name = TRUE,
  show_legend = FALSE
)
# 6) Heatmap
ht2 <- Heatmap(
  phi_mat,
  name = "Phi",
  col = phi_col,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = top_ha,
  left_annotation = left_ha,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8.5),
  row_names_gp = gpar(fontsize = 8.5),
  column_names_centered = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.4),
  row_title = NULL,
  heatmap_legend_param = list(
    title = "Phi",
    at = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)
  )
)
ht2
png("PanelA_G7_heatmap_clean6.png", width = 12, height = 10, units = "in", res = 600)
draw(ht2, heatmap_legend_side = "right", annotation_legend_side = "right", 
     merge_legends = TRUE, legend_grouping = "original")


dev.off()
