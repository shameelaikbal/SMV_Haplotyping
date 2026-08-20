#from the DE-selected MGs, identify haplotypes formed by crosshap, and the significant ones among them associated with the differential responses of soybean to soybean mosaic virus
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)
base_dir <- getwd()
de_dir   <- file.path(base_dir, "DE_selected_MGs/")

mg_stats <- readRDS("G7_mg_stats.rds")

newde_mg <- read.csv(
  "G7_DE_selected_MG.csv",
  stringsAsFactors = FALSE
)

de_mg_strings <- unique(newde_mg$MG_location)
length(de_mg_strings)
#ok now i need to pull out the ind var hap n fig of all the above 42 MGs only

crosshap_dir <- file.path(getwd(), "./G7_crosshap_results_DE_RNS/")
output_dir   <- file.path(getwd(), "./DE_selected_MGs/")

esc <- function(x) gsub("\\.", "\\\\.", as.character(x))

copy_for_one_mg <- function(target_mg_row) {
  row <- mg_stats[target_mg_row, ]
  
  mg_id_string <- row$MG_ID
  parts <- strsplit(mg_id_string, "_")[[1]]
  
  chr_num <- sub("^Gm", "", parts[1])
  pos_num <- parts[2]
  
  mgmin_val <- row$MGmin
  eps_val   <- row$epsilon
  
  pattern <- paste0(
    "^G7_RNS_viz_chr", sprintf("%02d", as.integer(chr_num)),
    "_pos", pos_num,
    "_mgmin", mgmin_val,
    "_eps", eps_val,
    ".*$"
  )
  
  matches <- list.files(
    path = crosshap_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  
  print(pattern)
  print(matches)
  
  if (length(matches) == 0) {
    warning(paste("No files matched for target_mg =", target_mg_row))
  } else {
    file.copy(matches, output_dir, overwrite = TRUE)
  }
}
#now connect it to my de_mg file
mg_rows_to_do <- unique(newde_mg$MG_ID)

for (target_mg in mg_rows_to_do) {
  copy_for_one_mg(target_mg_row = target_mg)
}
# now need to the hap presence absence
all_files <- list.files(de_dir, full.names = TRUE)

all_fisher_results <- list()

for (i in seq_along(de_mg_strings)) {
  
  target_MG <- de_mg_strings[i]
  mg_row_ex <- mg_stats[mg_stats$MG_ID == target_MG, ]
  
  if (nrow(mg_row_ex) == 0) next
  
  mg_id_string <- mg_row_ex$MG_ID[1]
  parts <- strsplit(mg_id_string, "_")[[1]]
  
  chr_num   <- sub("^Gm", "", parts[1])
  pos_num   <- parts[2]
  mgmin_val <- mg_row_ex$MGmin[1]
  eps_val   <- mg_row_ex$epsilon[1]
  
  pattern_hap_ex <- paste0(
    "G7_RNS_hap_file_chr", sprintf("%02d", as.integer(chr_num)),
    "_pos", pos_num,
    "_mgmin", mgmin_val,
    "_eps", eps_val
  )
  
  pattern_ind_ex <- paste0(
    "G7_RNS_ind_file_chr", sprintf("%02d", as.integer(chr_num)),
    "_pos", pos_num,
    "_mgmin", mgmin_val,
    "_eps", eps_val
  )
  
  hap_file_ex <- all_files[str_detect(all_files, pattern_hap_ex)]
  ind_file_ex <- all_files[str_detect(all_files, pattern_ind_ex)]
  
  if (length(hap_file_ex) == 0 || length(ind_file_ex) == 0) next
  
  # READ DATA
  hap_df <- read.csv(hap_file_ex[1], stringsAsFactors = FALSE, check.names = FALSE)
  ind_df <- read.csv(ind_file_ex[1], stringsAsFactors = FALSE, check.names = FALSE) 
  hap_df$hap <- trimws(as.character(hap_df$hap))
  ind_df$hap <- trimws(as.character(ind_df$hap))
  mg_col <- trimws(as.character(hap_df[[target_MG]]))
  # FIND ALL DE HAPS (MAY BE MULTIPLE)
  de_haps <- hap_df$hap[mg_col == "2"]
  cat("DE haplotypes found:", length(de_haps), "\n")
  if (length(de_haps) == 0) next  
  # RUN FISHER FOR EACH HAP SEPARATELY  
  for (hap_i in de_haps) {
    
    inds_with_hap <- unique(
      ind_df$Ind[ind_df$hap == hap_i]
    )
    
    pheno <- ind_df[, c("Ind", "Metadata")]  # Use Metadata column (N/R/S)
    pheno <- pheno[!duplicated(pheno$Ind), ]
    pheno$has_hap <- ifelse(pheno$Ind %in% inds_with_hap, "TRUE", "FALSE")
    tab <- table(pheno$has_hap, pheno$Metadata)  # Creates 2×3: FALSE/TRUE × N/R/S
    fisher_res <- fisher.test(tab)  # This does Freeman-Halton for 2×3
    hap_counts <- tab["TRUE", ]
    enriched_pheno <- names(hap_counts)[which.max(hap_counts)]
    enrichment_prop <- max(hap_counts) / sum(hap_counts)
    all_fisher_results[[paste(target_MG, hap_i, sep="__")]] <-
      data.frame(
        MG = target_MG,
        hap = hap_i,
        p_value = fisher_res$p.value,
        enriched_phenotype = enriched_pheno,
        enrichment_proportion = enrichment_prop
      )
  }
}
# SAVE OUTPUT
if (length(all_fisher_results) > 0) {
  
  out <- do.call(rbind, all_fisher_results)
  
  out$FDR <- p.adjust(out$p_value, method = "BH")
  
  write.csv(out, "G7_MG_hap_level_fisher_results_withenrichment.csv", row.names = FALSE)
  
  print(out)
}

fisher_summary <-  read.csv("G7_MG_hap_level_fisher_results_withenrichment.csv")
head(fisher_summary)
#pick out the significant ones but only for the pvalue 0.05
sig_fisher <- fisher_summary[fisher_summary$p_value < 0.05, ]
nrow(sig_fisher)
write.csv(sig_fisher, "significant_fisher_G7.csv", row.names = FALSE)
sig_haps <- unique(sig_fisher$hap)
#now take in only the unique haps cos there will be some repeating
sig_haps <- unique(sig_fisher$hap)
sig_haps
length(sig_haps)
head(sig_haps)

parse_hap <- function(hap) {
  parts <- strsplit(hap, "_")[[1]]
  
  chr <- gsub("Gm", "", parts[1])
  pos <- parts[2]
  
  list(chr = chr, pos = pos)
}

library(stringr)

all_mats <- list()

for (hap in sig_haps) {
  
  info <- parse_hap(hap)
  
  chr <- info$chr
  pos <- info$pos
  message("Processing: ", hap)
  # find matching file
  pattern <- paste0(
    "G7_RNS_ind_file_chr", chr,
    "_pos", pos
  )
  file <- list.files("DE_selected_MGs", full.names = TRUE)
  file <- file[str_detect(file, pattern)]
  if (length(file) == 0) {
    message("No file found for ", hap)
    next
  }
  df <- read.csv(file[1], stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- trimws(colnames(df))
  df$hap <- trimws(df$hap)
  # keep only target hap
  df_hap <- df[df$hap == hap, ]
  if (nrow(df_hap) == 0) {
    message("No carriers for ", hap)
    next
  }
  # build simple IND vector (presence = 1)
  mat <- data.frame(
    Ind = df_hap$Ind,
    hap = hap,
    value = 1
  )
  
  all_mats[[hap]] <- mat
}
long_df <- do.call(rbind, all_mats)
inds <- unique(long_df$Ind)
haps <- unique(long_df$hap)
hap_mat <- matrix(0,
                  nrow = length(inds),
                  ncol = length(haps),
                  dimnames = list(inds, haps))
for (i in 1:nrow(long_df)) {
  hap_mat[long_df$Ind[i], long_df$hap[i]] <- 1
}
hap_mat <- as.data.frame(hap_mat)
head(hap_mat)
#save the above objects
saveRDS(hap_mat, "hap_mat_G7.rds")
saveRDS(sig_fisher, "Sig_fisher_G7.rds")
saveRDS(sig_haps, "sig_haps_G7.rds")
