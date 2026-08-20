#local haplotyping using crosshap with automated parameter selection for the top 100 GWAS-SNPs associated with soybean mosaic virus resistance

library(crosshap)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(clustree)
library(vcfR)
library(dplyr)

#required inputs
#phenotype data
#top 100 GWAS-SNPs regions in a text file
#extract 200 kb region around each SNP - the region to be delimited depends on the LD decay for the plant species (vcf_file)
#calculate LD matrices for the SNP position using plink (ld_file)
#more documentation for running crosshap available here https://jacobimarsh.github.io/crosshap/

# Set base directory and file paths
base_dir <- 
setwd(base_dir)

# Create directory to save results
result_dir <- file.path(base_dir, "G7_crosshap_results_DE_RNS")
if (!dir.exists(result_dir)) dir.create(result_dir)

# Load phenotype
pheno_data <- read_pheno("SMV_G7_NRaszero_S1.txt")
colnames(pheno_data) <- c("Ind", "Pheno")
head(pheno_data)
nrow(pheno_data)

# Load metadata
metadata_G7 <- read_metadata("SMV_G7_common.txt")
colnames(metadata_G7) <- c("Ind", "Metadata")
head(metadata_G7)

# Read top 100 GWAS ranked SNPs
snp_data <- read.table("top100_G7_final.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(snp_data) <- c("chr_full", "pos")
head(snp_data)

# Extract chromosome number
snp_data$chr <- gsub(".*Gm", "", snp_data$chr_full)

# Parameter grid function for epsilon and mgmin
get_param_grids <- function(n_snps) {
  if (n_snps < 200) {
    mgmin_range <- c(5, 10, 15)
    epsilon_range <- c(0.2, 0.4, 0.6, 0.8, 1)
  } else if (n_snps >= 200 && n_snps < 800) {
    mgmin_range <- c(10, 18, 25)
    epsilon_range <- c(0.2, 0.4, 0.8, 1, 0.6)
  } else if (n_snps >= 800 && n_snps < 2000) {
    mgmin_range <- c(15, 25, 35)
    epsilon_range <- c(0.2, 0.4, 0.6, 0.8, 1)
  } else {
    mgmin_range <- c(20, 30, 40)
    epsilon_range <- c(0.2, 0.4, 0.6, 0.8, 1)
  }
  list(mgmin_range = mgmin_range, epsilon_range = epsilon_range)
}

# main crosshap loop - run haplotyping for all regions

for (i in 1:nrow(snp_data)) {
  chr <- snp_data$chr[i]
  pos <- snp_data$pos[i]
  
  cat("Processing region", i, "of", nrow(snp_data), "- Chr", chr, "Pos", pos, "\n")
  
  vcf_file <- sprintf("regions/G7_region_200kb_chr%s_%s.vcf", chr, pos)
  ld_file <- sprintf("ld/ld_chr%s_%s_clean.ld", chr, pos)
  
  # Read VCF with vcfR to count SNPs
  vcf_data <- tryCatch({
    read.vcfR(vcf_file, verbose = FALSE)
  }, error = function(e) {
    cat("  ERROR reading VCF:", e$message, "\n")
    cat("  Skipping region\n\n")
    return(NULL)
  })
  
  if(is.null(vcf_data)) {
    next
  }
  
  n_snps <- nrow(vcf_data@fix)
  cat("  SNPs in region:", n_snps, "\n")
  
  # Read vcf with crosshap 
  vcf_ch <- tryCatch({
    vcf_tmp <- read_vcf(vcf_file)
    
    # Try to add ID column
    vcf_tmp <- tryCatch({
      vcf_tmp %>% mutate(ID = paste0("Gm", chr, "_", POS))
    }, error = function(e) {
      cat("  ERROR: Cannot add ID column (POS issue):", e$message, "\n")
      cat("  Skipping region\n\n")
      return(NULL)
    })
    
    vcf_tmp
    
  }, error = function(e) {
    cat("  ERROR: crosshap read_vcf failed:", e$message, "\n")
    cat("  Skipping region\n\n")
    return(NULL)
  })
  
  # Skip if VCF reading failed
  if(is.null(vcf_ch)) {
    next
  }
  
  # Get parameter grids
  param_grids <- get_param_grids(n_snps)
  mgmin_range <- param_grids$mgmin_range
  epsilon_range <- param_grids$epsilon_range
  
  cat("  Testing", length(mgmin_range), "MGmin values:", paste(mgmin_range, collapse=", "), "\n")
  cat("  Testing", length(epsilon_range), "epsilon values:", paste(epsilon_range, collapse=", "), "\n")
  
  # Read LD file - ERROR HANDLING
  ld_data <- tryCatch({
    read_LD(ld_file, vcf = vcf_ch)
  }, error = function(e) {
    cat("  ERROR: crosshap read_LD failed:", e$message, "\n")
    cat("  Skipping region\n\n")
    return(NULL)
  })
  
  # Skip if LD reading failed
  if(is.null(ld_data)) {
    next
  }
  
  # Loop through parameter combinations
  for (mgmin in mgmin_range) {
    for (eps in epsilon_range) {
      
      haplo_result <- tryCatch(
        run_haplotyping(
          vcf = vcf_ch,
          LD = ld_data,
          pheno = pheno_data,
          metadata = metadata_G7,
          epsilon = eps,
          MGmin = mgmin,
          minHap = 3
        ),
        error = function(e) {
          message(sprintf("  Error at mgmin %d eps %.2f: %s", mgmin, eps, e$message))
          return(NULL)
        }
      )
      
   
    if (!is.null(haplo_result)) {
        hap_key <- sprintf("Haplotypes_MGmin%d_E%.1f", mgmin, eps)
        
        if (!is.null(haplo_result[[hap_key]])) {
          hap_viz <- crosshap_viz(HapObject = haplo_result, epsilon = eps)
          hap_data <- haplo_result[[hap_key]]
          
          # Calculate haplotype quality metrics
          # Basic metrics
          n_haplotypes <- length(unique(hap_data$Indfile$hap))
          n_samples <- nrow(hap_data$Indfile)
          
          # Get variants assigned to MGs (exclude unassigned SNPs with MG=0)
          mg_variants <- hap_data$Varfile %>%
            filter(!is.na(MGs) & MGs != 0 & MGs != "")
          
          n_MGs <- length(unique(mg_variants$MGs))
          
          # AAF quality metrics 
          if(nrow(mg_variants) > 0 & "AltAF" %in% colnames(mg_variants)) {
            mean_aaf <- mean(mg_variants$AltAF, na.rm = TRUE)
            usable <- sum(mg_variants$AltAF >= 0.05 & mg_variants$AltAF <= 0.95, na.rm = TRUE)
            pct_usable_aaf <- usable / nrow(mg_variants)
          } else {
            mean_aaf <- NA
            pct_usable_aaf <- 0
          }
          
          # Coverage: what % of total samples were assigned to haplotypes?
          coverage <- n_samples / length(unique(pheno_data$Ind))
          
          # Store metrics in Indfile
          hap_data$Indfile$n_haplotypes <- n_haplotypes
          hap_data$Indfile$n_MGs <- n_MGs
          hap_data$Indfile$coverage <- coverage
          hap_data$Indfile$mean_aaf <- mean_aaf
          hap_data$Indfile$pct_usable_aaf <- pct_usable_aaf
          
          # Add metadata       
          # Indfile metadata
          hap_data$Indfile$epsilon <- eps
          hap_data$Indfile$MGmin <- mgmin
          hap_data$Indfile$chr <- chr
          hap_data$Indfile$pos <- pos
          hap_data$Indfile$region_id <- paste0("Gm", chr, "_", pos)
          
          # Hapfile metadata
          hap_data$Hapfile$epsilon <- eps
          hap_data$Hapfile$MGmin <- mgmin
          hap_data$Hapfile$chr <- chr
          hap_data$Hapfile$pos <- pos
          
          # Varfile metadata
          hap_data$Varfile$epsilon <- eps
          hap_data$Varfile$MGmin <- mgmin
          hap_data$Varfile$chr <- chr
          hap_data$Varfile$gwas_pos <- pos
          
          # Create unique MG_ID in Var file
          # Format: Gm[chr]_[pos]_MG[number]
          hap_data$Varfile$MG_ID <- ifelse(
            !is.na(hap_data$Varfile$MGs) & hap_data$Varfile$MGs != 0 & hap_data$Varfile$MGs != "",
            paste0("Gm", chr, "_", pos, "_", hap_data$Varfile$MGs),
            NA
          )
          # Rename MG columns in hap file
          # Changes MG1, MG2... to Gm[chr]_[pos]_MG1, Gm[chr]_[pos]_MG2...
          mg_cols <- grep("^MG[0-9]+$", colnames(hap_data$Hapfile), value = TRUE)
          for (mg_col in mg_cols) {
            new_col_name <- paste0("Gm", chr, "_", pos, "_", mg_col)
            colnames(hap_data$Hapfile)[colnames(hap_data$Hapfile) == mg_col] <- new_col_name
          }
          
          # ====== SAVE FILES ======
          ind_file_name <- sprintf("G7_RNS_ind_file_chr%s_pos%s_mgmin%d_eps%.2f.csv", chr, pos, mgmin, eps)
          hap_file_name <- sprintf("G7_RNS_hap_file_chr%s_pos%s_mgmin%d_eps%.2f.csv", chr, pos, mgmin, eps)
          var_file_name <- sprintf("G7_RNS_var_file_chr%s_pos%s_mgmin%d_eps%.2f.csv", chr, pos, mgmin, eps)
          viz_file_name <- sprintf("G7_RNS_viz_chr%s_pos%s_mgmin%d_eps%.2f.jpg", chr, pos, mgmin, eps)
          
          write.csv(hap_data$Indfile, file.path(result_dir, ind_file_name), quote = FALSE, row.names = FALSE)
          write.csv(hap_data$Hapfile, file.path(result_dir, hap_file_name), quote = FALSE, row.names = FALSE)
          write.csv(hap_data$Varfile, file.path(result_dir, var_file_name), quote = FALSE, row.names = FALSE)
          ggsave(file.path(result_dir, viz_file_name), hap_viz, width=12, height=8)
          
        }
      }
      
    } # end epsilon loop
  } # end mgmin loop
  
  cat("  ✓ Region", i, "complete\n\n")
  
} # end region loop

# Rename hap labels with region prefix

cat("=== Adding region prefixes to haplotype labels ===\n\n")

ind_files <- list.files(result_dir, pattern = "G7_RNS_ind_file_chr.*csv", full.names = TRUE)
hap_files <- list.files(result_dir, pattern = "G7_RNS_hap_file_chr.*csv", full.names = TRUE)

cat("Found", length(ind_files), "ind files\n")
cat("Found", length(hap_files), "hap files\n\n")

for (f in seq_along(ind_files)) {
  ind_file <- ind_files[f]
  hap_file <- hap_files[f]
  
  # Extract chr and pos from filename
  chr <- str_match(ind_file, "chr(\\d+)")[,2]
  pos <- str_match(ind_file, "_pos(\\d+)")[,2]
  region_prefix <- paste0("Gm", chr, "_", pos, "_")
  
  # Read data
  ind_data <- read.csv(ind_file)
  head(ind_data)
  hap_data <- read.csv(hap_file)
  
  # Add region prefix to haplotype labels
  ind_data$hap <- paste0(region_prefix, ind_data$hap)
  hap_data$hap <- paste0(region_prefix, hap_data$hap)
  
  # Overwrite files with updated labels
  write.csv(ind_data, ind_file, quote = FALSE, row.names = FALSE)
  write.csv(hap_data, hap_file, quote = FALSE, row.names = FALSE)
}

cat("✓ Haplotype labels updated with region prefixes\n\n")
#automated parameter selection
# Helper function: Calculate phenotypic discrimination at HAPLOTYPE level
calc_pheno_discrimination <- function(ind_data, pheno_col = "Pheno") {
  # This measures: Do haplotypes separate R or N from S?
  # Mean phenotype per haplotype
  hap_means <- ind_data %>%
    group_by(hap) %>%
    summarise(
      hap_mean = mean(.data[[pheno_col]], na.rm = TRUE),
      hap_n = n(),
      .groups = "drop"
    )  
  # Overall mean
  grand_mean <- mean(ind_data[[pheno_col]], na.rm = TRUE)
  # Between-haplotype variance
  between_var <- sum(hap_means$hap_n * (hap_means$hap_mean - grand_mean)^2) / 
    sum(hap_means$hap_n)
  # Total variance
  total_var <- var(ind_data[[pheno_col]], na.rm = TRUE)
  # Eta-squared
  eta_sq <- between_var / total_var
  # Return 0 if undefined
  if(is.na(eta_sq) | is.nan(eta_sq) | is.infinite(eta_sq)) return(0)
  
  return(eta_sq)
}

# Read all ind files
cat("Loading all ind files...\n")
all_ind_data <- bind_rows(lapply(ind_files, read.csv))
cat("Loaded", nrow(all_ind_data), "rows from", length(ind_files), "files\n\n")

# Calculate phenotypic discrimination for each parameter combo
cat("Calculating phenotypic discrimination scores...\n")
cat("(This measures how well haplotypes separate Resistant from Susceptible)\n\n")

pheno_disc_scores <- all_ind_data %>%
  group_by(chr, pos, epsilon, MGmin) %>%
  group_modify(~ {
    data.frame(pheno_discrimination = calc_pheno_discrimination(.x))
  }) %>%
  ungroup()

cat("✓ Phenotypic discrimination calculated for all parameter combinations\n\n")

# Calculate all selection metrics
cat("Calculating selection metrics...\n")

region_metrics <- all_ind_data %>%
  group_by(chr, pos, epsilon, MGmin) %>%
  summarise(
    n_haplotypes = first(n_haplotypes),
    n_MGs = first(n_MGs),
    coverage = first(coverage),
    mean_aaf = first(mean_aaf),
    pct_usable_aaf = first(pct_usable_aaf),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  filter(n_samples >= 150) %>%  # Pre-filter: exclude combos with <150 samples
  # Add phenotypic discrimination
  left_join(pheno_disc_scores, by = c("chr", "pos", "epsilon", "MGmin")) %>%
  # Calculate stability scores (compare to next epsilon)
  group_by(chr, pos, MGmin) %>%
  arrange(epsilon) %>%
  mutate(
    # Stability: does haplotype count change much at next epsilon?
    next_n_hap = lead(n_haplotypes),
    stability = ifelse(
      !is.na(next_n_hap),
      1 - abs(n_haplotypes - next_n_hap) / pmax(n_haplotypes, next_n_hap, 1),
      1.0  # Highest epsilon, assume stable
    )
  ) %>%
  ungroup() %>%
  mutate  
    # 1. Phenotypic discrimination score 
    #    How well do haplotypes separate R from S?
    pheno_disc_score = ifelse(is.na(pheno_discrimination), 0, pheno_discrimination), 
    # 2. Haplotype diversity score (25%)
    #    Prefer 7-15 haplotypes (not too simple, not fragmented)
    hap_score = case_when(
      n_haplotypes < 3 ~ 0,
      n_haplotypes >= 3 & n_haplotypes <= 6 ~ 0.6,
      n_haplotypes >= 7 & n_haplotypes <= 15 ~ 1.0,
      n_haplotypes >= 16 & n_haplotypes <= 20 ~ 0.8,
      n_haplotypes > 20 ~ 0.5
    ), 
    # 3. Coverage score (20%)
    #    Higher % of samples assigned = better
    cov_score = coverage,  
    # 4. Stability score (10%)
    #    Prefer parameters where haplotypes are stable
    stab_score = stability,
    # 5. AAF quality score (5%)
    #    Prefer MGs with usable allele frequencies
    aaf_score = ifelse(is.na(pct_usable_aaf), 0, pct_usable_aaf),
    # 6. MG count score (5%)
    #    Prefer 2-8 MGs per region
    mg_score = case_when(
      n_MGs < 2 ~ 0.3,
      n_MGs >= 2 & n_MGs <= 8 ~ 1.0,
      n_MGs > 8 ~ 0.7
    ),
    
    selection_score = (pheno_disc_score * 0.35) +  # 35% - Phenotypic discrimination
      (hap_score * 0.25) +         # 25% - Haplotype diversity
      (cov_score * 0.20) +         # 20% - Coverage
      (stab_score * 0.10) +        # 10% - Stability
      (aaf_score * 0.05) +         #  5% - AAF quality
      (mg_score * 0.05)            #  5% - MG count
  )
# Select best parameter combo per region (highest selection score)
optimal_params_auto <- region_metrics %>%
  group_by(chr, pos) %>%
  slice_max(selection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(selection_method = "automated")
# Save all metrics
write.csv(region_metrics, file.path(result_dir, "all_parameter_metrics.csv"), 
          row.names = FALSE)
write.csv(optimal_params_auto, file.path(result_dir, "optimal_parameters_automated.csv"),
          row.names = FALSE)

cat("✓ Saved all_parameter_metrics.csv (all parameter combos tested)\n")
cat("✓ Saved optimal_parameters_automated.csv (best combo per region)\n\n")

# filter data to final selected parameters
# 1. Create file mapping with ZERO-PADDED chr (02 not 2)
cat("Creating file mappings from optimal parameters...\n")
optimal_files <- optimal_params_final %>%
  rowwise() %>%
  mutate(
    chr_padded = sprintf("%02d", chr),  # 2 → "02", 19 → "19"
    hap_file = file.path(result_dir, 
                         sprintf("G7_RNS_hap_file_chr%s_pos%d_mgmin%d_eps%.2f.csv", 
                                 chr_padded, pos, MGmin, epsilon)),
    ind_file = file.path(result_dir, 
                         sprintf("G7_RNS_ind_file_chr%s_pos%d_mgmin%d_eps%.2f.csv", 
                                 chr_padded, pos, MGmin, epsilon))
  ) %>%
  select(chr, chr_padded, pos, epsilon, MGmin, hap_file, ind_file) %>%
  ungroup()

print(head(optimal_files))  # Check filenames look right

# 2. Verify all files exist
cat("\nChecking files exist...\n")
existing_hap <- file.exists(optimal_files$hap_file)
existing_ind <- file.exists(optimal_files$ind_file)

cat("Hap files found:", sum(existing_hap), "/", nrow(optimal_files), "\n")
cat("Ind files found:", sum(existing_ind), "/", nrow(optimal_files), "\n")

if(any(!existing_hap)) {
  cat("⚠ Missing hap files:\n")
  print(optimal_files[!existing_hap, c("chr", "pos", "epsilon", "MGmin", "hap_file")])
}

if(any(!existing_ind)) {
  cat("⚠ Missing ind files:\n")
  print(optimal_files[!existing_ind, c("chr", "pos", "epsilon", "MGmin", "ind_file")])
}

cat("\n✓ File check complete\n")

# 3. Load ONLY selected files
cat("Loading", nrow(optimal_files), "selected hap/ind file pairs...\n")
hap_list <- lapply(optimal_files$hap_file, read.csv)
ind_list <- lapply(optimal_files$ind_file, read.csv)
# Bind all hap files into one big table
all_hap_data_filtered <- bind_rows(hap_list)
cat("✓ All hap files bound:", nrow(all_hap_data), "rows\n")
# Bind all ind files into one big table  
all_ind_data_filtered <- bind_rows(ind_list)
cat("✓ All ind files bound:", nrow(all_ind_data), "rows\n")


# CREATE MG MATRIX FOR DE OPTIMIZATIO
# Merge ind and hap data
cat("Merging individual and haplotype data...\n")
ind_hap_merged <- all_ind_data_filtered %>%
  left_join(all_hap_data_filtered, by = "hap", relationship = "many-to-one")

cat("Merged data:", nrow(ind_hap_merged), "rows\n\n")

# Extract MG columns (format: Gm[chr]_[pos]_MG[number])
mg_cols <- grep("^Gm\\d+_\\d+_MG\\d+$", colnames(ind_hap_merged), value = TRUE)
cat("Found", length(mg_cols), "unique MG columns across all regions\n\n")

if(length(mg_cols) == 0) {
  stop("ERROR: No MG columns found! Check column naming in hapfiles.")
}

# Collapse to one row per individual
# (Each individual appears once per region, we take max MG value across regions)
cat("Collapsing to one row per individual...\n")

mg_matrix <- ind_hap_merged %>%
  select(Ind, all_of(mg_cols)) %>%
  group_by(Ind) %>%
  summarise(across(everything(), ~max(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))

# Convert to matrix with individuals as rownames
mg_matrix_final <- mg_matrix %>%
  column_to_rownames("Ind") %>%
  as.matrix()

# Save matrix
saveRDS(mg_matrix_final, "G1_haplotype_matrix.rds")
cat("✓ Saved G1_haplotype_matrix.rds\n\n")
#create MG stats file
# Load data
mg_matrix_final <- readRDS("G7_haplotype_matrix.rds")

# Step 1: Create mg_stats from MG names
mg_stats <- data.frame(MG_ID = colnames(mg_matrix_final), stringsAsFactors = FALSE) %>%
  tidyr::extract(
    MG_ID,
    into = c("chromosome", "position", "mg_num"),
    regex = "Gm(\\d+)_(\\d+)_MG(\\d+)",
    remove = FALSE
  ) %>%
  mutate(
    chromosome = as.integer(chromosome),
    position = as.integer(position),
    mg_num = as.integer(mg_num)
  )

# Step 2: Add epsilon/MGmin from hap data
region_params <- all_hap_data_filtered %>%
  select(chr, pos, epsilon, MGmin) %>%
  distinct()

mg_stats <- mg_stats %>%
  left_join(region_params, by = c("chromosome" = "chr", "position" = "pos"))

# Step 3: Calculate AAF from matrix dosages
mg_stats$AAF <- colMeans(mg_matrix_final) / 2

# Step 4: Calculate phenotypic_diff
pheno_for_stats <- read.table("SMV_G7_NRaszero_S1.txt", header = TRUE)
colnames(pheno_for_stats) <- c("Ind", "Pheno")
pheno_for_stats <- pheno_for_stats %>% distinct(Ind, .keep_all = TRUE)
rownames(pheno_for_stats) <- pheno_for_stats$Ind

mg_stats$phenotypic_diff <- sapply(mg_stats$MG_ID, function(mg_id) {
  has_mg <- rownames(mg_matrix_final)[mg_matrix_final[, mg_id] >= 1]
  no_mg <- rownames(mg_matrix_final)[mg_matrix_final[, mg_id] == 0]
  
  if(length(has_mg) == 0 | length(no_mg) == 0) return(NA)
  
  mean_with <- mean(pheno_for_stats$Pheno[pheno_for_stats$Ind %in% has_mg], na.rm = TRUE)
  mean_without <- mean(pheno_for_stats$Pheno[pheno_for_stats$Ind %in% no_mg], na.rm = TRUE)
  
  return(mean_with - mean_without)
})

# Save
saveRDS(mg_stats, "G7_mg_stats.rds")
