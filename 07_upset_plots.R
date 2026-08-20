 library(ComplexUpset)
 library(ggplot2)
 library(dplyr)
 library(ggtext) 

#preprocessing - load phenotype data
pheno_df <- read.table("SMV_G7_common.txt", header = TRUE)
head(pheno_df)
colnames(pheno_df) = c("Ind", "Pheno")

phi_mat <-  readRDS("Phi_g7.rds")
hap_mat <- readRDS("hap_mat_G7.rds")

keep_haps <- colnames(phi_mat)[
  apply(phi_mat, 2, function(x) any(x > 0.2, na.rm = TRUE))
]

# keep only shared accessions
common_ids <- intersect(rownames(hap_mat), pheno_df$Ind)
hap_mat2 <- hap_mat[common_ids, , drop = FALSE]
head(hap_mat2)

pheno_df2 <- pheno_df[pheno_df$Ind %in% common_ids, , drop = FALSE]
pheno_df2 <- pheno_df2[match(rownames(hap_mat2), pheno_df2$Ind), ]
head(pheno_df2)

stopifnot(all(rownames(hap_mat2) == pheno_df2$Ind))

hap_mat_filt <- hap_mat2[, keep_haps, drop = FALSE]

stopifnot(all(rownames(hap_mat_filt) == pheno_df2$Ind))

upset_df <- as.data.frame(hap_mat_filt)

upset_df$Ind <- rownames(upset_df)
upset_df <- merge(upset_df, pheno_df2[, c("Ind", "Pheno")], by = "Ind", all.x = TRUE)


upset_df$status <- ifelse(upset_df$Pheno == "N", "SMV-G7 Necrotic", 
                          ifelse(upset_df$Pheno == "R", "SMV-G7 Resistant", 
                                 "SMV-G7 Susceptible"))

# 2. Delete the old column
upset_df$Pheno <- NULL
plot_haps <- keep_haps

 
 # start plotting the upset from here after all the preprocessing
 plot_df_clean <- upset_df
 colnames(plot_df_clean) [46] = "Phenotype"
 
 # Dynamic aggregation of the highly variable 10-11Mb interval on Chromosome 07
 chr07_all_cols <- colnames(plot_df_clean)[grep("^Gm07_10|^Gm07_11", colnames(plot_df_clean))]
 plot_df_clean$Gm07_10_11Mb <- ifelse(
   rowSums(plot_df_clean[, chr07_all_cols, drop = FALSE], na.rm = TRUE) > 0, 1, 0
 )
 
 # 10 final loci
 core_manuscript_haps <- c(
   "Gm07_10_11Mb", 
   "Gm10_6148847_E", 
   "Gm10_3474358_H",
   "Gm13_30293949_G",
   "Gm13_30293949_B",       
   "Gm14_2187298_G", 
   "Gm02_12365201_B",       
   "Gm06_48396005_A",       
   "Gm07_9176010_E",
   "Gm18_9299920_C"
 )
 
 # Enforce strict factor ordering for G7 status groups
 plot_df_clean$Phenotype <- factor(
   plot_df_clean$Phenotype,
   levels = c("SMV-G7 Resistant", "SMV-G7 Susceptible", "SMV-G7 Necrotic")
 )
 
 #colours
 pheno_cols <- c(
   "SMV-G7 Resistant"   = "#1a936f", 
   "SMV-G7 Susceptible" = "#f28f43", 
   "SMV-G7 Necrotic"    = "#0f4c5c"  
 )
 
 highlight_queries <- list(
   # 1. Susceptible Combo: Highlight ONLY the dots/lines in Orange
   upset_query(
     intersect = c("Gm07_10_11Mb", "Gm02_12365201_B", "Gm13_30293949_B"),
     color = "#f28f43", 
     only_components = c("intersections_matrix")
   ),
   # 2. Necrotic Combo: Highlight ONLY the dots/lines in dark blue
   upset_query(
     intersect = c("Gm06_48396005_A", "Gm14_2187298_G", "Gm13_30293949_G"),
     color = "#0f4c5c", 
     only_components = c("intersections_matrix")
   ),
   # 3. Resistant Combo: Highlight ONLY the dots/lines in Dark Blue
   upset_query(
     intersect = c("Gm07_10_11Mb", "Gm18_9299920_C", "Gm10_6148847_E", "Gm10_3474358_H"),
     color = "#1a936f",  
     only_components = c("intersections_matrix")
   )
 )
 
 # 4. Generate the final Master Figure
 p_master_figure4 <- upset(
   data = plot_df_clean,
   intersect = core_manuscript_haps,
   name = "Multi-locus haplotype combinations",
   min_size = 1,                       
   width_ratio = 0.25,
   min_degree = 3,
   guides = "over",                    # Keeps clean spacing and legend placement
   stripes = c("white", "grey95"),
   queries = highlight_queries,
   
   matrix = (
     intersection_matrix(
       geom = geom_point(size = 2.5)
     ) # Optional: re-append + scale_y_discrete(labels = coloured_labels) if using custom row text color
   ),
   
   base_annotations = list(
     "Intersection size" = intersection_size(
       counts = FALSE,
       position = "fill",              # Scales the top bars to a clean 100% proportion layout
       mapping = aes(fill = Phenotype)   
     ) + 
       scale_fill_manual(values = pheno_cols) +
       scale_y_continuous(labels = scales::percent_format()) +
       ylab("Proportion of Panel (%)")
   ),
   
   set_sizes = (
     upset_set_size(
       geom = geom_bar(aes(fill = Phenotype), position = "stack", width = 0.6),
       position = "left"
     )
     + scale_fill_manual(values = pheno_cols, guide = "none")
     + xlab("Set Size (Accession Count)")
   ),
   
   themes = upset_modify_themes(
     list(
       "Intersection size" = theme(
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.x = element_blank()
       ),
       "intersections_matrix" = theme(
         axis.title.y = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.text.y = element_markdown() 
       ),
       "set_sizes" = theme(
         axis.title.y = element_blank()  
       )
     )
   )
 )
 
 # Render clean figure
 print(p_master_figure4)
 
