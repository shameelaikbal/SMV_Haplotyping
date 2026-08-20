
# Soybean mosaic virus (SMV) resistance genomics

## Overview
Comprehensive genomics analysis identifying genome-wide haplotype combinations associated with differential responses of soybean to SMV strains, G1 and G7.

**Phenotype data:**
- SMV-G1 retrieved from USDA-GRIN https://npgsweb.ars-grin.gov/gringlobal/descriptordetail?id=51159
- SMV-G7 phenotype data retrieved from https://npgsweb.ars-grin.gov/gringlobal/descriptordetail?id=51165 

**Genotype data:**
- **Reference panel for genotype imputation**: [Valliyodan et al. (2021)](https://www.nature.com/articles/s41597-021-00834-w) 
- **Target VCF**: [Song et al. (2015)](https://www.soybase.org/tools/snp50k/) dataset

## Analysis Pipeline

### **1: GWAS**
- Association analysis for SMV-G1 and SMV-G7 using rMVP
- **Scripts**:```01_GWAS.R```

### **2: Imputation**
- Beagle 5.4 phasing and imputation
- Quality filtering (DR2 > 0.8, MAF > 0.03)
- **Scripts**: `02_imputation.sh`

### **3: Local haplotyping**
- Construct local haplotype structures around the top 100 GWAS-SNPs using crosshap
- Automated parameter selection for MGmin and Epsilon, based on SNP density for each region
- Create unified MG matrix for input to differential evolution
- **Scripts**: `03_local_haplotyping.R`
  
### **4. Differential evolution**
- Differential evolution algorithm to identify an optimal combination of MGs that maximises phenotypic prediction
- Fitness function based on 5-fold cross-validated Pearson correlation of observed and predicted phenotypes
- Obtain DE-selected MGs
- **Scripts**: `04_differential_evolution.R`

### **5. Significant haplotypes**
- From the DE-selected MGs, idenitfy haplotypes from crosshap outputs
- Conduct Fisher's exact/Freeman-Holtman extension of Fishers exact test to identify significant haplotypes associated with the differential responses of soybean to SMV
- **Scripts**: `05_significant_haplotypes.R`

### **6. Correlation matrix**
- Construct correlation matrix to show the pairwise relationships between the significant haplotypes
- Corr matrix using a signed Phi-based correlation matrix
- **Scripts**: `06_correlation_matrix.R`

### **7. Multi-locus haplotype combinations**
- From the haplotype pairs with correlation, form Upset plots to identify haplotype combinations
- **Scripts**: `07_upset_plots.R`

### **8. Variant effect prediction**
- Within haplotypes of interest, variant effect prediction using SnpEff

### **9. Candidate gene identification, functional annotation, and literature support**
- Candidate gene identification using genome browsers in SoyBase, Phytozome, SnpEff results
- Protein domain information retrieved from UniProt
- Literature review for gene families and candidate gene role in SMV

### **10.Results**
- Raw variant information, haplotype carriers, and visualisation of significant haplotypes associated with SMV-G1 and SMV-G7 available in `08_results/`
