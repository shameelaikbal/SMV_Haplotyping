#phasing and imputation of a target VCF using a reference panel
# Beagle 5.4 imputation pipeline with quality filtering

#phase target VCF and reference panel VCF
java -Xmx128g -jar beagle.06Aug24.a91.jar \
    gt=target.vcf \
    nthreads=16 \
    out=target_phased

  java -Xmx128g -jar beagle.06Aug24.a91.jar \
    gt=reference_panel.vcf \
    nthreads=16 \
    out=reference_panel_phased  

#impute phased target VCF using phased reference panel, done using an array

time srun -m block:block:block java -Xmx256g -jar beagle.06Aug24.a91.jar \
  ref=reference_panel_phased.vcf.gz\
  gt=target_phased.vcf.gz \
  nthreads=24 \
  out=target_imputed

#index the imputed file
bcftools index -t target_imputed.vcf.gz

#filter imputed for DR2 > 0.8
bcftools filter -i 'DR2>0.8'  target_imputed.vcf.gz -Oz -o target_imputed_filtered_for_DR0.8.vcf.gz

#index the DR2 filtered files
bcftools index -t target_imputed_filtered_for_DR0.8.vcf.gz

#if done separately for each chromosome, concat the DR2 filtered file
bcftools concat -a -O z -o target_imputed_filtered_for_DR0.8_merged.vcf.gz -f list_DR_filtered_files.txt

#filter for MAF 0.05
vcftools --xzvcf target_imputed_filtered_for_DR0.8_merged.vcf.gz --maf 0.05 --recode --out target_imputed_filtered_for_DR0.8_maf0.03_merged.vcf.gz
