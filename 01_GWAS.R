library (rMVP)

MVP.Data (fileVCF=".vcf",
          filePhe="phenotype.txt",
          fileKin=FALSE,
          filePC=FALSE,
          out="G7")

genotype <- attach.big.matrix("G7.geno.desc")
phenotype <- read.table("G7.phe",head=TRUE)
map <- read.table("G7.geno.map", head = TRUE)


fMVP <- MVP(
  phe=phenotype,
  geno=genotype,
  map=map,
  #K=Kinship,
  #CV.GLM=Covariates,     ##if you have additional covariates, please keep there open.
  #CV.MLM=Covariates,
  #CV.FarmCPU=Covariates,
  nPC.GLM=5,      ##if you have added PC into covariates, please keep there closed.
  nPC.MLM=5,
  nPC.FarmCPU=4,
  priority="memory",       ##for Kinship construction
  #ncpus=10,
  vc.method="BRENT",      ##only works for MLM
  maxLoop=10,
  method.bin="static",      ## "FaST-LMM", "static" (#only works for FarmCPU)
  #permutation.threshold=TRUE,
  #permutation.rep=100,
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"),
  out="G7_npc4")
