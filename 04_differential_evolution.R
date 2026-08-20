#Differential evolution haplotype optimisation for soybean mosaic virus resistance
library(data.table)
# Set random seed for reproducibility
set.seed(123)
#load SMV-G7 haplotype matrix from G7 GWAS and crosshap analysis
mg_haplotype_matrix <- readRDS("G7_haplotype_matrix.rds")
cat("  Dimensions:", nrow(mg_haplotype_matrix), "samples ×", 
    ncol(mg_haplotype_matrix), "marker groups\n")
head(mg_haplotype_matrix, n=5)
#check row names cos i think the above loaded the PI as row names
head(rownames(mg_haplotype_matrix), 5)
# "PI104708" "PI105579" "PI148260" "PI153292" "PI153321"
#convert to a dataframe
mg_haplotype_matrix <-  data.frame(PI = rownames(mg_haplotype_matrix),
                                   mg_haplotype_matrix,
                                   check.names = FALSE)
head(mg_haplotype_matrix, 5)
#now reset row names to default 1,2,3
rownames(mg_haplotype_matrix) <- NULL

# Load phenotype data for G7
pheno_data <- read.table("SMV_G7_NRaszero_S1.txt", header=FALSE, stringsAsFactors=FALSE)
#extract only MG columns in the haplotype matrix
mg_cols <- setdiff(colnames(mg_haplotype_matrix), "PI")
geno <-  as.matrix(mg_haplotype_matrix[, mg_cols])
rownames(geno) <-  mg_haplotype_matrix$PI
# Align samples between genotype and phenotype
samples <- intersect(rownames(geno), pheno_data$PI)
cat("  Common samples:", length(samples), "\n")
# Load G7 marker group statistics from crosshap analysis
cat("\nLoading G7 marker group statistics...\n")
mg_stats <- readRDS("G7_mg_stats.rds")
cat("  Marker groups in stats:", nrow(mg_stats), "\n")
#marker groups in G7 stats : 384 this is the total number of haplotypes for G7 before DE
#check number of columns in mg_haplotype_matrix
ncol(mg_haplotype_matrix)
# Check required columns
required_cols <- c("chromosome", "position", "AAF", "phenotypic_diff", "within_LD")
missing_cols <- setdiff(required_cols, colnames(mg_stats))
if(length(missing_cols) > 0) {
  stop(paste("ERROR: mg_stats missing required columns:", paste(missing_cols, collapse=", ")))
}
geno <- geno[samples, ]                    # 213 x 385
pheno_data <- pheno_data[match(samples, pheno_data$PI), ]
pheno <- pheno_data$pheno                  # Numeric phenotype vector 
# Verify everything matches
cat("geno dims:", nrow(geno), "x", ncol(geno), "\\n")      
cat("pheno length:", length(pheno), "\\n")                 
cat("mg_stats nrow:", nrow(mg_stats), "\\n")               

if(nrow(geno) != length(pheno) || ncol(geno) != nrow(mg_stats)) {
  stop("DIMENSIONS DON'T MATCH!")
}

# Prepare cross validation folds
n_samples <- length(pheno)
rndIndex <- sample(n_samples, replace=FALSE)
multiplier <- round(n_samples / 5, digits=0)

missindex1 <- rndIndex[1:multiplier]
missindex2 <- rndIndex[(multiplier+1):(2*multiplier)]
missindex3 <- rndIndex[(2*multiplier+1):(3*multiplier)]
missindex4 <- rndIndex[(3*multiplier+1):(4*multiplier)]
missindex5 <- rndIndex[(4*multiplier+1):n_samples]

All_missIndex <- list(missindex1, missindex2, missindex3, missindex4, missindex5)

#define fitness function
fitness <- function(allele) {
  # Extract number of MGs to use
  nMGs <- round(allele[1])
  # Constraints on number of MGs
  if(nMGs < minMGsize) nMGs <- minMGsize
  if(nMGs > numMG) nMGs <- numMG
  if(nMGs > 100) nMGs <- 100  # Hard maximum
  # Select top-ranked MGs based on allele weights
  randKey <- sort(allele[-1], index.return=TRUE, decreasing=TRUE)$ix
  selected_MGs <- randKey[1:nMGs]
  # Extract haplotype matrix for selected MGs
  X <- geno[, selected_MGs, drop=FALSE]
  # Cross-validated prediction accuracy
  trainAcc <- numeric(5)
  for(i in 1:5) {
    missindex <- All_missIndex[[i]]
    # Split into train and test
    train_X <- X[-missindex, , drop=FALSE]
    train_y <- pheno[-missindex]
    test_X <- X[missindex, , drop=FALSE]
    test_y <- pheno[missindex]
    
    # Fit model
    tryCatch({
      if(ncol(train_X) <= 10) {
        # Model with pairwise interactions (captures epistasis)
        fit <- lm(train_y ~ .^2, data=as.data.frame(train_X))
      } else {
        # Too many MGs - use additive model
        fit <- lm(train_y ~ ., data=as.data.frame(train_X))
      }
      
      # Predict on test set
      pred <- predict(fit, newdata=as.data.frame(test_X))
      
      # Calculate correlation
      trainAcc[i] <- cor(test_y, pred, use="complete.obs")
      
    }, error = function(e) {
      # Fallback to simple model if error
      fit <- lm(train_y ~ ., data=as.data.frame(train_X))
      pred <- predict(fit, newdata=as.data.frame(test_X))
      trainAcc[i] <- cor(test_y, pred, use="complete.obs")
    })
  }
  
  # Mean prediction accuracy
  prediction_score <- mean(trainAcc, na.rm=TRUE)
  
  if(is.na(prediction_score) || is.nan(prediction_score)) {
    prediction_score <- 0
  }
   
  fitness_value <-  prediction_score
  # Return: [nMGs, fitness]
  crit <- c(nMGs, fitness_value)
  
  return(crit)
}

# DE algorithm functions
# Challenge function (DE/current-to-gr-best/1 variant)
challengeDEc2gbest1 <- function(x) {
  
  # Select random indices for mutation
  index <- matrix(unlist(lapply(1:popsize, function(x) sample(c(1:popsize)[-x], 3))), 
                  popsize, 3, byrow=TRUE)
  index <- cbind(index, 0)
  
  # Create group and find best in group
  for(k in 1:popsize) {
    tmp <- 1:popsize
    tmp <- tmp[-index[k,]]
    groupindex <- sample(tmp, ceiling(0.15*popsize), FALSE)
    bestGindex <- groupindex[which(fit[groupindex] == max(fit[groupindex]))]
    index[k,4] <- bestGindex[1]
  }
  
  # Crossover probability matrix
  crtf <- matrix(CR > runif(allelesize*popsize), popsize, allelesize)
  
  # Find p-best for crossover
  p <<- ceiling((popsize/2) * (1 - ((x-1)/(numgen))))
  if(p < 5) p <- 5
  
  # Top-ranked vectors
  p_top_vectorIndex <- order(fit, decreasing=TRUE)[1:p]
  challenger <- pop[sample(p_top_vectorIndex, popsize, TRUE), ]
  
  # Mutation
  rndNum <- sample(1:popsize, popsize, FALSE)
  odd <- rndNum[1:round(0.6*popsize, digits=0)]
  even <- rndNum[(round(0.6*popsize, digits=0)+1):popsize]
  
  hold1 <- pop[index[odd,1],] - FR[odd] * (pop[index[odd,4],] - pop[index[odd,1],] + 
                                             pop[index[odd,2],] - pop[index[odd,3],])
  hold2 <- pop[index[even,1],] + FR[even] * (pop[index[even,4],] - pop[index[even,1],] + 
                                               pop[index[even,2],] - pop[index[even,3],])
  hold <- rbind(hold1, hold2)
  
  # Crossover
  challenger[crtf] <- hold[crtf]
  
  return(challenger)
}

# Run generation function
rungen <- function(x) {
  
  # Create challengers
  challenger <- challengeDEc2gbest1(x)
  
  # Calculate fitness of challengers
  findFitness <- apply(challenger, 1, fitness)
  challenger[,1] <- findFitness[1,]
  fitchal <- findFitness[2,]
  
  # Replace if better or equal
  index <- which(fitchal >= fit)
  fit[index] <<- fitchal[index]
  pop[index,] <<- challenger[index,]
  
  # Update FR (mutation factor) for next generation
  Fsuccess <- FR[index]
  if(length(Fsuccess) > 0) {
    meanPowFsuccess <- sum((Fsuccess / length(Fsuccess))^1.5)
    Wf <- 0.8 + 0.2 * runif(1, 0, 1)
    Fm <<- Wf * Fm + (1 - Wf) * meanPowFsuccess
  }
  
  FR <<- 1 - rcauchy(popsize, location=Fm, scale=0.1)
  
  # Ensure FR in valid range [0, 1]
  FR_index <- c(which(FR >= 1), which(FR <= 0))
  if(length(FR_index) > 0) {
    for(l in 1:length(FR_index)) {
      tmp <- FR[FR_index[l]]
      while(tmp >= 1 || tmp <= 0) {
        tmp <- 1 - rcauchy(1, location=Fm, scale=0.1)
      }
      FR[FR_index[l]] <<- tmp
    }
  }
  
  # Update CR (crossover rate) for next generation
  Csuccess <- CR[index]
  if(length(Csuccess) > 0) {
    meanPowCsuccess <- sum((Csuccess / length(Csuccess))^0.67)
    Wc <- 0.9 + 0.1 * runif(1, 0, 1)
    Cm <<- Wc * Cm + (1 - Wc) * meanPowCsuccess
  }
  
  CR <<- rnorm(popsize, mean=Cm, sd=0.1)
  
  # Restart if population converges prematurely
  if(100 * sd(fit) / mean(fit) < 1) {
    maxindex <- which(fit == max(fit))[1]
    tmp <- sample(-1000:1000, (popsize-1), FALSE)
    pop[-maxindex, 1] <<- pop[-maxindex, 1] + tmp
    fit <<- apply(pop, 1, fitness)[2,]
  }
  
  # Save progress periodically (every 20 generations)
  if(x %% 20 == 0) {
    bestsol <- pop[which(fit == max(fit))[1], ]
    nMG_best <- round(bestsol[1])
    
    if(nMG_best < minMGsize) nMG_best <- minMGsize
    if(nMG_best > numMG) nMG_best <- floor(numMG * 0.2)
    
    randKey <- sort(bestsol[-1], index.return=TRUE, decreasing=TRUE)$ix
    selected_best <- sort(randKey[1:nMG_best])
    
    saveRDS(selected_best, paste0(path, "/selected_MGs_gen", x, ".rds"))
    saveRDS(fit, paste0(path, "/fit_gen", x, ".rds"))
    
    cat(sprintf("[CHECKPOINT] Generation %d | Best fitness: %.4f | nMGs: %d\n", 
                x, max(fit), nMG_best))
  }
  
  # Print progress every 5 generations
  if(x %% 5 == 0) {
    bestsol <- pop[which(fit == max(fit))[1], ]
    cat(sprintf("Gen %3d | Fitness: %.4f | nMGs: %3d | Mean fit: %.4f\n", 
                x, max(fit), round(bestsol[1]), mean(fit)))
  }
  
  generation <<- generation + 1
  MaxFits[x] <<- max(fit)
}

#initialise DE parameters
cat("\n=== Initializing DE Parameters ===\n")

# Strain identifier
strain <- "G7"

# Create output directory
dirOut <- paste0("DE_output_", strain)
if(!dir.exists(dirOut)) {
  dir.create(dirOut)
  cat("Created output directory:", dirOut, "\n")
}
path <- dirOut

geno <-  geno

# DE parameters - OPTIMIZED FOR ~1 HOUR RUNTIME
numgen <- 200           # 200 generations
popsize <- 50           # Population size
numMG <- ncol(geno)
allelesize <- numMG + 1 # One for nMG, rest for MG weights
minMGsize <- 15         # Minimum MGs to select
target_size <- 50       # Target ~50 MGs

cat("DE Parameters:\n")
cat("  Generations:", numgen, "\n")
cat("  Population size:", popsize, "\n")
cat("  Total marker groups:", numMG, "\n")
cat("  Minimum MGs:", minMGsize, "\n")
cat("  Target MGs:", target_size, "\n")
cat("  Maximum MGs per chromosome:", 15, "\n")
cat("  Expected runtime: ~1 hour\n")
# Initialize population
generation <- 1
pop <- matrix(runif(popsize * allelesize), popsize, allelesize)
# Calculate initial fitness
cat("Calculating initial fitness...\n")
fit <- apply(pop, 1, fitness)[2,]
cat("Initial fitness calculated\n")
cat("  Max fitness:", round(max(fit), 4), "\n")
cat("  Mean fitness:", round(mean(fit), 4), "\n")
cat("  Min fitness:", round(min(fit), 4), "\n")
# Initialize tracking
MaxFits <- numeric(numgen)
MaxFits[1] <- max(fit)
# Initialize FR (mutation factor)
Fm <- 0.5
FR <- rcauchy(popsize, location=Fm, scale=0.1)
index_FR <- c(which(FR >= 1), which(FR <= 0))
if(length(index_FR) > 0) {
  for(l in 1:length(index_FR)) {
    tmp <- FR[index_FR[l]]
    while(tmp >= 1 || tmp <= 0) {
      tmp <- 1 - rcauchy(1, location=Fm, scale=0.1)
    }
    FR[index_FR[l]] <- tmp
  }
}
# Initialize CR (crossover rate)
Cm <- 0.6
CR <- rnorm(popsize, mean=Cm, sd=0.1)
# Run DE optimisation
lapply(1:numgen, rungen)
#extract and save final results
cat("=== Extracting Final Results ===\n")
# Get best solution
bestsol <- pop[which(fit == max(fit))[1], ]
nMG_final <- round(bestsol[1])
if(nMG_final < minMGsize) nMG_final <- minMGsize
if(nMG_final > numMG) nMG_final <- floor(numMG * 0.2)
# Get selected MG indices
randKey <- sort(bestsol[-1], index.return=TRUE, decreasing=TRUE)$ix
selected_MGs_final <- sort(randKey[1:nMG_final])
cat("Selected", nMG_final, "marker groups (target was ~", target_size, ")\n")
# Save final results
saveRDS(selected_MGs_final, paste0(path, "/selected_MGs_FINAL_", strain, ".rds"))
saveRDS(bestsol, paste0(path, "/best_solution_FINAL_", strain, ".rds"))
saveRDS(fit, paste0(path, "/fit_FINAL_", strain, ".rds"))
saveRDS(pop, paste0(path, "/pop_FINAL_", strain, ".rds"))
saveRDS(MaxFits, paste0(path, "/MaxFits_", strain, ".rds"))
