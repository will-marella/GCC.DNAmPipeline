# This script contains functions related to probe filtering and normalization

###############################################################################
#### Run Quantro Anova Test ####
###############################################################################

run_quantro_anovatest_across_variables <- function(mSet, phenotype_variables){
  results <- list()  # To store p-values for each phenotype variable
  
  for (variable in phenotype_variables) {
    # Extract the phenotype variable as a factor
    groupFactor <- pData(mSet)[[variable]]
    
    # Try-catch block to handle potential errors during quantro execution
    tryCatch({
      # Run quantro on the phenotype variable
      quantro_result <- quantro(mSet, groupFactor = groupFactor)
      
      # Extract the ANOVA results (summary table)
      anova_pval <- quantro_result@anova$`Pr(>F)`[1]
      
      # Store the p-value in the results list
      results[[variable]] <- anova_pval
      
    }, error = function(e) {
      cat(paste("Error in variable", variable, ":", e$message, "\n"))
      results[[variable]] <- NA  # Store NA for failed variables
    })
  }
  
  # Present results to the user
  anova_pvals_df <- data.frame(Phenotype = phenotype_variables, ANOVA_P_Value = unlist(results))
  cat("Statistically Significant ANOVA p-values suggest FunNorm is warranted\n")
  print(anova_pvals_df)
  
  return(anova_pvals_df)  # Return the data frame for further use
}


###############################################################################
#### Run Quantro Permutation Test ####
###############################################################################

run_quantro_permtest_across_variables <- function(mSet, phenotype_variables, num_permutations=1000){
  
  ##
  ## This will need to be synced up to doParallel on Ronin to use cores efficiently
  ## For many samples this is very time intensive, but it can be done in parallel
  ##
  
  results <- list()  # To store p-values for each phenotype variable
  
  for (variable in phenotype_variables) {
    # Extract the phenotype variable as a factor
    groupFactor <- pData(mSet)[[variable]]
    
    if (any(table(groupFactor) < 2)) {
      cat(paste("Skipping variable", variable, "- not enough samples in one or more groups\n"))
      results[[variable]] <- NA
      next
    }
    
    # Try-catch block to handle potential errors during quantro execution
    tryCatch({
      # Run quantro on the phenotype variable with permutation testing
      quantro_result <- quantro(mSet, groupFactor = groupFactor, B = num_permutations)
      
      # Extract the permutation-based p-value
      perm_pval <- quantro_result@quantroPvalPerm
      
      # Store the p-value in the results list
      results[[variable]] <- perm_pval
      
    }, error = function(e) {
      cat(paste("Error in variable", variable, ":", e$message, "\n"))
      results[[variable]] <- NA  # Store NA for failed variables
    })
    
  }
  
  # Present results to the user
  perm_pvals_df <- data.frame(Phenotype = phenotype_variables, perm_P_Value = unlist(results))
  cat("Statistically Significant Permutation Test p-values suggest FunNorm is warranted\n")
  print(perm_pvals_df)
  
  return(perm_pvals_df)  # Return the data frame for further use
}

###############################################################################
#### Run FunNorm ####
###############################################################################

?preprocessNoob

perform_Funnorm <- function(rgSet, sex_adjustment=FALSE, sex_column="Sex") {
  
  if (sex_adjustment == TRUE){grSet <- preprocessFunnorm(rgSet, sex = pData(rgSet)$sex_column)}
  if (sex_adjustment == FALSE){grSet <- preprocessFunnorm(rgSet)}
  
  return(grSet)
}

###############################################################################
#### Probe Filtering ####
###############################################################################

perform_probe_filtering <- function(grSet, rgSet, rgSet_extended, array_version){
  
  # Print the number of initial probes
  cat("Number of probes before filtering:", nrow(grSet), "\n")
  
  print("Entering Step 1")
  
  # 1. Remove SNP Probes
  grSet <- dropLociWithSnps(grSet, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)
  
  cat("Number of probes after filtering SNP probes:", nrow(grSet), "\n")
  
  print("Entering Step 2")
  
  # 2. Remove XY Probes
  # Load the EPIC manifest data
  data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  epicAnnotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  xyProbes <- rownames(epicAnnotation)[epicAnnotation$chr %in% c("chrX", "chrY")]
  
  # Remove XY probes
  probeNames <- featureNames(grSet)
  keep <- !(probeNames %in% xyProbes)
  grSet <- grSet[keep,]
  
  cat("Number of probes after filtering XY probes:", nrow(grSet), "\n")
  
  print("Entering Step 3")
  
  # 3. Remove Pidsley Cross-Reactive Probes
  if (array_version == "450k"){
    # These are the benton probes, the BOWTIE2 multi-map
    # Downloaded from sirselim github: HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt
    # That's benton's github
    benton <- "Data/bentonCR_450k.txt"
    benton_probes <- read.delim(benton, header = FALSE, sep = "\n")
    benton_probes <- as.character(benton_probes$V1)
    
    cat("# of benton probes", length(benton_probes), "\n")
    
    # These are the Chen et al. probes, which I also took from Benton's github
    # Downloaded as: 48639-non-specific-probes-Illumina450k.csv
    chen <- "Data/chenCR_450k.csv"
    chen_probes <- read.csv(chen)
    chen_probes <- as.character(chen_probes$TargetID)
    
    cat("# of chen probes", length(chen_probes), "\n")
    
    # Combine probes
    CR_probes <- unique(c(benton_probes, chen_probes))
    cat("# of unique probes", length(CR_probes), "\n")
  }
  
  if (array_version == "850k") {
    # This is: EPIC/13059_2016_1066_MOESM1_ESM.csv
    # From Pidsley et al. -- According to miles benton, mark chen
    pidsley_probes <- read.csv("Data/pidsleyCR_EPIC850k.csv", stringsAsFactors = FALSE, header = TRUE)
    pidsley_probes <- as.character(pidsley_probes[, 1])
    
    # This is: 1-s2.0-S221359601630071X-mmc2.txt
    # And: 1-s2.0-S221359601630071X-mmc3.txt
    # From Mccartney et al, This is how Mark Chen removed them
    # Includes CPG and nonCPG cr probes
    mccartney_CR <- readLines("Data/mccartney-mmc2_EPIC850k.txt")
    mccartney_CR <- as.character(mccartney_CR)
    mccartney_CR_nonCPG <- readLines("Data/mccartney-mmc3_EPIC850k.txt")
    mccartney_CR_nonCPG <- as.character(mccartney_CR_nonCPG)
    
    # CR probes will combine [id]
    CR_probes <- unique(c(pidsley_probes, mccartney_CR, mccartney_CR_nonCPG))
    
    cat("# of pidsley probes", length(pidsley_probes), "\n")
    cat("# of mccartney CPG probes", length(mccartney_CR), "\n")
    cat("# of mccartney non-CPG probes", length(mccartney_CR_nonCPG), "\n")
    cat("# of all unique probes", length(CR_probes), "\n")
    
  }
  
  probeNames <- featureNames(grSet)
  keep <- !(probeNames %in% CR_probes)
  grSet <- grSet[keep,]
  
  cat("Number of probes after filtering Pidsley Cross-Reactive Probes:", nrow(grSet), "\n")
  
  print("Entering Step 5")
  
  # 5. Remove probes with beadcount <3 in 5% of samples
  bc <- beadcount(rgSet_extended)
  bc_grSet <- bc[rownames(grSet), ]
  keep <- rowSums(is.na(bc_grSet) | bc_grSet < 3) <= 0.05 * ncol(bc_grSet)
  grSet <- grSet[keep,]
  
  cat("Number of probes after filtering by beadcount:", nrow(grSet), "\n")
  
  print("Entering Step 6")
  
  # 6. Remove probes with 1% of samples having a detection p-value > 0.01
  detP <- minfi::detectionP(rgSet)
  detP_grSet <- detP[rownames(grSet), ]
  keep <- rowSums(detP_grSet > 0.01) <= 0.01 * ncol(detP_grSet)
  grSet <- grSet[keep,]
  
  cat("Number of probes after filtering by detection p-values:", nrow(grSet), "\n")
  
  # Print the number of remaining probes
  cat("Number of probes after filtering:", nrow(grSet), "\n")
  
  return(grSet)
  
}



###############################################################################
#### Probe Type Normalization / BMIQ ####
###############################################################################

probe_type_normalization_BMIQ <- function(grSet, do_parallel = FALSE, n_cores = parallel::detectCores() - 1){

  # Pull off betas and remove NAs
  betas <-getBeta(grSet)
  betas <-na.omit(betas)
  
  # Pull off annotation data for probes present in RGset after removing rows with missing values
  norm_annot <-getAnnotation(na.omit(grSet))
  probe_type <-norm_annot[rownames(norm_annot) %in% rownames(betas), "Type"]
  
  # Change name so "I" in annotation is "1" and the rest are "2". 
  probe_type <-if_else(probe_type == "I", '1', '2')
  
  # Check the number of cores and decide whether to run in parallel or sequentially
  if ((n_cores > 1) && isTRUE(do_parallel)) {
    # Set up parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Ensure the cluster is stopped after function execution
    on.exit({
      stopCluster(cl)
      registerDoSEQ()
    }) 
    
    # Parallelize BMIQ normalization for each sample
    BMIQ_results_list <- foreach(col = colnames(betas), .combine = "cbind", 
                                 .packages = c("minfi", "wateRmelon")) %dopar% {
                                   tryCatch({
                                     result <- BMIQ(as.data.frame(betas[, col]), design.v = probe_type, nfit = 10000, plots=FALSE )$nbeta
                                     return(result)
                                   }, error = function(e) {
                                     message(paste("Error processing sample", col, ":", e$message))
                                     return(rep(NA, nrow(betas))) # Return NAs if there is an error
                                   })
                                 }  
    
    filtered_normalized_beta <- as.matrix(BMIQ_results_list)
    colnames(filtered_normalized_beta) <- colnames(betas)
    rownames(filtered_normalized_beta) <- rownames(betas)
    
    # Stop parallel processing
    stopCluster(cl)
    registerDoSEQ()
  } else {
    # Run sequentially if only one core is available
    filtered_normalized_beta <- do.call(cbind, lapply(colnames(betas), function(col) {
      tryCatch({
        BMIQ(as.data.frame(betas[, col]), design.v = probe_type, nfit = 10000, plots=FALSE )$nbeta
      }, error = function(e) {
        message(paste("Error processing sample", col, ":", e$message))
        return(rep(NA, nrow(betas))) # Return NAs if there is an error
      })
    }))
    
    # Assign names to the matrix
    colnames(filtered_normalized_beta) <- colnames(betas)
    rownames(filtered_normalized_beta) <- rownames(betas)
  }
  
  filtered_normalized_beta <- as.matrix(filtered_normalized_beta)
  
  # Density plot of pre- and post-normalization distributions
  par(mfrow=c(1,2))
  densityPlot(betas, main="Before BMIQ Normalization")
  densityPlot(filtered_normalized_beta, main="After BMIQ Normalization")
  
  # Check the correlation between samples before and after
  cor_before <- cor(betas)
  cor_after <- cor(filtered_normalized_beta)
  
  #Visualize the correlation difference
  par(mfrow=c(1,2))
  image(cor_before, main="Correlation Before BMIQ")
  image(cor_after, main="Correlation After BMIQ")
  
  return(filtered_normalized_beta)
  
}

?BMIQ

