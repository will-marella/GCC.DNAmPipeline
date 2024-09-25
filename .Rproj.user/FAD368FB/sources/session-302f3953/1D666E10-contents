# This script contains the functions required to get QC results
# For plotting functions, see plot_quality_control.R

###############################################################################
#### Run control metric QC ####
###############################################################################


control_metric_QC <- function(idat_dir){
  
  read_idats_fixed <- function(idat_dir) {
    # Save current working directory
    original_wd <- getwd()
    
    # Extract prefixes
    idat_files <- list.files(idat_dir, pattern = "\\.idat$", full.names = TRUE)
    prefixes <- unique(gsub("_(Grn|Red)\\.idat$", "", basename(idat_files)))
    
    # Change to the directory containing the idat files
    setwd(idat_dir)
    
    # Run the read_idats function
    data_for_control_metrics <- read_idats(prefixes)
    
    # Reset to original working directory
    setwd(original_wd)
    
    # Return the data
    return(data_for_control_metrics)
  }
  
  ## Using read_idats_fixed()
  # Read IDATS
  data_for_control_metrics <- read_idats_fixed(idat_dir)
  
  # Find control metrics with ewastools
  control_test_values = control_metrics(data_for_control_metrics)
  
  # Get sample names from one of the metrics
  sample_names <- data_for_control_metrics$meta$sample_id
  
  # Initialize empty data frames with sample names as columns and metric names as rows
  passed_QC <- data.frame(matrix(ncol = length(control_test_values), nrow = length(sample_names)))
  rownames(passed_QC) <- sample_names
  colnames(passed_QC) <- names(control_test_values)
  
  QC_result_values <- data.frame(matrix(ncol = length(control_test_values), nrow = length(sample_names)))
  rownames(QC_result_values) <- sample_names
  colnames(QC_result_values) <- names(control_test_values)
  
  # Loop through each metric to check pass/fail against threshold and calculate difference
  for (metric_name in names(control_test_values)){
    metric_value <- control_test_values[[metric_name]]
    threshold <- attr(metric_value, "threshold")
    
    # Check pass/fail status for each sample based on threshold
    pass_fail <- metric_value >= threshold
    # Calculate the difference from threshold
    differences <- metric_value - threshold
    
    # Add results to data frames
    passed_QC[, metric_name] <- pass_fail
    QC_result_values[, metric_name] <- metric_value
    
    # Set threshold as an attribute for each column
    attr(QC_result_values[, metric_name], "threshold") <- threshold
  }
  
  # Return both data frames as a list
  return(list(
    passed_QC = passed_QC,
    QC_result_values = QC_result_values
  ))
  
}

###############################################################################
#### Run log2 Intensity QC ####
###############################################################################


log2_intensity_QC <- function(passed_QC, QC_result_values, mSet, threshold = 10.5) {
  
  # Calculate log2 median intensities
  log2_intensities <- log2(cbind(
    Meth = apply(getMeth(mSet), 2, median),
    Unmeth = apply(getUnmeth(mSet), 2, median)
  ))
  
  # Create a data frame with the results
  log2_results <- data.frame(
    Sample = colnames(mSet),
    log2_Meth_intensity = log2_intensities[, "Meth"],
    log2_Unmeth_intensity = log2_intensities[, "Unmeth"],
    stringsAsFactors = FALSE
  )
  
  # Calculate pass/fail status
  log2_results$Passed_Log2Intensity_QC <- apply(log2_intensities >= threshold, 1, all)
  
  # Append results to passed_QC
  passed_QC$log2_Meth_intensity <- log2_results$Passed_Log2Intensity_QC
  passed_QC$log2_Unmeth_intensity <- log2_results$Passed_Log2Intensity_QC
  
  # Append results to QC_result_values
  QC_result_values$log2_Meth_intensity <- log2_results$log2_Meth_intensity
  QC_result_values$log2_Unmeth_intensity <- log2_results$log2_Unmeth_intensity
  
  # Set threshold as an attribute for the new columns
  attr(QC_result_values$log2_Meth_intensity, "threshold") <- threshold
  attr(QC_result_values$log2_Unmeth_intensity, "threshold") <- threshold
  
  return(list(
    passed_QC = passed_QC,
    QC_result_values = QC_result_values
  ))
}


###############################################################################
#### Calculate Detection P-Values ####
###############################################################################

avg_detection_p_QC <- function(passed_QC, QC_result_values, rgSet, threshold = 0.005) {
  
  # Calculate detection p-values
  detP <- minfi::detectionP(rgSet)
  avg_detP <- colMeans(detP)
  
  # Create a data frame with the results
  detP_results <- data.frame(
    Sample = names(avg_detP),
    avg_detection_p = avg_detP,
    stringsAsFactors = FALSE
  )
  
  # Calculate pass/fail status
  detP_results$Passed_AvgDetectionP_QC <- avg_detP <= threshold
  
  # Append results to passed_QC
  passed_QC$AvgDetectionP <- detP_results$Passed_AvgDetectionP_QC
  
  # Append results to QC_result_values
  QC_result_values$AvgDetectionP <- detP_results$avg_detection_p
  
  # Set threshold as an attribute for the new column
  attr(QC_result_values$AvgDetectionP, "threshold") <- threshold
  
  return(list(
    passed_QC = passed_QC,
    QC_result_values = QC_result_values
  ))
}

###############################################################################
#### Get BSCON Results ####
###############################################################################


bscon_QC <- function(passed_QC, QC_result_values, rgSet, threshold = 80) {
  # Perform bscon on rgSet
  bscon_result <- bscon(rgSet)
  
  # Create a data frame with the results
  bscon_results <- data.frame(
    Sample = names(bscon_result),
    bscon_value = bscon_result,
    stringsAsFactors = FALSE
  )
  
  # Calculate pass/fail status
  bscon_results$Passed_bscon_QC <- bscon_result >= threshold
  
  # Append results to passed_QC
  passed_QC$BSCon <- bscon_results$Passed_bscon_QC
  
  # Append results to QC_result_values
  QC_result_values$BSCon <- bscon_results$bscon_value
  
  # Set threshold as an attribute for the new column
  attr(QC_result_values$BSCon, "threshold") <- threshold
  
  return(list(
    passed_QC = passed_QC,
    QC_result_values = QC_result_values
  ))
}


###############################################################################
#### Run all non-FDR Sample QC ####
###############################################################################

get_sample_QC_results <- function(idat_dir, rgSet, mSet=NULL, log2intensity_threshold=10.5,
                                  detP_threshold=0.005, bscon_threshold=80){
  
  if (is.null(mSet)){mSet <- preprocessRaw(rgSet)}
  
  cat("Running control metric QC:\n")
  QC_results <- control_metric_QC(idat_dir)
  passed_QC <- QC_results$passed_QC
  QC_result_values <- QC_results$QC_result_values
  
  cat("Running log2 intensity QC:\n")
  QC_results <- log2_intensity_QC(passed_QC, QC_result_values, mSet, threshold=log2intensity_threshold)
  passed_QC <- QC_results$passed_QC
  QC_result_values <- QC_results$QC_result_values
  
  cat("Running detection p value QC:\n")
  QC_results <- avg_detection_p_QC(passed_QC, QC_result_values, rgSet, threshold=detP_threshold)
  passed_QC <- QC_results$passed_QC
  QC_result_values <- QC_results$QC_result_values
  
  cat("Running bisulfit conversion QC:\n")
  QC_results <- bscon_QC(passed_QC, QC_result_values, rgSet, threshold = bscon_threshold)
  
  return(QC_results)
  
}








###############################################################################
#### Get Global FDR results ####
###############################################################################

get_fdr_results <- function(any_set){
  
  glob_outliers_all <- data.frame(Sample = character(), Iteration = integer(), stringsAsFactors = FALSE)
  all_samples <- colnames(getBeta(any_set))
  iter <- 1
  pca_result_all <- NULL
  
  find_outliers <- function(any_set, iter){
    
    # Perform PCA and collect explained variance
    beta_values <- getBeta(any_set)
    beta_values <- beta_values[complete.cases(beta_values), ]
    beta_values_t <- t(beta_values)  # Transpose so rows are samples, columns are CpGs
    pca_result <- prcomp(beta_values_t, center = TRUE, scale. = TRUE)
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    
    if (is.null(pca_result_all)) {
      # Store the first PCA result for all samples
      pca_result_all <<- pca_result
    }
    
    # Calculate mean of PC1 and derive z-scores for samples (not CpGs)
    x = 1
    PC1 <- pca_result$x[, x]
    PC1_mean <- mean(PC1)
    squared_distances <- (PC1 - PC1_mean)^2
    z_scores <- scale(squared_distances)
    
    # Get sample names
    sample_names <- rownames(pca_result$x)  # Colnames correspond to samples after transposing
    
    # Convert into global FDRs
    p_values <- pnorm(z_scores, lower.tail = FALSE)
    adj_p_values <- p.adjust(p_values, method = "BH")
    
    # Find global FDR outliers
    glob_outliers <- sample_names[which(adj_p_values < 0.2)]
    
    # Output results
    cat("Iteration:", iter, "\n")
    cat("GlobFDR outliers:", paste(glob_outliers, collapse = ", "), "\n")
    
    # Append outliers with iteration number to the global data frame
    if(length(glob_outliers) > 0) {
      glob_outliers_all <<- rbind(glob_outliers_all, data.frame(Sample = glob_outliers, Iteration = iter))
    }
    
    # Remove outlier samples from the data
    outliers_to_remove <- unique(glob_outliers)
    any_set <- any_set[, !(colnames(getBeta(any_set)) %in% outliers_to_remove)]  # Remove samples (columns)
    
    # Return updated data and outliers found in this iteration
    return(list(any_set = any_set, glob_outliers = glob_outliers))
  }
  
  # Continue running the outlier detection until no more outliers are found
  while(TRUE){
    outliers_found <- find_outliers(any_set, iter)
    
    any_set <- outliers_found$any_set
    if (length(outliers_found$glob_outliers) == 0) {
      cat("No more outliers found. Stopping at iteration:", iter, "\n")
      break
    }
    
    iter <- iter + 1
  }
  
  # Return the lists of outliers, PCA results, and the final filtered dataset
  return(list(
    outliers = glob_outliers_all,
    filtered_set = any_set,
    pca_result_all = pca_result_all,
    all_samples = all_samples
  ))
}


# get_fdr_results <- function(grSet_filtered){
#   
#   glob_outliers_all <- data.frame(Sample = character(), Iteration = integer(), stringsAsFactors = FALSE)
#   all_samples <- colnames(getBeta(grSet_filtered))
#   iter <- 1
#   pca_result_all <- NULL
#   
#   find_outliers <- function(grSet_filtered, iter){
#     
#     # Perform PCA and collect explained variance
#     beta_values <- getBeta(grSet_filtered)
#     beta_values <- beta_values[complete.cases(beta_values), ]
#     beta_values_t <- t(beta_values)  # Transpose so rows are samples, columns are CpGs
#     pca_result <- prcomp(beta_values_t, center = TRUE, scale. = TRUE)
#     var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
#     
#     if (is.null(pca_result_all)) {
#       # Store the first PCA result for all samples
#       pca_result_all <<- pca_result
#     }
#     
#     # Calculate mean of PC1 and derive z-scores for samples (not CpGs)
#     x = 1
#     PC1 <- pca_result$x[, x]
#     PC1_mean <- mean(PC1)
#     squared_distances <- (PC1 - PC1_mean)^2
#     z_scores <- scale(squared_distances)
#     
#     # Get sample names
#     sample_names <- rownames(pca_result$x)  # Colnames correspond to samples after transposing
#     
#     # Convert into global FDRs
#     p_values <- pnorm(z_scores, lower.tail = FALSE)
#     adj_p_values <- p.adjust(p_values, method = "BH")
#     
#     # Find global FDR outliers
#     glob_outliers <- sample_names[which(adj_p_values < 0.2)]
#     
#     # Output results
#     cat("Iteration:", iter, "\n")
#     cat("GlobFDR outliers:", paste(glob_outliers, collapse = ", "), "\n")
#     
#     # Append outliers with iteration number to the global data frame
#     if(length(glob_outliers) > 0) {
#       glob_outliers_all <<- rbind(glob_outliers_all, data.frame(Sample = glob_outliers, Iteration = iter))
#     }
#     
#     # Remove outlier samples from the data
#     outliers_to_remove <- unique(glob_outliers)
#     grSet_filtered <- grSet_filtered[, !(colnames(getBeta(grSet_filtered)) %in% outliers_to_remove)]  # Remove samples (columns)
#     
#     # Return updated data and outliers found in this iteration
#     return(list(grSet_filtered = grSet_filtered, glob_outliers = glob_outliers))
#   }
#   
#   # Continue running the outlier detection until no more outliers are found
#   while(TRUE){
#     outliers_found <- find_outliers(grSet_filtered, iter)
#     
#     grSet_filtered <- outliers_found$grSet_filtered
#     if (length(outliers_found$glob_outliers) == 0) {
#       cat("No more outliers found. Stopping at iteration:", iter, "\n")
#       break
#     }
#     
#     iter <- iter + 1
#   }
#   
#   # Return the lists of outliers, PCA results, and the final filtered dataset
#   return(list(
#     glob_outliers_all = glob_outliers_all,
#     final_grSet_filtered = grSet_filtered,
#     pca_result_all = pca_result_all,
#     all_samples = all_samples
#   ))
# }

###############################################################################
#### FDR with local FDR... doesn't work well on this data####
###############################################################################


##
## KEEPING THIS IN CASE I WANT TO USE LOCFDR, BUT RN IT SEEMS IFFY
## It didn't work with this data, but it might work with others
##

# library(locfdr)
# plot_fdr <- function(grSet_filtered){
#   
#   loc_outliers_all <- list()
#   glob_outliers_all <- list()
#   iter <- 1
#   
#   find_outliers <- function(grSet_filtered, iter){
#     
#     # Perform PCA and collect explained variance
#     beta_values <- getBeta(grSet_filtered)
#     beta_values <- beta_values[complete.cases(beta_values), ]
#     beta_values_t <- t(beta_values)  # Transpose so rows are samples, columns are CpGs
#     pca_result <- prcomp(beta_values_t, center = TRUE, scale. = TRUE)
#     var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
#     
#     # Calculate mean of PC1 and derive z-scores for samples (not CpGs)
#     x = 1
#     PC1 <- pca_result$x[, x]
#     PC1_mean <- mean(PC1)
#     squared_distances <- (PC1 - PC1_mean)^2
#     z_scores <- scale(squared_distances)
#     
#     # Convert into local FDRs
#     locfdr_results <- locfdr(z_scores)
#     locfdr_values <- locfdr_results$fdr  # Extract locfdr values
#     
#     # Get sample names
#     sample_names <- colnames(beta_values)  # Colnames correspond to samples after transposing
#     
#     # Find local FDR outliers based on threshold
#     loc_outliers <- sample_names[which(locfdr_values < 0.2)]
#     
#     # Convert into global FDRs
#     p_values <- pnorm(z_scores, lower.tail = FALSE)
#     adj_p_values <- p.adjust(p_values, method = "BH")
#     
#     # Find global FDR outliers
#     glob_outliers <- sample_names[which(adj_p_values < 0.2)]
#     
#     # Output results
#     cat("Iteration:", iter, "\n")
#     cat("LocFDR outliers:", loc_outliers, "\n")
#     cat("GlobFDR outliers:", glob_outliers, "\n")
#     
#     # Append outliers with iteration number
#     loc_outliers_all[[iter]] <- list(outliers = loc_outliers, iteration = iter)
#     glob_outliers_all[[iter]] <- list(outliers = glob_outliers, iteration = iter)
#     
#     # Remove outlier samples from the data
#     outliers_to_remove <- unique(c(loc_outliers, glob_outliers))
#     grSet_filtered <- grSet_filtered[, !(colnames(getBeta(grSet_filtered)) %in% outliers_to_remove)]  # Remove samples (columns)
#     
#     # Return updated data and outliers found in this iteration
#     return(list(grSet_filtered = grSet_filtered, loc_outliers = loc_outliers, glob_outliers = glob_outliers))
#   }
#   
#   # Continue running the outlier detection until no more outliers are found
#   while(TRUE){
#     outliers_found <- find_outliers(grSet_filtered, iter)
#     
#     grSet_filtered <- outliers_found$grSet_filtered
#     if (length(outliers_found$loc_outliers) == 0 && length(outliers_found$glob_outliers) == 0) {
#       cat("No more outliers found. Stopping at iteration:", iter, "\n")
#       break
#     }
#     
#     iter <- iter + 1
#   }
#   
#   # Return the lists of outliers and the final filtered dataset
#   return(list(loc_outliers_all = loc_outliers_all, glob_outliers_all = glob_outliers_all, final_grSet_filtered = grSet_filtered))
# }
# 
# # Call the function
# result <- plot_fdr(grSet_filtered)
