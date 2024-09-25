# This script is for getting control probe PC results

###############################################################################
#### Get control probe PCA results ####
###############################################################################

control_probe_PC <- function(rgSet, sample_identifier_column, mSet=NULL, variance_threshold = 0.9){
  
  # Get control probes
  control_probes <- getControlAddress(rgSet)
  num_samples <- length(rgSet$sample_identifier_column)
  
  
  
  # Estimate beta values
  red_intensities <- getRed(rgSet)[control_probes,]
  green_intensities <- getGreen(rgSet)[control_probes,]
  beta_values <- red_intensities / (red_intensities + green_intensities + 100)
  beta_values <- na.omit(beta_values)
  
  # Run PCA
  # pc <- pca(t(beta_values), nPcs = num_samples, method = "ppca")
  pca_result <- prcomp(t(beta_values), center = TRUE, scale. = TRUE)
  
  # Calculate the proportion of variance explained by each principal component
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  # Get explained variance
  # sd_pca <- sDev(pc)
  # sd_pca <- sd_pca^2 / sum(sd_pca^2)
  
  # Find the first principal component where cumulative variance reaches 90%
  cutoff <- which(cumsum(var_explained) >= variance_threshold)[1]
  
  # # Get PC scores corresponding to the top components explaining 90% variance
  # pca_df <- as.data.frame(scores(pc)[, 1:cutoff])
  # 
  # # Add sample names to the dataframe
  # pca_df$Sample <- rownames(pca_df)
  # 
  # # Reorder so that sample name is the first column
  # pca_df <- pca_df[, c(ncol(pca_df), 1:(ncol(pca_df) - 1))]
  
  # Get PC scores corresponding to the top components explaining variance_threshold variance
  control_probe_PCA_results <- as.data.frame(pca_result$x[, 1:cutoff])
  
  # Add sample names to the dataframe
  control_probe_PCA_results$Sample <- rownames(control_probe_PCA_results)
  
  # Reorder so that sample name is the first column
  control_probe_PCA_results <- control_probe_PCA_results[, c(ncol(control_probe_PCA_results), 1:(ncol(control_probe_PCA_results) - 1))]
  
  return(control_probe_PCA_results)
}