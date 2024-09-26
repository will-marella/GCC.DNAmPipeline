# This script is for functions that perform or assist general visualizations
# Namely, beta distribution plots and PCA plots

###############################################################################
#### Plot Beta Distribution ####
###############################################################################


plot_beta_distribution <- function(any_set_or_beta_matrix, sample_identifier_column, color_by=NULL, pheno_data=NULL, n_sites = 100000) {
  
  # Error checking for input type
  if (!(is.matrix(any_set_or_beta_matrix) ||
        inherits(any_set_or_beta_matrix, c("RGChannelSet", "RGChannelSetExtended", "MethylSet", "GenomicRatioSet", "GenomicMethylSet")))) {
    stop("Input must be a matrix or a recognized set object (RGChannelSet, MethylSet, GenomicRatioSet, GenomicMethylSet).")
  }
  
  if (is.matrix(any_set_or_beta_matrix)){
    raw_betas <- any_set_or_beta_matrix
    plot_type <- "Betas"
  }
  else{
    raw_betas <- getBeta(any_set_or_beta_matrix)
    plot_type <- "Set"
  }
  
  # Select a certain # of sites for plotting feasibility
  selected_sites <- sample(1:nrow(raw_betas), n_sites)
  beta_sample <- raw_betas[selected_sites, ]
  
  if (!is.data.frame(pheno_data)) {
    pheno_data <- as.data.frame(pheno_data)
  }
  
  # Ensure phenotype data is available for matrix input
  if (plot_type == "Betas" && is.null(pheno_data)) {
    stop("Phenotype data must be provided when the input is a matrix.")
  }
  
  # Check if sample identifiers match
  if (!all(colnames(raw_betas) == rownames(pheno_data))) {
    stop("Sample identifiers in beta values do not match those in phenotype data.")
  }
  
  # Prepare data for plotting
  beta_long <- beta_sample %>%
    as.data.frame() %>%
    rownames_to_column("CpG") %>%
    pivot_longer(-CpG, names_to = sample_identifier_column, values_to = "beta_value") %>%
    filter(!is.na(beta_value))
  
  # Merge with phenotype data
  beta_long <- beta_long %>%
    left_join(pheno_data, by = sample_identifier_column)
  
  # Set up color aesthetics
  if (!is.null(color_by) && color_by %in% colnames(pheno_data)) {
    color_aes <- aes(color = !!sym(color_by))
    title_suffix <- paste("colored by", color_by)
    legend_position = "right"
  } else {
    color_aes <- aes(color = !!sym(sample_identifier_column))
    title_suffix <- "colored by sample"
    legend_position = "none"
  }
  
  # Generate the plot
  plot <- ggplot(beta_long, aes(x = beta_value, group = !!sym(sample_identifier_column))) +
    geom_density(linewidth = 0.5, alpha = 0.5, color_aes, show.legend = TRUE) +
    theme_gray() +
    theme(legend.position = legend_position) +
    labs(
      title = paste("Beta Value Distribution (", n_sites, " random sites),", title_suffix),
      x = "Beta Value",
      y = "Density"
    ) +
    scale_x_continuous(limits = c(0, 1))
  
  print(plot)
  return(plot)
}


###############################################################################
#### Get PCA Results ####
###############################################################################

get_PCA_results <- function(any_set){
  
  # Extract beta values
  beta_values <- getBeta(any_set)
  beta_values <- beta_values[complete.cases(beta_values), ]
  
  # Transpose the data for PCA (samples should be rows)
  beta_values_t <- t(beta_values)
  pca_result <- prcomp(beta_values_t, center = TRUE, scale. = TRUE)
  cat("Number of PCs calculated:", ncol(pca_result$x), "\n")
  
  # Calculate variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  return(list(pca_result = pca_result, var_explained = var_explained))
  
}

###############################################################################
#### Plot PCA Results ####
###############################################################################

plot_PCA <- function(pca_data, any_set, x=1, y=2, color_by=NULL, label_by=NULL){
  
  # Unpack PCA results
  pca_result <- pca_data$pca_result
  var_explained <- pca_data$var_explained
  
  # Create PCA dataframe
  pca_data <- data.frame(
    PCx = pca_result$x[, x],
    PCy = pca_result$x[, y],
    Sample = rownames(pca_result$x)
  )
  
  # If phenotype data is available, merge it to the plotting data frame
  if (ncol(pData(any_set)) > 0) {
    pca_data <- merge(pca_data, pData(any_set), by.x = "Sample", by.y = "row.names", all.x = TRUE)
  }
  
  # Create PCA plots -- Dependent on color_by and label_by
  if (!is.null(color_by) && !is.null(label_by)){
    
    pca_plot <- ggplot(pca_data, aes(x = PCx, y = PCy,
                                     color = .data[[color_by]],
                                     shape = .data[[label_by]],
                                     label = .data[[label_by]])) +
      geom_point(alpha = 0, show.legend = TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      labs(
        x = paste0("PC", x, " (", round(var_explained[x] * 100, 2), "%)"),
        y = paste0("PC", y, " (", round(var_explained[y] * 100, 2), "%)"),
        title = paste("PCA of Methylation Data: PC", x, " vs PC", y)
      ) +
      theme_gray() + 
      scale_color_manual(values = palette, name = color_by) +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (!is.null(color_by) && is.null(label_by)){
    pca_plot <- ggplot(pca_data, aes(x = PCx, y = PCy,
                                     color = .data[[color_by]])) +
      geom_point(alpha = 1, show.legend = TRUE) +
      labs(
        x = paste0("PC", x, " (", round(var_explained[x] * 100, 2), "%)"),
        y = paste0("PC", y, " (", round(var_explained[y] * 100, 2), "%)"),
        title = paste("PCA of Methylation Data: PC", x, " vs PC", y)
      ) +
      theme_gray() + 
      scale_color_manual(values = palette, name = color_by) +
      scale_shape_discrete(name = "Sample") + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && !is.null(label_by)){
    pca_plot <- ggplot(pca_data, aes(x = PCx, y = PCy,
                                     shape = .data[[label_by]],
                                     label = .data[[label_by]])) +
      geom_point(alpha = 0, show.legend = TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      labs(
        x = paste0("PC", x, " (", round(var_explained[x] * 100, 2), "%)"),
        y = paste0("PC", y, " (", round(var_explained[y] * 100, 2), "%)"),
        title = paste("PCA of Methylation Data: PC", x, " vs PC", y)
      ) +
      theme_gray() + 
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && is.null(label_by)){
    pca_plot <- ggplot(pca_data, aes(x = PCx, y = PCy)) + 
      geom_point(alpha = 1, show.legend = TRUE) +
      labs(
        x = paste0("PC", x, " (", round(var_explained[x] * 100, 2), "%)"),
        y = paste0("PC", y, " (", round(var_explained[y] * 100, 2), "%)"),
        title = paste("PCA of Methylation Data: PC", x, " vs PC", y)
      ) +
      theme_gray() + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  print(pca_plot)
  return(pca_plot)
  summary(pca_result)
  
}