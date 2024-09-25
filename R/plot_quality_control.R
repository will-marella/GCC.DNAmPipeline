# This script contains the functions required to plot quality control results

###############################################################################
#### Plot log2 Intensity ####
###############################################################################

# Function to plot median intensities for each sample
plot_log2_intensity_results <- function(QC_result_values, passed_QC, any_set, sample_identifier_column, threshold=10.5, color_by=NULL, label_by=NULL){
  
  # Extract log2 intensity data from QC_result_values
  log2_meth_results <- data.frame(
    Sample = rownames(QC_result_values),
    log2_Meth = QC_result_values$log2_Meth_intensity,
    log2_Unmeth = QC_result_values$log2_Unmeth_intensity,
    stringsAsFactors = FALSE
  )
  
  # Add the pass/fail status from passed_QC
  log2_meth_results$Passed_Log2Intensity_QC <- passed_QC$log2_Meth_intensity
  
  pheno_data <- pData(any_set)
  
  # Merge QC_results with pData(any_set)
  plot_data <- merge(log2_meth_results, pheno_data, by.x = "Sample", by.y = sample_identifier_column, all.x = TRUE)
  
  if (!is.null(color_by) && !is.null(label_by)){
    p <- ggplot(plot_data, aes(x = log2_Meth, y = log2_Unmeth, 
                               color = .data[[color_by]], 
                               shape = .data[[label_by]],  # Add shape aesthetic
                               label = .data[[label_by]])) +
      geom_point(alpha = 0, show.legend = TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
      labs(x = "Median Methylated Median Intensity (log2)", 
           y = "Median Unmethylated Median Intensity (log2)", 
           title = paste("Median Intensities for Each Sample (Threshold =", threshold, ")")) +
      theme_gray() +
      coord_cartesian(xlim = c(min(plot_data$log2_Meth, threshold) - 0.1, 
                               max(plot_data$log2_Meth, threshold) + 0.1),
                      ylim = c(min(plot_data$log2_Unmeth, threshold) - 0.1, 
                               max(plot_data$log2_Unmeth, threshold) + 0.1)) +
      scale_color_manual(values = palette, name = color_by) +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (!is.null(color_by) && is.null(label_by)){
    p <- ggplot(plot_data, aes(x = log2_Meth, y = log2_Unmeth, 
                               color = .data[[color_by]])) +
      geom_point(alpha = 1, show.legend = TRUE) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
      labs(x = "Median Methylated Median Intensity (log2)", 
           y = "Median Unmethylated Median Intensity (log2)", 
           title = paste("Median Intensities for Each Sample (Threshold =", threshold, ")")) +
      theme_gray() +
      coord_cartesian(xlim = c(min(plot_data$log2_Meth, threshold) - 0.1, 
                               max(plot_data$log2_Meth, threshold) + 0.1),
                      ylim = c(min(plot_data$log2_Unmeth, threshold) - 0.1, 
                               max(plot_data$log2_Unmeth, threshold) + 0.1)) +
      scale_color_manual(values = palette, name = color_by) +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && !is.null(label_by)){
    p <- ggplot(plot_data, aes(x = log2_Meth, y = log2_Unmeth, 
                               color = .data[["Passed_Log2Intensity_QC"]], 
                               shape = .data[[label_by]],  # Add shape aesthetic
                               label = .data[[label_by]])) +
      geom_point(alpha = 0, show.legend = TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
      labs(x = "Median Methylated Median Intensity (log2)", 
           y = "Median Unmethylated Median Intensity (log2)", 
           title = paste("Median Intensities for Each Sample (Threshold =", threshold, ")")) +
      theme_gray() +
      coord_cartesian(xlim = c(min(plot_data$log2_Meth, threshold) - 0.1, 
                               max(plot_data$log2_Meth, threshold) + 0.1),
                      ylim = c(min(plot_data$log2_Unmeth, threshold) - 0.1, 
                               max(plot_data$log2_Unmeth, threshold) + 0.1)) +
      scale_color_manual(values = palette, name = "Passed") +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && is.null(label_by)){
    p <- ggplot(plot_data, aes(x = log2_Meth, y = log2_Unmeth, 
                               color = .data[["Passed_Log2Intensity_QC"]])) +
      geom_point(alpha = 1, show.legend = TRUE) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
      labs(x = "Median Methylated Median Intensity (log2)", 
           y = "Median Unmethylated Median Intensity (log2)", 
           title = paste("Median Intensities for Each Sample (Threshold =", threshold, ")")) +
      theme_gray() +
      coord_cartesian(xlim = c(min(plot_data$log2_Meth, threshold) - 0.1, 
                               max(plot_data$log2_Meth, threshold) + 0.1),
                      ylim = c(min(plot_data$log2_Unmeth, threshold) - 0.1, 
                               max(plot_data$log2_Unmeth, threshold) + 0.1)) +
      scale_color_manual(values = palette, name = "Passed") +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  
  print(p)
  return(p)
}

###############################################################################
#### Plot Detection P-Values ####
###############################################################################

plot_avg_detection_p_values <- function(QC_result_values, passed_QC, any_set, sample_identifier_column, threshold=0.005, color_by=NULL, label_by=NULL){
  
  # Extract detP data from QC_result_values
  detP_df <- data.frame(
    Sample = rownames(QC_result_values),
    AveragePValue = QC_result_values$AvgDetectionP,
    stringsAsFactors = FALSE
  )
  
  # Add the pass/fail status from passed_QC
  detP_df$Passed_DetectionP <- passed_QC$AvgDetectionP

  pheno_data <- pData(any_set)
  detP_df <- merge(detP_df, pheno_data, by.x = "Sample", by.y = sample_identifier_column, all.x = TRUE)
  
  if (!is.null(color_by) && !is.null(label_by)) {
    p <- ggplot(detP_df, aes(x = Sample, y = AveragePValue,
                             color = .data[[color_by]],
                             shape = .data[[label_by]],
                             label = .data[[label_by]])) +
      geom_point(alpha=0, size = 3, show.legend=TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Samples", y = "Average Detection P-value", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = color_by) +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (!is.null(color_by) && is.null(label_by)) {
    p <- ggplot(detP_df, aes(x = Sample, y = AveragePValue,
                             color = .data[[color_by]])) +
      geom_point(alpha = 1, size = 2.8, stroke = 0.5, color = "black", show.legend = FALSE) +  # Add black outline
      geom_point(alpha=1, size = 2, show.legend=TRUE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Samples", y = "Average Detection P-value", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = color_by) +
      scale_shape_discrete(name = "Sample") + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && !is.null(label_by)) {
    p <- ggplot(detP_df, aes(x = Sample, y = AveragePValue,
                             shape = .data[[label_by]],
                             label = .data[[label_by]],
                             color = .data[["Passed_DetectionP"]])) +
      geom_point(alpha=0, size = 3, show.legend=TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Samples", y = "Average Detection P-value", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = "Passed") +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && is.null(label_by)) {
    p <- ggplot(detP_df, aes(x = Sample, y = AveragePValue, color = .data[["Passed_DetectionP"]])) +
      geom_point(alpha = 1, size = 2.8, stroke = 0.5, color = "black", show.legend = FALSE) +  # Add black outline
      geom_point(alpha=1, size = 2, show.legend=TRUE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Samples", y = "Average Detection P-value", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = "Passed") +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      ) 
  }
  
  print(p)
  return(p)
  
}

###############################################################################
#### Plot BSCON Results ####
###############################################################################

plot_bscon <- function(QC_result_values, passed_QC, any_set, sample_identifier_column, color_by=NULL, label_by=NULL, threshold=80){
  
  # Extract detP data from QC_result_values
  bscon_df <- data.frame(
    Sample = rownames(QC_result_values),
    bscon_value = QC_result_values$BSCon,
    stringsAsFactors = FALSE
  )
  
  # Add the pass/fail status from passed_QC
  bscon_df$Passed_BSCon <- passed_QC$BSCon
  
  
  pheno_data <- pData(any_set)
  
  # Merge QC results with phenotype data
  bscon_plot_data <- merge(pheno_data, bscon_df, by.x = sample_identifier_column, by.y = "Sample", all.x = TRUE)
  
  if (!is.null(color_by) && !is.null(label_by)) {
    p <- ggplot(bscon_plot_data, aes(x = Sample, y = bscon_value,
                                     color = .data[[color_by]],
                                     shape = .data[[label_by]],
                                     label = .data[[label_by]])) +
      geom_point(alpha=0, size = 3, show.legend=TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Sample", y = "Bisulfite Conversion Score", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = color_by) +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (!is.null(color_by) && is.null(label_by)) {
    p <- ggplot(bscon_plot_data, aes(x = Sample, y = bscon_value,
                                     color = .data[[color_by]])) +
      geom_point(alpha = 1, size = 2.8, stroke = 0.5, color = "black", show.legend = FALSE) +  # Add black outline
      geom_point(alpha=1, size = 2, show.legend=TRUE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Sample", y = "Bisulfite Conversion Score", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = color_by) +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && !is.null(label_by)) {
    p <- ggplot(bscon_plot_data, aes(x = Sample, y = bscon_value,
                                     color = .data[["Passed_BSCon"]],
                                     shape = .data[[label_by]],
                                     label = .data[[label_by]])) +
      geom_point(alpha=0, size = 3, show.legend=TRUE) +
      geom_text(size = 3, show.legend = FALSE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Sample", y = "Bisulfite Conversion Score", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_shape_discrete(name = label_by) + 
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  if (is.null(color_by) && is.null(label_by)) {
    p <- ggplot(bscon_plot_data, aes(x = Sample, y = bscon_value,
                                     color = .data[["Passed_BSCon"]])) +
      geom_point(alpha = 1, size = 2.8, stroke = 0.5, color = "black", show.legend = FALSE) +  # Add black outline
      geom_point(alpha=1, size = 2, show.legend=TRUE) +
      theme_gray() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
      labs(x = "Sample", y = "Bisulfite Conversion Score", 
           title = "Average Detection P-values per Sample") +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
      scale_color_manual(values = palette, name = "Passed_BSCon") +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
  }
  
  print(p)
  return(p)
  
}


###############################################################################
#### Plot Sex Cross-Check ####
###############################################################################

plot_sex_cross_check <- function(any_set, sample_identifier_column, palette){
  
  beta_values <- getBeta(any_set)
  pheno_data <- pData(any_set)
  pheno_data <- as.data.frame(pheno_data)
  
  sex_probes <- c("cg03554089", "cg12653510", "cg05533223")
  sex_betas <- beta_values[rownames(beta_values) %in% sex_probes, , drop = FALSE]
  
  sex_data <- sex_betas %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(sample_identifier_column)
  
  sex_cross_check_data <- pheno_data %>%
    left_join(sex_data, by = sample_identifier_column)
  
  simple_boundary <- function(data) {
    gmm_fit <- Mclust(na.omit(data), G=2, verbose=FALSE)
    mean(gmm_fit$parameters$mean)
  }
  
  cg03554089_boundary <- simple_boundary(sex_cross_check_data$cg03554089)
  cg05533223_boundary <- simple_boundary(sex_cross_check_data$cg05533223)
  cg12653510_boundary <- simple_boundary(sex_cross_check_data$cg12653510)
  
  # Determine y-axis limits
  y_min <- min(sex_cross_check_data[, sex_probes], na.rm = TRUE)
  y_max <- max(sex_cross_check_data[, sex_probes], na.rm = TRUE)
  
  create_plot <- function(data, y_var, title, show.legend = FALSE, show.y.axis = FALSE, boundary) {
    p <- ggplot(data, aes(x = Sample, y = .data[[y_var]],
                          color = Sex)) +
      geom_point(alpha = 1, size = 2.8, stroke = 0.5, color = "black", show.legend = FALSE) +  # Add black outline
      geom_point(alpha=1, size = 2, show.legend=show.legend) +
      geom_hline(yintercept = boundary, color = "red", linetype = "dashed") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        axis.line.x = element_line(color = "black", linewidth = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.length.y = unit(0.2, "cm")
      ) +
      labs(title = title) +
      scale_color_manual(values = palette, name = "Sex") +
      scale_y_continuous(limits = c(y_min, y_max)) +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
        shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
      )
    
    if (!show.y.axis) {
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    } else {
      p <- p + theme(axis.line.y = element_line(color = "black", size = 0.5))
    }
    
    return(p)
  }
  
  p1 <- create_plot(sex_cross_check_data, "cg03554089", "cg03554089", show.legend = FALSE, show.y.axis = TRUE, cg03554089_boundary)
  p2 <- create_plot(sex_cross_check_data, "cg05533223", "cg05533223", show.legend = FALSE, show.y.axis = FALSE, cg05533223_boundary)
  p3 <- create_plot(sex_cross_check_data, "cg12653510", "cg12653510", show.legend = TRUE, show.y.axis = FALSE, cg12653510_boundary)
  
  # Combine plots using patchwork
  combined_plot <- p1 + p2 + p3 + 
    plot_layout(ncol = 3, widths = c(1.2, 1, 1.2)) +
    plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))
  
  # Add y-label only to the left plot
  combined_plot[[1]] <- combined_plot[[1]] + 
    theme(axis.title.y = element_text(size = 10, margin = margin(t = 10))) +
    ylab("Beta Value")
  
  # Add x-label only to the middle plot
  combined_plot[[2]] <- combined_plot[[2]] + 
    theme(axis.title.x = element_text(size = 10, margin = margin(t = 10))) +
    xlab(sample_identifier_column)
  
  print(combined_plot)
  return(combined_plot)
  
}

###############################################################################
#### Plot Global FDR Results ####
###############################################################################

# Plotting function to overlay outliers on PC1 and PC2
plot_fdr_results <- function(fdr_results){
  
  # Extract necessary data from the results
  pca_result <- fdr_results$pca_result_all
  glob_outliers_all <- fdr_results$glob_outliers_all
  all_samples <- fdr_results$all_samples
  
  # Prepare data for plotting
  pca_df <- data.frame(PC1 = pca_result$x[, 1],
                       PC2 = pca_result$x[, 2],
                       Sample = rownames(pca_result$x),
                       stringsAsFactors = FALSE)
  
  # Add a column for group (inliers or outliers by iteration)
  pca_df$Group <- "Inliers"
  
  # If there are any outliers, map them to their respective iterations
  if (nrow(glob_outliers_all) > 0) {
    # Find the earliest iteration each sample was flagged as an outlier
    outlier_info <- glob_outliers_all %>%
      group_by(Sample) %>%
      summarise(Iteration = min(Iteration)) %>%
      ungroup()
    
    # Merge with PCA data
    pca_df <- pca_df %>%
      left_join(outlier_info, by = "Sample") %>%
      mutate(Group = ifelse(is.na(Iteration), "Inliers", paste("Outliers, iteration", Iteration)))
  }
  
  # Plot PC1 vs PC2 with groups
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(alpha = 1, size = 2.8, stroke = 0.5, color = "black", show.legend = FALSE) +  # Add black outline
    geom_point(alpha=1, size = 2, show.legend=TRUE) +
    labs(title = "FDR Outliers by Iteration", x = "PC1", y = "PC2") +
    theme_gray() +
    scale_color_manual(values = palette) +
    guides(
      color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)),
      shape = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1))
    )
  
  print(p)
  return(p)
}