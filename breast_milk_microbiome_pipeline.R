#!/usr/bin/env Rscript

#' Breast Milk Microbiome Community Analysis Pipeline
#' 
#' A comprehensive pipeline for analyzing breast milk microbiome data using
#' Dirichlet Multinomial Mixtures (DMM) to identify community types.
#'
#' This pipeline includes:
#' 1. Data preprocessing and validation
#' 2. DMM modeling
#' 3. Model validation
#' 4. Statistical analysis
#' 5. Visualization
#' 6. Community prediction for new samples
#'
#' @author Claude
#' @version 1.1.0
#' @date February 2025

# =====================================================================
# INSTALL AND LOAD PACKAGES
# =====================================================================

# List of required packages with minimum versions
required_packages <- c(
  "BiocManager", "tidyverse", "vegan", "reshape2", "ggplot2", 
  "plotly", "parallel", "doParallel", "foreach"
)

# Install basic packages if needed
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Install Bioconductor packages if needed
bioc_packages <- c("phyloseq", "DirichletMultinomial", "DESeq2")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# =====================================================================
# DATA LOADING AND PROCESSING
# =====================================================================

#' Create phyloseq object from OTU table, taxonomy, and metadata
#' 
#' @param otu_path Path to OTU count table (CSV with taxa as rows, samples as columns)
#' @param tax_path Path to taxonomy table (CSV with taxa as rows)
#' @param meta_path Path to sample metadata (CSV with samples as rows)
#' @return A phyloseq object
create_phyloseq <- function(otu_path, tax_path, meta_path) {
  # Check that files exist
  if (!file.exists(otu_path)) stop("OTU file does not exist:", otu_path)
  if (!file.exists(tax_path)) stop("Taxonomy file does not exist:", tax_path)
  if (!file.exists(meta_path)) stop("Metadata file does not exist:", meta_path)
  
  # Read files
  otu_data <- read.csv(otu_path, row.names = 1, check.names = FALSE)
  tax_data <- read.csv(tax_path, row.names = 1, check.names = FALSE)
  meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
  
  # Check data consistency
  missing_taxa <- setdiff(rownames(otu_data), rownames(tax_data))
  if (length(missing_taxa) > 0) {
    warning(paste("Removing", length(missing_taxa), "taxa missing from taxonomy table"))
    otu_data <- otu_data[rownames(otu_data) %in% rownames(tax_data), ]
  }
  
  missing_samples <- setdiff(colnames(otu_data), rownames(meta_data))
  if (length(missing_samples) > 0) {
    warning(paste("Removing", length(missing_samples), "samples missing from metadata"))
    otu_data <- otu_data[, colnames(otu_data) %in% rownames(meta_data)]
  }
  
  # Convert to matrices
  otu_matrix <- as.matrix(otu_data)
  tax_matrix <- as.matrix(tax_data[rownames(otu_data), ])
  
  # Create phyloseq components
  OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  TAX <- tax_table(tax_matrix)
  META <- sample_data(meta_data[colnames(otu_data), ])
  
  # Create phyloseq object
  physeq <- phyloseq(OTU, TAX, META)
  
  message(sprintf("Created phyloseq object with %d taxa and %d samples", 
                  ntaxa(physeq), nsamples(physeq)))
  
  return(physeq)
}

#' Preprocess microbiome data by filtering and normalizing
#' 
#' @param physeq Phyloseq object
#' @param min_prevalence Minimum prevalence (fraction of samples) for a taxon to be kept
#' @param min_abundance Minimum relative abundance for a taxon to be kept
#' @param rarefy_depth Sequencing depth for rarefaction (NULL for no rarefaction)
#' @param normalize Method for normalization: "relative", "log", or NULL
#' @return Processed phyloseq object
preprocess_data <- function(physeq, min_prevalence = 0.1, min_abundance = 0.001,
                            rarefy_depth = NULL, normalize = "relative") {
  # Track initial dimensions
  initial_taxa <- ntaxa(physeq)
  initial_samples <- nsamples(physeq)
  
  # Filter by prevalence
  if (min_prevalence > 0) {
    prevalence <- apply(otu_table(physeq) > 0, 1, mean)
    physeq <- prune_taxa(prevalence >= min_prevalence, physeq)
    message(sprintf("Removed %d taxa with prevalence < %s", 
                    initial_taxa - ntaxa(physeq), min_prevalence))
  }
  
  # Filter by abundance
  if (min_abundance > 0) {
    total_counts <- sum(otu_table(physeq))
    rel_abundance <- rowSums(otu_table(physeq)) / total_counts
    physeq <- prune_taxa(rel_abundance >= min_abundance, physeq)
    message(sprintf("Removed %d taxa with relative abundance < %s",
                    initial_taxa - ntaxa(physeq), min_abundance))
  }
  
  # Rarefy if requested
  if (!is.null(rarefy_depth)) {
    sample_counts <- sample_sums(physeq)
    if (any(sample_counts < rarefy_depth)) {
      warning(sprintf("%d samples have fewer than %d counts and will be removed",
                      sum(sample_counts < rarefy_depth), rarefy_depth))
    }
    physeq <- rarefy_even_depth(physeq, sample.size = rarefy_depth, rngseed = 123)
    message(sprintf("Rarefied to %d counts per sample", rarefy_depth))
  }
  
  # Normalize if requested
  if (!is.null(normalize)) {
    if (normalize == "relative") {
      physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
      message("Converted to relative abundances")
    } else if (normalize == "log") {
      physeq <- transform_sample_counts(physeq, function(x) log1p(x))
      message("Applied log(x+1) transformation")
    } else {
      warning("Unknown normalization method:", normalize)
    }
  }
  
  message(sprintf("Final dataset: %d taxa across %d samples", 
                  ntaxa(physeq), nsamples(physeq)))
  
  return(physeq)
}

# =====================================================================
# DMM MODELING
# =====================================================================

#' Fit Dirichlet Multinomial Mixture models for different numbers of communities
#' 
#' @param physeq Phyloseq object
#' @param k_range Range of community numbers to test
#' @param n_cores Number of CPU cores to use (parallel)
#' @param seed Random seed for reproducibility
#' @return List with models and fit statistics
fit_dmm_models <- function(physeq, k_range = 1:7, n_cores = 1, seed = 123) {
  # Extract count table
  counts <- as(otu_table(physeq), "matrix")
  
  # Set up parallel processing if requested
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    message(sprintf("Using %d cores for parallel processing", n_cores))
    
    # Fit models in parallel
    models <- foreach(k = k_range, .packages = "DirichletMultinomial") %dopar% {
      set.seed(seed + k)
      message(sprintf("Fitting model with k=%d communities", k))
      DirichletMultinomial::dmn(counts, k, verbose = FALSE)
    }
    
    stopCluster(cl)
  } else {
    # Sequential processing
    models <- lapply(k_range, function(k) {
      set.seed(seed + k)
      message(sprintf("Fitting model with k=%d communities", k))
      DirichletMultinomial::dmn(counts, k, verbose = TRUE)
    })
  }
  
  names(models) <- paste0("k", k_range)
  
  # Calculate fit statistics
  fit_stats <- data.frame(
    k = k_range,
    laplace = sapply(models, DirichletMultinomial::laplace),
    AIC = sapply(models, DirichletMultinomial::AIC),
    BIC = sapply(models, DirichletMultinomial::BIC)
  )
  
  # Find optimal k by different criteria
  best_k <- list(
    laplace = fit_stats$k[which.min(fit_stats$laplace)],
    AIC = fit_stats$k[which.min(fit_stats$AIC)],
    BIC = fit_stats$k[which.min(fit_stats$BIC)]
  )
  
  message("Model fitting complete")
  message(sprintf("Best k by Laplace: %d, by AIC: %d, by BIC: %d",
                  best_k$laplace, best_k$AIC, best_k$BIC))
  
  return(list(
    models = models,
    fit_stats = fit_stats,
    best_k = best_k
  ))
}

#' Assign samples to communities based on optimal DMM model
#' 
#' @param dmm_fit Results from fit_dmm_models
#' @param physeq Original phyloseq object
#' @param criterion Criterion for model selection: "laplace", "AIC", or "BIC"
#' @return Phyloseq object with community assignments and optimal model
assign_communities <- function(dmm_fit, physeq, criterion = "laplace") {
  # Validate criterion
  if (!criterion %in% c("laplace", "AIC", "BIC")) {
    stop("criterion must be one of 'laplace', 'AIC', or 'BIC'")
  }
  
  # Select k based on criterion
  best_k <- dmm_fit$best_k[[criterion]]
  k_idx <- which(dmm_fit$fit_stats$k == best_k)
  optimal_model <- dmm_fit$models[[k_idx]]
  
  message(sprintf("Using optimal model with k=%d (criterion: %s)", 
                  best_k, criterion))
  
  # Get mixture proportions and community assignments
  mixtures <- mixture(optimal_model)
  assignments <- apply(mixtures, 1, which.max)
  probabilities <- apply(mixtures, 1, max)
  
  # Create assignment data frame
  assignment_df <- data.frame(
    Sample = names(assignments),
    Community = factor(assignments),
    Probability = probabilities
  )
  
  # Add assignments to phyloseq sample data
  sample_data(physeq)$Community <- factor(
    assignment_df$Community[match(sample_names(physeq), assignment_df$Sample)]
  )
  
  sample_data(physeq)$CommunityProb <- 
    assignment_df$Probability[match(sample_names(physeq), assignment_df$Sample)]
  
  return(list(
    physeq = physeq,
    model = optimal_model,
    assignments = assignment_df
  ))
}

#' Extract community characteristic taxa
#' 
#' @param model Optimal DMM model
#' @param physeq Phyloseq object with taxonomy
#' @param tax_level Taxonomic level to use (e.g., "Genus")
#' @param top_n Number of top taxa to report per community
#' @return Data frame with characteristic taxa for each community
get_characteristic_taxa <- function(model, physeq, tax_level = "Genus", top_n = 20) {
  # Extract fitted parameters
  fit_params <- fitted(model)
  
  # Get contributions of each taxon to each community
  contributions <- fit_params$pi
  
  # Create results data frame
  results <- data.frame()
  
  # For each community
  for (i in 1:ncol(contributions)) {
    # Get top taxa
    top_idx <- order(contributions[, i], decreasing = TRUE)[1:min(top_n, nrow(contributions))]
    
    # Get taxonomy data
    taxa_info <- tax_table(physeq)[rownames(contributions)[top_idx], ]
    
    # Combine into data frame
    community_df <- data.frame(
      Community = i,
      Taxon = rownames(contributions)[top_idx],
      Contribution = contributions[top_idx, i],
      taxa_info[, tax_level],
      stringsAsFactors = FALSE
    )
    
    colnames(community_df)[4] <- tax_level
    
    # Add to results
    results <- rbind(results, community_df)
  }
  
  return(results)
}

# =====================================================================
# MODEL VALIDATION
# =====================================================================

#' Perform cross-validation for DMM model
#' 
#' @param physeq Phyloseq object
#' @param k Number of communities
#' @param folds Number of CV folds
#' @param seed Random seed
#' @return Cross-validation results
cross_validate_model <- function(physeq, k, folds = 5, seed = 123) {
  # Extract counts
  counts <- as(otu_table(physeq), "matrix")
  
  # Set seed
  set.seed(seed)
  
  # Create folds
  sample_idx <- sample(rep(1:folds, length.out = ncol(counts)))
  
  # Results containers
  fold_results <- data.frame()
  
  # For each fold
  for (fold in 1:folds) {
    message(sprintf("Processing fold %d/%d", fold, folds))
    
    # Split data
    train_idx <- sample_idx != fold
    test_idx <- sample_idx == fold
    
    train_counts <- counts[, train_idx]
    test_counts <- counts[, test_idx]
    
    # Train model
    fold_model <- DirichletMultinomial::dmn(train_counts, k, verbose = FALSE)
    
    # Test on held-out data
    pred_mix <- predict(fold_model, test_counts)
    assignments <- apply(pred_mix, 1, which.max)
    probabilities <- apply(pred_mix, 1, max)
    
    # Store results
    fold_df <- data.frame(
      Fold = fold,
      Sample = colnames(test_counts),
      Community = assignments,
      Probability = probabilities
    )
    
    fold_results <- rbind(fold_results, fold_df)
  }
  
  # Create summary
  fold_summary <- aggregate(Probability ~ Fold, data = fold_results, FUN = mean)
  
  message("Cross-validation complete")
  message(sprintf("Mean assignment probability: %.2f", mean(fold_results$Probability)))
  
  return(list(
    sample_results = fold_results,
    fold_summary = fold_summary
  ))
}

#' Assess stability of community assignments through bootstrapping
#' 
#' @param physeq Phyloseq object
#' @param k Number of communities
#' @param n_bootstrap Number of bootstrap samples
#' @param seed Random seed
#' @return Stability assessment
assess_stability <- function(physeq, k, n_bootstrap = 10, seed = 123) {
  # Extract counts
  counts <- as(otu_table(physeq), "matrix")
  
  # Set seed
  set.seed(seed)
  
  # Train full model
  full_model <- DirichletMultinomial::dmn(counts, k, verbose = FALSE)
  full_assignments <- apply(mixture(full_model), 1, which.max)
  
  # Results container
  stability_df <- data.frame()
  
  # For each bootstrap iteration
  for (b in 1:n_bootstrap) {
    message(sprintf("Bootstrap iteration %d/%d", b, n_bootstrap))
    
    # Create bootstrap sample
    boot_idx <- sample(1:ncol(counts), replace = TRUE)
    boot_counts <- counts[, boot_idx]
    
    # Train model
    boot_model <- DirichletMultinomial::dmn(boot_counts, k, verbose = FALSE)
    
    # Get assignments
    boot_mix <- mixture(boot_model)
    boot_assignments <- apply(boot_mix, 1, which.max)
    
    # Calculate agreement with full model
    # This is simplified - in practice need to handle label switching
    overlap_table <- table(
      full_assignments,
      boot_assignments[match(names(full_assignments), names(boot_assignments))]
    )
    agreement <- sum(diag(overlap_table)) / sum(overlap_table)
    
    # Store results
    stability_df <- rbind(stability_df, data.frame(
      Bootstrap = b,
      Agreement = agreement
    ))
  }
  
  message("Stability assessment complete")
  message(sprintf("Mean bootstrap agreement: %.2f", mean(stability_df$Agreement)))
  
  return(stability_df)
}

# =====================================================================
# VISUALIZATION
# =====================================================================

#' Create a heatmap of community composition
#' 
#' @param physeq Phyloseq object with community assignments
#' @param tax_level Taxonomic level to display
#' @param top_n Number of top taxa to display
#' @return ggplot object
plot_community_heatmap <- function(physeq, tax_level = "Genus", top_n = 15) {
  # Check for community assignments
  if (!"Community" %in% sample_variables(physeq)) {
    stop("Community assignments not found in phyloseq object")
  }
  
  # Aggregate at taxonomic level
  physeq_glom <- tax_glom(physeq, taxrank = tax_level)
  
  # Convert to relative abundance
  physeq_rel <- transform_sample_counts(physeq_glom, function(x) x / sum(x))
  
  # Melt to data frame
  melted <- psmelt(physeq_rel)
  
  # Calculate mean abundances by community
  community_means <- aggregate(Abundance ~ Community + get(tax_level), 
                               data = melted, FUN = mean)
  colnames(community_means)[2] <- tax_level
  
  # Select top taxa
  top_taxa <- community_means %>%
    group_by(!!as.name(tax_level)) %>%
    summarise(MaxMean = max(Abundance)) %>%
    arrange(desc(MaxMean)) %>%
    head(top_n) %>%
    pull(!!as.name(tax_level))
  
  # Filter data to top taxa
  plot_data <- community_means[community_means[[tax_level]] %in% top_taxa, ]
  
  # Create heatmap
  p <- ggplot(plot_data, aes(x = Community, y = .data[[tax_level]], fill = Abundance)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkblue", 
                        name = "Mean\nRelative\nAbundance") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          panel.grid = element_blank()) +
    labs(title = "Community Composition",
         x = "Community",
         y = tax_level)
  
  return(p)
}

#' Plot model fit statistics
#' 
#' @param dmm_fit Results from fit_dmm_models
#' @return ggplot object
plot_fit_stats <- function(dmm_fit) {
  # Reshape fit statistics to long format
  fit_long <- reshape2::melt(dmm_fit$fit_stats, 
                             id.vars = "k",
                             variable.name = "Criterion",
                             value.name = "Value")
  
  # Create plot
  p <- ggplot(fit_long, aes(x = k, y = Value, color = Criterion)) +
    geom_point(size = 3) +
    geom_line() +
    facet_wrap(~Criterion, scales = "free_y") +
    theme_minimal() +
    labs(title = "DMM Model Fit Statistics",
         x = "Number of Communities (k)",
         y = "Value (lower is better)")
  
  return(p)
}

#' Create an ordination plot colored by community
#' 
#' @param physeq Phyloseq object with community assignments
#' @param method Ordination method: "PCoA", "NMDS", or "PCA"
#' @param distance Distance measure for PCoA/NMDS: "bray", "jaccard", etc.
#' @return ggplot object
plot_ordination_communities <- function(physeq, method = "PCoA", distance = "bray") {
  # Check for community assignments
  if (!"Community" %in% sample_variables(physeq)) {
    stop("Community assignments not found in phyloseq object")
  }
  
  # Create ordination
  if (method == "PCA") {
    ord <- ordinate(physeq, method = "RDA", distance = "euclidean")
  } else {
    ord <- ordinate(physeq, method = method, distance = distance)
  }
  
  # Plot
  p <- plot_ordination(physeq, ord, color = "Community") +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    labs(title = sprintf("%s of Breast Milk Microbiomes", method))
  
  return(p)
}

# =====================================================================
# PREDICTION FOR NEW SAMPLES
# =====================================================================

#' Predict community assignments for new samples
#' 
#' @param model DMM model
#' @param new_counts Matrix of counts for new samples (taxa as rows)
#' @return Data frame with predicted communities
predict_communities <- function(model, new_counts) {
  # Check input
  if (!inherits(model, "DMN")) {
    stop("model must be a DMN object")
  }
  
  if (!is.matrix(new_counts)) {
    stop("new_counts must be a matrix with taxa as rows")
  }
  
  # Align taxa with model
  model_taxa <- rownames(model@count)
  missing_taxa <- setdiff(model_taxa, rownames(new_counts))
  
  if (length(missing_taxa) > 0) {
    warning(sprintf("%d taxa in model are missing from new data - adding zeros", 
                    length(missing_taxa)))
    
    # Add missing taxa as zeros
    zero_counts <- matrix(0, nrow = length(missing_taxa), ncol = ncol(new_counts))
    rownames(zero_counts) <- missing_taxa
    new_counts <- rbind(new_counts, zero_counts)
  }
  
  # Ensure same order as model
  new_counts <- new_counts[model_taxa, ]
  
  # Predict
  mixtures <- predict(model, new_counts)
  assignments <- apply(mixtures, 1, which.max)
  probabilities <- apply(mixtures, 1, max)
  
  # Create results
  results <- data.frame(
    Sample = colnames(new_counts),
    Community = assignments,
    Probability = probabilities
  )
  
  return(results)
}

#' Save model and results to file
#' 
#' @param model DMM model to save
#' @param community_chars Community characteristics
#' @param dir Output directory
#' @param prefix File name prefix
save_model <- function(model, community_chars, dir = ".", prefix = "breast_milk_dmm") {
  # Create directory if needed
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  # Save model
  model_file <- file.path(dir, paste0(prefix, "_model.rds"))
  saveRDS(model, model_file)
  
  # Save community characteristics
  chars_file <- file.path(dir, paste0(prefix, "_community_chars.csv"))
  write.csv(community_chars, chars_file, row.names = FALSE)
  
  message("Saved model to: ", model_file)
  message("Saved community characteristics to: ", chars_file)
}

# =====================================================================
# COMPLETE WORKFLOW FUNCTION
# =====================================================================

#' Run complete breast milk microbiome analysis workflow
#' 
#' @param otu_path Path to OTU table
#' @param tax_path Path to taxonomy table
#' @param meta_path Path to metadata
#' @param output_dir Output directory
#' @param k_range Range of communities to test
#' @param n_cores Number of CPU cores to use
#' @return List with analysis results
analyze_breast_milk <- function(otu_path, tax_path, meta_path, output_dir = "results",
                                k_range = 1:7, n_cores = 1) {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Start logging
  log_file <- file.path(output_dir, "analysis_log.txt")
  sink(log_file)
  cat("Breast Milk Microbiome Analysis\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Record session info
  cat("Session Info:\n")
  print(sessionInfo())
  cat("\n")
  
  # Load data
  cat("Loading data...\n")
  physeq <- create_phyloseq(otu_path, tax_path, meta_path)
  
  # Preprocess
  cat("\nPreprocessing data...\n")
  physeq_proc <- preprocess_data(physeq, min_prevalence = 0.1, min_abundance = 0.001)
  
  # Fit DMM models
  cat("\nFitting DMM models...\n")
  dmm_fit <- fit_dmm_models(physeq_proc, k_range = k_range, n_cores = n_cores)
  
  # Select optimal model and assign communities
  cat("\nAssigning communities...\n")
  community_results <- assign_communities(dmm_fit, physeq_proc, criterion = "laplace")
  
  # Extract community characteristics
  cat("\nExtracting community characteristics...\n")
  chars <- get_characteristic_taxa(community_results$model, physeq_proc)
  
  # Cross-validate
  cat("\nPerforming cross-validation...\n")
  cv_results <- cross_validate_model(physeq_proc, k = dmm_fit$best_k$laplace)
  
  # Save results
  cat("\nSaving results...\n")
  save_model(community_results$model, chars, dir = output_dir)
  
  # Create plots
  cat("\nCreating visualizations...\n")
  p1 <- plot_fit_stats(dmm_fit)
  ggsave(file.path(output_dir, "fit_statistics.pdf"), p1, width = 8, height = 6)
  
  p2 <- plot_community_heatmap(community_results$physeq)
  ggsave(file.path(output_dir, "community_heatmap.pdf"), p2, width = 10, height = 8)
  
  p3 <- plot_ordination_communities(community_results$physeq)
  ggsave(file.path(output_dir, "ordination.pdf"), p3, width = 8, height = 6)
  
  # Stop logging
  sink()
  
  cat("Analysis complete! Results saved to:", output_dir, "\n")
  
  # Return results
  return(list(
    physeq = physeq,
    physeq_processed = physeq_proc,
    dmm_fit = dmm_fit,
    community_results = community_results,
    community_chars = chars,
    cv_results = cv_results
  ))
}

# =====================================================================
# USAGE EXAMPLE
# =====================================================================

# Example usage
if (FALSE) {
  # Define paths to data files
  otu_file <- "data/breast_milk_otu_table.csv"
  tax_file <- "data/breast_milk_taxonomy.csv"
  meta_file <- "data/breast_milk_metadata.csv"
  
  # Run analysis
  results <- analyze_breast_milk(
    otu_path = otu_file,
    tax_path = tax_file,
    meta_path = meta_file,
    output_dir = "breast_milk_analysis",
    k_range = 1:10,
    n_cores = 4
  )
  
  # Predict communities for new samples
  new_otu_file <- "data/new_samples_otu.csv"
  new_counts <- as.matrix(read.csv(new_otu_file, row.names = 1))
  
  predictions <- predict_communities(
    model = results$community_results$model,
    new_counts = new_counts
  )
  
  # Save predictions
  write.csv(predictions, "breast_milk_analysis/new_sample_predictions.csv", row.names = FALSE)
}