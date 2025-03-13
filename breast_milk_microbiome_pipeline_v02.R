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
#' @version 1.3.0
#' @date February 2025

# =====================================================================
# INSTALL AND LOAD PACKAGES
# =====================================================================

required_packages <- c("BiocManager", "tidyverse", "vegan", "reshape2", "ggplot2", 
                       "plotly", "parallel", "doParallel", "foreach", "clue")

installed_packages <- local({
  ip <- installed.packages()
  ip[, "Package"]
})

needed <- setdiff(required_packages, installed_packages)
if (length(needed) > 0) {
  install.packages(needed)
}
lapply(required_packages, library, character.only = TRUE)

bioc_packages <- c("phyloseq", "DirichletMultinomial", "DESeq2")
needed_bioc <- setdiff(bioc_packages, installed_packages)
if (length(needed_bioc) > 0) {
  BiocManager::install(needed_bioc)
}
lapply(bioc_packages, library, character.only = TRUE)

# =====================================================================
# DATA LOADING AND PROCESSING
# =====================================================================

create_phyloseq <- function(otu_path, tax_path, meta_path) {
  otu_data <- read.csv(otu_path, row.names = 1, check.names = FALSE)
  tax_data <- read.csv(tax_path, row.names = 1, check.names = FALSE)
  meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
  
  otu_matrix <- as.matrix(otu_data)
  tax_matrix <- as.matrix(tax_data[rownames(otu_data), ])
  
  OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  TAX <- tax_table(tax_matrix)
  META <- sample_data(meta_data[colnames(otu_data), ])
  
  physeq <- phyloseq(OTU, TAX, META)
  return(physeq)
}

# =====================================================================
# IMPROVED PARALLEL PROCESSING
# =====================================================================

initialize_cluster <- function(n_cores) {
  if (!is.numeric(n_cores) || n_cores < 1) {
    stop("n_cores must be a positive integer")
  }
  
  if (n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(n_cores, type = "PSOCK")
    } else {
      cl <- makeCluster(n_cores)
    }
    registerDoParallel(cl)
    return(cl)
  }
  return(NULL)
}

stop_cluster <- function(cl) {
  if (!is.null(cl)) stopCluster(cl)
}

# =====================================================================
# DMM MODELING
# =====================================================================

fit_dmm_models <- function(physeq, k_range = 1:7, n_cores = 1) {
  counts <- as(otu_table(physeq), "matrix")
  cl <- initialize_cluster(n_cores)
  
  models <- lapply(k_range, function(k) {
    set.seed(123 + k)
    DirichletMultinomial::dmn(counts, k, verbose = FALSE)
  })
  
  stop_cluster(cl)
  names(models) <- paste0("k", k_range)
  return(models)
}

# =====================================================================
# COMMUNITY ASSIGNMENT WITH THRESHOLD
# =====================================================================

assign_communities <- function(dmm_fit, physeq, criterion = "laplace", threshold = 0.6) {
  best_k <- dmm_fit[[criterion]]
  optimal_model <- dmm_fit[[paste0("k", best_k)]]
  
  mixtures <- mixture(optimal_model)
  assignments <- apply(mixtures, 1, which.max)
  probabilities <- apply(mixtures, 1, max)
  assignments[probabilities < threshold] <- NA
  
  assignment_df <- data.frame(
    Sample = rownames(mixtures),
    Community = factor(assignments),
    Probability = probabilities
  )
  
  sample_data(physeq)$Community <- factor(
    assignment_df$Community[match(sample_names(physeq), assignment_df$Sample)]
  )
  
  return(physeq)
}

# =====================================================================
# VISUALIZATION
# =====================================================================

plot_community_heatmap <- function(physeq, tax_level = "Genus", top_n = 15) {
  physeq_glom <- tax_glom(physeq, taxrank = tax_level)
  physeq_rel <- transform_sample_counts(physeq_glom, function(x) scale(x))
  melted <- psmelt(physeq_rel)
  
  top_taxa <- melted %>%
    group_by(.data[[tax_level]]) %>%
    summarise(MaxMean = max(Abundance)) %>%
    arrange(desc(MaxMean)) %>%
    head(top_n) %>%
    pull(.data[[tax_level]])
  
  plot_data <- melted[melted[[tax_level]] %in% top_taxa, ]
  
  p <- ggplot(plot_data, aes(x = Community, y = .data[[tax_level]], fill = Abundance)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkblue", name = "Z-score Scaled Abundance") +
    theme_minimal() +
    labs(title = "Community Composition", x = "Community", y = tax_level)
  
  return(p)
}
