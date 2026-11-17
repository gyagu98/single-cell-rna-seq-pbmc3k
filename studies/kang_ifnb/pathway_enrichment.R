# Pre-ranked GSEA and Hallmark pathway enrichment utilities

# Calculate running enrichment score (Kolmogorov-Smirnov style)
calc_gsea_score <- function(ranked_genes, gene_set, p_weight = 1.0) {
  if (length(ranked_genes) == 0 || length(gene_set) == 0) {
    stop("Ranked genes and gene set must be non-empty")
  }
  
  in_set <- names(ranked_genes) %in% gene_set
  n_hit <- sum(in_set)
  n_miss <- length(ranked_genes) - n_hit
  
  if (n_hit == 0) {
    return(list(es = 0, running_es = rep(0, length(ranked_genes)), lead_genes = character(0)))
  }
  
  weights <- abs(ranked_genes) ^ p_weight
  norm_hit <- sum(weights[in_set])
  if (norm_hit == 0) norm_hit <- 1
  
  p_hit <- ifelse(in_set, weights / norm_hit, 0)
  p_miss <- ifelse(!in_set, 1 / n_miss, 0)
  
  running_es <- cumsum(p_hit - p_miss)
  max_es <- max(running_es)
  min_es <- min(running_es)
  
  es <- if (abs(max_es) >= abs(min_es)) max_es else min_es
  
  # Identify leading edge genes
  peak_idx <- which.max(running_es)
  lead_genes <- names(ranked_genes)[which(in_set & seq_along(ranked_genes) <= peak_idx)]
  
  list(es = es, running_es = running_es, lead_genes = lead_genes, tag_fraction = n_hit / length(gene_set))
}

# Run GSEA across multiple gene sets
run_hallmark_gsea <- function(ranked_genes, gene_sets, min_size = 5, max_size = 500) {
  if (is.null(names(ranked_genes))) stop("ranked_genes must be a named numeric vector")
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  results <- list()
  for (set_name in names(gene_sets)) {
    gset <- gene_sets[[set_name]]
    overlap <- intersect(names(ranked_genes), gset)
    if (length(overlap) < min_size || length(overlap) > max_size) next
    
    score <- calc_gsea_score(ranked_genes, overlap)
    results[[set_name]] <- data.frame(
      pathway = set_name,
      set_size = length(gset),
      overlap = length(overlap),
      enrichment_score = score$es,
      tag_fraction = score$tag_fraction,
      n_leading_edge = length(score$lead_genes),
      leading_edge = paste(head(score$lead_genes, 10), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(results) == 0) {
    return(data.frame(pathway = character(0), enrichment_score = numeric(0), stringsAsFactors = FALSE))
  }
  
  res_df <- do.call(rbind, results)
  res_df <- res_df[order(-res_df$enrichment_score), ]
  rownames(res_df) <- NULL
  res_df
}
