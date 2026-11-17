# STRING protein interaction network utilities

# Parse STRING API JSON or tabular results into graph nodes and edges
parse_string_interactions <- function(interactions_df, min_score = 0.700) {
  if (nrow(interactions_df) == 0) {
    return(list(
      edges = data.frame(source = character(0), target = character(0), score = numeric(0), stringsAsFactors = FALSE),
      nodes = data.frame(gene = character(0), degree = integer(0), stringsAsFactors = FALSE)
    ))
  }
  
  # Ensure columns exist
  pA_col <- if ("preferredName_A" %in% names(interactions_df)) "preferredName_A" else "source"
  pB_col <- if ("preferredName_B" %in% names(interactions_df)) "preferredName_B" else "target"
  score_col <- if ("score" %in% names(interactions_df)) "score" else "combined_score"
  
  df <- interactions_df[interactions_df[[score_col]] >= min_score, ]
  if (nrow(df) == 0) {
    return(list(
      edges = data.frame(source = character(0), target = character(0), score = numeric(0), stringsAsFactors = FALSE),
      nodes = data.frame(gene = character(0), degree = integer(0), stringsAsFactors = FALSE)
    ))
  }
  
  # Canonical undirected edges (sort endpoints)
  src <- as.character(df[[pA_col]])
  tgt <- as.character(df[[pB_col]])
  sc <- as.numeric(df[[score_col]])
  
  ordered_edges <- t(apply(cbind(src, tgt), 1, sort))
  edge_df <- data.frame(
    source = ordered_edges[, 1],
    target = ordered_edges[, 2],
    score = sc,
    stringsAsFactors = FALSE
  )
  
  # Deduplicate edges
  edge_df <- unique(edge_df)
  edge_df <- edge_df[edge_df$source != edge_df$target, ] # Remove self-loops
  
  # Compute node degree
  all_nodes <- unique(c(edge_df$source, edge_df$target))
  degrees <- sapply(all_nodes, function(n) {
    sum(edge_df$source == n | edge_df$target == n)
  })
  
  node_df <- data.frame(
    gene = all_nodes,
    degree = as.integer(degrees),
    stringsAsFactors = FALSE
  )
  node_df <- node_df[order(-node_df$degree), ]
  rownames(node_df) <- NULL
  
  list(edges = edge_df, nodes = node_df)
}

# Assign functional clusters to interferon effectors
annotate_interferon_effectors <- function(node_df) {
  effector_clusters <- list(
    "ISGylation cascade" = c("ISG15", "USP18", "HERC5", "UBE2L6"),
    "Viral RNA recognition (IFIT)" = c("IFIT1", "IFIT2", "IFIT3", "IFIT5"),
    "2'-5'-OAS / RNase L pathway" = c("OAS1", "OAS2", "OAS3", "OASL"),
    "Dynamin-like GTPase antiviral trap" = c("MX1", "MX2"),
    "Viperin nucleotide synthesis" = c("RSAD2", "CMPK2"),
    "Pro-inflammatory chemokine" = c("CXCL10", "CXCL11", "CXCL9", "CCL8"),
    "Cytosolic RNA sensor" = c("DDX58", "IFIH1"),
    "Translation blockade" = c("EIF2AK2"),
    "JAK-STAT transcription factor" = c("STAT1", "STAT2", "IRF7", "IRF9"),
    "Viral tethering / release inhibitor" = c("BST2")
  )
  
  cluster_map <- character(nrow(node_df))
  for (i in seq_len(nrow(node_df))) {
    g <- node_df$gene[i]
    matched <- "Antiviral effector"
    for (cname in names(effector_clusters)) {
      if (g %in% effector_clusters[[cname]]) {
        matched <- cname
        break
      }
    }
    cluster_map[i] <- matched
  }
  
  node_df$functional_cluster <- cluster_map
  node_df
}
