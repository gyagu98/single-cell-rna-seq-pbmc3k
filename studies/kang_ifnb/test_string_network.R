source("studies/kang_ifnb/string_network.R")

# Test 1: Parse mock interactions and verify degree and edge deduplication
mock_interactions <- data.frame(
  preferredName_A = c("ISG15", "USP18", "IFIT1", "IFIT2", "ISG15", "LOW_CONF"),
  preferredName_B = c("USP18", "ISG15", "IFIT2", "IFIT1", "IFIT1", "GENE_X"),
  score = c(0.999, 0.999, 0.998, 0.998, 0.850, 0.400),
  stringsAsFactors = FALSE
)

net <- parse_string_interactions(mock_interactions, min_score = 0.700)

# LOW_CONF (<0.7) should be excluded
stopifnot(!("LOW_CONF" %in% net$nodes$gene))
# Duplicate undirected edges should be deduplicated
stopifnot(nrow(net$edges) == 3) # (ISG15-USP18, IFIT1-IFIT2, IFIT1-ISG15)

# Verify degrees
stopifnot(subset(net$nodes, gene == "ISG15")$degree == 2)
stopifnot(subset(net$nodes, gene == "IFIT1")$degree == 2)
stopifnot(subset(net$nodes, gene == "USP18")$degree == 1)

# Test 2: Effector functional clustering annotation
annotated <- annotate_interferon_effectors(net$nodes)
stopifnot("functional_cluster" %in% names(annotated),
          subset(annotated, gene == "ISG15")$functional_cluster == "ISGylation cascade",
          subset(annotated, gene == "IFIT1")$functional_cluster == "Viral RNA recognition (IFIT)")

# Test 3: Empty inputs
empty_net <- parse_string_interactions(data.frame(), min_score = 0.700)
stopifnot(nrow(empty_net$edges) == 0, nrow(empty_net$nodes) == 0)

cat("STRING network parsing and effector annotation tests passed.\n")
