source("studies/kang_ifnb/pathway_enrichment.R")

# Test 1: Validate running score with enriched gene set
set.seed(42)
all_genes <- paste0("gene_", 1:100)
stats <- seq(5, -5, length.out = 100)
names(stats) <- all_genes

# Strongly enriched set at top
target_set <- paste0("gene_", 1:15)
res <- calc_gsea_score(stats, target_set)
stopifnot(res$es > 0.8, length(res$lead_genes) == 15)

# Negative enriched set at bottom
bot_set <- paste0("gene_", 85:100)
res_bot <- calc_gsea_score(stats, bot_set)
stopifnot(res_bot$es < -0.8)

# Test 2: Multi-pathway runner
gene_sets <- list(
  HALLMARK_INTERFERON_ALPHA = paste0("gene_", 1:12),
  HALLMARK_INTERFERON_GAMMA = paste0("gene_", 5:20),
  UNRELATED_PATHWAY = paste0("gene_", 45:55)
)

gsea_res <- run_hallmark_gsea(stats, gene_sets, min_size = 5)
stopifnot(nrow(gsea_res) == 3,
          gsea_res$pathway[1] == "HALLMARK_INTERFERON_ALPHA",
          gsea_res$enrichment_score[1] > 0.7)

# Test 3: Edge cases (empty sets, size filters)
fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
fails(calc_gsea_score(numeric(0), target_set))
fails(calc_gsea_score(stats, character(0)))

empty_res <- run_hallmark_gsea(stats, list(TINY = c("gene_1")), min_size = 5)
stopifnot(nrow(empty_res) == 0)

cat("Pathway enrichment and GSEA score calculation tests passed.\n")
