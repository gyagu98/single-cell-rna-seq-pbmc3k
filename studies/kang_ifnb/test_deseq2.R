source("studies/kang_ifnb/pseudobulk.R")
set.seed(42)
meta <- expand.grid(cell=1:25, donor=paste0("d", 1:4),
                    condition=c("control", "stim"), stringsAsFactors=FALSE)
meta$cell_type <- "monocyte"
rownames(meta) <- paste0("c", seq_len(nrow(meta)))
counts <- vapply(seq_len(nrow(meta)), function(i) {
  mu <- rep(5, 500)
  if (meta$condition[i] == "stim") mu[1:30] <- 80
  rnbinom(500, mu=mu, size=5)
}, numeric(500))
dimnames(counts) <- list(paste0("gene",1:500), rownames(meta))
pb <- paired_pseudobulk(counts, meta, "monocyte", "control", "stim")
fit <- fit_paired_deseq2(pb, "control", "stim")
stopifnot(median(fit$results$log2FoldChange[fit$results$gene %in% paste0("gene",1:30)])>2,
          any(fit$results$padj<0.05, na.rm=TRUE))
cat("Paired DESeq2 synthetic-effect integration passed.\n")
