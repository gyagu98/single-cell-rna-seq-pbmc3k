source("studies/kang_ifnb/pseudobulk.R")
meta <- expand.grid(cell=1:3, donor=paste0("d",1:4),
                    condition=c("control", "stim"), stringsAsFactors=FALSE)
meta$cell_type <- "monocyte"
rownames(meta) <- paste0("c", seq_len(nrow(meta)))
x <- matrix(seq_len(2*nrow(meta)), nrow=2,
            dimnames=list(c("g1", "g2"), rownames(meta)))
pb <- paired_pseudobulk(x, meta, "monocyte", "control", "stim", min_cells=3)
stopifnot(ncol(pb$counts)==8, identical(colnames(pb$counts), rownames(pb$metadata)))
for (i in seq_len(ncol(pb$counts))) {
  idx <- meta$donor==pb$metadata$donor[i] & meta$condition==pb$metadata$condition[i]
  stopifnot(identical(unname(pb$counts[,i]), unname(rowSums(x[,idx]))))
}
fails <- function(expr) stopifnot(inherits(try(force(expr), silent=TRUE), "try-error"))
fails(paired_pseudobulk(x, meta[nrow(meta):1, ], "monocyte", "control", "stim", 3))
bad <- x; bad[1,1] <- -1
fails(paired_pseudobulk(bad, meta, "monocyte", "control", "stim", 3))
badmeta <- meta; badmeta$donor[1] <- NA
fails(paired_pseudobulk(x, badmeta, "monocyte", "control", "stim", 3))
# Remove one condition for d1: retain the remaining three complete pairs.
keep <- !(meta$donor=="d1" & meta$condition=="stim")
paired <- paired_pseudobulk(x[,keep], meta[keep,], "monocyte", "control", "stim", 3)
stopifnot(ncol(paired$counts)==6, !any(paired$metadata$donor=="d1"))
fails(paired_pseudobulk(x[,keep], meta[keep,], "monocyte", "control", "stim", 3, 4))
cat("Pseudobulk aggregation, pairing and invalid-input tests passed.\n")
