# Aggregate one cell type using independent donor-condition pairs.
# Counts: genes x cells; metadata rows must match cell IDs exactly.
paired_pseudobulk <- function(counts, metadata, cell_type, reference, comparison,
                             min_cells=20L, min_pairs=3L) {
  required <- c("donor", "condition", "cell_type")
  if (!all(required %in% names(metadata)) ||
      is.null(colnames(counts)) || is.null(rownames(counts)) ||
      anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts)) ||
      !identical(colnames(counts), rownames(metadata)))
    stop("Counts require unique gene/cell IDs and exactly aligned metadata")
  if (anyNA(metadata[, required]) ||
      any(vapply(metadata[, required], function(x) any(!nzchar(as.character(x))), logical(1))))
    stop("Missing donor, condition or cell-type metadata")
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0) ||
      any(values != round(values))) stop("Raw nonnegative integer counts required")
  if (length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 1 ||
      min_cells != as.integer(min_cells) || length(min_pairs) != 1L ||
      !is.finite(min_pairs) || min_pairs < 3 || min_pairs != as.integer(min_pairs))
    stop("Require positive integer min_cells and at least three donor pairs")
  if (identical(reference, comparison) ||
      !setequal(unique(as.character(metadata$condition)), c(reference, comparison)))
    stop("Exactly two conditions matching the explicit contrast are required")
  select <- which(as.character(metadata$cell_type) == cell_type)
  if (!length(select)) stop("Requested cell type is absent")
  m <- metadata[select, , drop=FALSE]
  m$donor <- as.character(m$donor)
  m$condition <- factor(as.character(m$condition), levels=c(reference, comparison))
  audit <- as.data.frame(table(donor=m$donor, condition=m$condition),
                         stringsAsFactors=FALSE)
  names(audit)[3] <- "n_cells"
  complete <- tapply(audit$n_cells >= min_cells, audit$donor, all)
  donors <- names(complete)[complete]
  audit$retained <- as.character(audit$donor) %in% donors
  if (length(donors) < min_pairs)
    stop("Insufficient complete donor pairs after the cell-count threshold")
  sample_metadata <- audit[audit$retained, c("donor", "condition", "n_cells")]
  rownames(sample_metadata) <- paste0("pb", seq_len(nrow(sample_metadata)))
  summed <- vapply(seq_len(nrow(sample_metadata)), function(i) {
    idx <- select[m$donor == sample_metadata$donor[i] &
                  m$condition == sample_metadata$condition[i]]
    if (inherits(counts, "Matrix")) Matrix::rowSums(counts[, idx, drop=FALSE])
    else rowSums(counts[, idx, drop=FALSE])
  }, numeric(nrow(counts)))
  rownames(summed) <- rownames(counts)
  colnames(summed) <- rownames(sample_metadata)
  sample_metadata$donor <- factor(sample_metadata$donor)
  sample_metadata$condition <- factor(sample_metadata$condition,
                                      levels=c(reference, comparison))
  design <- model.matrix(~ donor + condition, data=sample_metadata)
  if (qr(design)$rank != ncol(design)) stop("Paired design is not full rank")
  list(counts=summed, metadata=sample_metadata, audit=audit)
}

fit_paired_deseq2 <- function(pb, reference, comparison) {
  if (!requireNamespace("DESeq2", quietly=TRUE)) stop("Install DESeq2")
  if (!identical(colnames(pb$counts), rownames(pb$metadata)))
    stop("Pseudobulk samples are misaligned")
  dds <- DESeq2::DESeqDataSetFromMatrix(pb$counts, pb$metadata, ~ donor + condition)
  keep <- rowSums(pb$counts >= 10) >= 3L
  if (sum(keep) < 2L) stop("Insufficient genes after count filtering")
  dds <- DESeq2::DESeq(dds[keep, ])
  res <- DESeq2::results(dds, contrast=c("condition", comparison, reference), alpha=0.05)
  list(dds=dds, results=data.frame(gene=rownames(res), as.data.frame(res), row.names=NULL))
}
