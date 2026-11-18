#!/usr/bin/env Rscript
# Rebuild the Kang et al. IFN-beta study outputs from the cached or downloaded
# muscData/ExperimentHub SingleCellExperiment.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(script_file), mustWork = TRUE)
usage <- function() {
  cat("Usage: Rscript studies/kang_ifnb/run_study.R [options]\n\n",
      "Options:\n",
      "  --output-dir PATH   Results directory (default: studies/kang_ifnb/results)\n",
      "  --cache-dir PATH    RDS cache directory (default: studies/kang_ifnb/data)\n",
      "  --refresh-data      Fetch the ExperimentHub/muscData object again\n",
      "  --help              Show this help\n\n",
      "The first run downloads the public Kang18_8vs8 object when no cache exists.\n",
      sep = "")
}

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) { usage(); quit(status = 0) }
parse_args <- function(x) {
  opts <- list(output_dir = file.path(script_dir, "results"),
               cache_dir = file.path(script_dir, "data"), refresh_data = FALSE)
  i <- 1L
  while (i <= length(x)) {
    if (x[i] == "--refresh-data") { opts$refresh_data <- TRUE; i <- i + 1L
    } else if (x[i] %in% c("--output-dir", "--cache-dir")) {
      if (i == length(x)) stop("Missing path after ", x[i], call. = FALSE)
      key <- if (x[i] == "--output-dir") "output_dir" else "cache_dir"
      opts[[key]] <- x[i + 1L]; i <- i + 2L
    } else stop("Unknown option: ", x[i], " (use --help)", call. = FALSE)
  }
  opts
}
opts <- parse_args(args)
out_dir <- normalizePath(opts$output_dir, mustWork = FALSE)
cache_dir <- normalizePath(opts$cache_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

required_cran <- c("jsonlite", "digest", "base64enc", "msigdbr")
required_bioc <- c("SummarizedExperiment", "SingleCellExperiment", "ExperimentHub",
                   "muscData", "DESeq2", "clusterProfiler")
missing <- c(required_cran, required_bioc)[!vapply(c(required_cran, required_bioc),
  requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
  ". See the one-time installation instructions in studies/kang_ifnb/README.md.", call. = FALSE)

source(file.path(script_dir, "pseudobulk.R"), local = TRUE)
source(file.path(script_dir, "string_network.R"), local = TRUE)

cache_file <- file.path(cache_dir, "Kang18_8vs8.rds")
if (opts$refresh_data && file.exists(cache_file)) unlink(cache_file)
if (!file.exists(cache_file)) {
  message("Downloading/loading Kang18_8vs8 from Bioconductor ExperimentHub...")
  sce <- tryCatch(muscData::Kang18_8vs8(), error = function(e) {
    message("muscData accessor unavailable; trying ExperimentHub EH2259: ", conditionMessage(e))
    eh <- ExperimentHub::ExperimentHub()
    eh[["EH2259"]]
  })
  saveRDS(sce, cache_file, compress = FALSE)
} else {
  message("Using cached dataset: ", cache_file)
  sce <- readRDS(cache_file)
}
if (!inherits(sce, "SingleCellExperiment")) stop("Loaded data is not a SingleCellExperiment.", call. = FALSE)
if (!"counts" %in% SummarizedExperiment::assayNames(sce)) stop("Dataset has no raw 'counts' assay.", call. = FALSE)
counts <- SummarizedExperiment::assay(sce, "counts")
cd <- as.data.frame(SummarizedExperiment::colData(sce))
needed <- c("multiplets", "cell", "ind", "stim")
if (!all(needed %in% names(cd))) stop("Expected metadata columns missing: ", paste(setdiff(needed, names(cd)), collapse = ", "), call. = FALSE)
keep <- !is.na(cd$multiplets) & as.character(cd$multiplets) == "singlet" & !is.na(cd$cell) & !is.na(cd$ind) & !is.na(cd$stim)
cd <- cd[keep, , drop = FALSE]
counts <- counts[, keep, drop = FALSE]
meta <- data.frame(donor = as.character(cd$ind), condition = as.character(cd$stim),
                   cell_type = as.character(cd$cell), row.names = colnames(counts),
                   stringsAsFactors = FALSE)
if (!all(c("ctrl", "stim") %in% unique(meta$condition))) stop("Expected ctrl and stim conditions in dataset.", call. = FALSE)

# Clear only known generated result files, so a failed run cannot leave stale
# files that look like outputs from this run.
cell_types <- sort(unique(meta$cell_type))
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
known_outputs <- c("celltype_response_summary.csv", "hallmark_gsea_results.csv",
  "gsea_hallmark_dotplot.png", "gsea_hallmark_dotplot.svg", "string_network_edges.csv",
  "string_network_nodes.csv", "string_network_raw.tsv", "string_network.png",
  "string_network.svg", "kang_interferon_report.html", "study_manifest.json", "manifest.json",
  paste0("deseq2_", safe_name(cell_types), ".csv"))
for (f in known_outputs) if (file.exists(file.path(out_dir, f))) unlink(file.path(out_dir, f))
# Prior versions used a slightly different filename sanitizer; remove only
# generated DESeq2 CSVs so legacy files cannot be mistaken for this run.
old_de_files <- list.files(out_dir, pattern = "^deseq2_.*\\.csv$", full.names = TRUE)
if (length(old_de_files)) unlink(old_de_files)

message("Running paired pseudobulk and DESeq2 for ", length(cell_types), " annotated cell types...")
summary_rows <- list(); de_tables <- list(); skipped <- character()
for (ct in cell_types) {
  message("  - ", ct)
  # The shared helper expects cells x genes. Preserve integer UMI values.
  fit <- tryCatch({
    pb_r <- paired_pseudobulk(counts = counts, metadata = meta,
      cell_type = ct, reference = "ctrl", comparison = "stim", min_cells = 20L, min_pairs = 3L)
    fit_paired_deseq2(pb_r, reference = "ctrl", comparison = "stim")
  }, error = function(e) e)
  if (inherits(fit, "error")) { skipped <- c(skipped, paste(ct, conditionMessage(fit), sep = ": ")); next }
  res <- fit$results
  res$gene <- as.character(res$gene)
  res <- res[order(res$padj, na.last = TRUE), , drop = FALSE]
  filename <- paste0("deseq2_", safe_name(ct), ".csv")
  utils::write.csv(res, file.path(out_dir, filename), row.names = FALSE, na = "")
  de_tables[[ct]] <- res
  summary_rows[[ct]] <- data.frame(cell_type = ct,
    complete_donor_pairs = length(unique(pb_r$metadata$donor)), genes_tested = nrow(res),
    significant_degs = sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1),
    upregulated = sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 1),
    downregulated = sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange <= -1),
    stringsAsFactors = FALSE)
}
if (!length(de_tables)) stop("No cell types produced a valid paired DESeq2 result. Check metadata and cell thresholds.", call. = FALSE)
summary_df <- do.call(rbind, summary_rows)
# Add canonical interferon effectors to the cross-lineage summary.
for (gene in c("ISG15", "IFIT1", "MX1", "OAS1")) {
  summary_df[[paste0(gene, "_log2FC")]] <- vapply(summary_df$cell_type, function(ct) {
    d <- de_tables[[ct]]; value <- d$log2FoldChange[match(gene, d$gene)]
    if (!length(value) || is.na(value)) NA_real_ else value
  }, numeric(1))
}
utils::write.csv(summary_df, file.path(out_dir, "celltype_response_summary.csv"), row.names = FALSE, na = "")

# Hallmark GSEA on CD14+ monocytes from the same DESeq2 run.
primary <- "CD14+ Monocytes"
if (!primary %in% names(de_tables)) stop("Primary GSEA cell type 'CD14+ Monocytes' was not analyzable.", call. = FALSE)
d <- de_tables[[primary]]
ranks <- d$stat
names(ranks) <- d$gene
valid <- is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))
ranks <- ranks[valid]
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)
msig_args <- names(formals(msigdbr::msigdbr))
msig <- if ("collection" %in% msig_args) msigdbr::msigdbr(species = "Homo sapiens", collection = "H") else msigdbr::msigdbr(species = "Homo sapiens", category = "H")
term2gene <- unique(msig[, c("gs_name", "gene_symbol")])
set.seed(42)
gsea <- clusterProfiler::GSEA(geneList = ranks, TERM2GENE = term2gene,
  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH",
  seed = TRUE, verbose = FALSE)
gsea_df <- if (is.null(gsea)) data.frame() else as.data.frame(gsea)
if (nrow(gsea_df)) {
  if ("core_enrichment" %in% names(gsea_df)) gsea_df$leading_edge <- gsea_df$core_enrichment
  utils::write.csv(gsea_df, file.path(out_dir, "hallmark_gsea_results.csv"), row.names = FALSE, na = "")
} else {
  gsea_df <- data.frame(ID = character(), Description = character(), setSize = integer(),
                        enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
                        p.adjust = numeric(), leading_edge = character())
  utils::write.csv(gsea_df, file.path(out_dir, "hallmark_gsea_results.csv"), row.names = FALSE)
}

plot_gsea <- function(path, device = c("png", "svg")) {
  device <- match.arg(device)
  if (device == "png") grDevices::png(path, width = 1600, height = 1000, res = 160) else grDevices::svg(path, width = 10, height = 6.25)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!nrow(gsea_df)) { plot.new(); title("Hallmark GSEA: no pathways returned"); return(invisible(NULL)) }
  idx <- order(abs(gsea_df$NES), decreasing = TRUE)[seq_len(min(15, nrow(gsea_df)))]
  z <- gsea_df[idx, , drop = FALSE]; z <- z[order(z$NES), , drop = FALSE]
  q <- pmax(as.numeric(z$p.adjust), 1e-300)
  colors <- grDevices::colorRampPalette(c("#2b6cb0", "#f6ad55", "#c53030"))(100)
  col_index <- pmax(1, pmin(100, round((z$NES - min(z$NES)) / max(diff(range(z$NES)), 1e-8) * 99) + 1))
  par(mar = c(5, 19, 4, 2)); plot(z$NES, seq_len(nrow(z)), pch = 19, cex = 1 + 2 * pmin(1, -log10(q) / 10),
       col = colors[col_index], yaxt = "n", ylab = "", xlab = "Normalized enrichment score (NES)",
       main = "Kang IFN-beta: Hallmark GSEA")
  axis(2, at = seq_len(nrow(z)), labels = gsub("^HALLMARK_", "", z$Description), las = 2, cex.axis = .78)
  abline(v = 0, lty = 2, col = "grey50")
}
plot_gsea(file.path(out_dir, "gsea_hallmark_dotplot.png"), "png")
plot_gsea(file.path(out_dir, "gsea_hallmark_dotplot.svg"), "svg")

# Query and preserve the raw human STRING response, then create deterministic
# filtered tables and circular-layout figures from that captured response.
effectors <- c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "USP18", "HERC5", "BST2", "DDX58", "IFIH1", "CXCL10", "CXCL11", "STAT1", "EIF2AK2")
string_url <- paste0("https://string-db.org/api/tsv/network?identifiers=", utils::URLencode(paste(effectors, collapse = "\r"), reserved = TRUE), "&species=9606&required_score=700")
string_raw <- tryCatch(utils::read.delim(string_url, sep = "\t", header = TRUE, check.names = FALSE, quote = "", stringsAsFactors = FALSE),
                       error = function(e) stop("STRING request failed: ", conditionMessage(e), call. = FALSE))
utils::write.table(string_raw, file.path(out_dir, "string_network_raw.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (nrow(string_raw)) {
  score <- as.numeric(string_raw$score)
  finite_scores <- score[is.finite(score)]
  if (length(finite_scores) && max(finite_scores) > 1) score <- score / 1000
  interactions <- data.frame(preferredName_A = string_raw$preferredName_A,
    preferredName_B = string_raw$preferredName_B, score = score, stringsAsFactors = FALSE)
} else interactions <- data.frame(preferredName_A = character(), preferredName_B = character(), score = numeric())
net <- parse_string_interactions(interactions, min_score = 0.700)
all_nodes <- sort(unique(c(effectors, as.character(net$nodes$gene))))
node_df <- data.frame(gene = all_nodes, degree = vapply(all_nodes, function(g) as.integer(sum(net$edges$source == g | net$edges$target == g)), integer(1)), stringsAsFactors = FALSE)
node_df <- annotate_interferon_effectors(node_df)
utils::write.csv(net$edges, file.path(out_dir, "string_network_edges.csv"), row.names = FALSE)
utils::write.csv(node_df, file.path(out_dir, "string_network_nodes.csv"), row.names = FALSE)
plot_network <- function(path, device = c("png", "svg")) {
  device <- match.arg(device)
  if (device == "png") grDevices::png(path, width = 1400, height = 1200, res = 160) else grDevices::svg(path, width = 8.75, height = 7.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  n <- nrow(node_df); theta <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n + 1L)[seq_len(n)]
  xy <- cbind(cos(theta), sin(theta)); rownames(xy) <- node_df$gene
  plot(NA, xlim = c(-1.35, 1.35), ylim = c(-1.2, 1.2), axes = FALSE, xlab = "", ylab = "", asp = 1,
       main = "STRING human interaction network (score >= 0.700)")
  if (nrow(net$edges)) for (i in seq_len(nrow(net$edges))) {
    a <- xy[net$edges$source[i], ]; b <- xy[net$edges$target[i], ]
    if (all(is.finite(c(a, b)))) segments(a[1], a[2], b[1], b[2], col = grDevices::adjustcolor("#718096", alpha.f = .45), lwd = .6 + 2 * net$edges$score[i])
  }
  colors <- ifelse(grepl("ISGylation", node_df$functional_cluster), "#d53f8c",
    ifelse(grepl("IFIT", node_df$functional_cluster), "#3182ce",
    ifelse(grepl("OAS", node_df$functional_cluster), "#38a169", "#dd6b20")))
  symbols <- ifelse(node_df$degree > 0, 21, 21)
  points(xy[,1], xy[,2], pch = symbols, bg = colors, col = "white", cex = 1.4 + .12 * node_df$degree)
  text(xy[,1] * 1.13, xy[,2] * 1.13, labels = node_df$gene, cex = .62, xpd = TRUE)
  mtext(paste(nrow(net$edges), "edges;", sum(node_df$degree > 0), "connected nodes; STRING response saved as TSV"), side = 1, line = 1, cex = .75)
}
plot_network(file.path(out_dir, "string_network.png"), "png")
plot_network(file.path(out_dir, "string_network.svg"), "svg")

# Self-contained HTML with embedded figures and tables.
escape_html <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE); x <- gsub("\"", "&quot;", x, fixed = TRUE); x }
table_html <- function(df, n = 25L) {
  if (!nrow(df)) return("<p>No rows were produced.</p>")
  df <- head(df, n); head_html <- paste0("<tr>", paste0("<th>", escape_html(names(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(row) paste0("<tr>", paste0("<td>", escape_html(row), "</td>", collapse = ""), "</tr>"))
  paste0("<div class='table-wrap'><table>", head_html, paste(rows, collapse = ""), "</table></div>")
}
embed_png <- function(file) paste0("data:image/png;base64,", paste(base64enc::base64encode(readBin(file, "raw", n = file.info(file)$size), linewidth = 0), collapse = ""))
report <- paste0("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>Kang et al. 2018 IFN-beta study</title><style>body{font:16px/1.55 system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#243044}h1,h2{color:#152238}.notice{background:#fff4d6;border-left:5px solid #d69e2e;padding:1rem}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.9rem}td,th{border-bottom:1px solid #dce3eb;padding:.45rem;text-align:left}th{background:#edf2f7}img{max-width:100%;height:auto;border:1px solid #e2e8f0}</style><body>",
  "<h1>Kang et al. (2018): PBMC interferon-beta response</h1><p>Rebuilt from the cached/downloaded muscData <code>Kang18_8vs8</code> dataset (GEO <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583'>GSE96583</a>; ExperimentHub EH2259). Paired donor pseudobulk counts were tested with DESeq2 using <code>~ donor + condition</code>.</p>",
  "<div class='notice'><strong>Research prototype:</strong> this reproducible secondary analysis is for research and demonstration, not clinical decision-making. Pathway interpretation and STRING associations are descriptive, not proof of causal mechanisms. STRING results come from an external live service; the exact response used is saved alongside the outputs.</div>",
  "<h2>Cell-type differential expression</h2>", table_html(summary_df),
  "<h2>Hallmark GSEA</h2><p>Pre-ranked by the CD14+ monocyte DESeq2 Wald statistic; MSigDB Hallmark gene sets.</p><img alt='Hallmark GSEA dot plot' src='", embed_png(file.path(out_dir, "gsea_hallmark_dotplot.png")), "'>",
  "<h2>STRING protein interaction network</h2><p>Human taxon 9606; required interaction score 700/1000 (0.700).</p><img alt='STRING interaction network' src='", embed_png(file.path(out_dir, "string_network.png")), "'>",
  "<h2>Top enriched pathways</h2>", table_html(if (nrow(gsea_df)) gsea_df[order(gsea_df$p.adjust), intersect(c("Description", "NES", "pvalue", "p.adjust", "leading_edge"), names(gsea_df)), drop = FALSE] else gsea_df),
  "<h2>Reproducibility</h2><p>Run command: <code>Rscript studies/kang_ifnb/run_study.R</code>. The adjacent machine-readable manifest records runtime package versions, parameters, random seed, dataset-cache checksum, output paths, and SHA-256 hashes. Output tables and figures are generated from the same run.</p><p>Reference: Kang et al., Nature Biotechnology (2018), DOI <a href='https://doi.org/10.1038/nbt.4042'>10.1038/nbt.4042</a>.</p></body></html>")
report_file <- file.path(out_dir, "kang_interferon_report.html")
writeLines(report, report_file, useBytes = TRUE)

# The cached source object is fingerprinted; manifests do not hash themselves.
output_files <- file.path(out_dir, known_outputs)
output_files <- output_files[file.exists(output_files)]
artifact_descriptions <- c("celltype_response_summary.csv" = "Per-cell-type donor-paired DE summary and canonical ISG effects",
  "hallmark_gsea_results.csv" = "ClusterProfiler GSEA results for MSigDB Hallmark sets",
  "gsea_hallmark_dotplot.png" = "Hallmark GSEA dot plot", "gsea_hallmark_dotplot.svg" = "Vector Hallmark GSEA dot plot",
  "string_network_raw.tsv" = "Captured raw STRING API response for this run", "string_network_edges.csv" = "Filtered high-confidence STRING edges",
  "string_network_nodes.csv" = "STRING network nodes, degrees, and functional labels", "string_network.png" = "STRING network figure",
  "string_network.svg" = "Vector STRING network figure", "kang_interferon_report.html" = "Self-contained HTML report")
package_names <- unique(c(required_cran, required_bioc))
versions <- setNames(vapply(package_names, function(p) as.character(utils::packageVersion(p)), character(1)), package_names)
manifest <- list(manifest_version = "2.0.0", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  study = list(geo_accession = "GSE96583", experimenthub_id = "EH2259", object = "muscData::Kang18_8vs8", organism = "Homo sapiens", taxid = 9606),
  input_data = list(cache_file = normalizePath(cache_file, mustWork = TRUE), sha256 = digest::digest(file = cache_file, algo = "sha256"),
                    cells_used = ncol(counts), features = nrow(counts)),
  software = list(R = R.version.string, Bioconductor = if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else NA_character_, packages = as.list(versions)),
  parameters = list(conditions = c("ctrl", "stim"), design = "~ donor + condition", contrast = c("condition", "stim", "ctrl"),
    minimum_cells_per_donor_condition = 20L, minimum_complete_pairs = 3L, DEG_adjusted_p_threshold = 0.05,
    DEG_absolute_log2_fold_change_threshold = 1, GSEA_collection = "MSigDB Hallmark", GSEA_min_size = 15L,
    GSEA_max_size = 500L, GSEA_seed = 42L, STRING_taxid = 9606L, STRING_required_score = 700L, STRING_score_threshold = 0.7),
  skipped_cell_types = as.list(skipped),
  output_files = list())
manifest$output_files <- lapply(output_files, function(f) {
  desc <- artifact_descriptions[basename(f)]
  list(path = basename(f), format = tools::file_ext(f), size_bytes = unname(file.info(f)$size),
       sha256 = digest::digest(file = f, algo = "sha256"),
       description = if (length(desc) && !is.na(desc)) unname(desc) else paste("Generated", tools::file_ext(f), "artifact"))
})
jsonlite::write_json(manifest, file.path(out_dir, "study_manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
file.copy(file.path(out_dir, "study_manifest.json"), file.path(out_dir, "manifest.json"), overwrite = TRUE)
message("Completed. Results and self-contained HTML report: ", normalizePath(out_dir, mustWork = TRUE))
if (length(skipped)) message("Skipped cell types (insufficient pairs or model issue): ", paste(skipped, collapse = " | "))
