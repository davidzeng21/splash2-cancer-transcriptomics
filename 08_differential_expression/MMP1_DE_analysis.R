#!/usr/bin/env Rscript
# MMP1 Differential Expression Analysis: CPTAC-3 LUAD (Tumor vs. Normal)
# Uses TCGAbiolinks to locate/download GDC files, then reads MMP1 directly
# (no GDCprepare — avoids loading all 60K genes just for one gene)

# ── 0. Install dependencies if needed ─────────────────────────────────────────
cran_repo  <- "https://cloud.r-project.org"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = cran_repo)

all_pkgs <- c("TCGAbiolinks", "data.table", "ggplot2", "dplyr")

missing_pkgs <- all_pkgs[!all_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install_ok <- tryCatch({
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE); TRUE
  }, error = function(e) {
    message("BiocManager install failed (", conditionMessage(e), "); using fallback.")
    FALSE
  })
  if (!install_ok) {
    bioc_ver   <- tryCatch(as.character(BiocManager::version()), error = function(e) "3.22")
    bioc_repos <- c(
      CRAN     = cran_repo,
      BioCsoft = paste0("https://bioconductor.org/packages/", bioc_ver, "/bioc")
    )
    still_missing <- missing_pkgs[!missing_pkgs %in% rownames(installed.packages())]
    if (length(still_missing) > 0)
      install.packages(still_missing, repos = bioc_repos)
  }
}

# ── 1. Load libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ── 2. Parse sample IDs from the input manifest ────────────────────────────────
input_file <- "input_CPTAC_107.txt"
if (!file.exists(input_file))
  stop("Input file not found: ", input_file,
       "\nPlease run this script from the repo root directory.")

manifest <- readLines(input_file)
manifest <- manifest[nzchar(trimws(manifest))]   # drop blank lines

# Each line: "<patient_id>-tumor" or "<patient_id>-normal"
# Patient IDs are the text before the first space (the FASTQ path follows)
sample_labels <- sub("\\s.*", "", manifest)       # keep only ID column
patient_ids   <- unique(sub("-(tumor|normal)$", "", sample_labels))

message("Parsed ", length(patient_ids), " unique patient IDs from ", input_file)

# ── 3. GDC query ──────────────────────────────────────────────────────────────
# NOTE: In CPTAC-3, GDC sample barcodes are CPT0* format, which is DIFFERENT
# from the case submitter IDs (C3L-* / C3N-*) in our input file.
# Using barcode= with C3L-* IDs causes "None of the barcodes were matched".
# Fix: query without barcode filter, then subset the manifest by cases_submitter_id.
message("\nQuerying GDC for all CPTAC-3 RNA-seq files ...")

query <- GDCquery(
  project       = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Filter manifest to only our 107 patients
manifest_df <- getResults(query)
message("Total files returned by GDC query: ", nrow(manifest_df))
message("Available columns: ", paste(colnames(manifest_df), collapse = ", "))

# cases_submitter_id contains the C3L-* / C3N-* IDs we need
case_col <- intersect(
  c("cases_submitter_id", "cases.submitter_id", "case_submitter_id", "submitter_id"),
  colnames(manifest_df)
)[1]

if (is.na(case_col))
  stop("Cannot identify case submitter ID column. Columns available: ",
       paste(colnames(manifest_df), collapse = ", "))

# cases_submitter_id may be a list-column in some TCGAbiolinks versions; flatten
if (is.list(manifest_df[[case_col]])) {
  manifest_df[[case_col]] <- vapply(manifest_df[[case_col]], `[[`, character(1), 1)
}

filtered_manifest <- manifest_df[manifest_df[[case_col]] %in% patient_ids, ]
n_matched <- length(unique(filtered_manifest[[case_col]]))
message("Matched ", nrow(filtered_manifest), " files for ", n_matched,
        " / ", length(patient_ids), " patients")

if (nrow(filtered_manifest) == 0)
  stop("No files matched any of the 107 patient IDs.\n",
       "Head of '", case_col, "' column: ",
       paste(head(manifest_df[[case_col]], 5), collapse = ", "))

query$results[[1]] <- filtered_manifest

# ── 4. Download (skips files already present) ──────────────────────────────────
message("Downloading count files (this may take a while on first run) ...")
GDCdownload(query, method = "api", files.per.chunk = 10)

# ── 5. Read MMP1 counts directly from the downloaded files ────────────────────
# Cache: if results/mmp1_counts_cache.rds exists, skip file reading entirely.
CACHE_FILE <- "results/mmp1_counts_cache.rds"
# GDCprepare loads all ~60 K genes — wasteful for one gene.
# Instead, scan each STAR .tsv directly: grab the MMP1 row and total counts
# for library-size (CPM) normalization.  This takes ~seconds.
MMP1_ENSEMBL <- "ENSG00000196611"   # version suffix stripped for robustness

if (file.exists(CACHE_FILE)) {
  message("Loading cached count data from ", CACHE_FILE)
  cache       <- readRDS(CACHE_FILE)
  mf          <- cache$mf
  mmp1_ensembl <- cache$mmp1_ensembl
} else {

# Build a map: file_id -> (cases.submitter_id, sample_type) from the manifest
mf <- getResults(query)
# Flatten list-columns, then take only the first semicolon-delimited token.
# GDC manifest sometimes stores multi-aliquot entries as "Primary Tumor;Primary Tumor"
flatten1 <- function(x) {
  v <- if (is.list(x)) vapply(x, function(e) e[1], character(1)) else as.character(x)
  sub(";.*", "", v)   # keep only the first semicolon-delimited token
}
mf$patient     <- flatten1(mf$cases.submitter_id)
mf$sample_type <- flatten1(mf$sample_type)

# Keep Primary Tumor and Solid Tissue Normal only
mf <- mf[mf$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]

# Locate each downloaded file on disk by building a filename → path index.
# This avoids relying on the manifest 'id' column matching the directory UUID.
gdc_root <- file.path("GDCdata", "CPTAC-3",
                      "Transcriptome_Profiling",
                      "Gene_Expression_Quantification")

all_tsvs  <- list.files(gdc_root, pattern = "\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
tsv_index <- setNames(all_tsvs, basename(all_tsvs))
message("TSV files found on disk: ", length(tsv_index))

read_mmp1_counts <- function(file_name) {
  tsv <- tsv_index[file_name]
  if (is.na(tsv) || !file.exists(tsv)) return(c(mmp1 = NA_real_, total = NA_real_))
  dt <- data.table::fread(tsv, sep = "\t", skip = 1L, header = TRUE,
                          select = c("gene_id", "gene_name", "unstranded"),
                          showProgress = FALSE)
  dt <- dt[!startsWith(dt$gene_id, "N_"), ]
  mmp1_row <- dt[startsWith(dt$gene_id, MMP1_ENSEMBL), ]
  mmp1_cnt <- if (nrow(mmp1_row) > 0) mmp1_row$unstranded[1] else NA_real_
  c(mmp1 = as.numeric(mmp1_cnt), total = sum(dt$unstranded, na.rm = TRUE))
}

message("Reading MMP1 counts from ", nrow(mf), " files ...")
counts_mat <- sapply(mf$file_name, read_mmp1_counts)
# counts_mat is 2 x n (rows: mmp1, total)
mf$mmp1_count <- counts_mat["mmp1", ]
mf$lib_size   <- counts_mat["total", ]
mf <- mf[!is.na(mf$mmp1_count), ]

# CPM normalization: counts / library_size * 1e6
mf$mmp1_cpm    <- mf$mmp1_count / mf$lib_size * 1e6
mf$log2_mmp1   <- log2(mf$mmp1_cpm + 1)

# Identify MMP1 Ensembl ID (with version) from first available file
mmp1_ensembl <- {
  tsv1 <- tsv_index[mf$file_name[1]]
  dt1  <- data.table::fread(tsv1, sep = "\t", skip = 1L, header = TRUE,
                             select = c("gene_id", "gene_name"),
                             showProgress = FALSE)
  row  <- dt1[startsWith(dt1$gene_id, MMP1_ENSEMBL), ]
  if (nrow(row)) row$gene_id[1] else MMP1_ENSEMBL
}
message("MMP1 Ensembl ID: ", mmp1_ensembl)
message("Samples read: ", nrow(mf))

# Save cache for fast plot-only reruns
dir.create("results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(mf = mf, mmp1_ensembl = mmp1_ensembl), CACHE_FILE)
message("Cached count data to ", CACHE_FILE)

} # end else (no cache)

# ── 6. Keep matched tumor/normal pairs ────────────────────────────────────────
pair_counts <- mf %>%
  group_by(patient) %>%
  summarise(
    has_tumor  = any(sample_type == "Primary Tumor"),
    has_normal = any(sample_type == "Solid Tissue Normal"),
    .groups = "drop"
  ) %>%
  filter(has_tumor & has_normal)

complete_patients <- pair_counts$patient
mf <- mf[mf$patient %in% complete_patients, ]
message("Complete matched pairs: ", length(complete_patients),
        "  (", sum(mf$sample_type == "Primary Tumor"), " tumor + ",
        sum(mf$sample_type == "Solid Tissue Normal"), " normal)")

# ── 7. Paired Wilcoxon signed-rank test ───────────────────────────────────────
tumor_df  <- mf[mf$sample_type == "Primary Tumor", ]
normal_df <- mf[mf$sample_type == "Solid Tissue Normal", ]

common_pats    <- intersect(tumor_df$patient, normal_df$patient)
tumor_ordered  <- tumor_df$log2_mmp1[match(common_pats, tumor_df$patient)]
normal_ordered <- normal_df$log2_mmp1[match(common_pats, normal_df$patient)]

wt <- wilcox.test(tumor_ordered, normal_ordered, paired = TRUE, exact = FALSE)

# log2FC from medians of CPM values
med_tumor  <- median(tumor_df$mmp1_cpm)
med_normal <- median(normal_df$mmp1_cpm)
log2fc     <- log2((med_tumor + 1) / (med_normal + 1))

# Vectors needed for plotting
sample_type_vec <- mf$sample_type
patient_vec     <- mf$patient
log2_mmp1       <- mf$log2_mmp1

mmp1_res <- data.frame(
  gene_name              = "MMP1",
  ensembl_id             = mmp1_ensembl,
  log2FC_tumor_vs_normal = round(log2fc, 3),
  median_tumor_cpm       = round(med_tumor, 1),
  median_normal_cpm      = round(med_normal, 1),
  wilcox_statistic       = wt$statistic,
  p_value                = wt$p.value,
  p_value_BH             = p.adjust(wt$p.value, method = "BH", n = 1),
  n_pairs                = length(common_pats),
  test                   = "paired Wilcoxon signed-rank"
)

message("\n── MMP1 DE Result ──────────────────────────────────────────────────────")
print(t(mmp1_res))

# ── 8. Save result ─────────────────────────────────────────────────────────────
dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(mmp1_res, file = "results/MMP1_DE_result.csv", row.names = FALSE)
message("Saved: results/MMP1_DE_result.csv")

# ── 9. Paired boxplot ──────────────────────────────────────────────────────────
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

plot_df <- data.frame(
  patient    = patient_vec,
  group      = factor(sample_type_vec,
                      levels = c("Solid Tissue Normal", "Primary Tumor")),
  expression = log2_mmp1
)

p_label <- if (mmp1_res$p_value < 0.001) "p < 0.001" else
             sprintf("p = %.3f", mmp1_res$p_value)

p_box <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.4) +
  geom_point(aes(color = group), position = position_jitter(width = 0.05),
             size = 1.5, alpha = 0.7) +
  geom_line(aes(group = patient), alpha = 0.25, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(values  = c("Solid Tissue Normal" = "#4393C3",
                                "Primary Tumor"       = "#D6604D")) +
  scale_color_manual(values = c("Solid Tissue Normal" = "#2166AC",
                                "Primary Tumor"       = "#B2182B")) +
  labs(
    title    = "MMP1 expression: Tumor vs. Normal",
    subtitle = sprintf("CPTAC-3 LUAD  |  log2FC = %s  |  %s  (paired Wilcoxon)",
                       round(log2fc, 2), p_label),
    x        = NULL,
    y        = "log2(CPM + 1)",
    caption  = paste0("n = ", length(common_pats), " matched pairs")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position  = "none",
        plot.title        = element_text(face = "bold"),
        plot.subtitle     = element_text(size = 11),
        axis.text.x       = element_text(size = 11),
        plot.margin       = margin(t = 10, r = 20, b = 10, l = 10))

ggsave("plots/MMP1_boxplot.pdf", p_box, width = 6.5, height = 5.5)
ggsave("plots/MMP1_boxplot.png", p_box, width = 6.5, height = 5.5, dpi = 200)
message("Saved: plots/MMP1_boxplot.pdf")

message("\nAnalysis complete.")
