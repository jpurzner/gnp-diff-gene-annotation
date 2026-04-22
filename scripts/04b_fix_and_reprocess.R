suppressMessages({library(tidyverse)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. RE-PARSE ALL BATCHES WITH LENIENT PARSING
# =============================================================================
cat("Re-parsing all batches with tolerance for extra columns...\n")

batch_files <- list.files("agent_batches", pattern="batch_[0-9]+_results\\.tsv$", full.names=TRUE)

parse_batch_lenient <- function(f) {
  lines <- readLines(f, warn=FALSE)
  if (length(lines) < 2) return(NULL)

  # Parse each line manually
  parsed <- lapply(lines, function(l) strsplit(l, "\t", fixed=TRUE)[[1]])

  # Expected 17 cols
  expected_cols <- c("gene","protein_class","molecular_action","key_substrates",
    "subcellular","tissues","cell_types","pathways","ko_phenotype","ko_severity",
    "neural_phenotype","cerebellar_evidence","cancer_types","cancer_role",
    "neuro_disease","differentiation_role","one_line_summary")

  # For rows with >17 fields, merge extras into the last field
  data_rows <- parsed[-1]
  fixed <- lapply(data_rows, function(r) {
    if (length(r) == 17) return(r)
    if (length(r) > 17) {
      # Merge overflow into field 17 (summary)
      return(c(r[1:16], paste(r[17:length(r)], collapse=" ")))
    }
    if (length(r) < 17) {
      return(c(r, rep(NA, 17-length(r))))
    }
  })

  df <- as.data.frame(do.call(rbind, fixed), stringsAsFactors=FALSE)
  colnames(df) <- expected_cols
  df
}

all_parsed <- lapply(batch_files, parse_batch_lenient)
all_parsed <- all_parsed[!sapply(all_parsed, is.null)]

merged <- do.call(rbind, all_parsed)
cat("Total rows parsed:", nrow(merged), "\n")
merged <- merged[!duplicated(merged$gene) & !is.na(merged$gene) & merged$gene != "", ]
cat("Unique genes:", nrow(merged), "\n")

# Check missing genes
master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))
diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude))

all_diff <- diff_genes$mgi_symbol
missing <- setdiff(all_diff, merged$gene)
cat("Genes still missing:", length(missing), "\n")
if (length(missing) > 0) cat("  ", paste(head(missing, 30), collapse=", "), "\n")

# Save re-parsed
write.csv(merged, "agent_extracted_annotations.csv", row.names=FALSE)
cat("Saved: agent_extracted_annotations.csv (re-parsed)\n")

# =============================================================================
# 2. IDENTIFY GENES NEEDING RE-PROCESSING
# =============================================================================
# Count fields filled
fields <- c("protein_class","molecular_action","key_substrates","subcellular","tissues",
  "cell_types","pathways","ko_phenotype","ko_severity","neural_phenotype",
  "cerebellar_evidence","cancer_types","cancer_role","neuro_disease",
  "differentiation_role","one_line_summary")

is_filled <- function(r) {
  sum(!is.na(r) & r != "" & r != "NA" & r != "unknown" & r != "none" & r != "none known" &
      r != "other" & r != "no")
}
merged$n_filled <- apply(merged[, fields], 1, is_filled)

# Load abstracts
compiled <- readRDS("diff_gene_tiered_compiled.rds")
merged_with_abs <- merged %>%
  dplyr::left_join(compiled %>% dplyr::select(symbol, n_papers),
    by=c("gene"="symbol"))

# Genes with >= 5 abstracts but < 6 fields filled
to_reprocess <- merged_with_abs %>%
  dplyr::filter(n_papers >= 5, n_filled < 6) %>%
  dplyr::pull(gene) %>%
  unique()

# Add any truly missing genes that have abstracts
missing_with_abs <- compiled %>%
  dplyr::filter(symbol %in% missing, n_papers >= 5) %>%
  dplyr::pull(symbol)

to_reprocess <- unique(c(to_reprocess, missing_with_abs))
cat("\nGenes needing reprocessing (>=5 papers but poor annotation):", length(to_reprocess), "\n")
cat("  Gene list:", paste(to_reprocess, collapse=", "), "\n")

# =============================================================================
# 3. CREATE REPROCESS BATCHES
# =============================================================================
if (length(to_reprocess) > 0) {
  cat("\nCreating reprocess batch files...\n")

  # Load all data sources
  uniprot <- read.csv("diff_gene_uniprot_functions.csv", stringsAsFactors=FALSE)
  pheno <- read.csv("diff_gene_mgi_phenotypes.csv", stringsAsFactors=FALSE)

  gene_data <- diff_genes %>%
    dplyr::select(ensembl_gene_id, mgi_symbol, description, H3K27me3_bin, t50_cat) %>%
    dplyr::distinct(mgi_symbol, .keep_all=TRUE) %>%
    dplyr::left_join(uniprot %>% dplyr::select(mgi_symbol, uniprot_function), by="mgi_symbol") %>%
    dplyr::left_join(compiled %>% dplyr::select(symbol, combined_text, titles, n_papers, n_tier1),
      by=c("mgi_symbol"="symbol")) %>%
    dplyr::left_join(pheno %>% dplyr::select(mgi_symbol, phenotype_terms), by="mgi_symbol")

  batch <- gene_data %>% dplyr::filter(mgi_symbol %in% to_reprocess)

  # Truncate abstracts helper — MORE context for reprocess (more papers, more chars)
  truncate_abstracts <- function(combined_text, max_papers=15, max_chars_per=600) {
    if (is.na(combined_text) || nchar(combined_text) < 10) return("")
    papers <- strsplit(combined_text, " \\|\\|\\| ")[[1]]
    papers <- head(papers, max_papers)
    papers <- sapply(papers, function(p) {
      if (nchar(p) > max_chars_per) paste0(substr(p, 1, max_chars_per-3), "...") else p
    })
    paste(papers, collapse="\n---\n")
  }

  # Split into batches of ~20 genes
  n_batches <- max(1, ceiling(nrow(batch) / 20))
  if (n_batches == 1) {
    batch_idx <- list(`1` = 1:nrow(batch))
  } else {
    batch_idx <- split(1:nrow(batch), cut(1:nrow(batch), n_batches, labels=FALSE))
  }

  for (b in seq_along(batch_idx)) {
    sub <- batch[batch_idx[[b]], ]
    lines <- character()
    lines <- c(lines, paste0("# REPROCESS Batch ", b, " — ", nrow(sub), " genes"))
    lines <- c(lines, paste0("# These genes had >=5 abstracts but <6 filled fields in initial pass"))
    lines <- c(lines, "")

    for (i in 1:nrow(sub)) {
      g <- sub[i, ]
      lines <- c(lines, paste0("========== GENE: ", g$mgi_symbol, " =========="))
      lines <- c(lines, paste0("Description: ", ifelse(is.na(g$description), "NA", g$description)))
      lines <- c(lines, paste0("H3K27me3: ", g$H3K27me3_bin, " | Timing: ", g$t50_cat))
      lines <- c(lines, "")

      if (!is.na(g$uniprot_function) && nchar(g$uniprot_function) > 10) {
        lines <- c(lines, "--- UniProt Function ---")
        lines <- c(lines, substr(g$uniprot_function, 1, 1500))
        lines <- c(lines, "")
      }

      if (!is.na(g$phenotype_terms) && nchar(g$phenotype_terms) > 5) {
        lines <- c(lines, "--- MGI Knockout Phenotypes ---")
        pt <- g$phenotype_terms
        if (nchar(pt) > 800) pt <- paste0(substr(pt, 1, 797), "...")
        lines <- c(lines, pt)
        lines <- c(lines, "")
      }

      if (!is.na(g$combined_text) && nchar(g$combined_text) > 50) {
        lines <- c(lines, paste0("--- Abstracts (", g$n_papers, " papers, ", g$n_tier1, " tier1) ---"))
        abs_text <- truncate_abstracts(g$combined_text, max_papers=15, max_chars_per=600)
        lines <- c(lines, abs_text)
        lines <- c(lines, "")
      }
      lines <- c(lines, "")
    }

    fname <- sprintf("agent_batches/reprocess_batch_%02d.txt", b)
    writeLines(lines, fname)
    cat(sprintf("  Wrote: %s (%d genes)\n", fname, nrow(sub)))
  }

  # Save list
  writeLines(to_reprocess, "genes_to_reprocess.txt")
}

cat("\n=== DONE ===\n")
