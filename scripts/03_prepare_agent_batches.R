suppressMessages({library(tidyverse)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

cat("Loading data...\n")

compiled <- readRDS("diff_gene_tiered_compiled.rds")
uniprot <- read.csv("diff_gene_uniprot_functions.csv", stringsAsFactors=FALSE)
umap_df <- read.csv("gene_umap_k10_coordinates.csv", stringsAsFactors=FALSE)
pheno <- read.csv("diff_gene_mgi_phenotypes.csv", stringsAsFactors=FALSE)

master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude)) %>%
  dplyr::select(ensembl_gene_id, mgi_symbol, description, H3K27me3_bin, t50_cat) %>%
  dplyr::distinct(mgi_symbol, .keep_all=TRUE)

# Merge everything
gene_data <- diff_genes %>%
  dplyr::left_join(umap_df %>% dplyr::select(gene, cluster), by=c("mgi_symbol"="gene")) %>%
  dplyr::left_join(uniprot %>% dplyr::select(mgi_symbol, uniprot_function), by="mgi_symbol") %>%
  dplyr::left_join(compiled %>% dplyr::select(symbol, combined_text, titles, n_papers, n_tier1),
    by=c("mgi_symbol"="symbol")) %>%
  dplyr::left_join(pheno %>% dplyr::select(mgi_symbol, phenotype_terms),
    by="mgi_symbol")

# Replace NA cluster with 0
gene_data$cluster[is.na(gene_data$cluster)] <- 0

cat("Genes:", nrow(gene_data), "\n")
cat("Cluster distribution:\n")
print(table(gene_data$cluster))

# Build per-gene text bundles
# For each gene: gene name, description, UniProt, top abstract excerpts, MGI phenotypes
dir.create("agent_batches", showWarnings=FALSE)

# Split abstracts into individual papers (separated by |||)
# Keep only the first ~500 chars of each abstract to stay within limits
truncate_abstracts <- function(combined_text, max_papers=10, max_chars_per=500) {
  if (is.na(combined_text) || nchar(combined_text) < 10) return("")
  papers <- strsplit(combined_text, " \\|\\|\\| ")[[1]]
  papers <- head(papers, max_papers)
  papers <- sapply(papers, function(p) {
    if (nchar(p) > max_chars_per) paste0(substr(p, 1, max_chars_per-3), "...") else p
  })
  paste(papers, collapse="\n---\n")
}

# Group genes by cluster for batching
# Target ~40 genes per batch, split large clusters
batches <- list()
batch_id <- 1

for (cl in sort(unique(gene_data$cluster))) {
  cl_genes <- gene_data %>% dplyr::filter(cluster == cl)
  n_cl <- nrow(cl_genes)

  if (n_cl <= 45) {
    batches[[batch_id]] <- cl_genes
    batch_id <- batch_id + 1
  } else {
    # Split into sub-batches
    n_sub <- ceiling(n_cl / 40)
    idx <- split(1:n_cl, cut(1:n_cl, n_sub, labels=FALSE))
    for (sub in idx) {
      batches[[batch_id]] <- cl_genes[sub, ]
      batch_id <- batch_id + 1
    }
  }
}

cat("\nTotal batches:", length(batches), "\n")
for (b in seq_along(batches)) {
  cat(sprintf("  Batch %d: %d genes (cluster %s)\n", b, nrow(batches[[b]]),
    paste(unique(batches[[b]]$cluster), collapse=",")))
}

# Write each batch as a text file the agent can read
for (b in seq_along(batches)) {
  batch <- batches[[b]]

  lines <- character()
  lines <- c(lines, paste0("# Batch ", b, " — ", nrow(batch), " genes"))
  lines <- c(lines, paste0("# Clusters: ", paste(unique(batch$cluster), collapse=", ")))
  lines <- c(lines, "")

  for (i in 1:nrow(batch)) {
    g <- batch[i, ]
    lines <- c(lines, paste0("========== GENE: ", g$mgi_symbol, " =========="))
    lines <- c(lines, paste0("Description: ", ifelse(is.na(g$description), "NA", g$description)))
    lines <- c(lines, paste0("H3K27me3: ", g$H3K27me3_bin, " | Timing: ", g$t50_cat))
    lines <- c(lines, "")

    if (!is.na(g$uniprot_function) && nchar(g$uniprot_function) > 10) {
      lines <- c(lines, "--- UniProt Function ---")
      up_text <- substr(g$uniprot_function, 1, 1000)
      lines <- c(lines, up_text)
      lines <- c(lines, "")
    }

    if (!is.na(g$phenotype_terms) && nchar(g$phenotype_terms) > 5) {
      lines <- c(lines, "--- MGI Knockout Phenotypes ---")
      # Truncate long phenotype lists
      pt <- g$phenotype_terms
      if (nchar(pt) > 500) pt <- paste0(substr(pt, 1, 497), "...")
      lines <- c(lines, pt)
      lines <- c(lines, "")
    }

    if (!is.na(g$combined_text) && nchar(g$combined_text) > 50) {
      lines <- c(lines, paste0("--- Abstracts (", g$n_papers, " papers, ", g$n_tier1, " tier1) ---"))
      abs_text <- truncate_abstracts(g$combined_text, max_papers=8, max_chars_per=400)
      lines <- c(lines, abs_text)
      lines <- c(lines, "")
    }

    lines <- c(lines, "")
  }

  fname <- sprintf("agent_batches/batch_%02d.txt", b)
  writeLines(lines, fname)
}

cat("\nWrote", length(batches), "batch files to agent_batches/\n")

# Save batch metadata
batch_meta <- data.frame(
  batch = seq_along(batches),
  n_genes = sapply(batches, nrow),
  clusters = sapply(batches, function(b) paste(unique(b$cluster), collapse=",")),
  genes = sapply(batches, function(b) paste(b$mgi_symbol, collapse=",")),
  stringsAsFactors = FALSE
)
write.csv(batch_meta, "agent_batches/batch_metadata.csv", row.names=FALSE)

cat("\n=== DONE ===\n")
