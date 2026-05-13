suppressMessages({library(tidyverse); library(org.Mm.eg.db); library(AnnotationDbi)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

cat("=== PREP EMBRYONIC BATCHES ===\n")
emb <- read.csv("gene_list_embryonic.csv", stringsAsFactors=FALSE)
compiled <- readRDS("tiered_compiled_EMB.rds")

mouse_up <- read.table("/tmp/uniprot_mouse_functions.tsv", sep="\t", header=TRUE,
  stringsAsFactors=FALSE, quote="", fill=TRUE, comment.char="")
mouse_up$primary_gene <- sapply(strsplit(mouse_up$Gene.Names, " "), function(x) x[1])
mouse_up$func_clean <- gsub("\\{ECO:[^\\}]+\\}", "", mouse_up$Function..CC.)
mouse_up$func_clean <- gsub("^FUNCTION: ", "", mouse_up$func_clean)
mouse_up$func_clean <- gsub("\\s+", " ", trimws(mouse_up$func_clean))

pheno <- read.csv("diff_gene_mgi_phenotypes.csv", stringsAsFactors=FALSE)

# 782 genes / 22 per batch = 36 batches
n_per_batch <- 22
emb$idx <- 1:nrow(emb)
emb$batch <- ceiling(emb$idx / n_per_batch)
n_batches <- max(emb$batch)

dir.create("agent_batches_emb", showWarnings=FALSE)

truncate_abstracts <- function(combined_text, max_papers=10, max_chars_per=500) {
  if (is.na(combined_text) || nchar(combined_text) < 10) return("")
  papers <- strsplit(combined_text, " \\|\\|\\| ")[[1]]
  papers <- head(papers, max_papers)
  papers <- sapply(papers, function(p) if (nchar(p) > max_chars_per) paste0(substr(p, 1, max_chars_per-3), "...") else p)
  paste(papers, collapse="\n---\n")
}

for (b in 1:n_batches) {
  sub <- emb %>% dplyr::filter(batch == b)
  lines <- c(paste0("# EMBRYONIC Batch ", b, " — ", nrow(sub), " genes"), "")

  for (i in 1:nrow(sub)) {
    g <- sub[i, ]; sym <- g$mgi_symbol
    lines <- c(lines, paste0("========== GENE: ", sym, " =========="))
    lines <- c(lines, paste0("Description: ", ifelse(is.na(g$description), "NA", g$description)))
    lines <- c(lines, paste0("max_cat: ", g$max_cat, " | t50: ", g$t50_cat,
                              " | H3K27me3: ", g$H3K27me3_bin, " | avg_fc: ", round(g$avg_fc, 2)))
    lines <- c(lines, "")

    up <- mouse_up$func_clean[toupper(mouse_up$primary_gene) == toupper(sym)]
    up <- up[!is.na(up)]
    if (length(up) > 0 && nchar(up[1]) > 10) {
      lines <- c(lines, "--- UniProt Function ---")
      lines <- c(lines, substr(up[1], 1, 1000)); lines <- c(lines, "")
    }

    p_row <- pheno %>% dplyr::filter(mgi_symbol == sym)
    if (nrow(p_row) > 0 && !is.na(p_row$phenotype_terms[1])) {
      lines <- c(lines, "--- MGI Knockout Phenotypes ---")
      pt <- p_row$phenotype_terms[1]
      if (nchar(pt) > 500) pt <- paste0(substr(pt, 1, 497), "...")
      lines <- c(lines, pt); lines <- c(lines, "")
    }

    c_row <- compiled %>% dplyr::filter(symbol == sym)
    if (nrow(c_row) > 0 && !is.na(c_row$combined_text[1])) {
      lines <- c(lines, paste0("--- Abstracts (", c_row$n_papers[1], " papers, ",
        c_row$n_tier1[1], " tier1) ---"))
      lines <- c(lines, truncate_abstracts(c_row$combined_text[1])); lines <- c(lines, "")
    }
    lines <- c(lines, "")
  }

  fname <- sprintf("agent_batches_emb/batch_%02d.txt", b)
  writeLines(lines, fname)
}
cat("Wrote", n_batches, "batch files to agent_batches_emb/\n")
