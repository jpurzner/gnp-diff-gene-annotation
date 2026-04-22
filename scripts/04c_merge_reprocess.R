suppressMessages({library(tidyverse); library(openxlsx)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

cat("Loading data...\n")
ann <- read.csv("agent_extracted_annotations.csv", stringsAsFactors=FALSE)
cat("Current annotations:", nrow(ann), "\n")

# Load reprocess results (lenient parse)
parse_batch_lenient <- function(f) {
  lines <- readLines(f, warn=FALSE)
  if (length(lines) < 2) return(NULL)
  parsed <- lapply(lines, function(l) strsplit(l, "\t", fixed=TRUE)[[1]])
  expected_cols <- c("gene","protein_class","molecular_action","key_substrates",
    "subcellular","tissues","cell_types","pathways","ko_phenotype","ko_severity",
    "neural_phenotype","cerebellar_evidence","cancer_types","cancer_role",
    "neuro_disease","differentiation_role","one_line_summary")
  data_rows <- parsed[-1]
  fixed <- lapply(data_rows, function(r) {
    if (length(r) == 17) return(r)
    if (length(r) > 17) return(c(r[1:16], paste(r[17:length(r)], collapse=" ")))
    if (length(r) < 17) return(c(r, rep(NA, 17-length(r))))
  })
  df <- as.data.frame(do.call(rbind, fixed), stringsAsFactors=FALSE)
  colnames(df) <- expected_cols
  df
}

reprocess_files <- list.files("agent_batches", pattern="reprocess_batch_[0-9]+_results\\.tsv$", full.names=TRUE)
cat("Reprocess files:", length(reprocess_files), "\n")

reprocess <- do.call(rbind, lapply(reprocess_files, parse_batch_lenient))
cat("Reprocess rows:", nrow(reprocess), "\n")
cat("Genes updated:\n  ", paste(reprocess$gene, collapse=", "), "\n")

# Replace old entries with reprocess entries
ann_updated <- ann %>% dplyr::filter(!gene %in% reprocess$gene)
ann_updated <- rbind(ann_updated, reprocess[, colnames(ann_updated)[colnames(ann_updated) %in% colnames(reprocess)]])

# Check structure
cat("\nFinal annotations:", nrow(ann_updated), "\n")

# Count filled fields
fields <- c("protein_class","molecular_action","key_substrates","subcellular","tissues",
  "cell_types","pathways","ko_phenotype","ko_severity","neural_phenotype",
  "cerebellar_evidence","cancer_types","cancer_role","neuro_disease",
  "differentiation_role","one_line_summary")
ann_updated$n_filled <- apply(ann_updated[, fields], 1, function(r) {
  sum(!is.na(r) & r != "" & r != "NA" & r != "unknown" & r != "none" & r != "none known" &
      r != "other" & r != "no")
})

cat("\nCoverage after reprocess:\n")
cat("  Median filled fields:", median(ann_updated$n_filled), "\n")
cat("  Mean filled fields:", round(mean(ann_updated$n_filled), 1), "\n")
for (field in fields) {
  n <- sum(!is.na(ann_updated[[field]]) & ann_updated[[field]] != "" & ann_updated[[field]] != "NA" &
    !ann_updated[[field]] %in% c("unknown","none","none known","other","no"))
  cat(sprintf("  %-25s %3d / %d (%2.0f%%)\n", field, n, nrow(ann_updated), 100*n/nrow(ann_updated)))
}

# Save updated
write.csv(ann_updated %>% dplyr::select(-n_filled), "agent_extracted_annotations.csv", row.names=FALSE)
cat("\nSaved: agent_extracted_annotations.csv (updated)\n")

# =============================================================================
# REBUILD EXCEL with updated annotations
# =============================================================================
cat("\nRebuilding Excel...\n")

# Load other annotation sources
master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude))

base_ann <- read.csv("diff_gene_full_annotations.csv", stringsAsFactors=FALSE)
umap_df <- read.csv("gene_umap_k10_coordinates.csv", stringsAsFactors=FALSE)
tiered <- read.csv("diff_gene_tiered_abstract_summaries.csv", stringsAsFactors=FALSE)
uniprot <- read.csv("diff_gene_uniprot_functions.csv", stringsAsFactors=FALSE)

combined <- base_ann %>%
  dplyr::left_join(umap_df %>% dplyr::select(gene, cluster) %>% dplyr::rename(cluster_k10=cluster),
    by=c("mgi_symbol"="gene")) %>%
  dplyr::mutate(cluster = ifelse(is.na(cluster), cluster_k10, cluster)) %>%
  dplyr::left_join(uniprot %>% dplyr::select(mgi_symbol, uniprot_protein, uniprot_function), by="mgi_symbol") %>%
  dplyr::left_join(tiered %>% dplyr::select(symbol, n_papers, n_tier1), by=c("mgi_symbol"="symbol")) %>%
  dplyr::left_join(ann_updated %>% dplyr::select(-n_filled), by=c("mgi_symbol"="gene"))

wb <- createWorkbook()
hs <- createStyle(textDecoration="bold", fgFill="#4472C4", fontColour="white",
  halign="center", border="Bottom", wrapText=TRUE)
ws <- createStyle(wrapText=TRUE, valign="top")

# --- Sheet 1: Gene annotations ---
addWorksheet(wb, "Gene Annotations")

s1 <- combined %>%
  dplyr::select(
    mgi_symbol, GENENAME, protein_class, molecular_action, one_line_summary,
    uniprot_function, ensembl_gene_id, cluster, t50_cat, H3K27me3_bin, H3K27me3,
    key_substrates, subcellular, tissues, cell_types, pathways,
    ko_phenotype, ko_severity, neural_phenotype, cerebellar_evidence,
    cancer_types, cancer_role, neuro_disease, differentiation_role,
    n_pubmed, n_papers, n_tier1,
    n_phenotypes, has_neural, has_lethality, phenotype_terms
  ) %>%
  dplyr::arrange(cluster, mgi_symbol)

colnames(s1) <- c(
  "Gene","Gene Name","Protein Class","Molecular Action","One-Line Summary",
  "UniProt Function","Ensembl ID","Cluster","Timing","H3K27me3","H3K27me3 Level",
  "Substrates/Partners","Subcellular","Tissues","Cell Types","Pathways",
  "KO Phenotype","KO Severity","Neural Phenotype","Cerebellar Evidence",
  "Cancer Types","Cancer Role","Neuro Disease","Differentiation Role",
  "PubMed Refs","Tiered Papers","Tier1 Papers",
  "MGI Phenotypes","Neural Pheno (MGI)","Lethality (MGI)","MGI Phenotype Terms"
)

writeData(wb, "Gene Annotations", s1)
addStyle(wb, "Gene Annotations", hs, rows=1, cols=1:ncol(s1))
addStyle(wb, "Gene Annotations", ws, rows=2:(nrow(s1)+1),
  cols=c(3,4,5,6,12,13,14,15,16,17,19,20,21,23,31), gridExpand=TRUE, stack=TRUE)
setColWidths(wb, "Gene Annotations", cols=1:ncol(s1),
  widths=c(12,30,18,35,50,45,18,6,6,6,6,25,20,20,20,20,35,12,25,20,20,15,20,15,8,8,8,8,6,6,45))
freezePane(wb, "Gene Annotations", firstRow=TRUE, firstCol=TRUE)

saveWorkbook(wb, "Supplementary_Table_Diff_Genes.xlsx", overwrite=TRUE)
cat("Saved: Supplementary_Table_Diff_Genes.xlsx\n")

cat("\n=== DONE ===\n")
