suppressMessages({library(tidyverse); library(openxlsx)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. MERGE ALL BATCH RESULTS
# =============================================================================
cat("Merging batch results...\n")

batch_files <- list.files("agent_batches", pattern = "batch_[0-9]+_results\\.tsv$", full.names = TRUE)
cat("Batch files found:", length(batch_files), "\n")

all_results <- lapply(batch_files, function(f) {
  tryCatch(read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE,
    quote="", fill=TRUE, comment.char=""),
    error = function(e) { cat("ERROR reading", f, ":", e$message, "\n"); NULL })
})

# Remove NULLs
all_results <- all_results[!sapply(all_results, is.null)]
cat("Successfully read:", length(all_results), "batches\n")

# Standardize column names across batches
expected_cols <- c("gene","protein_class","molecular_action","key_substrates",
  "subcellular","tissues","cell_types","pathways","ko_phenotype","ko_severity",
  "neural_phenotype","cerebellar_evidence","cancer_types","cancer_role",
  "neuro_disease","differentiation_role","one_line_summary")

for (i in seq_along(all_results)) {
  # Ensure all expected columns exist
  missing <- setdiff(expected_cols, colnames(all_results[[i]]))
  if (length(missing) > 0) {
    for (m in missing) all_results[[i]][[m]] <- NA
    cat("  Batch", i, "— added missing cols:", paste(missing, collapse=", "), "\n")
  }
  # Keep only expected columns in expected order
  all_results[[i]] <- all_results[[i]][, expected_cols]
}

merged <- do.call(rbind, all_results)
cat("Total genes merged:", nrow(merged), "\n")
cat("Unique genes:", n_distinct(merged$gene), "\n")

# Remove exact duplicates
merged <- merged[!duplicated(merged$gene), ]
cat("After dedup:", nrow(merged), "\n")

# =============================================================================
# 2. QUALITY CHECK
# =============================================================================
cat("\n=== QUALITY CHECK ===\n")
for (col in expected_cols) {
  n_filled <- sum(!is.na(merged[[col]]) & merged[[col]] != "" & merged[[col]] != "NA")
  cat(sprintf("  %-25s %3d / %d (%2.0f%%)\n", col, n_filled, nrow(merged), 100*n_filled/nrow(merged)))
}

cat("\n=== PROTEIN CLASS DISTRIBUTION ===\n")
pc <- sort(table(merged$protein_class), decreasing=TRUE)
for (cls in names(pc)) cat(sprintf("  %-25s %3d\n", cls, pc[cls]))

cat("\n=== CANCER ROLE ===\n")
print(table(merged$cancer_role, useNA="ifany"))

cat("\n=== KO SEVERITY ===\n")
print(table(merged$ko_severity, useNA="ifany"))

cat("\n=== CEREBELLAR EVIDENCE ===\n")
ce <- merged$cerebellar_evidence
ce_clean <- ifelse(is.na(ce) | ce == "NA" | ce == "" | ce == "no" | ce == "No", "no", "yes")
cat("  Yes:", sum(ce_clean=="yes"), " No:", sum(ce_clean=="no"), "\n")

# =============================================================================
# 3. MERGE WITH EXISTING ANNOTATIONS
# =============================================================================
cat("\nMerging with existing annotation data...\n")

ann <- read.csv("diff_gene_full_annotations.csv", stringsAsFactors=FALSE)
umap_df <- read.csv("gene_umap_k10_coordinates.csv", stringsAsFactors=FALSE)
tiered <- read.csv("diff_gene_tiered_abstract_summaries.csv", stringsAsFactors=FALSE)
uniprot <- read.csv("diff_gene_uniprot_functions.csv", stringsAsFactors=FALSE)

ann <- ann %>%
  dplyr::left_join(umap_df %>% dplyr::select(gene, cluster) %>% dplyr::rename(cluster_k10=cluster),
    by=c("mgi_symbol"="gene")) %>%
  dplyr::mutate(cluster = ifelse(is.na(cluster), cluster_k10, cluster)) %>%
  dplyr::left_join(uniprot %>% dplyr::select(mgi_symbol, uniprot_protein, uniprot_function), by="mgi_symbol") %>%
  dplyr::left_join(tiered %>% dplyr::select(symbol, n_papers, n_tier1, n_tier2, functional_keywords),
    by=c("mgi_symbol"="symbol")) %>%
  dplyr::left_join(merged, by=c("mgi_symbol"="gene"))

cat("Final table:", nrow(ann), "genes x", ncol(ann), "columns\n")

# =============================================================================
# 4. BUILD EXCEL
# =============================================================================
cat("\nBuilding Excel...\n")

wb <- createWorkbook()
hs <- createStyle(textDecoration="bold", fgFill="#4472C4", fontColour="white",
  halign="center", border="Bottom", wrapText=TRUE)
ws <- createStyle(wrapText=TRUE, valign="top")

# --- Sheet 1: Full gene table ---
addWorksheet(wb, "Gene Annotations")

s1 <- ann %>%
  dplyr::select(
    mgi_symbol, GENENAME, protein_class, molecular_action, one_line_summary,
    uniprot_function,
    ensembl_gene_id, cluster, t50_cat, H3K27me3_bin, H3K27me3,
    key_substrates, subcellular, tissues, cell_types, pathways,
    ko_phenotype, ko_severity, neural_phenotype, cerebellar_evidence,
    cancer_types, cancer_role, neuro_disease, differentiation_role,
    n_pubmed, n_papers, n_tier1,
    n_phenotypes, has_neural, has_lethality,
    functional_keywords, phenotype_terms
  ) %>%
  dplyr::arrange(cluster, mgi_symbol)

colnames(s1) <- c(
  "Gene","Gene Name","Protein Class","Molecular Action","One-Line Summary",
  "UniProt Function",
  "Ensembl ID","Cluster","Timing","H3K27me3","H3K27me3 Level",
  "Substrates/Partners","Subcellular","Tissues","Cell Types","Pathways",
  "KO Phenotype","KO Severity","Neural Phenotype","Cerebellar Evidence",
  "Cancer Types","Cancer Role","Neuro Disease","Differentiation Role",
  "PubMed Refs","Tiered Papers","Tier1 Papers",
  "MGI Phenotypes","Neural Pheno (MGI)","Lethality (MGI)",
  "Abstract Keywords","MGI Phenotype Terms"
)

writeData(wb, "Gene Annotations", s1)
addStyle(wb, "Gene Annotations", hs, rows=1, cols=1:ncol(s1))
addStyle(wb, "Gene Annotations", ws, rows=2:(nrow(s1)+1),
  cols=c(3,4,5,6,12,13,14,15,16,17,19,20,21,23,31,32), gridExpand=TRUE, stack=TRUE)
setColWidths(wb, "Gene Annotations", cols=1:ncol(s1),
  widths=c(12,30,18,35,50,45,18,6,6,6,6,25,20,20,20,20,35,12,25,20,20,15,20,15,8,8,8,8,6,6,25,45))
freezePane(wb, "Gene Annotations", firstRow=TRUE, firstCol=TRUE)

# --- Sheet 2: Cluster Summaries ---
addWorksheet(wb, "Cluster Summaries")

cs <- data.frame(
  Cluster = 0:10,
  N = sapply(0:10, function(cl) sum(ann$cluster==cl, na.rm=TRUE)),
  Pct_K27 = sapply(0:10, function(cl) {
    sub <- ann %>% dplyr::filter(cluster==cl)
    if(nrow(sub)==0) return(NA)
    round(100*mean(sub$H3K27me3_bin==1, na.rm=TRUE))
  }),
  Top_classes = sapply(0:10, function(cl) {
    sub <- ann %>% dplyr::filter(cluster==cl, !is.na(protein_class), protein_class!="NA")
    if(nrow(sub)==0) return("NA")
    tc <- sort(table(sub$protein_class), decreasing=TRUE)
    paste(paste0(names(tc)[1:min(3,length(tc))], " (", tc[1:min(3,length(tc))], ")"), collapse="; ")
  }),
  Pct_cerebellar = sapply(0:10, function(cl) {
    sub <- ann %>% dplyr::filter(cluster==cl)
    if(nrow(sub)==0) return(NA)
    ce <- sub$cerebellar_evidence
    round(100*mean(!is.na(ce) & ce != "NA" & ce != "" & ce != "no" & ce != "No", na.rm=TRUE))
  }),
  Pct_neural_pheno = sapply(0:10, function(cl) {
    sub <- ann %>% dplyr::filter(cluster==cl)
    if(nrow(sub)==0) return(NA)
    round(100*mean(sub$has_neural, na.rm=TRUE))
  }),
  stringsAsFactors=FALSE
)

writeData(wb, "Cluster Summaries", cs)
addStyle(wb, "Cluster Summaries", hs, rows=1, cols=1:ncol(cs))
setColWidths(wb, "Cluster Summaries", cols=1:ncol(cs), widths=c(8,6,8,50,12,12))

# --- Sheet 3: Summary Stats ---
addWorksheet(wb, "Summary Stats")

stats <- data.frame(
  Metric = c(
    "Total genes", "",
    "--- Agent Extraction Coverage ---",
    "Protein class assigned", "Molecular action described", "One-line summary",
    "KO phenotype", "Neural phenotype", "Cerebellar evidence",
    "Cancer types identified", "Cancer role assigned",
    "Neurological disease", "Differentiation role", "",
    "--- Protein Class Breakdown ---",
    names(pc)),
  Value = c(
    nrow(merged), NA, NA,
    sum(!is.na(merged$protein_class) & merged$protein_class != "NA" & merged$protein_class != "other"),
    sum(!is.na(merged$molecular_action) & merged$molecular_action != "NA"),
    sum(!is.na(merged$one_line_summary) & merged$one_line_summary != "NA"),
    sum(!is.na(merged$ko_phenotype) & merged$ko_phenotype != "NA"),
    sum(!is.na(merged$neural_phenotype) & merged$neural_phenotype != "NA"),
    sum(ce_clean == "yes"),
    sum(!is.na(merged$cancer_types) & merged$cancer_types != "NA" & merged$cancer_types != "none"),
    sum(!is.na(merged$cancer_role) & !merged$cancer_role %in% c("NA","none","unknown")),
    sum(!is.na(merged$neuro_disease) & merged$neuro_disease != "NA"),
    sum(!is.na(merged$differentiation_role) & !merged$differentiation_role %in% c("NA","none known","none")),
    NA, NA, as.integer(pc)),
  Pct = c(
    "100%", "", "",
    paste0(round(100*sum(!is.na(merged$protein_class) & merged$protein_class != "NA" & merged$protein_class != "other")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$molecular_action) & merged$molecular_action != "NA")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$one_line_summary) & merged$one_line_summary != "NA")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$ko_phenotype) & merged$ko_phenotype != "NA")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$neural_phenotype) & merged$neural_phenotype != "NA")/nrow(merged),1),"%"),
    paste0(round(100*sum(ce_clean=="yes")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$cancer_types) & merged$cancer_types != "NA" & merged$cancer_types != "none")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$cancer_role) & !merged$cancer_role %in% c("NA","none","unknown"))/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$neuro_disease) & merged$neuro_disease != "NA")/nrow(merged),1),"%"),
    paste0(round(100*sum(!is.na(merged$differentiation_role) & !merged$differentiation_role %in% c("NA","none known","none"))/nrow(merged),1),"%"),
    "", "", paste0(round(100*as.integer(pc)/nrow(merged),1),"%")),
  stringsAsFactors=FALSE
)

writeData(wb, "Summary Stats", stats)
addStyle(wb, "Summary Stats", hs, rows=1, cols=1:3)
setColWidths(wb, "Summary Stats", cols=1:3, widths=c(35,12,12))

# Save
saveWorkbook(wb, "Supplementary_Table_Diff_Genes.xlsx", overwrite=TRUE)
cat("\nSaved: Supplementary_Table_Diff_Genes.xlsx\n")
cat("  Sheet 1: Gene Annotations (", nrow(s1), "genes x", ncol(s1), "cols)\n")
cat("  Sheet 2: Cluster Summaries\n")
cat("  Sheet 3: Summary Stats\n")

# Also save merged TSV for reference
write.csv(merged, "agent_extracted_annotations.csv", row.names=FALSE)
cat("Saved: agent_extracted_annotations.csv\n")
