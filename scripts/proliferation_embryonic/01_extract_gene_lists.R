suppressMessages({library(tidyverse); library(openxlsx)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

m <- read.table("gnp_MB_histone_k_toc_t50.txt", header=TRUE, sep="\t",
  stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(m) <- gsub("\"","",colnames(m))

# =============================================================================
# 1. EMBRYONIC GENES (783) — peak at E15 + dynamic + >2 log2FC
# =============================================================================
embryonic <- m %>%
  dplyr::filter(
    max_cat == "E15",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC",
    !grepl("_exclude", max_cat_exclude)) %>%
  dplyr::select(ensembl_gene_id, mgi_symbol, hg_ensembl_gene_id, hgnc_symbol,
    description, max_cat, t50_cat, group_1, group_2,
    avg_max, avg_min, avg_fc, avg_t50,
    H3K27me3_bin, H3K27me3, H3K4me3_bin, H3K4me3) %>%
  dplyr::distinct(mgi_symbol, .keep_all=TRUE) %>%
  dplyr::arrange(mgi_symbol)
cat("Embryonic (E15 + Dynamic + >2 log2FC):", nrow(embryonic), "\n")

# =============================================================================
# 2. PROLIFERATION GENES — Cell cycle + Early/Mid t50
# =============================================================================
proliferation <- m %>%
  dplyr::filter(
    group_1 == "Cell cycle",
    t50_cat %in% c("Early", "Mid"),
    !grepl("_exclude", max_cat_exclude)) %>%
  dplyr::select(ensembl_gene_id, mgi_symbol, hg_ensembl_gene_id, hgnc_symbol,
    description, max_cat, t50_cat, group_1, group_2,
    avg_max, avg_min, avg_fc, avg_t50,
    H3K27me3_bin, H3K27me3, H3K4me3_bin, H3K4me3) %>%
  dplyr::distinct(mgi_symbol, .keep_all=TRUE) %>%
  dplyr::arrange(mgi_symbol)
cat("Proliferation (Cell cycle + Early/Mid t50):", nrow(proliferation), "\n")

# Also try with cat_avg_fc != "0-1" filter (more conservative)
prolif_dyn <- m %>%
  dplyr::filter(
    group_1 == "Cell cycle",
    t50_cat %in% c("Early", "Mid"),
    cat_avg_fc != "0-1 log2FC",
    !grepl("_exclude", max_cat_exclude)) %>%
  dplyr::distinct(mgi_symbol, .keep_all=TRUE)
cat("Proliferation (+ avg_fc>=1):", nrow(prolif_dyn), "\n")

# =============================================================================
# 3. DIFFERENTIATION GENES (786) — for reference
# =============================================================================
diff_genes <- m %>%
  dplyr::filter(
    max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC",
    t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude)) %>%
  dplyr::distinct(mgi_symbol, .keep_all=TRUE)
cat("Differentiation (for ref):", nrow(diff_genes), "\n")

# =============================================================================
# 4. SAVE
# =============================================================================
write.csv(embryonic,    "gene_list_embryonic.csv",    row.names=FALSE)
write.csv(proliferation,"gene_list_proliferation.csv",row.names=FALSE)
cat("\nSaved: gene_list_embryonic.csv (", nrow(embryonic), "genes)\n")
cat("Saved: gene_list_proliferation.csv (", nrow(proliferation), "genes)\n")

# Combined Excel
wb <- createWorkbook()
hs <- createStyle(textDecoration="bold", fgFill="#4472C4", fontColour="white",
  halign="center", border="Bottom")
addWorksheet(wb, "Embryonic (n=783)")
writeData(wb, "Embryonic (n=783)", embryonic)
addStyle(wb, "Embryonic (n=783)", hs, rows=1, cols=1:ncol(embryonic))
setColWidths(wb, "Embryonic (n=783)", cols=1:ncol(embryonic),
  widths=c(20,12,20,12,40,8,8,40,15,8,8,8,8,6,8,6,8))
freezePane(wb, "Embryonic (n=783)", firstRow=TRUE)

addWorksheet(wb, paste0("Proliferation (n=", nrow(proliferation), ")"))
writeData(wb, paste0("Proliferation (n=", nrow(proliferation), ")"), proliferation)
addStyle(wb, paste0("Proliferation (n=", nrow(proliferation), ")"), hs,
  rows=1, cols=1:ncol(proliferation))
setColWidths(wb, paste0("Proliferation (n=", nrow(proliferation), ")"),
  cols=1:ncol(proliferation),
  widths=c(20,12,20,12,40,8,8,40,15,8,8,8,8,6,8,6,8))
freezePane(wb, paste0("Proliferation (n=", nrow(proliferation), ")"), firstRow=TRUE)

# Definitions sheet
defs <- data.frame(
  Group = c("Differentiation", "Embryonic", "Proliferation"),
  N = c(nrow(diff_genes), nrow(embryonic), nrow(proliferation)),
  Filter = c(
    "max_cat=='P14_P56' & group_1=='Dynamic GN lineage' & group_2=='>2 log2FC' & t50_cat in {Mid,Late} & no _exclude",
    "max_cat=='E15' & group_1=='Dynamic GN lineage' & group_2=='>2 log2FC' & no _exclude",
    "group_1=='Cell cycle' & t50_cat in {Early,Mid} & no _exclude"),
  stringsAsFactors=FALSE)
addWorksheet(wb, "Definitions")
writeData(wb, "Definitions", defs)
addStyle(wb, "Definitions", hs, rows=1, cols=1:ncol(defs))
setColWidths(wb, "Definitions", cols=1:ncol(defs), widths=c(18, 8, 100))

saveWorkbook(wb, "Gene_lists_emb_prolif_diff.xlsx", overwrite=TRUE)
cat("\nSaved: Gene_lists_emb_prolif_diff.xlsx (3 sheets)\n")
