# Usage: Rscript run_merge_and_excel.R PROLIF
#    or: Rscript run_merge_and_excel.R EMB
suppressMessages({library(tidyverse); library(openxlsx)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

args <- commandArgs(trailingOnly=TRUE)
SET <- args[1]  # "PROLIF" or "EMB"
if (!SET %in% c("PROLIF", "EMB")) stop("Set must be PROLIF or EMB")

cat("=== MERGE + EXCEL FOR", SET, "===\n")

# Set-specific paths
if (SET == "PROLIF") {
  gene_list <- read.csv("gene_list_proliferation.csv", stringsAsFactors=FALSE)
  batch_dir <- "agent_batches_prolif"
  pmid_assignments <- readRDS("tiered_pmid_assignments_PROLIF.rds")
  abstracts_raw <- readRDS("tiered_abstracts_raw_PROLIF.rds")
  compiled <- readRDS("tiered_compiled_PROLIF.rds")
  out_xlsx <- "Supplementary_Table_Proliferation_Genes.xlsx"
  out_csv <- "agent_extracted_annotations_PROLIF.csv"
  set_label <- "Proliferation"
  expected_n <- 438
} else {
  gene_list <- read.csv("gene_list_embryonic.csv", stringsAsFactors=FALSE)
  batch_dir <- "agent_batches_emb"
  pmid_assignments <- readRDS("tiered_pmid_assignments_EMB.rds")
  abstracts_raw <- readRDS("tiered_abstracts_raw_EMB.rds")
  compiled <- readRDS("tiered_compiled_EMB.rds")
  out_xlsx <- "Supplementary_Table_Embryonic_Genes.xlsx"
  out_csv <- "agent_extracted_annotations_EMB.csv"
  set_label <- "Embryonic"
  expected_n <- 782
}

# =============================================================================
# 1. PARSE BATCH RESULTS WITH LENIENT TSV PARSING
# =============================================================================
batch_files <- list.files(batch_dir, pattern="batch_[0-9]+_results\\.tsv$", full.names=TRUE)
cat("Batch files:", length(batch_files), "\n")

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

all_results <- lapply(batch_files, parse_batch_lenient)
all_results <- all_results[!sapply(all_results, is.null)]
merged <- do.call(rbind, all_results)
merged <- merged[!duplicated(merged$gene) & !is.na(merged$gene) & merged$gene != "", ]
cat("Merged genes:", nrow(merged), "/", expected_n, "\n")

write.csv(merged, out_csv, row.names=FALSE)

# =============================================================================
# 2. BUILD CITATIONS PER GENE
# =============================================================================
build_citations <- function(sym) {
  pdf <- pmid_assignments[[sym]]
  if (is.null(pdf) || nrow(pdf) == 0) {
    return(list(pmids_all="", urls_all="", tier1_pmids="",
      tier1_titles="", n_pmids=0, n_tier1=0))
  }
  all_pmids <- pdf$pmid
  tier1_pmids <- pdf$pmid[pdf$tier == grep("Tier1", unique(pdf$tier), value=TRUE)[1]]
  urls_all <- paste0("https://pubmed.ncbi.nlm.nih.gov/", all_pmids, "/")
  top3_titles <- character()
  for (pmid in head(tier1_pmids, 3)) {
    art <- abstracts_raw[[pmid]]
    if (!is.null(art) && !is.null(art$title)) {
      year <- ifelse(is.null(art$year) || is.na(art$year), "", paste0(" (", art$year, ")"))
      top3_titles <- c(top3_titles, paste0("[PMID:", pmid, "] ", art$title, year))
    }
  }
  list(pmids_all=paste(all_pmids, collapse="; "),
    urls_all=paste(urls_all, collapse="; "),
    tier1_pmids=paste(tier1_pmids, collapse="; "),
    tier1_titles=paste(top3_titles, collapse=" | "),
    n_pmids=length(all_pmids), n_tier1=length(tier1_pmids))
}
cite_list <- lapply(merged$gene, build_citations)
cite_df <- data.frame(
  gene = merged$gene,
  n_pmids_total = sapply(cite_list, function(x) x$n_pmids),
  n_pmids_tier1 = sapply(cite_list, function(x) x$n_tier1),
  tier1_pmids = sapply(cite_list, function(x) x$tier1_pmids),
  pubmed_urls = sapply(cite_list, function(x) x$urls_all),
  tier1_top_titles = sapply(cite_list, function(x) x$tier1_titles),
  stringsAsFactors=FALSE)

# =============================================================================
# 3. UMAP CLUSTER (k=10) using agent features
# =============================================================================
cat("\nClustering with agent features...\n")
suppressMessages({library(uwot); library(irlba); library(cluster)})

# Build a feature matrix from agent columns
ann <- merged
agent <- merged
pc_levels <- c("ion channel","transporter","GPCR","receptor","kinase","phosphatase",
  "transcription factor","enzyme","structural protein","adhesion molecule",
  "scaffold/adaptor","secreted factor","GTPase/regulator","ubiquitin ligase",
  "protease","motor protein","channel regulator","ECM protein","other")
agent$protein_class[is.na(agent$protein_class) | !agent$protein_class %in% pc_levels] <- "other"
pc_mat <- model.matrix(~ 0 + factor(agent$protein_class, levels=pc_levels))
colnames(pc_mat) <- paste0("PC:", pc_levels)
rownames(pc_mat) <- agent$gene

multi_label_mat <- function(col_data, term_list, prefix) {
  m <- matrix(0, nrow=length(col_data), ncol=length(term_list))
  colnames(m) <- paste0(prefix, term_list)
  for (i in seq_along(col_data)) {
    if (!is.na(col_data[i]) && col_data[i] != "NA") {
      for (j in seq_along(term_list)) {
        if (grepl(term_list[j], col_data[i], ignore.case=TRUE)) m[i,j] <- 1
      }
    }
  }
  m
}
subcel_mat <- multi_label_mat(agent$subcellular,
  c("plasma membrane","synapse","nucleus","cytoplasm","mitochondria","ER","Golgi",
    "axon","dendrite","growth cone","extracellular"), "Loc:")
tissue_mat <- multi_label_mat(agent$tissues,
  c("brain","cerebellum","cortex","hippocampus","spinal cord","retina","heart",
    "muscle","kidney","liver","lung","bone","pancreas","intestine"), "Tiss:")
ct_mat <- multi_label_mat(agent$cell_types,
  c("granule neuron","Purkinje","interneuron","astrocyte","oligodendrocyte",
    "neuron","neural progenitor","microglia","endothelial","epithelial",
    "macrophage","cardiomyocyte","fibroblast","stem cell"), "Cell:")
pw_mat <- multi_label_mat(agent$pathways,
  c("Wnt","Hedgehog","Notch","TGF","MAPK","PI3K","cAMP","Calcium","NF-kB",
    "JAK/STAT","Hippo"), "PW:")
rownames(subcel_mat) <- rownames(tissue_mat) <- rownames(ct_mat) <- rownames(pw_mat) <- agent$gene

mat_combined <- cbind(pc_mat, subcel_mat, tissue_mat, ct_mat, pw_mat)
col_var <- apply(mat_combined, 2, var)
mat_combined <- mat_combined[, col_var > 0]

set.seed(42)
n_pcs <- min(50, ncol(mat_combined)-1, nrow(mat_combined)-1)
pca <- prcomp_irlba(mat_combined, n=n_pcs, center=TRUE, scale.=FALSE)
umap_res <- umap(pca$x, n_neighbors=15, min_dist=0.15, metric="cosine",
  n_components=2, n_epochs=500, spread=2.0)
km <- kmeans(umap_res, centers=10, nstart=50, iter.max=300)
sil <- round(mean(silhouette(km$cluster, dist(umap_res))[,3]), 3)
cat("Silhouette:", sil, "\n")

umap_df <- data.frame(gene=rownames(mat_combined),
  UMAP1=umap_res[,1], UMAP2=umap_res[,2], cluster=km$cluster, stringsAsFactors=FALSE)

# =============================================================================
# 4. MERGE EVERYTHING + UMAP coords + citations + master annotations
# =============================================================================
gene_meta <- gene_list %>%
  dplyr::select(mgi_symbol, hgnc_symbol, description,
    max_cat, t50_cat, H3K27me3_bin, H3K27me3, H3K4me3_bin)

combined <- gene_meta %>%
  dplyr::left_join(umap_df %>% dplyr::select(gene, cluster, UMAP1, UMAP2),
    by=c("mgi_symbol"="gene")) %>%
  dplyr::left_join(merged, by=c("mgi_symbol"="gene")) %>%
  dplyr::left_join(cite_df, by=c("mgi_symbol"="gene")) %>%
  dplyr::left_join(compiled %>% dplyr::select(symbol, n_papers, n_tier1, n_tier2),
    by=c("mgi_symbol"="symbol"))

# Save merged data + UMAP coords for downstream figures
write.csv(combined, paste0("combined_annotations_", SET, ".csv"), row.names=FALSE)
write.csv(umap_df, paste0("umap_coords_", SET, ".csv"), row.names=FALSE)

# =============================================================================
# 5. EXCEL
# =============================================================================
wb <- createWorkbook()
hs <- createStyle(textDecoration="bold", fgFill="#4472C4", fontColour="white",
  halign="center", border="Bottom", wrapText=TRUE)
ws <- createStyle(wrapText=TRUE, valign="top")
hyperStyle <- createStyle(fontColour="#0563C1", textDecoration="underline",
  wrapText=TRUE, valign="top")

# Sheet 1: Gene Annotations
addWorksheet(wb, "Gene Annotations")
s1 <- combined %>%
  dplyr::select(mgi_symbol, hgnc_symbol, description, protein_class,
    molecular_action, one_line_summary,
    cluster, t50_cat, max_cat, H3K27me3_bin,
    n_pmids_total, n_pmids_tier1, tier1_top_titles, tier1_pmids, pubmed_urls,
    key_substrates, subcellular, tissues, cell_types, pathways,
    ko_phenotype, ko_severity, neural_phenotype, cerebellar_evidence,
    cancer_types, cancer_role, neuro_disease, differentiation_role) %>%
  dplyr::arrange(cluster, mgi_symbol)
colnames(s1) <- c("Gene","Human Ortholog","Description","Protein Class",
  "Molecular Action","One-Line Summary",
  "Cluster","Timing","max_cat","H3K27me3",
  "PubMed Refs (total)","Tier1 Refs","Top Tier1 Titles","Tier1 PMIDs","All PubMed URLs",
  "Substrates/Partners","Subcellular","Tissues","Cell Types","Pathways",
  "KO Phenotype","KO Severity","Neural Phenotype","Cerebellar Evidence",
  "Cancer Types","Cancer Role","Neuro Disease","Differentiation Role")
writeData(wb, "Gene Annotations", s1)
addStyle(wb, "Gene Annotations", hs, rows=1, cols=1:ncol(s1))
addStyle(wb, "Gene Annotations", ws, rows=2:(nrow(s1)+1),
  cols=c(3,4,5,6,13,14,15,16,17,18,19,20,21,23,24,28), gridExpand=TRUE, stack=TRUE)
addStyle(wb, "Gene Annotations", hyperStyle, rows=2:(nrow(s1)+1), cols=15,
  gridExpand=TRUE, stack=TRUE)
setColWidths(wb, "Gene Annotations", cols=1:ncol(s1),
  widths=c(12,12,30,18,35,50,8,8,8,8,10,10,60,25,50,20,20,20,20,20,30,12,25,20,20,15,20,15))
freezePane(wb, "Gene Annotations", firstRow=TRUE, firstCol=TRUE)
setRowHeights(wb, "Gene Annotations", rows=2:(nrow(s1)+1), heights=80)

# Sheet 2: Cluster Summaries
addWorksheet(wb, "Cluster Summaries")
cs <- data.frame(
  Cluster = 1:10,
  N = sapply(1:10, function(cl) sum(combined$cluster==cl, na.rm=TRUE)),
  Pct_K27 = sapply(1:10, function(cl) {
    sub <- combined %>% dplyr::filter(cluster==cl)
    if (nrow(sub)==0) return(NA)
    round(100*mean(sub$H3K27me3_bin==1, na.rm=TRUE))
  }),
  Top_classes = sapply(1:10, function(cl) {
    sub <- combined %>% dplyr::filter(cluster==cl, !is.na(protein_class), protein_class != "NA")
    if (nrow(sub)==0) return("NA")
    tc <- sort(table(sub$protein_class), decreasing=TRUE)
    paste(paste0(names(tc)[1:min(3,length(tc))], " (", tc[1:min(3,length(tc))], ")"), collapse="; ")
  }),
  stringsAsFactors=FALSE)
writeData(wb, "Cluster Summaries", cs)
addStyle(wb, "Cluster Summaries", hs, rows=1, cols=1:ncol(cs))
setColWidths(wb, "Cluster Summaries", cols=1:ncol(cs), widths=c(8,6,8,55))

# Sheet 3: PubMed References (focused)
addWorksheet(wb, "PubMed References")
ref_sheet <- combined %>%
  dplyr::select(mgi_symbol, description, cluster, n_pmids_total, n_pmids_tier1,
    tier1_top_titles, tier1_pmids, pubmed_urls) %>%
  dplyr::arrange(cluster, mgi_symbol)
colnames(ref_sheet) <- c("Gene","Description","Cluster","PMIDs (total)",
  "Tier1 PMIDs (count)","Top Tier1 Paper Titles","Tier1 PMID list","All PubMed URLs")
writeData(wb, "PubMed References", ref_sheet)
addStyle(wb, "PubMed References", hs, rows=1, cols=1:ncol(ref_sheet))
addStyle(wb, "PubMed References", ws, rows=2:(nrow(ref_sheet)+1),
  cols=c(2,6,7,8), gridExpand=TRUE, stack=TRUE)
addStyle(wb, "PubMed References", hyperStyle, rows=2:(nrow(ref_sheet)+1),
  cols=8, gridExpand=TRUE, stack=TRUE)
setColWidths(wb, "PubMed References", cols=1:ncol(ref_sheet),
  widths=c(12,30,8,10,10,70,35,65))
setRowHeights(wb, "PubMed References", rows=2:(nrow(ref_sheet)+1), heights=80)

saveWorkbook(wb, out_xlsx, overwrite=TRUE)
cat("Saved:", out_xlsx, "\n")
