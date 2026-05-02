suppressMessages({library(tidyverse); library(openxlsx)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. LOAD ALL DATA SOURCES
# =============================================================================
cat("Loading data...\n")

ann <- read.csv("agent_extracted_annotations.csv", stringsAsFactors=FALSE)
base_ann <- read.csv("diff_gene_full_annotations.csv", stringsAsFactors=FALSE)
umap_df <- read.csv("gene_umap_agent_k10_coordinates.csv", stringsAsFactors=FALSE)
uniprot <- read.csv("diff_gene_uniprot_functions.csv", stringsAsFactors=FALSE)
tiered_summary <- read.csv("diff_gene_tiered_abstract_summaries.csv", stringsAsFactors=FALSE)

# PubMed source data
pmid_per_gene <- readRDS("tiered_pmid_assignments.rds")
abstracts_raw <- readRDS("tiered_abstracts_raw.rds")

cat("Genes with PMID assignments:", length(pmid_per_gene), "\n")
cat("Total abstracts indexed:", length(abstracts_raw), "\n")

# =============================================================================
# 2. BUILD CITATION COLUMNS PER GENE
# =============================================================================
cat("\nBuilding citation columns...\n")

# Per gene: tier1 PMIDs (most relevant: cerebellum/MB/differentiation)
# Per gene: all PMIDs as a single semicolon-separated list of URLs
# Per gene: top 3 tier1 paper titles for quick context

build_citations <- function(sym) {
  pdf <- pmid_per_gene[[sym]]
  if (is.null(pdf) || nrow(pdf) == 0) {
    return(list(pmids_all="", urls_all="", tier1_pmids="",
      tier1_titles="", n_pmids=0, n_tier1=0))
  }

  all_pmids <- pdf$pmid
  tier1_pmids <- pdf$pmid[pdf$tier == "Tier1_cerebellum_MB_diff"]

  urls_all <- paste0("https://pubmed.ncbi.nlm.nih.gov/", all_pmids, "/")

  # Top 3 tier1 paper titles
  top3_titles <- character()
  for (pmid in head(tier1_pmids, 3)) {
    art <- abstracts_raw[[pmid]]
    if (!is.null(art) && !is.null(art$title)) {
      title <- art$title
      year <- ifelse(is.null(art$year) || is.na(art$year), "", paste0(" (", art$year, ")"))
      top3_titles <- c(top3_titles, paste0("[PMID:", pmid, "] ", title, year))
    }
  }

  list(
    pmids_all = paste(all_pmids, collapse="; "),
    urls_all = paste(urls_all, collapse="; "),
    tier1_pmids = paste(tier1_pmids, collapse="; "),
    tier1_titles = paste(top3_titles, collapse=" | "),
    n_pmids = length(all_pmids),
    n_tier1 = length(tier1_pmids)
  )
}

# Build for all genes in the agent table
cite_list <- lapply(ann$gene, build_citations)
cite_df <- data.frame(
  gene = ann$gene,
  n_pmids_total = sapply(cite_list, function(x) x$n_pmids),
  n_pmids_tier1 = sapply(cite_list, function(x) x$n_tier1),
  tier1_pmids = sapply(cite_list, function(x) x$tier1_pmids),
  pubmed_urls = sapply(cite_list, function(x) x$urls_all),
  tier1_top_titles = sapply(cite_list, function(x) x$tier1_titles),
  stringsAsFactors = FALSE
)

cat("Genes with at least 1 PMID:", sum(cite_df$n_pmids_total > 0), "\n")
cat("Genes with at least 1 tier1 PMID:", sum(cite_df$n_pmids_tier1 > 0), "\n")

# =============================================================================
# 3. MERGE EVERYTHING
# =============================================================================
combined <- base_ann %>%
  dplyr::left_join(umap_df %>% dplyr::select(gene, cluster) %>%
    dplyr::rename(cluster_k10=cluster), by=c("mgi_symbol"="gene")) %>%
  dplyr::mutate(cluster = ifelse(is.na(cluster), cluster_k10, cluster)) %>%
  dplyr::left_join(uniprot %>% dplyr::select(mgi_symbol, uniprot_protein, uniprot_function),
    by="mgi_symbol") %>%
  dplyr::left_join(tiered_summary %>% dplyr::select(symbol, n_papers, n_tier1, n_tier2),
    by=c("mgi_symbol"="symbol")) %>%
  dplyr::left_join(ann, by=c("mgi_symbol"="gene")) %>%
  dplyr::left_join(cite_df, by=c("mgi_symbol"="gene"))

# =============================================================================
# 4. BUILD EXCEL WITH HYPERLINKS
# =============================================================================
cat("\nBuilding Excel with hyperlinks...\n")
wb <- createWorkbook()
hs <- createStyle(textDecoration="bold", fgFill="#4472C4", fontColour="white",
  halign="center", border="Bottom", wrapText=TRUE)
ws <- createStyle(wrapText=TRUE, valign="top")
hyperStyle <- createStyle(fontColour="#0563C1", textDecoration="underline",
  wrapText=TRUE, valign="top")

# --- Sheet 1: Gene Annotations with citations ---
addWorksheet(wb, "Gene Annotations")

s1 <- combined %>%
  dplyr::select(
    mgi_symbol, GENENAME, protein_class, molecular_action, one_line_summary,
    uniprot_function,
    n_pmids_total, n_pmids_tier1, tier1_top_titles, tier1_pmids, pubmed_urls,
    ensembl_gene_id, cluster, t50_cat, H3K27me3_bin, H3K27me3,
    key_substrates, subcellular, tissues, cell_types, pathways,
    ko_phenotype, ko_severity, neural_phenotype, cerebellar_evidence,
    cancer_types, cancer_role, neuro_disease, differentiation_role,
    n_pubmed, n_papers, n_phenotypes, has_neural, has_lethality, phenotype_terms
  ) %>%
  dplyr::arrange(cluster, mgi_symbol)

colnames(s1) <- c(
  "Gene","Gene Name","Protein Class","Molecular Action","One-Line Summary",
  "UniProt Function",
  "PubMed Refs (total)","Tier1 Refs (cereb/MB/diff)",
  "Top Tier1 Paper Titles","Tier1 PMIDs","All PubMed URLs",
  "Ensembl ID","Cluster","Timing","H3K27me3","H3K27me3 Level",
  "Substrates/Partners","Subcellular","Tissues","Cell Types","Pathways",
  "KO Phenotype","KO Severity","Neural Phenotype","Cerebellar Evidence",
  "Cancer Types","Cancer Role","Neuro Disease","Differentiation Role",
  "PubMed Refs (org.Mm.eg.db)","Tiered Papers","MGI Phenotypes","Neural Pheno (MGI)","Lethality (MGI)","MGI Phenotype Terms"
)

writeData(wb, "Gene Annotations", s1)
addStyle(wb, "Gene Annotations", hs, rows=1, cols=1:ncol(s1))
addStyle(wb, "Gene Annotations", ws, rows=2:(nrow(s1)+1),
  cols=c(3,4,5,6,9,10,11,17,18,19,20,21,23,24,28,35), gridExpand=TRUE, stack=TRUE)

# Style URL columns differently (blue, underlined)
addStyle(wb, "Gene Annotations", hyperStyle,
  rows=2:(nrow(s1)+1), cols=c(11), gridExpand=TRUE, stack=TRUE)

setColWidths(wb, "Gene Annotations", cols=1:ncol(s1),
  widths=c(12, 30, 18, 35, 50, 45,
           10, 10, 60, 25, 50,
           18, 6, 6, 6, 6,
           25, 20, 20, 20, 20,
           35, 12, 25, 20, 20, 15, 20, 15,
           8, 8, 8, 6, 6, 45))
freezePane(wb, "Gene Annotations", firstRow=TRUE, firstCol=TRUE)
setRowHeights(wb, "Gene Annotations", rows=2:(nrow(s1)+1), heights=80)

# --- Sheet 2: Cluster summaries (carry over) ---
addWorksheet(wb, "Cluster Summaries")
cs <- data.frame(
  Cluster = 0:10,
  N = sapply(0:10, function(cl) sum(combined$cluster==cl, na.rm=TRUE)),
  Pct_K27 = sapply(0:10, function(cl) {
    sub <- combined %>% dplyr::filter(cluster==cl)
    if(nrow(sub)==0) return(NA)
    round(100*mean(sub$H3K27me3_bin==1, na.rm=TRUE))
  }),
  Top_classes = sapply(0:10, function(cl) {
    sub <- combined %>% dplyr::filter(cluster==cl, !is.na(protein_class), protein_class!="NA")
    if(nrow(sub)==0) return("NA")
    tc <- sort(table(sub$protein_class), decreasing=TRUE)
    paste(paste0(names(tc)[1:min(3,length(tc))], " (", tc[1:min(3,length(tc))], ")"), collapse="; ")
  }),
  Pct_with_MB_evidence = sapply(0:10, function(cl) {
    sub <- combined %>% dplyr::filter(cluster==cl)
    if(nrow(sub)==0) return(NA)
    round(100 * sum(grepl("medulloblastoma", sub$cancer_types, ignore.case=TRUE), na.rm=TRUE) / nrow(sub))
  }),
  stringsAsFactors=FALSE
)
writeData(wb, "Cluster Summaries", cs)
addStyle(wb, "Cluster Summaries", hs, rows=1, cols=1:ncol(cs))
setColWidths(wb, "Cluster Summaries", cols=1:ncol(cs), widths=c(8,6,8,55,15))

# --- Sheet 3: PubMed reference inventory ---
# One row per gene with full PMID list + tier breakdown
addWorksheet(wb, "PubMed References")

ref_sheet <- combined %>%
  dplyr::select(mgi_symbol, GENENAME, cluster,
    n_pmids_total, n_pmids_tier1,
    tier1_top_titles, tier1_pmids, pubmed_urls) %>%
  dplyr::arrange(cluster, mgi_symbol)

colnames(ref_sheet) <- c("Gene","Gene Name","Cluster",
  "PMIDs (total)","Tier1 PMIDs (count)",
  "Top Tier1 Paper Titles","Tier1 PMID list","All PubMed URLs")

writeData(wb, "PubMed References", ref_sheet)
addStyle(wb, "PubMed References", hs, rows=1, cols=1:ncol(ref_sheet))
addStyle(wb, "PubMed References", ws, rows=2:(nrow(ref_sheet)+1),
  cols=c(2,6,7,8), gridExpand=TRUE, stack=TRUE)
addStyle(wb, "PubMed References", hyperStyle,
  rows=2:(nrow(ref_sheet)+1), cols=8, gridExpand=TRUE, stack=TRUE)
setColWidths(wb, "PubMed References", cols=1:ncol(ref_sheet),
  widths=c(12, 30, 8, 10, 10, 70, 35, 65))
setRowHeights(wb, "PubMed References", rows=2:(nrow(ref_sheet)+1), heights=80)
freezePane(wb, "PubMed References", firstRow=TRUE)

# --- Sheet 4: Tier definitions ---
addWorksheet(wb, "Tier Definitions")
tier_def <- data.frame(
  Tier = c("Tier 1", "Tier 2", "Tier 3", "Tier 4", "Tier 5", "Tier 6"),
  Search_terms = c(
    "cerebellum, cerebellar, granule neuron, granule cell, medulloblastoma, differentiation",
    "brain, neuron, neural, CNS, cortex, hippocampus, synapse, synaptic, axon, dendrite",
    "brain tumor, glioma, glioblastoma, neuroblastoma, brain cancer, neuro-oncology",
    "cancer, tumor, oncogene, tumor suppressor, carcinoma, malignant, metastasis",
    "development, embryo, progenitor, stem cell, lineage, morphogenesis, patterning",
    "(no filter — any paper containing the gene name)"),
  Description = c(
    "Most relevant: directly studies cerebellar GNP biology or medulloblastoma",
    "Brain / neuron / synapse — relevant to neural function generally",
    "Brain cancer / glioma / neuroblastoma — relevant to CNS cancers",
    "General cancer literature",
    "Development & differentiation literature",
    "Catch-all if no higher-tier papers found"),
  stringsAsFactors=FALSE
)
writeData(wb, "Tier Definitions", tier_def)
addStyle(wb, "Tier Definitions", hs, rows=1, cols=1:ncol(tier_def))
setColWidths(wb, "Tier Definitions", cols=1:ncol(tier_def), widths=c(8, 70, 60))

saveWorkbook(wb, "Supplementary_Table_Diff_Genes.xlsx", overwrite=TRUE)
cat("\nSaved: Supplementary_Table_Diff_Genes.xlsx\n")
cat("  Sheet 1: Gene Annotations (", nrow(s1), "x", ncol(s1), ")\n")
cat("  Sheet 2: Cluster Summaries\n")
cat("  Sheet 3: PubMed References (", nrow(ref_sheet), "rows)\n")
cat("  Sheet 4: Tier Definitions\n")
