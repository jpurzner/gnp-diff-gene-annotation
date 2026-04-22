suppressMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(patchwork)
  library(ggrepel)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# PAPER COLORS
# =============================================================================
col_ezh2_low  <- "orange2"
col_ezh2_high <- "navajowhite4"

# =============================================================================
# 1. LOAD DATA + DEFINE DIFF GENES
# =============================================================================
cat("Loading data...\n")

master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "\"")
colnames(master) <- gsub("\"", "", colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC",
    t50_cat %in% c("Mid", "Late"),
    !grepl("_exclude", max_cat_exclude))

cat("Differentiation genes:", nrow(diff_genes), "\n")
cat("  H3K27me3+:", sum(diff_genes$H3K27me3_bin == 1), "\n")
cat("  H3K27me3-:", sum(diff_genes$H3K27me3_bin == 0), "\n")

# =============================================================================
# 2. CONVERT TO ENTREZ IDs
# =============================================================================
cat("\nConverting gene IDs...\n")

gene_entrez <- bitr(diff_genes$ensembl_gene_id,
  fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Mapped to Entrez:", nrow(gene_entrez), "/", nrow(diff_genes), "\n")

# Also get symbol mapping for later
gene_symbols <- bitr(diff_genes$ensembl_gene_id,
  fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Mm.eg.db)

# Merge H3K27me3 status
gene_symbols <- gene_symbols %>%
  dplyr::left_join(diff_genes %>% dplyr::select(ensembl_gene_id, H3K27me3_bin, t50_cat),
    by = c("ENSEMBL" = "ensembl_gene_id"))

# Universe: all genes in master table
universe_entrez <- bitr(master$ensembl_gene_id,
  fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Universe genes:", nrow(universe_entrez), "\n")

# =============================================================================
# 3. GO ENRICHMENT (BP, MF, CC)
# =============================================================================
cat("\nRunning GO enrichment...\n")

ego_BP <- enrichGO(gene = gene_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
  readable = TRUE)

ego_MF <- enrichGO(gene = gene_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db, ont = "MF",
  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
  readable = TRUE)

ego_CC <- enrichGO(gene = gene_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db, ont = "CC",
  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
  readable = TRUE)

cat("BP terms:", nrow(ego_BP@result %>% dplyr::filter(p.adjust < 0.05)), "\n")
cat("MF terms:", nrow(ego_MF@result %>% dplyr::filter(p.adjust < 0.05)), "\n")
cat("CC terms:", nrow(ego_CC@result %>% dplyr::filter(p.adjust < 0.05)), "\n")

# Simplify to remove redundant terms
ego_BP_simp <- simplify(ego_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_MF_simp <- simplify(ego_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

cat("BP simplified:", nrow(ego_BP_simp@result %>% dplyr::filter(p.adjust < 0.05)), "\n")
cat("MF simplified:", nrow(ego_MF_simp@result %>% dplyr::filter(p.adjust < 0.05)), "\n")

# =============================================================================
# 4. KEGG PATHWAY ENRICHMENT
# =============================================================================
cat("\nRunning KEGG enrichment...\n")

kegg_res <- enrichKEGG(gene = gene_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  organism = "mmu", pvalueCutoff = 0.05)

cat("KEGG pathways:", nrow(kegg_res@result %>% dplyr::filter(p.adjust < 0.05)), "\n")

# Make readable
kegg_readable <- setReadable(kegg_res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# =============================================================================
# 5. REACTOME PATHWAY ENRICHMENT
# =============================================================================
cat("Running Reactome enrichment...\n")

reactome_res <- enrichPathway(gene = gene_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  organism = "mouse", pvalueCutoff = 0.05, readable = TRUE)

cat("Reactome pathways:", nrow(reactome_res@result %>% dplyr::filter(p.adjust < 0.05)), "\n")

# =============================================================================
# 6. FIGURE 1: DOTPLOT — TOP GO BP TERMS
# =============================================================================
cat("\nGenerating figures...\n")

p_bp <- dotplot(ego_BP_simp, showCategory = 20, font.size = 12,
  title = "GO Biological Process") +
  theme(axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 15, face = "bold"))

ggsave("Fig_pathway_GO_BP.pdf", p_bp, width = 10, height = 9)
cat("Saved: Fig_pathway_GO_BP.pdf\n")

# =============================================================================
# 7. FIGURE 2: DOTPLOT — TOP GO MF TERMS
# =============================================================================
p_mf <- dotplot(ego_MF_simp, showCategory = 20, font.size = 12,
  title = "GO Molecular Function") +
  theme(axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 15, face = "bold"))

ggsave("Fig_pathway_GO_MF.pdf", p_mf, width = 10, height = 9)
cat("Saved: Fig_pathway_GO_MF.pdf\n")

# =============================================================================
# 8. FIGURE 3: KEGG DOTPLOT
# =============================================================================
if (nrow(kegg_readable@result %>% dplyr::filter(p.adjust < 0.05)) > 0) {
  p_kegg <- dotplot(kegg_readable, showCategory = 15, font.size = 12,
    title = "KEGG Pathways") +
    theme(axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 15, face = "bold"))
  ggsave("Fig_pathway_KEGG.pdf", p_kegg, width = 10, height = 7)
  cat("Saved: Fig_pathway_KEGG.pdf\n")
}

# =============================================================================
# 9. FIGURE 4: REACTOME DOTPLOT
# =============================================================================
if (nrow(reactome_res@result %>% dplyr::filter(p.adjust < 0.05)) > 0) {
  p_react <- dotplot(reactome_res, showCategory = 15, font.size = 12,
    title = "Reactome Pathways") +
    theme(axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 15, face = "bold"))
  ggsave("Fig_pathway_Reactome.pdf", p_react, width = 10, height = 7)
  cat("Saved: Fig_pathway_Reactome.pdf\n")
}

# =============================================================================
# 10. FIGURE 5: CNETPLOT — BP (gene-concept network)
# =============================================================================
p_cnet <- cnetplot(ego_BP_simp, showCategory = 10,
  categorySize = "pvalue", foldChange = NULL,
  cex.params = list(category_label = 1.2, gene_label = 0.7)) +
  ggtitle("GO BP gene-concept network") +
  theme(plot.title = element_text(size = 15, face = "bold"))

ggsave("Fig_pathway_cnetplot_BP.pdf", p_cnet, width = 16, height = 14)
cat("Saved: Fig_pathway_cnetplot_BP.pdf\n")

# =============================================================================
# 11. STRATIFY BY H3K27me3 STATUS
# =============================================================================
cat("\nStratifying by H3K27me3 status...\n")

k27_pos_entrez <- gene_symbols %>%
  dplyr::filter(H3K27me3_bin == 1) %>% dplyr::pull(ENTREZID) %>% unique()
k27_neg_entrez <- gene_symbols %>%
  dplyr::filter(H3K27me3_bin == 0) %>% dplyr::pull(ENTREZID) %>% unique()

cat("H3K27me3+ genes with Entrez:", length(k27_pos_entrez), "\n")
cat("H3K27me3- genes with Entrez:", length(k27_neg_entrez), "\n")

# Compare GO enrichment: K27+ vs K27-
ego_k27pos <- enrichGO(gene = k27_pos_entrez,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
  readable = TRUE)

ego_k27neg <- enrichGO(gene = k27_neg_entrez,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
  readable = TRUE)

ego_k27pos_simp <- simplify(ego_k27pos, cutoff = 0.7)
ego_k27neg_simp <- simplify(ego_k27neg, cutoff = 0.7)

cat("K27+ BP terms:", nrow(ego_k27pos_simp@result %>% dplyr::filter(p.adjust < 0.05)), "\n")
cat("K27- BP terms:", nrow(ego_k27neg_simp@result %>% dplyr::filter(p.adjust < 0.05)), "\n")

# =============================================================================
# 12. FIGURE 6: SIDE-BY-SIDE K27+ vs K27- DOTPLOTS
# =============================================================================
p_k27pos <- dotplot(ego_k27pos_simp, showCategory = 15, font.size = 11,
  title = "H3K27me3-bound differentiation genes") +
  theme(axis.text.y = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold"))

p_k27neg <- dotplot(ego_k27neg_simp, showCategory = 15, font.size = 11,
  title = "H3K27me3-unbound differentiation genes") +
  theme(axis.text.y = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold"))

p_k27_compare <- p_k27pos + p_k27neg +
  plot_annotation(title = "GO Biological Process: H3K27me3-bound vs unbound differentiation genes",
    theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("Fig_pathway_K27_comparison.pdf", p_k27_compare, width = 20, height = 10)
cat("Saved: Fig_pathway_K27_comparison.pdf\n")

# =============================================================================
# 13. COMBINED SUMMARY FIGURE
# =============================================================================
cat("Generating combined summary figure...\n")

combined_pathway <- (p_bp + p_mf) / (p_k27pos + p_k27neg) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("Fig_pathway_combined.pdf", combined_pathway, width = 22, height = 20)
cat("Saved: Fig_pathway_combined.pdf\n")

# =============================================================================
# 14. PRINT TOP RESULTS SUMMARY
# =============================================================================
cat("\n=== TOP GO BIOLOGICAL PROCESS TERMS ===\n")
top_bp <- ego_BP_simp@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, Count, p.adjust, GeneRatio) %>%
  dplyr::slice_head(n = 20)
print(as.data.frame(top_bp), row.names = FALSE)

cat("\n=== TOP GO MOLECULAR FUNCTION TERMS ===\n")
top_mf <- ego_MF_simp@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, Count, p.adjust, GeneRatio) %>%
  dplyr::slice_head(n = 15)
print(as.data.frame(top_mf), row.names = FALSE)

cat("\n=== TOP KEGG PATHWAYS ===\n")
top_kegg <- kegg_readable@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, Count, p.adjust, GeneRatio) %>%
  dplyr::slice_head(n = 15)
print(as.data.frame(top_kegg), row.names = FALSE)

cat("\n=== TOP REACTOME PATHWAYS ===\n")
top_reactome <- reactome_res@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, Count, p.adjust, GeneRatio) %>%
  dplyr::slice_head(n = 15)
print(as.data.frame(top_reactome), row.names = FALSE)

cat("\n=== FUNCTIONAL CATEGORIES SUMMARY ===\n")

# Count genes in key functional categories from BP
key_terms_bp <- ego_BP@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(Description, Count, geneID)

# Synaptic/ion channel
synapse_genes <- key_terms_bp %>%
  dplyr::filter(grepl("synap|neurotransmitter|ion trans|channel|potassium|calcium|sodium", Description, ignore.case = TRUE)) %>%
  dplyr::pull(geneID) %>% paste(collapse = "/") %>% strsplit("/") %>% unlist() %>% unique()
cat("Synaptic signaling / ion transport genes:", length(synapse_genes), "\n")

# Morphogenesis/adhesion
morpho_genes <- key_terms_bp %>%
  dplyr::filter(grepl("morpho|adhesion|migration|axon|dendrit", Description, ignore.case = TRUE)) %>%
  dplyr::pull(geneID) %>% paste(collapse = "/") %>% strsplit("/") %>% unlist() %>% unique()
cat("Morphogenesis / cell adhesion genes:", length(morpho_genes), "\n")

# Signaling
signal_genes <- key_terms_bp %>%
  dplyr::filter(grepl("signal|cAMP|GPCR|G protein|Wnt|Notch|kinase", Description, ignore.case = TRUE)) %>%
  dplyr::pull(geneID) %>% paste(collapse = "/") %>% strsplit("/") %>% unlist() %>% unique()
cat("Signaling pathway genes:", length(signal_genes), "\n")

# Save results tables
write.table(ego_BP@result %>% dplyr::filter(p.adjust < 0.05),
  "pathway_GO_BP_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ego_MF@result %>% dplyr::filter(p.adjust < 0.05),
  "pathway_GO_MF_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(kegg_readable@result %>% dplyr::filter(p.adjust < 0.05),
  "pathway_KEGG_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(reactome_res@result %>% dplyr::filter(p.adjust < 0.05),
  "pathway_Reactome_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nSaved: pathway_GO_BP_results.tsv, pathway_GO_MF_results.tsv,")
cat("\n        pathway_KEGG_results.tsv, pathway_Reactome_results.tsv\n")

cat("\n=== ALL OUTPUTS GENERATED SUCCESSFULLY ===\n")
