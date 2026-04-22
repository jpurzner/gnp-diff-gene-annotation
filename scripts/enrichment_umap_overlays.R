suppressMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(patchwork)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

umap_df <- read.csv("gene_umap_v3_coordinates.csv", stringsAsFactors=FALSE)
cat("UMAP genes:", nrow(umap_df), "\n")

master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude))

gene_ids <- bitr(diff_genes$ensembl_gene_id, fromType="ENSEMBL",
  toType=c("ENTREZID","SYMBOL"), OrgDb=org.Mm.eg.db)
gene_ids <- gene_ids[!duplicated(gene_ids$ENSEMBL),]
entrez <- unique(gene_ids$ENTREZID)

# Enrichment
cat("Running enrichment...\n")
ego_bp <- enrichGO(gene=entrez, OrgDb=org.Mm.eg.db, ont="BP",
  pAdjustMethod="BH", pvalueCutoff=0.01, readable=TRUE, minGSSize=10, maxGSSize=500)
ego_mf <- enrichGO(gene=entrez, OrgDb=org.Mm.eg.db, ont="MF",
  pAdjustMethod="BH", pvalueCutoff=0.01, readable=TRUE, minGSSize=10, maxGSSize=500)
ekegg <- enrichKEGG(gene=entrez, organism="mmu", pAdjustMethod="BH",
  pvalueCutoff=0.05, minGSSize=10, maxGSSize=500)

ego_bp_simp <- simplify(ego_bp, cutoff=0.6)
ego_mf_simp <- simplify(ego_mf, cutoff=0.6)

bp_top <- ego_bp_simp@result %>% filter(p.adjust<0.01) %>% arrange(p.adjust) %>%
  slice_head(n=10) %>% mutate(source="GO:BP")
mf_top <- ego_mf_simp@result %>% filter(p.adjust<0.01) %>% arrange(p.adjust) %>%
  slice_head(n=6) %>% mutate(source="GO:MF")
kegg_top <- ekegg@result %>% filter(p.adjust<0.05) %>% arrange(p.adjust) %>%
  slice_head(n=4) %>% mutate(source="KEGG")

# Convert KEGG entrez to symbols
entrez_to_sym <- setNames(gene_ids$SYMBOL, gene_ids$ENTREZID)
kegg_top$geneID <- sapply(kegg_top$geneID, function(x) {
  ids <- strsplit(x, "/")[[1]]
  syms <- entrez_to_sym[ids]
  paste(syms[!is.na(syms)], collapse="/")
})

# Clean KEGG descriptions
kegg_top$Description <- gsub(" - Mus musculus.*$", "", kegg_top$Description)

all_top <- rbind(
  bp_top %>% dplyr::select(ID, Description, p.adjust, Count, geneID, source),
  mf_top %>% dplyr::select(ID, Description, p.adjust, Count, geneID, source),
  kegg_top %>% dplyr::select(ID, Description, p.adjust, Count, geneID, source)
)

cat("Total overlay terms:", nrow(all_top), "\n")
for (i in 1:nrow(all_top)) {
  cat(sprintf("  [%s] %s — %d genes, FDR=%.1e\n",
    all_top$source[i], all_top$Description[i], all_top$Count[i], all_top$p.adjust[i]))
}

# =============================================
# UMAP overlays
# =============================================
cat("\nGenerating UMAP overlays...\n")

overlay_plots <- list()
for (i in 1:nrow(all_top)) {
  term_desc <- all_top$Description[i]
  term_genes_raw <- strsplit(all_top$geneID[i], "/")[[1]]
  in_term <- umap_df$gene[toupper(umap_df$gene) %in% toupper(term_genes_raw)]
  n_in <- length(in_term)

  umap_tmp <- umap_df %>% mutate(highlight = ifelse(gene %in% in_term, "Yes", "No"))

  desc_short <- ifelse(nchar(term_desc) > 42,
    paste0(substr(term_desc, 1, 39), "..."), term_desc)

  p <- ggplot(umap_tmp, aes(x=UMAP1, y=UMAP2)) +
    geom_point(data=umap_tmp %>% filter(highlight=="No"),
      color="grey88", size=0.5, alpha=0.3) +
    geom_point(data=umap_tmp %>% filter(highlight=="Yes"),
      color="#D6604D", size=2, alpha=0.75) +
    labs(title=desc_short,
      subtitle=paste0(n_in, " genes | ", all_top$source[i],
        " | FDR=", format(all_top$p.adjust[i], digits=2))) +
    theme_void(base_size=10) +
    theme(plot.title=element_text(size=10, face="bold", hjust=0.5),
      plot.subtitle=element_text(size=8, hjust=0.5, color="grey40"),
      plot.margin=margin(3,3,3,3,"pt"))

  overlay_plots[[i]] <- p
}

n_plots <- length(overlay_plots)
n_cols <- 5
n_rows <- ceiling(n_plots / n_cols)

p_grid <- wrap_plots(overlay_plots, ncol=n_cols) +
  plot_annotation(
    title="Enriched pathways in GNP differentiation genes mapped onto functional UMAP",
    subtitle="Red = member genes of each enriched term (clusterProfiler ORA, simplified, FDR-ranked)",
    theme=theme(
      plot.title=element_text(size=16, face="bold", hjust=0.5),
      plot.subtitle=element_text(size=12, hjust=0.5, color="grey40")))

ggsave("Fig_gene_umap_v3_term_overlays.pdf", p_grid, width=22, height=n_rows*4.5)
cat("Saved: Fig_gene_umap_v3_term_overlays.pdf\n")

# Dotplots
pdf("Fig_enrichment_dotplot_BP.pdf", width=10, height=10)
print(dotplot(ego_bp_simp, showCategory=20, title="GO Biological Process") +
  theme(axis.text.y=element_text(size=10)))
dev.off()
cat("Saved: Fig_enrichment_dotplot_BP.pdf\n")

pdf("Fig_enrichment_dotplot_MF.pdf", width=10, height=8)
print(dotplot(ego_mf_simp, showCategory=15, title="GO Molecular Function") +
  theme(axis.text.y=element_text(size=10)))
dev.off()
cat("Saved: Fig_enrichment_dotplot_MF.pdf\n")

write.csv(all_top, "enrichment_top_terms.csv", row.names=FALSE)
cat("\n=== DONE ===\n")
