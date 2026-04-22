suppressMessages({
  library(tidyverse)
  library(patchwork)
  library(cluster)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data...\n")

umap_df <- read.csv("gene_umap_v3_coordinates.csv", stringsAsFactors=FALSE)
umap_ann <- umap_df %>% dplyr::filter(cluster > 0)
umap_mat <- as.matrix(umap_ann[, c("UMAP1","UMAP2")])
umap_unann <- umap_df %>% dplyr::filter(cluster == 0)
mat_ann <- readRDS("gene_term_matrix_v3.rds")

# Gene ID mapping
master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude))

gene_ids <- AnnotationDbi::select(org.Mm.eg.db,
  keys = diff_genes$ensembl_gene_id, keytype = "ENSEMBL",
  columns = c("ENTREZID","SYMBOL"))
gene_ids <- gene_ids[!duplicated(gene_ids$ENSEMBL),]
# Map symbol to ensembl for later
sym_to_entrez <- setNames(gene_ids$ENTREZID, gene_ids$SYMBOL)

# =============================================================================
# 2. K=10 CLUSTERING
# =============================================================================
cat("Clustering with k=10...\n")
set.seed(42)
km10 <- kmeans(umap_mat, centers = 10, nstart = 50, iter.max = 300)
sil_val <- round(mean(silhouette(km10$cluster, dist(umap_mat))[, 3]), 3)
cat("Silhouette:", sil_val, "\n")

umap_ann$cluster <- km10$cluster

# Clean term name helper
clean_term <- function(t) {
  t <- gsub("^(GO_BP:|GO_MF:|GO_CC:|MSig_[^:]+:|KEGG:|Reactome:|TF_family:|PFAM:|Desc:|Feature:|DescWord:)", "", t)
  repeat {
    t_old <- t
    t <- gsub("^(GOBP[_ ]|GOMF[_ ]|GOCC[_ ])", "", t)
    t <- gsub("^(GO:MF:GOMF[_ ]?|GO:BP:GOBP[_ ]?|GO:CC:GOCC[_ ]?)", "", t)
    t <- gsub("^(HALLMARK[_ ]|HP[_ ]|HPO:HP[_ ])", "", t)
    t <- gsub("^(REACTOME[_ ]|CP:REACTOME:REACTOME[_ ])", "", t)
    t <- gsub("^(WP[_ ]|KEGG[_ ]?LEGACY[_ ]?|PID[_ ]|BIOCARTA[_ ])", "", t)
    t <- gsub("^(TFT:GTRD:|TFT:TFT[_ ]?LEGACY:)", "", t)
    t <- gsub("^(CGN[_ ]|MANNO[_ ]|GAUTAM[_ ]|DESCARTES[_ ])", "", t)
    t <- gsub("^(PF[0-9]+)$", "Pfam \\1", t)
    if (t == t_old) break
  }
  t <- gsub("_", " ", t)
  t <- gsub("  +", " ", trimws(t))
  if (nchar(t) > 0) t <- paste0(toupper(substr(t,1,1)), tolower(substr(t,2,nchar(t))))
  t
}

generic_roots <- tolower(c(
  "molecular function", "cellular component", "biological process",
  "binding", "catalytic activity", "protein binding",
  "cell", "cell part", "organelle", "intracellular",
  "cytoplasm", "membrane", "nucleus",
  "cellular process", "metabolic process",
  "multicellular organism development", "system development",
  "anatomical structure development", "regulation of biological process",
  "ion binding", "metal ion binding", "cation binding",
  "transferase activity", "hydrolase activity",
  "signaling", "signal transduction", "cellular response to stimulus",
  "response to stimulus", "developmental process", "multicellular organismal process"
))

# =============================================================================
# 3. PER-CLUSTER ENRICHMENT VIA clusterProfiler (proper ORA with background)
# =============================================================================
cat("\nRunning per-cluster ORA enrichment...\n")

all_entrez <- na.omit(unique(sym_to_entrez[umap_ann$gene]))

cluster_ora_results <- list()

for (cl in 1:10) {
  cl_genes <- umap_ann$gene[umap_ann$cluster == cl]
  cl_entrez <- na.omit(sym_to_entrez[cl_genes])
  n_cl <- length(cl_genes)

  cat(sprintf("  Cluster %d: %d genes, %d with entrez...\n", cl, n_cl, length(cl_entrez)))

  # GO BP
  ego_bp <- tryCatch(enrichGO(gene = cl_entrez, universe = all_entrez,
    OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, readable = TRUE, minGSSize = 5, maxGSSize = 500),
    error = function(e) NULL)

  # GO MF
  ego_mf <- tryCatch(enrichGO(gene = cl_entrez, universe = all_entrez,
    OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, readable = TRUE, minGSSize = 5, maxGSSize = 500),
    error = function(e) NULL)

  # GO CC
  ego_cc <- tryCatch(enrichGO(gene = cl_entrez, universe = all_entrez,
    OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, readable = TRUE, minGSSize = 5, maxGSSize = 500),
    error = function(e) NULL)

  # KEGG
  ekegg <- tryCatch(enrichKEGG(gene = cl_entrez, universe = all_entrez,
    organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 0.05,
    minGSSize = 5, maxGSSize = 500),
    error = function(e) NULL)

  # Combine results
  combine_res <- function(res, source_name) {
    if (is.null(res) || nrow(res@result %>% dplyr::filter(p.adjust < 0.05)) == 0) return(NULL)
    res@result %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::mutate(source = source_name) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, source)
  }

  combined <- rbind(
    combine_res(ego_bp, "GO:BP"),
    combine_res(ego_mf, "GO:MF"),
    combine_res(ego_cc, "GO:CC"),
    combine_res(ekegg, "KEGG")
  )

  if (!is.null(combined) && nrow(combined) > 0) {
    # Clean KEGG descriptions
    combined$Description <- gsub(" - Mus musculus.*$", "", combined$Description)

    # Parse GeneRatio to numeric
    combined <- combined %>%
      dplyr::mutate(
        gene_ratio_num = sapply(GeneRatio, function(x) {
          parts <- as.numeric(strsplit(x, "/")[[1]])
          parts[1] / parts[2]
        }),
        neg_log10_fdr = -log10(p.adjust)
      )

    cluster_ora_results[[cl]] <- combined
    cat(sprintf("    -> %d sig terms (BP:%d, MF:%d, CC:%d, KEGG:%d)\n",
      nrow(combined),
      sum(combined$source == "GO:BP"),
      sum(combined$source == "GO:MF"),
      sum(combined$source == "GO:CC"),
      sum(combined$source == "KEGG")))
  } else {
    cat("    -> 0 sig terms\n")
    cluster_ora_results[[cl]] <- NULL
  }
}

# =============================================================================
# 4. GENERATE PER-CLUSTER DOTPLOTS
# =============================================================================
cat("\nGenerating per-cluster dotplots...\n")

source_colors <- c("GO:BP" = "#66C2A5", "GO:MF" = "#FC8D62",
                    "GO:CC" = "#8DA0CB", "KEGG" = "#E78AC3")

for (cl in 1:10) {
  res <- cluster_ora_results[[cl]]
  if (is.null(res) || nrow(res) == 0) {
    cat(sprintf("  Cluster %d: no sig terms, skipping dotplot\n", cl))
    next
  }

  # Simplify: remove very generic terms
  res <- res %>% dplyr::filter(!tolower(Description) %in% generic_roots)

  # Take top terms per source, then combine
  top_res <- res %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(source, p.adjust)

  # Truncate long descriptions
  top_res$Description_short <- ifelse(nchar(top_res$Description) > 50,
    paste0(substr(top_res$Description, 1, 47), "..."), top_res$Description)

  # Order by source then FDR
  top_res$Description_short <- factor(top_res$Description_short,
    levels = rev(top_res$Description_short))

  n_cl <- sum(umap_ann$cluster == cl)
  k27_pct <- round(100 * mean(umap_ann$H3K27me3_bin[umap_ann$cluster == cl] == 1, na.rm = TRUE))

  p_dot <- ggplot(top_res, aes(x = gene_ratio_num, y = Description_short)) +
    geom_point(aes(size = Count, color = neg_log10_fdr), alpha = 0.8) +
    scale_color_gradient(low = "grey70", high = "firebrick3",
      name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(2, 8), name = "Gene count") +
    facet_grid(source ~ ., scales = "free_y", space = "free_y",
      switch = "y") +
    labs(x = "Gene ratio", y = NULL,
      title = paste0("Cluster ", cl, " (n=", n_cl, ", ", k27_pct, "% H3K27me3+)")) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  h <- max(4, 0.35 * nrow(top_res) + 2)
  fname <- sprintf("Fig_cluster%d_dotplot.pdf", cl)
  ggsave(fname, p_dot, width = 10, height = h)
  cat(sprintf("  Saved: %s\n", fname))
}

# =============================================================================
# 5. REGENERATE UMAP WITH K=10
# =============================================================================
cat("\nGenerating k=10 UMAP figure...\n")

# Label clusters using matrix-based Fisher enrichment (same as v3)
label_terms <- grep("^(GO_|MSig_|KEGG|Reactome|TF_|PFAM)", colnames(mat_ann), value = TRUE)
n_total <- nrow(mat_ann)

centroids <- umap_ann %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = median(UMAP1), y = median(UMAP2), n = n(), .groups = "drop")

cluster_short_labels <- character(10)
for (cl in 1:10) {
  cl_genes <- umap_ann$gene[umap_ann$cluster == cl]
  cl_genes_in_mat <- intersect(cl_genes, rownames(mat_ann))
  n_cl <- length(cl_genes_in_mat)

  if (n_cl < 3) { cluster_short_labels[cl] <- "Small"; next }

  cl_mat <- mat_ann[cl_genes_in_mat, label_terms, drop = FALSE]
  rest_genes <- setdiff(rownames(mat_ann), cl_genes_in_mat)
  in_cl <- colSums(cl_mat)
  in_rest <- colSums(mat_ann[rest_genes, label_terms, drop = FALSE])

  candidates <- which(in_cl >= 3)
  if (length(candidates) == 0) { cluster_short_labels[cl] <- "Diverse"; next }

  fisher_p <- sapply(candidates, function(j) {
    fisher.test(matrix(c(in_cl[j], n_cl - in_cl[j],
      in_rest[j], (n_total - n_cl) - in_rest[j]), nrow = 2),
      alternative = "greater")$p.value
  })

  res <- data.frame(
    term_raw = label_terms[candidates], fisher_p = fisher_p,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(fdr = p.adjust(fisher_p, method = "BH"),
      term_clean = sapply(term_raw, clean_term),
      term_lower = tolower(term_clean)) %>%
    dplyr::filter(fdr < 0.05, !term_lower %in% generic_roots) %>%
    dplyr::filter(!duplicated(term_lower)) %>%
    dplyr::arrange(fisher_p) %>%
    dplyr::slice_head(n = 2)

  if (nrow(res) == 0) { cluster_short_labels[cl] <- "Diverse"; next }

  short <- paste(res$term_clean, collapse = " / ")
  if (nchar(short) > 50) short <- paste0(substr(short, 1, 47), "...")
  cluster_short_labels[cl] <- short
}

centroids$label <- cluster_short_labels[centroids$cluster]

# Plot center for term placement
plot_cx <- mean(range(umap_ann$UMAP1))
plot_cy <- mean(range(umap_ann$UMAP2))

# Build enriched term annotations radiating outward
term_annot_list <- list()
for (cl in 1:10) {
  ct <- centroids[centroids$cluster == cl, ]
  res <- cluster_ora_results[[cl]]
  if (is.null(res) || nrow(res) == 0) next

  # Top 4 non-generic terms
  top_terms <- res %>%
    dplyr::filter(!tolower(Description) %in% generic_roots) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 4)

  if (nrow(top_terms) == 0) next

  dx <- ct$x - plot_cx
  dy <- ct$y - plot_cy
  dist_fc <- sqrt(dx^2 + dy^2)
  if (dist_fc < 0.01) { dx <- 1; dy <- 0; dist_fc <- 1 }
  dx_n <- dx / dist_fc
  dy_n <- dy / dist_fc

  for (j in 1:nrow(top_terms)) {
    perp_offset <- (j - (nrow(top_terms) + 1) / 2) * 0.9
    desc_short <- ifelse(nchar(top_terms$Description[j]) > 35,
      paste0(substr(top_terms$Description[j], 1, 32), "..."),
      top_terms$Description[j])
    term_annot_list[[paste0(cl,"_",j)]] <- data.frame(
      cluster = cl,
      term = paste0(desc_short, " (", top_terms$Count[j], ")"),
      x_anchor = ct$x, y_anchor = ct$y,
      x_label = ct$x + dx_n * 4.0 - dy_n * perp_offset,
      y_label = ct$y + dy_n * 4.0 + dx_n * perp_offset,
      stringsAsFactors = FALSE)
  }
}
term_annot <- do.call(rbind, term_annot_list)

# Colors
pal10 <- c(scales::hue_pal()(10))
names(pal10) <- as.character(1:10)

umap_ann$cluster_f <- factor(umap_ann$cluster)

p_umap <- ggplot(umap_ann, aes(x = UMAP1, y = UMAP2, color = cluster_f)) +
  geom_point(size = 1.6, alpha = 0.55) +
  {if (nrow(umap_unann) > 0)
    geom_point(data = umap_unann, color = "grey60", size = 0.8, alpha = 0.3,
      inherit.aes = FALSE, aes(x = UMAP1, y = UMAP2))} +
  geom_label(data = centroids,
    aes(x = x, y = y, label = paste0(cluster, " (n=", n, ")")),
    size = 3.5, fontface = "bold", fill = alpha("white", 0.8),
    label.size = 0.3, color = "black", inherit.aes = FALSE) +
  geom_segment(data = term_annot,
    aes(x = x_anchor, y = y_anchor, xend = x_label, yend = y_label),
    color = "grey50", linewidth = 0.15, alpha = 0.4, inherit.aes = FALSE) +
  geom_text(data = term_annot,
    aes(x = x_label, y = y_label, label = term),
    size = 2.6, color = "grey20", hjust = 0.5, vjust = 0.5,
    fontface = "italic", inherit.aes = FALSE) +
  scale_color_manual(values = pal10, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(title = "Functional landscape of GNP differentiation genes (k=10)") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.margin = margin(20, 90, 20, 90, "pt"))

ggsave("Fig_gene_umap_k10_clusters.pdf", p_umap, width = 20, height = 16)
cat("Saved: Fig_gene_umap_k10_clusters.pdf\n")

# =============================================================================
# 6. CLUSTER SUMMARY TABLE
# =============================================================================
cat("\n=== CLUSTER SUMMARIES ===\n")

for (cl in 1:10) {
  cl_genes <- umap_ann$gene[umap_ann$cluster == cl]
  n_cl <- length(cl_genes)
  k27_pct <- round(100 * mean(umap_ann$H3K27me3_bin[umap_ann$cluster == cl] == 1, na.rm = TRUE))

  res <- cluster_ora_results[[cl]]
  if (!is.null(res) && nrow(res) > 0) {
    top5 <- res %>%
      dplyr::filter(!tolower(Description) %in% generic_roots) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 5)
    terms_str <- paste(top5$Description, collapse = "; ")
  } else {
    terms_str <- "No significant enrichment"
  }

  genes_str <- paste(head(sort(cl_genes), 8), collapse = ", ")

  cat(sprintf("\nCluster %d: %s\n  N=%d, %d%% H3K27me3+\n  Top terms: %s\n  Genes: %s\n",
    cl, cluster_short_labels[cl], n_cl, k27_pct, terms_str, genes_str))
}

# Save updated coordinates
umap_full <- rbind(
  umap_ann %>% dplyr::select(gene, UMAP1, UMAP2, cluster, n_annotations, H3K27me3_bin, t50_cat),
  umap_unann %>% dplyr::select(gene, UMAP1, UMAP2, cluster, n_annotations, H3K27me3_bin, t50_cat)
)
write.csv(umap_full, "gene_umap_k10_coordinates.csv", row.names = FALSE)

cat("\n=== ALL OUTPUTS GENERATED SUCCESSFULLY ===\n")
