suppressMessages({
  library(tidyverse)
  library(org.Mm.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(msigdbr)
  library(uwot)
  library(cluster)
  library(ggrepel)
  library(patchwork)
  library(irlba)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. LOAD DATA
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

gene_map <- bitr(diff_genes$ensembl_gene_id,
  fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Mm.eg.db)
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]
gene_map <- gene_map %>%
  dplyr::left_join(diff_genes %>%
    dplyr::select(ensembl_gene_id, H3K27me3_bin, mgi_symbol, t50_cat, description),
    by = c("ENSEMBL" = "ensembl_gene_id"))

all_symbols <- unique(gene_map$SYMBOL)
all_entrez <- unique(gene_map$ENTREZID)
entrez_to_sym <- setNames(gene_map$SYMBOL, gene_map$ENTREZID)
sym_to_entrez <- setNames(gene_map$ENTREZID, gene_map$SYMBOL)

cat("Diff genes with IDs:", length(all_symbols), "\n")

# =============================================================================
# 2. SOURCE 1: ALL GO ANNOTATIONS (not just enriched — every annotation per gene)
# =============================================================================
cat("\n--- Source 1: Direct GO annotations per gene ---\n")

# Get all GO terms for each diff gene directly from org.Mm.eg.db
go_all <- AnnotationDbi::select(org.Mm.eg.db, keys = all_entrez,
  columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")
go_all <- go_all %>%
  dplyr::filter(!is.na(GO), ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  dplyr::mutate(gene = entrez_to_sym[ENTREZID]) %>%
  dplyr::filter(!is.na(gene))

# Get GO term names
go_terms_used <- unique(go_all$GO)
go_names <- AnnotationDbi::select(GO.db, keys = go_terms_used, columns = "TERM", keytype = "GOID")
go_all <- go_all %>%
  dplyr::left_join(go_names, by = c("GO" = "GOID")) %>%
  dplyr::mutate(term = paste0("GO_", ONTOLOGY, ":", TERM)) %>%
  dplyr::filter(!is.na(TERM)) %>%
  dplyr::select(gene, term) %>%
  dplyr::distinct()

cat("  GO annotations:", nrow(go_all), " unique terms:", length(unique(go_all$term)), "\n")

# =============================================================================
# 3. SOURCE 2: MSigDB — ALL collections
# =============================================================================
cat("\n--- Source 2: MSigDB gene sets ---\n")

msig_collections <- list(
  c("H", NA),           # Hallmark
  c("C2", "CP:KEGG_LEGACY"),
  c("C2", "CP:REACTOME"),
  c("C2", "CP:WIKIPATHWAYS"),
  c("C2", "CP:PID"),
  c("C2", "CP:BIOCARTA"),
  c("C3", "TFT:GTRD"),       # TF targets
  c("C3", "TFT:TFT_LEGACY"),
  c("C4", "CGN"),             # Cancer gene neighborhoods
  c("C5", "GO:BP"),
  c("C5", "GO:MF"),
  c("C5", "GO:CC"),
  c("C5", "HPO"),             # Human phenotype ontology
  c("C8", NA)                 # Cell type signatures
)

msig_long_list <- list()
for (coll in msig_collections) {
  cat_name <- coll[1]; subcat <- coll[2]
  label <- if (is.na(subcat)) cat_name else paste0(cat_name, ":", subcat)

  gs <- tryCatch({
    if (is.na(subcat)) {
      msigdbr(species = "Mus musculus", category = cat_name)
    } else {
      msigdbr(species = "Mus musculus", category = cat_name, subcategory = subcat)
    }
  }, error = function(e) NULL)

  if (is.null(gs) || nrow(gs) == 0) {
    cat("  ", label, ": skipped (no data)\n")
    next
  }

  # Direct membership — no enrichment filter, just which diff genes belong to which sets
  gs_filt <- gs %>%
    dplyr::filter(as.character(entrez_gene) %in% all_entrez) %>%
    dplyr::mutate(gene = entrez_to_sym[as.character(entrez_gene)],
      term = paste0("MSig_", gsub(" ", "_", label), ":", gs_name)) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::select(gene, term) %>%
    dplyr::distinct()

  cat("  ", label, ":", length(unique(gs_filt$term)), "sets,",
      length(unique(gs_filt$gene)), "genes\n")
  msig_long_list[[label]] <- gs_filt
}

msig_long <- do.call(rbind, msig_long_list)
cat("  Total MSigDB pairs:", nrow(msig_long), "\n")

# =============================================================================
# 4. SOURCE 3: TF annotation from project files
# =============================================================================
cat("\n--- Source 3: TF family annotation ---\n")

tf_anno <- read.csv("~/Dropbox/GNP_MB_integrate/annotation/MM_TF.csv", stringsAsFactors = FALSE)
tf_diff <- tf_anno %>%
  dplyr::filter(Symbol %in% all_symbols) %>%
  dplyr::mutate(term = paste0("TF_family:", Family)) %>%
  dplyr::select(gene = Symbol, term) %>%
  dplyr::distinct()

cat("  TF genes:", length(unique(tf_diff$gene)), " families:", length(unique(tf_diff$term)), "\n")

# Also add a generic "is_TF" annotation
tf_generic <- data.frame(gene = unique(tf_diff$gene), term = "Feature:Transcription_factor",
  stringsAsFactors = FALSE)

# =============================================================================
# 5. SOURCE 4: Protein domains (InterPro via org.Mm.eg.db)
# =============================================================================
cat("\n--- Source 4: InterPro domains ---\n")

# Use PFAM as proxy for protein domains
pfam_data <- tryCatch({
  AnnotationDbi::select(org.Mm.eg.db, keys = all_entrez,
    columns = "PFAM", keytype = "ENTREZID")
}, error = function(e) NULL)

if (!is.null(pfam_data)) {
  pfam_long <- pfam_data %>%
    dplyr::filter(!is.na(PFAM)) %>%
    dplyr::mutate(gene = entrez_to_sym[ENTREZID],
      term = paste0("PFAM:", PFAM)) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::select(gene, term) %>%
    dplyr::distinct()
  cat("  PFAM domains:", length(unique(pfam_long$term)), " covering",
      length(unique(pfam_long$gene)), "genes\n")
} else {
  pfam_long <- data.frame(gene = character(), term = character())
  cat("  PFAM: not available\n")
}

# =============================================================================
# 6. SOURCE 5: Description-based keywords
# =============================================================================
cat("\n--- Source 5: Gene description keywords ---\n")

desc_map <- setNames(tolower(gene_map$description), gene_map$SYMBOL)

keyword_rules <- list(
  "Desc:ion_channel" = "ion channel|voltage.gated|ligand.gated channel",
  "Desc:receptor" = "receptor(?!.*olfact)",
  "Desc:kinase" = "kinase",
  "Desc:phosphatase" = "phosphatase",
  "Desc:GTPase" = "gtpase|gtp.binding|ras.related|rho family|rab family",
  "Desc:transporter" = "transporter|solute carrier|carrier family",
  "Desc:cell_adhesion" = "adhesion|cadherin|integrin|selectin",
  "Desc:cytoskeleton" = "actin|tubulin|myosin|dynein|kinesin|spectrin|filament|cytoskelet",
  "Desc:ECM" = "extracellular matrix|collagen|laminin|proteoglycan|fibronectin",
  "Desc:enzyme" = "synthase|transferase|dehydrogenase|oxidase|reductase|hydrolase|esterase|ligase|protease|peptidase|lyase|dismutase|epimerase|isomerase",
  "Desc:ubiquitin" = "ubiquitin|proteasom|sumo",
  "Desc:vesicle_traffic" = "vesicl|endosom|lysosom|golgi|exocyt|trafficking",
  "Desc:RNA_binding" = "rna.binding|splicing|ribosom|translat|hnrnp|rbm|rbfox",
  "Desc:lipid" = "lipid|fatty acid|sphingo|phospholipid|ceramide|sterol|cholest",
  "Desc:apoptosis" = "apoptosis|death|caspase|bcl",
  "Desc:immune" = "immun|inflam|cytokine|interleukin|complement",
  "Desc:calcium_binding" = "calcium.bind|calmodulin|ef.hand|annexin|s100",
  "Desc:zinc_finger" = "zinc finger",
  "Desc:homeobox" = "homeobox|homeodomain",
  "Desc:bHLH" = "helix.loop.helix|bhlh",
  "Desc:GPCR" = "g protein.coupled|gpcr|7 transmembrane",
  "Desc:growth_factor" = "growth factor",
  "Desc:scaffold_adaptor" = "scaffold|adaptor|pdz|sh3|sh2 domain",
  "Desc:nuclear_receptor" = "nuclear receptor|hormone receptor",
  "Desc:secreted" = "secreted|extracellular",
  "Desc:predicted_uncharacterized" = "riken|predicted gene|expressed sequence|hypothetical|uncharacterized|cdna"
)

desc_long_list <- list()
for (kw_name in names(keyword_rules)) {
  pattern <- keyword_rules[[kw_name]]
  matched <- names(desc_map)[grepl(pattern, desc_map, perl = TRUE)]
  if (length(matched) > 0) {
    desc_long_list[[kw_name]] <- data.frame(gene = matched, term = kw_name,
      stringsAsFactors = FALSE)
  }
}
desc_long <- do.call(rbind, desc_long_list)
cat("  Keyword annotations:", nrow(desc_long), " covering",
    length(unique(desc_long$gene)), "genes\n")

# =============================================================================
# 7. SOURCE 6: H3K27me3 and timing as features
# =============================================================================
cat("\n--- Source 6: Epigenetic and timing features ---\n")

epi_features <- gene_map %>%
  dplyr::select(SYMBOL, H3K27me3_bin, t50_cat) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

epi_long <- rbind(
  epi_features %>% dplyr::filter(H3K27me3_bin == 1) %>%
    dplyr::mutate(term = "Feature:H3K27me3_bound") %>%
    dplyr::select(gene = SYMBOL, term),
  epi_features %>% dplyr::mutate(term = paste0("Feature:t50_", t50_cat)) %>%
    dplyr::select(gene = SYMBOL, term)
)
cat("  Epigenetic features:", nrow(epi_long), "\n")

# =============================================================================
# 8. COMBINE ALL SOURCES
# =============================================================================
cat("\n=== COMBINING ALL SOURCES ===\n")

all_long <- rbind(go_all, msig_long, tf_diff, tf_generic, pfam_long, desc_long, epi_long) %>%
  dplyr::filter(gene %in% all_symbols) %>%
  dplyr::distinct()

cat("Total gene-term pairs:", nrow(all_long), "\n")
cat("Unique terms:", length(unique(all_long$term)), "\n")
cat("Genes covered:", length(unique(all_long$gene)), "/", length(all_symbols), "\n")
cat("Genes still with 0 terms:", sum(!all_symbols %in% all_long$gene), "\n")

# =============================================================================
# 9. PRUNE TERMS — remove too rare, too common, or uninformative
# =============================================================================
cat("\nPruning terms...\n")

term_stats <- all_long %>%
  dplyr::group_by(term) %>%
  dplyr::summarise(n_genes = n(), .groups = "drop")

# Remove terms annotating <3 genes (no discriminative power) or >400 genes (too generic)
# Also remove very generic GO terms
generic_go <- c(
  "GO_BP:biological_process", "GO_BP:cellular process",
  "GO_BP:metabolic process", "GO_MF:binding",
  "GO_MF:catalytic activity", "GO_MF:protein binding",
  "GO_CC:cell", "GO_CC:cell part", "GO_CC:organelle",
  "GO_CC:intracellular", "GO_CC:cytoplasm",
  "GO_CC:membrane", "GO_CC:nucleus"
)

term_stats$keep <- term_stats$n_genes >= 3 & term_stats$n_genes <= 400 &
  !term_stats$term %in% generic_go

# Also prune GO terms that are essentially duplicates of parent terms
# (keep more specific ones by preferring those with fewer genes)
cat("  Before pruning:", nrow(term_stats), "terms\n")
cat("  After pruning:", sum(term_stats$keep), "terms\n")
cat("  Removed <3 genes:", sum(term_stats$n_genes < 3), "\n")
cat("  Removed >400 genes:", sum(term_stats$n_genes > 400), "\n")

good_terms <- term_stats %>% dplyr::filter(keep) %>% dplyr::pull(term)
all_long_filt <- all_long %>% dplyr::filter(term %in% good_terms)

cat("Filtered gene-term pairs:", nrow(all_long_filt), "\n")
cat("Genes covered after filter:", length(unique(all_long_filt$gene)), "/", length(all_symbols), "\n")

# Check coverage
still_missing <- setdiff(all_symbols, unique(all_long_filt$gene))
if (length(still_missing) > 0) {
  cat("\nStill missing genes:", length(still_missing), "\n")
  cat("Adding description-based fallback annotations...\n")

  # For any remaining gene, add a description keyword if possible
  for (g in still_missing) {
    d <- desc_map[g]
    if (!is.na(d) && nchar(d) > 0) {
      # Extract first meaningful word from description
      words <- strsplit(d, "[ ,;]+")[[1]]
      words <- words[nchar(words) > 3 & !words %in% c("gene", "protein", "family", "member",
        "containing", "like", "type", "domain", "with", "encoding", "predicted")]
      if (length(words) > 0) {
        fallback_term <- paste0("DescWord:", paste(head(words, 3), collapse = "_"))
        all_long_filt <- rbind(all_long_filt,
          data.frame(gene = g, term = fallback_term, stringsAsFactors = FALSE))
      }
    }
  }

  # Re-check
  still_missing2 <- setdiff(all_symbols, unique(all_long_filt$gene))
  cat("After fallback:", length(still_missing2), "genes still missing\n")
  if (length(still_missing2) > 0) {
    cat("Truly unannotated:", paste(head(still_missing2, 10), collapse = ", "), "\n")
  }
}

# =============================================================================
# 10. BUILD MATRIX
# =============================================================================
cat("\nBuilding gene x term matrix...\n")

gene_term_wide <- all_long_filt %>%
  dplyr::mutate(val = 1) %>%
  tidyr::pivot_wider(names_from = term, values_from = val, values_fill = 0)

# Add truly missing genes with all zeros
missing_final <- setdiff(all_symbols, gene_term_wide$gene)
if (length(missing_final) > 0) {
  zero_df <- data.frame(gene = missing_final, stringsAsFactors = FALSE)
  for (col in colnames(gene_term_wide)[-1]) zero_df[[col]] <- 0
  gene_term_wide <- rbind(gene_term_wide, zero_df)
}

gene_names <- gene_term_wide$gene
mat <- as.matrix(gene_term_wide[, -1])
rownames(mat) <- gene_names

n_annot <- rowSums(mat)
cat("Matrix:", nrow(mat), "genes x", ncol(mat), "terms\n")
cat("Genes with 0 annotations:", sum(n_annot == 0), "\n")
cat("Median annotations per gene:", median(n_annot), "\n")
cat("Mean annotations per gene:", round(mean(n_annot), 1), "\n")
cat("Quantiles:\n")
print(quantile(n_annot, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)))

# =============================================================================
# 11. PCA THEN UMAP
# =============================================================================
cat("\nRunning PCA (", ncol(mat), "terms -> top PCs)...\n")

# Separate annotated vs unannotated
annotated_idx <- n_annot > 0
mat_ann <- mat[annotated_idx, ]
unannotated_genes <- rownames(mat)[!annotated_idx]

# Save matrix for reuse
saveRDS(mat_ann, "gene_term_matrix_v3.rds")
cat("Saved: gene_term_matrix_v3.rds\n")

cat("Annotated genes:", nrow(mat_ann), "  Unannotated:", length(unannotated_genes), "\n")

# Remove zero-variance columns for PCA
col_var <- apply(mat_ann, 2, var)
mat_ann_nz <- mat_ann[, col_var > 0]
cat("Non-zero-variance terms:", ncol(mat_ann_nz), "\n")

# Truncated SVD / PCA
n_pcs <- min(100, ncol(mat_ann_nz) - 1, nrow(mat_ann_nz) - 1)
cat("Computing top", n_pcs, "PCs via irlba...\n")

set.seed(42)
pca_res <- irlba::prcomp_irlba(mat_ann_nz, n = n_pcs, center = TRUE, scale. = FALSE)
pca_coords <- pca_res$x

var_explained <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
cat("Variance explained by top PCs:\n")
cat("  PC10:", round(var_explained[10]*100, 1), "%\n")
cat("  PC30:", round(var_explained[30]*100, 1), "%\n")
cat("  PC50:", round(var_explained[min(50, n_pcs)]*100, 1), "%\n")
cat("  PC100:", round(var_explained[min(100, n_pcs)]*100, 1), "%\n")

# Use enough PCs to capture ~80% variance
n_use <- min(which(var_explained >= 0.80))
if (is.na(n_use) || length(n_use) == 0) n_use <- n_pcs
cat("Using", n_use, "PCs (", round(var_explained[n_use]*100, 1), "% variance)\n")

pca_use <- pca_coords[, 1:n_use]

# UMAP on PCA coordinates
cat("Running UMAP on PCA coordinates...\n")
set.seed(42)
umap_res <- umap(pca_use, n_neighbors = 15, min_dist = 0.15, metric = "cosine",
  n_components = 2, n_epochs = 1000, spread = 2.0)

# =============================================================================
# 12. CLUSTERING
# =============================================================================
cat("Clustering...\n")

sil_scores <- sapply(6:18, function(k) {
  set.seed(42)
  km <- kmeans(umap_res, centers = k, nstart = 30, iter.max = 200)
  mean(silhouette(km$cluster, dist(umap_res))[, 3])
})
names(sil_scores) <- 6:18
cat("Silhouette scores:\n")
print(round(sil_scores, 3))
best_k <- as.integer(names(which.max(sil_scores)))
cat("Best k:", best_k, "\n")

use_k <- max(best_k, 8)
cat("Using k =", use_k, "\n")

set.seed(42)
km_final <- kmeans(umap_res, centers = use_k, nstart = 50, iter.max = 300)

# Build data frame
umap_df <- data.frame(
  gene = rownames(mat_ann),
  UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2],
  cluster = km_final$cluster,
  n_annotations = rowSums(mat_ann),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(gene_map %>% dplyr::select(SYMBOL, H3K27me3_bin, t50_cat) %>% distinct(),
    by = c("gene" = "SYMBOL"))

# Add unannotated at periphery
if (length(unannotated_genes) > 0) {
  xr <- range(umap_df$UMAP1); yr <- range(umap_df$UMAP2)
  set.seed(42)
  umap_df <- rbind(umap_df, data.frame(
    gene = unannotated_genes,
    UMAP1 = runif(length(unannotated_genes), xr[2]+1, xr[2]+3),
    UMAP2 = runif(length(unannotated_genes), yr[1], yr[2]),
    cluster = 0, n_annotations = 0,
    H3K27me3_bin = gene_map$H3K27me3_bin[match(unannotated_genes, gene_map$SYMBOL)],
    t50_cat = gene_map$t50_cat[match(unannotated_genes, gene_map$SYMBOL)],
    stringsAsFactors = FALSE))
}

cat("\nCluster sizes:\n")
print(sort(table(umap_df$cluster[umap_df$cluster > 0]), decreasing = TRUE))
if (sum(umap_df$cluster == 0) > 0) cat("Unannotated (cluster 0):", sum(umap_df$cluster == 0), "\n")

# =============================================================================
# 13. LABEL CLUSTERS — Fisher's exact enrichment ranking
# =============================================================================
cat("\nLabeling clusters via Fisher's exact enrichment...\n")

# Only use biologically meaningful terms for labeling
label_terms <- grep("^(GO_|MSig_|KEGG|Reactome|TF_|PFAM)", colnames(mat_ann), value = TRUE)

n_total <- nrow(mat_ann)

# Clean term name helper — aggressively strip all prefixes
clean_term <- function(t) {
  # 1. Strip the source prefix we added (GO_BP:, MSig_...:, etc)
  t <- gsub("^(GO_BP:|GO_MF:|GO_CC:|MSig_[^:]+:|KEGG:|Reactome:|TF_family:|PFAM:|Desc:|Feature:|DescWord:)", "", t)
  # 2. Now strip any remaining MSigDB-internal prefixes — iteratively
  # These come from MSigDB gene set names like "GOMF_TRANSPORTER_ACTIVITY"
  repeat {
    t_old <- t
    t <- gsub("^(GOBP[_ ]|GOMF[_ ]|GOCC[_ ])", "", t)
    t <- gsub("^(GO:MF:GOMF[_ ]?|GO:BP:GOBP[_ ]?|GO:CC:GOCC[_ ]?)", "", t)
    t <- gsub("^(GO:MF:[_ ]?|GO:BP:[_ ]?|GO:CC:[_ ]?)", "", t)
    t <- gsub("^(HALLMARK[_ ]|HP[_ ]|HPO:HP[_ ])", "", t)
    t <- gsub("^(REACTOME[_ ]|CP:REACTOME:REACTOME[_ ])", "", t)
    t <- gsub("^(WP[_ ]|KEGG[_ ]?LEGACY[_ ]?|PID[_ ]|BIOCARTA[_ ])", "", t)
    t <- gsub("^(TFT:GTRD:|TFT:TFT[_ ]?LEGACY:)", "", t)
    t <- gsub("^(CGN[_ ]|MANNO[_ ]|GAUTAM[_ ]|DESCARTES[_ ])", "", t)
    t <- gsub("^(PF[0-9]+)$", "Pfam \\1", t)  # keep Pfam IDs readable
    if (t == t_old) break
  }
  t <- gsub("_", " ", t)
  t <- gsub("  +", " ", trimws(t))
  # Sentence case
  if (nchar(t) > 0) {
    t <- paste0(toupper(substr(t, 1, 1)), tolower(substr(t, 2, nchar(t))))
  }
  t
}

# Generic root GO terms to exclude from labeling
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

cluster_labels <- list()
cluster_enriched_terms <- list()

for (cl in sort(unique(umap_df$cluster[umap_df$cluster > 0]))) {
  cl_genes <- umap_df$gene[umap_df$cluster == cl]
  cl_mat <- mat_ann[cl_genes, label_terms, drop = FALSE]
  n_cl <- length(cl_genes)

  # Fisher's exact test for each term: enrichment in cluster vs rest
  in_cl <- colSums(cl_mat)
  in_rest <- colSums(mat_ann[setdiff(rownames(mat_ann), cl_genes), label_terms, drop = FALSE])

  fisher_results <- data.frame(
    term_raw = label_terms,
    in_cluster = in_cl,
    in_rest = in_rest,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(in_cluster >= 3) %>%  # must have >=3 genes
    dplyr::rowwise() %>%
    dplyr::mutate(
      fisher_p = fisher.test(matrix(c(in_cluster, n_cl - in_cluster,
        in_rest, (n_total - n_cl) - in_rest), nrow = 2),
        alternative = "greater")$p.value,
      odds_ratio = fisher.test(matrix(c(in_cluster, n_cl - in_cluster,
        in_rest, (n_total - n_cl) - in_rest), nrow = 2),
        alternative = "greater")$estimate,
      pct_cluster = round(100 * in_cluster / n_cl, 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      fdr = p.adjust(fisher_p, method = "BH"),
      term_clean = sapply(term_raw, clean_term)
    ) %>%
    dplyr::filter(fdr < 0.05) %>%
    dplyr::arrange(fisher_p)

  # Remove generic root terms
  fisher_results$term_lower <- tolower(fisher_results$term_clean)
  fisher_results <- fisher_results %>%
    dplyr::filter(!term_lower %in% generic_roots)

  # Deduplicate near-identical term names (keep most significant)
  fisher_results <- fisher_results[!duplicated(fisher_results$term_lower), ]

  # Top terms for this cluster
  if (nrow(fisher_results) == 0) {
    # No specific enrichment — label as diverse/miscellaneous
    top_n <- 0
    top_terms <- character(0)
    short <- "Diverse / ECM / cell surface"
    fisher_results <- data.frame(term_clean = character(), pct_cluster = numeric(),
      fdr = numeric(), stringsAsFactors = FALSE)
  } else {
    top_n <- min(6, nrow(fisher_results))
    top_terms <- fisher_results$term_clean[1:top_n]
    short <- paste(head(top_terms, 2), collapse = " / ")
    if (nchar(short) > 55) short <- paste0(substr(short, 1, 52), "...")
  }

  k27_pct <- round(100 * mean(umap_df$H3K27me3_bin[umap_df$cluster == cl] == 1, na.rm = TRUE))

  cluster_labels[[as.character(cl)]] <- data.frame(
    cluster = cl, n = n_cl, label = short, k27_pct = k27_pct,
    stringsAsFactors = FALSE)

  # Store enriched terms for annotation around the plot
  cluster_enriched_terms[[as.character(cl)]] <- fisher_results %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::select(term_clean, pct_cluster, fdr)

  cat(sprintf("  Cluster %d (%d genes, %d%% K27+): %s\n", cl, n_cl, k27_pct, short))
  for (tt in 1:min(top_n, 4)) {
    cat(sprintf("    - %s (%s%%, FDR=%.1e)\n",
      fisher_results$term_clean[tt], fisher_results$pct_cluster[tt], fisher_results$fdr[tt]))
  }
}

label_df <- do.call(rbind, cluster_labels)

# Add "Unannotated" label for cluster 0
if (any(umap_df$cluster == 0)) {
  n_unann <- sum(umap_df$cluster == 0)
  k27_unann <- round(100 * mean(umap_df$H3K27me3_bin[umap_df$cluster == 0] == 1, na.rm = TRUE))
  label_df <- rbind(label_df, data.frame(cluster = 0, n = n_unann,
    label = "Unannotated", k27_pct = k27_unann, stringsAsFactors = FALSE))
}

umap_df <- umap_df %>%
  dplyr::left_join(label_df %>% dplyr::select(cluster, label), by = "cluster")

centroids <- umap_df %>%
  dplyr::filter(cluster > 0) %>%
  dplyr::group_by(cluster, label) %>%
  dplyr::summarise(x = median(UMAP1), y = median(UMAP2), n = n(), .groups = "drop")

# =============================================================================
# 14. BUILD TERM ANNOTATION DATA FRAME (placed around cluster periphery)
# =============================================================================
cat("\nBuilding term annotations for figure...\n")

# For each cluster, place enriched terms around the cluster centroid
# radiating outward from center of the whole plot
plot_center_x <- mean(range(umap_df$UMAP1[umap_df$cluster > 0]))
plot_center_y <- mean(range(umap_df$UMAP2[umap_df$cluster > 0]))

term_annot_list <- list()
for (cl in sort(unique(umap_df$cluster[umap_df$cluster > 0]))) {
  ct <- centroids[centroids$cluster == cl, ]
  terms_df <- cluster_enriched_terms[[as.character(cl)]]
  if (is.null(terms_df) || nrow(terms_df) == 0) next

  # Direction from plot center to cluster centroid
  dx <- ct$x - plot_center_x
  dy <- ct$y - plot_center_y
  dist_from_center <- sqrt(dx^2 + dy^2)
  if (dist_from_center < 0.01) { dx <- 1; dy <- 0 }

  # Normalize and push labels outward
  dx_n <- dx / sqrt(dx^2 + dy^2)
  dy_n <- dy / sqrt(dx^2 + dy^2)

  n_terms <- nrow(terms_df)
  offset_base <- 3.5  # distance from centroid

  for (j in 1:n_terms) {
    # Spread terms along perpendicular axis
    perp_offset <- (j - (n_terms + 1) / 2) * 0.8
    term_annot_list[[paste0(cl, "_", j)]] <- data.frame(
      cluster = cl,
      term = terms_df$term_clean[j],
      pct = terms_df$pct_cluster[j],
      x_anchor = ct$x,
      y_anchor = ct$y,
      x_label = ct$x + dx_n * offset_base - dy_n * perp_offset,
      y_label = ct$y + dy_n * offset_base + dx_n * perp_offset,
      stringsAsFactors = FALSE
    )
  }
}

term_annot <- do.call(rbind, term_annot_list)
# Truncate long terms
term_annot$term_short <- ifelse(nchar(term_annot$term) > 40,
  paste0(substr(term_annot$term, 1, 37), "..."), term_annot$term)
term_annot$term_label <- paste0(term_annot$term_short, " (", term_annot$pct, "%)")

# =============================================================================
# 15. FIGURES — clean, no axes, enriched terms around clusters
# =============================================================================
cat("\nGenerating figures...\n")

n_actual <- max(umap_df$cluster)
cluster_pal <- c(scales::hue_pal()(n_actual), "grey60")
names(cluster_pal) <- c(as.character(1:n_actual), "0")

umap_df$cluster_f <- factor(umap_df$cluster)

# Main figure — no axis labels, terms radiating outward
p_main <- ggplot(umap_df %>% dplyr::filter(cluster > 0),
    aes(x = UMAP1, y = UMAP2, color = cluster_f)) +
  geom_point(size = 1.8, alpha = 0.55) +
  {if (any(umap_df$cluster == 0))
    geom_point(data = umap_df %>% dplyr::filter(cluster == 0),
      color = "grey60", size = 1, alpha = 0.3)} +
  # Cluster number labels at centroids
  geom_label(data = centroids,
    aes(x = x, y = y, label = paste0(cluster, " (n=", n, ")")),
    size = 4, fontface = "bold", fill = alpha("white", 0.8),
    label.size = 0.3, color = "black", inherit.aes = FALSE) +
  # Enriched term labels radiating outward
  geom_segment(data = term_annot,
    aes(x = x_anchor, y = y_anchor, xend = x_label, yend = y_label),
    color = "grey50", linewidth = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  geom_text(data = term_annot,
    aes(x = x_label, y = y_label, label = term_label),
    size = 2.8, color = "grey20", hjust = 0.5, vjust = 0.5,
    fontface = "italic", inherit.aes = FALSE) +
  scale_color_manual(values = cluster_pal, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(title = "Functional landscape of GNP differentiation genes",
       subtitle = paste0(nrow(umap_df), " genes | ", ncol(mat),
         " annotations | PCA -> UMAP | Fisher's exact enrichment per cluster")) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5),
        plot.margin = margin(20, 80, 20, 80, "pt"))

ggsave("Fig_gene_umap_v3_clusters.pdf", p_main, width = 18, height = 14)
cat("Saved: Fig_gene_umap_v3_clusters.pdf\n")

# K27 overlay — also clean, no axes
p_k27 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,
    color = factor(H3K27me3_bin, labels = c("Unbound", "H3K27me3-bound")))) +
  geom_point(size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("Unbound" = "grey75", "H3K27me3-bound" = "#4393C3"),
    name = NULL, na.value = "grey90") +
  labs(title = "H3K27me3 binding status") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.position = c(0.85, 0.1), legend.text = element_text(size = 12))

ggsave("Fig_gene_umap_v3_k27.pdf", p_k27, width = 10, height = 8)
cat("Saved: Fig_gene_umap_v3_k27.pdf\n")

# Timing overlay
p_timing <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = t50_cat)) +
  geom_point(size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("Mid" = "#FDB863", "Late" = "#B2182B"),
    name = "Timing", na.value = "grey90") +
  labs(title = "Differentiation timing (t50)") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.position = c(0.85, 0.1), legend.text = element_text(size = 12))

ggsave("Fig_gene_umap_v3_timing.pdf", p_timing, width = 10, height = 8)
cat("Saved: Fig_gene_umap_v3_timing.pdf\n")

# Combined multi-panel
p_combined <- (p_main) / (p_k27 | p_timing) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("Fig_gene_umap_v3_combined.pdf", p_combined, width = 18, height = 20)
cat("Saved: Fig_gene_umap_v3_combined.pdf\n")

# =============================================================================
# 16. SUMMARY
# =============================================================================
cat("\n=== CLUSTER SUMMARY ===\n")
for (i in 1:nrow(label_df)) {
  cl <- label_df$cluster[i]
  cl_genes <- sort(umap_df$gene[umap_df$cluster == cl])
  terms_for_cl <- cluster_enriched_terms[[as.character(cl)]]
  terms_str <- if (!is.null(terms_for_cl)) {
    paste(terms_for_cl$term_clean, collapse = "; ")
  } else { "N/A" }
  cat(sprintf("\nCluster %d: %s\n  N=%d, %d%% H3K27me3+\n  Top enriched: %s\n  Genes: %s\n",
    cl, label_df$label[i], label_df$n[i], label_df$k27_pct[i],
    terms_str, paste(head(cl_genes, 8), collapse = ", ")))
}

write.csv(umap_df, "gene_umap_v3_coordinates.csv", row.names = FALSE)
write.csv(label_df, "gene_umap_v3_cluster_labels.csv", row.names = FALSE)

cat("\n=== ALL OUTPUTS GENERATED SUCCESSFULLY ===\n")
