suppressMessages({
  library(tidyverse)
  library(uwot)
  library(irlba)
  library(cluster)
  library(patchwork)
  library(ggrepel)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. LOAD EXISTING MATRIX + AGENT RESULTS
# =============================================================================
cat("Loading data...\n")

mat_ann <- readRDS("gene_term_matrix_v3.rds")  # 773 x 8828
agent <- read.csv("agent_extracted_annotations.csv", stringsAsFactors=FALSE)

cat("Database matrix:", nrow(mat_ann), "x", ncol(mat_ann), "\n")
cat("Agent annotations:", nrow(agent), "\n")

# Gene metadata
master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude)) %>%
  dplyr::select(mgi_symbol, H3K27me3_bin, t50_cat) %>%
  dplyr::distinct(mgi_symbol, .keep_all=TRUE)

# =============================================================================
# 2. BUILD AGENT FEATURE MATRIX
# =============================================================================
cat("\nBuilding agent feature matrix...\n")

# --- Protein class (one-hot) ---
agent$protein_class[is.na(agent$protein_class) | agent$protein_class == "NA"] <- "other"
pc_levels <- c("ion channel","transporter","GPCR","receptor","kinase","phosphatase",
  "transcription factor","enzyme","structural protein","adhesion molecule",
  "scaffold/adaptor","secreted factor","GTPase/regulator","ubiquitin ligase",
  "protease","motor protein","channel regulator","ECM protein","other")
agent$protein_class[!agent$protein_class %in% pc_levels] <- "other"

pc_mat <- model.matrix(~ 0 + factor(agent$protein_class, levels=pc_levels))
colnames(pc_mat) <- paste0("AgPC:", pc_levels)
rownames(pc_mat) <- agent$gene

# --- Subcellular (multi-label) ---
subcel_terms <- c("plasma membrane","synapse","nucleus","cytoplasm","mitochondria",
  "ER","Golgi","axon","dendrite","growth cone","extracellular")
subcel_mat <- matrix(0, nrow=nrow(agent), ncol=length(subcel_terms))
colnames(subcel_mat) <- paste0("AgLoc:", subcel_terms)
rownames(subcel_mat) <- agent$gene
for (i in 1:nrow(agent)) {
  if (!is.na(agent$subcellular[i]) && agent$subcellular[i] != "NA") {
    for (j in seq_along(subcel_terms)) {
      if (grepl(subcel_terms[j], agent$subcellular[i], ignore.case=TRUE))
        subcel_mat[i,j] <- 1
    }
  }
}

# --- Tissues (multi-label) ---
tissue_terms <- c("brain","cerebellum","cortex","hippocampus","spinal cord","retina",
  "heart","muscle","kidney","liver","lung","bone","pancreas","intestine")
tissue_mat <- matrix(0, nrow=nrow(agent), ncol=length(tissue_terms))
colnames(tissue_mat) <- paste0("AgTiss:", tissue_terms)
rownames(tissue_mat) <- agent$gene
for (i in 1:nrow(agent)) {
  if (!is.na(agent$tissues[i]) && agent$tissues[i] != "NA") {
    for (j in seq_along(tissue_terms)) {
      if (grepl(tissue_terms[j], agent$tissues[i], ignore.case=TRUE))
        tissue_mat[i,j] <- 1
    }
  }
}

# --- Cell types (multi-label) ---
ct_terms <- c("granule neuron","Purkinje","interneuron","astrocyte","oligodendrocyte",
  "Schwann","neuron","neural progenitor","microglia","endothelial","epithelial",
  "macrophage","cardiomyocyte","fibroblast")
ct_mat <- matrix(0, nrow=nrow(agent), ncol=length(ct_terms))
colnames(ct_mat) <- paste0("AgCell:", ct_terms)
rownames(ct_mat) <- agent$gene
for (i in 1:nrow(agent)) {
  if (!is.na(agent$cell_types[i]) && agent$cell_types[i] != "NA") {
    for (j in seq_along(ct_terms)) {
      if (grepl(ct_terms[j], agent$cell_types[i], ignore.case=TRUE))
        ct_mat[i,j] <- 1
    }
  }
}

# --- Pathways (multi-label) ---
pw_terms <- c("Wnt","Hedgehog","Notch","TGF-beta/BMP","MAPK/ERK","PI3K/AKT/mTOR",
  "cAMP/PKA","Calcium","NF-kB","JAK/STAT","Hippo/YAP")
pw_mat <- matrix(0, nrow=nrow(agent), ncol=length(pw_terms))
colnames(pw_mat) <- paste0("AgPW:", pw_terms)
rownames(pw_mat) <- agent$gene
for (i in 1:nrow(agent)) {
  if (!is.na(agent$pathways[i]) && agent$pathways[i] != "NA") {
    for (j in seq_along(pw_terms)) {
      if (grepl(pw_terms[j], agent$pathways[i], ignore.case=TRUE, fixed=FALSE))
        pw_mat[i,j] <- 1
    }
  }
}

# --- KO severity (one-hot) ---
sev_levels <- c("embryonic lethal","perinatal lethal","postnatal lethal","subviable","viable","unknown")
agent$ko_sev_clean <- agent$ko_severity
agent$ko_sev_clean[is.na(agent$ko_sev_clean) | agent$ko_sev_clean == "NA"] <- "unknown"
agent$ko_sev_clean[!agent$ko_sev_clean %in% sev_levels] <- "unknown"
sev_mat <- model.matrix(~ 0 + factor(agent$ko_sev_clean, levels=sev_levels))
colnames(sev_mat) <- paste0("AgKO:", sev_levels)
rownames(sev_mat) <- agent$gene

# --- Cancer role (one-hot) ---
cr_levels <- c("oncogene","tumor suppressor","biomarker only","none","unknown")
agent$cr_clean <- agent$cancer_role
agent$cr_clean[is.na(agent$cr_clean) | agent$cr_clean == "NA"] <- "unknown"
agent$cr_clean[!agent$cr_clean %in% cr_levels] <- "unknown"
cr_mat <- model.matrix(~ 0 + factor(agent$cr_clean, levels=cr_levels))
colnames(cr_mat) <- paste0("AgCR:", cr_levels)
rownames(cr_mat) <- agent$gene

# --- Differentiation role (one-hot) ---
dr_levels <- c("promotes","inhibits","required for","marker of","associated","none known")
agent$dr_clean <- agent$differentiation_role
agent$dr_clean[is.na(agent$dr_clean) | agent$dr_clean == "NA"] <- "none known"
agent$dr_clean[!agent$dr_clean %in% dr_levels] <- "none known"
dr_mat <- model.matrix(~ 0 + factor(agent$dr_clean, levels=dr_levels))
colnames(dr_mat) <- paste0("AgDiff:", dr_levels)
rownames(dr_mat) <- agent$gene

# --- Cerebellar evidence (binary) ---
ce_mat <- matrix(0, nrow=nrow(agent), ncol=1)
colnames(ce_mat) <- "AgCereb:evidence"
rownames(ce_mat) <- agent$gene
ce <- agent$cerebellar_evidence
ce_mat[,1] <- as.integer(!is.na(ce) & ce != "NA" & ce != "" & ce != "no" & ce != "No")

# --- Neural phenotype (binary) ---
np_mat <- matrix(0, nrow=nrow(agent), ncol=1)
colnames(np_mat) <- "AgNeuro:has_phenotype"
rownames(np_mat) <- agent$gene
np <- agent$neural_phenotype
np_mat[,1] <- as.integer(!is.na(np) & np != "NA" & np != "" & np != "none" & np != "None")

# Combine all agent features
agent_full <- cbind(pc_mat, subcel_mat, tissue_mat, ct_mat, pw_mat, sev_mat, cr_mat, dr_mat, ce_mat, np_mat)
cat("Agent feature matrix:", nrow(agent_full), "x", ncol(agent_full), "\n")

# =============================================================================
# 3. MERGE WITH DATABASE MATRIX — UNION OF GENES (pad zeros for missing)
# =============================================================================
cat("\nMerging matrices (union of genes)...\n")

# Union of all genes
all_gene_set <- union(rownames(mat_ann), rownames(agent_full))
cat("Union of genes:", length(all_gene_set), "\n")
cat("  In both db + agent:", length(intersect(rownames(mat_ann), rownames(agent_full))), "\n")
cat("  Only in db (will pad agent cols with 0):", length(setdiff(rownames(mat_ann), rownames(agent_full))), "\n")
cat("  Only in agent (will pad db cols with 0):", length(setdiff(rownames(agent_full), rownames(mat_ann))), "\n")

# Pad mat_ann with zero rows for genes not present
missing_from_db <- setdiff(all_gene_set, rownames(mat_ann))
if (length(missing_from_db) > 0) {
  pad_db <- matrix(0, nrow = length(missing_from_db), ncol = ncol(mat_ann))
  rownames(pad_db) <- missing_from_db
  colnames(pad_db) <- colnames(mat_ann)
  mat_ann_full <- rbind(mat_ann, pad_db)
} else {
  mat_ann_full <- mat_ann
}

# Pad agent_full with zero rows for genes not present
missing_from_agent <- setdiff(all_gene_set, rownames(agent_full))
if (length(missing_from_agent) > 0) {
  pad_agent <- matrix(0, nrow = length(missing_from_agent), ncol = ncol(agent_full))
  rownames(pad_agent) <- missing_from_agent
  colnames(pad_agent) <- colnames(agent_full)
  agent_full_padded <- rbind(agent_full, pad_agent)
} else {
  agent_full_padded <- agent_full
}

# Align rows and combine
mat_combined <- cbind(mat_ann_full[all_gene_set, ], agent_full_padded[all_gene_set, ])
cat("Combined matrix:", nrow(mat_combined), "x", ncol(mat_combined), "\n")

# Remove zero-variance columns
col_var <- apply(mat_combined, 2, var)
mat_combined <- mat_combined[, col_var > 0]
cat("After zero-var removal:", ncol(mat_combined), "\n")

# Genes with near-zero total features may still be placed poorly
n_features_per_gene <- rowSums(mat_combined)
cat("Genes with 0 total features:", sum(n_features_per_gene == 0), "\n")
cat("Median features per gene:", median(n_features_per_gene), "\n")

# =============================================================================
# 4. PCA + UMAP
# =============================================================================
cat("\nRunning PCA...\n")
n_pcs <- min(100, ncol(mat_combined)-1, nrow(mat_combined)-1)
pca <- prcomp_irlba(mat_combined, n=n_pcs, center=TRUE, scale.=FALSE)
var_exp <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
n_use <- min(which(var_exp >= 0.80))
if (is.na(n_use)) n_use <- n_pcs
cat("Using", n_use, "PCs (", round(var_exp[n_use]*100,1), "% variance)\n")

cat("Running UMAP...\n")
set.seed(42)
umap_res <- umap(pca$x[, 1:n_use], n_neighbors=15, min_dist=0.15,
  metric="cosine", n_components=2, n_epochs=1000, spread=2.0)

# =============================================================================
# 5. CLUSTERING (k=13)
# =============================================================================
K_CLUSTERS <- 10
SHOW_TITLE <- FALSE  # Set TRUE to include title/subtitle in main figure
cat("Clustering k=", K_CLUSTERS, "...\n")
set.seed(42)
km <- kmeans(umap_res, centers=K_CLUSTERS, nstart=50, iter.max=300)
sil <- round(mean(silhouette(km$cluster, dist(umap_res))[,3]), 3)
cat("Silhouette:", sil, "\n")

umap_df <- data.frame(
  gene=rownames(mat_combined), UMAP1=umap_res[,1], UMAP2=umap_res[,2],
  cluster=km$cluster, stringsAsFactors=FALSE
) %>%
  dplyr::left_join(diff_genes, by=c("gene"="mgi_symbol")) %>%
  dplyr::left_join(agent %>% dplyr::select(gene, protein_class, one_line_summary), by="gene")

# Add any truly unannotated genes (0 features in combined matrix) at periphery
all_syms <- diff_genes$mgi_symbol
missing <- setdiff(all_syms, rownames(mat_combined))
if (length(missing) > 0) {
  xr <- range(umap_df$UMAP1); yr <- range(umap_df$UMAP2)
  set.seed(42)
  umap_df <- rbind(umap_df, data.frame(
    gene=missing, UMAP1=runif(length(missing), xr[2]+1, xr[2]+3),
    UMAP2=runif(length(missing), yr[1], yr[2]),
    cluster=0,
    H3K27me3_bin=diff_genes$H3K27me3_bin[match(missing, diff_genes$mgi_symbol)],
    t50_cat=diff_genes$t50_cat[match(missing, diff_genes$mgi_symbol)],
    protein_class=NA, one_line_summary=NA,
    stringsAsFactors=FALSE))
}

cat("\nCluster sizes:\n")
print(sort(table(umap_df$cluster[umap_df$cluster>0]), decreasing=TRUE))

# =============================================================================
# 6. LABEL CLUSTERS — using agent protein_class + Fisher enrichment
# =============================================================================
cat("\nLabeling clusters...\n")

# Use combined matrix for Fisher enrichment
all_col_sums_c <- colSums(mat_combined)
n_total <- nrow(mat_combined)

clean_term <- function(t) {
  t <- gsub("^(GO_BP:|GO_MF:|GO_CC:|MSig_[^:]+:|KEGG:|Reactome:|TF_family:|PFAM:|Desc:|Feature:|DescWord:)", "", t)
  t <- gsub("^(Ag[A-Z]+:)", "", t)  # Strip agent prefixes for display
  repeat {
    t_old <- t
    t <- gsub("^(GOBP[_ ]|GOMF[_ ]|GOCC[_ ])", "", t)
    t <- gsub("^(GO:MF:GOMF[_ ]?|GO:BP:GOBP[_ ]?|GO:CC:GOCC[_ ]?)", "", t)
    t <- gsub("^(HALLMARK[_ ]|HP[_ ]|HPO:HP[_ ])", "", t)
    t <- gsub("^(REACTOME[_ ]|CP:REACTOME:REACTOME[_ ])", "", t)
    t <- gsub("^(WP[_ ]|KEGG[_ ]?LEGACY[_ ]?|PID[_ ]|BIOCARTA[_ ])", "", t)
    t <- gsub("^(TFT:GTRD:|TFT:TFT[_ ]?LEGACY:)", "", t)
    t <- gsub("^(CGN[_ ]|MANNO[_ ]|GAUTAM[_ ]|DESCARTES[_ ])", "", t)
    t <- gsub("^(CP:WIKIPATHWAYS:WP[_ ]?|CP:PID:PID[_ ]?|CP:BIOCARTA:BIOCARTA[_ ]?|CP:KEGG_LEGACY:KEGG[_ ]?)", "", t)
    if (t == t_old) break
  }
  t <- gsub("_", " ", t); t <- gsub("  +", " ", trimws(t))
  if (nchar(t) > 0) t <- paste0(toupper(substr(t,1,1)), tolower(substr(t,2,nchar(t))))
  t
}

generic_roots <- tolower(c(
  # Generic GO roots
  "molecular function","cellular component","biological process",
  "binding","catalytic activity","protein binding","cell","cell part",
  "organelle","intracellular","cytoplasm","membrane","nucleus",
  "cellular process","metabolic process","system development",
  "anatomical structure development","regulation of biological process",
  "ion binding","metal ion binding","cation binding","transferase activity",
  "hydrolase activity","signaling","signal transduction",
  "response to stimulus","developmental process","multicellular organismal process",
  "extracellular region","extracellular space",
  # Agent classification labels (not biological terms)
  "other","unknown","none","none known","viable","subviable",
  "embryonic lethal","perinatal lethal","postnatal lethal",
  "biomarker only","oncogene","tumor suppressor",
  "promotes","inhibits","required for","marker of","associated",
  "has phenotype","evidence",
  # Meaningless / overly broad
  "predicted uncharacterized","positive regulation","negative regulation",
  "regulation","broad","protein","domain","family","activity",
  "development","differentiation","process","function"
))

# Agent-derived feature columns (one-hot / categorical) — don't use as labels
is_agent_feature_col <- function(tname) {
  grepl("^(AgPC:|AgLoc:|AgTiss:|AgCell:|AgPW:|AgKO:|AgCR:|AgDiff:|AgCereb:|AgNeuro:)",
    tname)
}

# Pfam IDs (PF00001 etc) — not informative without description (vectorized)
is_pfam_id_only <- function(tname) {
  cleaned <- gsub("^PFAM:", "", tname)
  grepl("^PF[0-9]+$", cleaned) | grepl("^Pf[0-9]+$", cleaned)
}

centroids <- umap_df %>% dplyr::filter(cluster>0) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x=median(UMAP1), y=median(UMAP2), n=n(), .groups="drop")

cluster_labels <- character(K_CLUSTERS)
cluster_term_lists <- list()

for (cl in 1:K_CLUSTERS) {
  cl_genes <- umap_df$gene[umap_df$cluster==cl]
  cl_in <- intersect(cl_genes, rownames(mat_combined))
  n_cl <- length(cl_in)
  if (n_cl < 3) { cluster_labels[cl] <- "Small"; next }

  # Top protein classes in this cluster
  cl_agent <- agent %>% dplyr::filter(gene %in% cl_genes)
  top_pc <- sort(table(cl_agent$protein_class), decreasing=TRUE)
  pc_label <- paste(paste0(names(top_pc)[1:min(2,length(top_pc))], " (",
    top_pc[1:min(2,length(top_pc))], ")"), collapse=", ")

  # Fisher enrichment for top terms
  cl_mat <- mat_combined[cl_in, , drop=FALSE]
  in_cl <- colSums(cl_mat)
  cands <- which(in_cl >= 3)

  if (length(cands) > 0) {
    fisher_p <- numeric(length(cands))
    for (idx in seq_along(cands)) {
      j <- cands[idx]
      a <- in_cl[j]; b <- n_cl-a
      c_v <- all_col_sums_c[j]-a; d <- (n_total-n_cl)-c_v
      fisher_p[idx] <- fisher.test(matrix(c(a,b,c_v,d), nrow=2), alternative="greater")$p.value
    }
    res <- data.frame(term=names(cands), fisher_p=fisher_p, in_cl=in_cl[cands],
      pct=round(100*in_cl[cands]/n_cl,1), stringsAsFactors=FALSE) %>%
      dplyr::mutate(fdr=p.adjust(fisher_p, method="BH"),
        term_clean=sapply(term, clean_term), term_lower=tolower(term_clean),
        is_ag_feat=is_agent_feature_col(term),
        is_pfam_only=is_pfam_id_only(term),
        n_words=sapply(term_clean, function(x) length(strsplit(x," ")[[1]]))) %>%
      dplyr::filter(fdr < 0.05,
        !term_lower %in% generic_roots,
        !is_ag_feat,           # Exclude agent feature columns from labels
        !is_pfam_only,         # Exclude bare Pfam IDs
        n_words >= 2,          # Require >=2 words (avoid "Viable", "Other")
        nchar(term_clean) >= 6) %>%
      dplyr::filter(!duplicated(term_lower)) %>%
      dplyr::arrange(fisher_p)
    cluster_term_lists[[cl]] <- res %>% dplyr::slice_head(n=4)
  } else {
    cluster_term_lists[[cl]] <- data.frame()
  }

  # Build label from top enriched + protein class
  if (nrow(cluster_term_lists[[cl]]) > 0) {
    top2 <- head(cluster_term_lists[[cl]]$term_clean, 2)
    short <- paste(top2, collapse=" / ")
  } else {
    short <- pc_label
  }
  if (nchar(short) > 50) short <- paste0(substr(short, 1, 47), "...")

  k27 <- round(100*mean(umap_df$H3K27me3_bin[umap_df$cluster==cl]==1, na.rm=TRUE))
  cluster_labels[cl] <- short
  cat(sprintf("  Cl %d (%d genes, %d%% K27+): %s | Classes: %s\n", cl, n_cl, k27, short, pc_label))
}

centroids$label <- cluster_labels[centroids$cluster]

# =============================================================================
# 7. FIGURES
# =============================================================================
cat("\nGenerating figures...\n")

# Term annotations — ONE BOX PER CLUSTER with terms stacked newline-separated,
# placed outward from centroid toward plot edge
plot_cx <- mean(range(umap_df$UMAP1[umap_df$cluster>0]))
plot_cy <- mean(range(umap_df$UMAP2[umap_df$cluster>0]))

term_annot_list <- list()
for (cl in 1:K_CLUSTERS) {
  ct <- centroids[centroids$cluster==cl,]
  terms_df <- cluster_term_lists[[cl]]
  if (is.null(terms_df) || nrow(terms_df)==0) next

  # Stack up to 4 terms in a single multi-line string
  lines <- character()
  for (j in 1:min(4, nrow(terms_df))) {
    d <- terms_df$term_clean[j]
    if (nchar(d) > 30) d <- paste0(substr(d, 1, 27), "...")
    lines <- c(lines, paste0(d, " (", terms_df$pct[j], "%)"))
  }
  label_txt <- paste(lines, collapse="\n")

  # Offset outward from plot center
  dx <- ct$x - plot_cx; dy <- ct$y - plot_cy
  dist_fc <- sqrt(dx^2 + dy^2)
  if (dist_fc < 0.01) { dx <- 1; dy <- 0; dist_fc <- 1 }
  dx_n <- dx / dist_fc; dy_n <- dy / dist_fc
  offset <- 3.2

  term_annot_list[[as.character(cl)]] <- data.frame(
    cluster = cl, label = label_txt,
    x_anchor = ct$x, y_anchor = ct$y,
    x_label = ct$x + dx_n * offset,
    y_label = ct$y + dy_n * offset,
    stringsAsFactors = FALSE)
}
term_annot <- do.call(rbind, term_annot_list)

pal10 <- scales::hue_pal()(K_CLUSTERS)
names(pal10) <- as.character(1:K_CLUSTERS)
umap_df$cluster_f <- factor(umap_df$cluster)

# Main UMAP — single box per cluster, no radiating lines
p_main <- ggplot(umap_df %>% dplyr::filter(cluster>0),
    aes(x=UMAP1, y=UMAP2, color=cluster_f)) +
  geom_point(size=1.6, alpha=0.55) +
  {if(any(umap_df$cluster==0))
    geom_point(data=umap_df %>% filter(cluster==0), color="grey60",
      size=0.8, alpha=0.3, inherit.aes=FALSE, aes(x=UMAP1,y=UMAP2))} +
  # Cluster number at centroid
  geom_label(data=centroids,
    aes(x=x, y=y, label=paste0(cluster," (n=",n,")")),
    size=4, fontface="bold", fill=alpha("white",0.85),
    label.size=0.3, color="black", inherit.aes=FALSE) +
  # Single term box per cluster, terms stacked newline-separated
  geom_label(data=term_annot,
    aes(x=x_label, y=y_label, label=label),
    size=2.8, color="grey15", hjust=0.5, vjust=0.5,
    fontface="italic", fill=alpha("white", 0.9),
    label.size=0.2, label.padding=unit(0.25, "lines"),
    inherit.aes=FALSE, lineheight=0.95) +
  scale_color_manual(values=pal10, guide="none") +
  coord_cartesian(clip="off") +
  {if (SHOW_TITLE)
    labs(title="Functional landscape of GNP differentiation genes",
      subtitle=paste0(sum(umap_df$cluster>0), " genes | ", ncol(mat_combined),
        " features (database + agent-extracted) | k=", K_CLUSTERS, " | sil=", sil))
   else labs(title=NULL, subtitle=NULL)} +
  theme_void(base_size=14) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5),
    plot.subtitle=element_text(size=11, face="italic", hjust=0.5),
    plot.margin=margin(20, 90, 20, 90, "pt"))

ggsave("Fig_gene_umap_agent_k10.pdf", p_main, width=20, height=16)
cat("Saved: Fig_gene_umap_agent_k10.pdf\n")

# H3K27me3 overlay
p_k27 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2,
    color=factor(H3K27me3_bin, labels=c("Unbound","H3K27me3-bound")))) +
  geom_point(size=1.5, alpha=0.5) +
  scale_color_manual(values=c("Unbound"="grey75","H3K27me3-bound"="#4393C3"),
    name=NULL, na.value="grey90") +
  labs(title="H3K27me3 binding status") +
  theme_void(base_size=14) +
  theme(plot.title=element_text(size=15, face="bold", hjust=0.5),
    legend.position=c(0.85,0.1), legend.text=element_text(size=12))

ggsave("Fig_gene_umap_agent_k27.pdf", p_k27, width=10, height=8)
cat("Saved: Fig_gene_umap_agent_k27.pdf\n")

# Protein class overlay
pc_cols <- c("ion channel"="#E41A1C","transporter"="#377EB8","GPCR"="#4DAF4A",
  "receptor"="#984EA3","kinase"="#FF7F00","transcription factor"="#A65628",
  "enzyme"="#F781BF","structural protein"="#999999","adhesion molecule"="#66C2A5",
  "scaffold/adaptor"="#FC8D62","secreted factor"="#8DA0CB",
  "GTPase/regulator"="#E78AC3","channel regulator"="#A6D854",
  "ECM protein"="#FFD92F","ubiquitin ligase"="#E5C494",
  "protease"="#B3B3B3","phosphatase"="#FDDAEC","motor protein"="#CBD5E8","other"="grey80")

p_pc <- ggplot(umap_df %>% dplyr::filter(cluster>0, !is.na(protein_class)),
    aes(x=UMAP1, y=UMAP2, color=protein_class)) +
  geom_point(size=1.5, alpha=0.6) +
  scale_color_manual(values=pc_cols, name="Protein class") +
  labs(title="Protein class (agent-extracted)") +
  theme_void(base_size=14) +
  theme(plot.title=element_text(size=15, face="bold", hjust=0.5),
    legend.position="right", legend.text=element_text(size=9)) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1), ncol=1))

ggsave("Fig_gene_umap_agent_proteinclass.pdf", p_pc, width=14, height=10)
cat("Saved: Fig_gene_umap_agent_proteinclass.pdf\n")

# Save coordinates
write.csv(umap_df, "gene_umap_agent_k10_coordinates.csv", row.names=FALSE)
cat("Saved: gene_umap_agent_k10_coordinates.csv\n")

cat("\n=== DONE ===\n")
