# Usage: Rscript run_umap_figure.R PROLIF
suppressMessages({
  library(tidyverse); library(uwot); library(irlba); library(cluster)
  library(patchwork); library(ggrepel); library(ggforce); library(scales)
})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

args <- commandArgs(trailingOnly=TRUE)
SET <- args[1]
if (!SET %in% c("PROLIF","EMB")) stop("Set must be PROLIF or EMB")

if (SET == "PROLIF") {
  combined <- read.csv("combined_annotations_PROLIF.csv", stringsAsFactors=FALSE)
  set_name <- "Proliferation"
  out_pdf <- "Fig_umap_PROLIF.pdf"
} else {
  combined <- read.csv("combined_annotations_EMB.csv", stringsAsFactors=FALSE)
  set_name <- "Embryonic"
  out_pdf <- "Fig_umap_EMB.pdf"
}

cat("=== UMAP FIGURE for", set_name, "===\n")

K_CLUSTERS <- 10
N_EXAMPLES <- 3

# UMAP coords + cluster
umap_df <- combined %>%
  dplyr::filter(!is.na(UMAP1), !is.na(cluster), cluster > 0) %>%
  dplyr::mutate(cluster_f = factor(cluster))

# Centroids
centroids <- umap_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x=median(UMAP1), y=median(UMAP2), n=n(), .groups="drop")

# Per-cluster top enriched terms — use the protein_class + cell_types + pathways
# columns from the agent data as a simple "vocabulary" for labels
cluster_term_lists <- list()
for (cl in 1:K_CLUSTERS) {
  sub <- umap_df %>% dplyr::filter(cluster == cl)
  if (nrow(sub) < 3) { cluster_term_lists[[cl]] <- data.frame(); next }
  pcs <- sort(table(sub$protein_class), decreasing=TRUE)
  pcs <- pcs[!names(pcs) %in% c("other","NA","")]
  cell_terms <- unlist(strsplit(paste(sub$cell_types, collapse="; "), "; "))
  cell_terms <- cell_terms[!is.na(cell_terms) & cell_terms != "" & cell_terms != "NA"]
  cells <- sort(table(cell_terms), decreasing=TRUE)
  pw_terms <- unlist(strsplit(paste(sub$pathways, collapse="; "), "; "))
  pw_terms <- pw_terms[!is.na(pw_terms) & pw_terms != "" & pw_terms != "NA"]
  pws <- sort(table(pw_terms), decreasing=TRUE)

  top <- character()
  if (length(pcs) > 0) top <- c(top, paste0(names(pcs)[1:min(2,length(pcs))],
    " (", round(100*pcs[1:min(2,length(pcs))]/nrow(sub)), "%)"))
  if (length(cells) > 0) top <- c(top, paste0(names(cells)[1:min(1,length(cells))],
    " cells (", round(100*cells[1:min(1,length(cells))]/nrow(sub)), "%)"))
  if (length(pws) > 0) top <- c(top, paste0(names(pws)[1:min(1,length(pws))],
    " (", round(100*pws[1:min(1,length(pws))]/nrow(sub)), "%)"))
  cluster_term_lists[[cl]] <- data.frame(label=paste(top, collapse="\n"),
    stringsAsFactors=FALSE)
}

# Hull annotation data
hull_annot <- centroids %>%
  dplyr::mutate(
    cluster_id = paste0("Cl", cluster, " (n=", n, ")"),
    label = sapply(cluster, function(cl) {
      tl <- cluster_term_lists[[cl]]
      if (nrow(tl) == 0) "Mixed / Other" else tl$label[1]
    }))

umap_df <- umap_df %>%
  dplyr::left_join(hull_annot %>% dplyr::select(cluster, cluster_id, label), by="cluster")

# Dominant cluster fill
pal10 <- scales::hue_pal()(K_CLUSTERS)
names(pal10) <- as.character(1:K_CLUSTERS)

# Protein class palette
pc_cols <- c("ion channel"="#E41A1C","transporter"="#377EB8","GPCR"="#4DAF4A",
  "receptor"="#984EA3","kinase"="#FF7F00","transcription factor"="#A65628",
  "enzyme"="#F781BF","structural protein"="#999999","adhesion molecule"="#66C2A5",
  "scaffold/adaptor"="#FC8D62","secreted factor"="#8DA0CB",
  "GTPase/regulator"="#E78AC3","channel regulator"="#A6D854",
  "ECM protein"="#FFD92F","ubiquitin ligase"="#E5C494",
  "protease"="#B3B3B3","phosphatase"="#FDDAEC","motor protein"="#CBD5E8",
  "other"="grey80")

umap_df$protein_class[is.na(umap_df$protein_class) | umap_df$protein_class == "NA"] <- "other"

# Example genes per cluster
is_predicted <- function(g) grepl("^(Gm[0-9]+|[0-9]{4}[A-Z][0-9]+Rik|[0-9]+[A-Z]{1,2}[0-9]+Rik|BC[0-9]+|AI[0-9]+|LOC[0-9]+)$", g)
example_list <- list()
for (cl in 1:K_CLUSTERS) {
  sub <- umap_df %>% dplyr::filter(cluster == cl)
  ct <- centroids[centroids$cluster == cl, ]
  sub$dist <- sqrt((sub$UMAP1 - ct$x)^2 + (sub$UMAP2 - ct$y)^2)
  sub$pred <- sapply(sub$mgi_symbol, is_predicted)
  picked <- sub %>% dplyr::filter(!pred,
    !is.na(one_line_summary), nchar(one_line_summary) > 30) %>%
    dplyr::arrange(dist) %>% dplyr::slice_head(n=N_EXAMPLES)
  if (nrow(picked) < N_EXAMPLES) {
    more <- sub %>% dplyr::filter(!mgi_symbol %in% picked$mgi_symbol) %>%
      dplyr::arrange(dist) %>% dplyr::slice_head(n=N_EXAMPLES - nrow(picked))
    picked <- rbind(picked, more)
  }
  example_list[[cl]] <- picked %>% dplyr::select(mgi_symbol, UMAP1, UMAP2, cluster)
}
examples_df <- do.call(rbind, example_list)
colnames(examples_df)[1] <- "gene"

# Side padding so labels go left/right
x_range <- range(umap_df$UMAP1); y_range <- range(umap_df$UMAP2)
x_pad <- 0.55 * diff(x_range); y_pad <- 0.10 * diff(y_range)

p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  ggforce::geom_mark_hull(
    aes(fill=cluster_f, group=cluster, label=cluster_id, description=label),
    colour="grey40", linewidth=0.5, alpha=0.22,
    expand=unit(4,"mm"), radius=unit(4,"mm"), concavity=2,
    label.fontsize=c(28, 22),
    label.fontface=c("bold","italic"),
    label.margin=margin(4,5,4,5,"mm"),
    label.colour="grey10",
    label.buffer=unit(8,"mm"),
    label.width=unit(140,"mm"),
    label.minwidth=unit(110,"mm"),
    con.colour="grey40", con.size=0.7,
    con.cap=unit(2,"mm"), con.type="elbow") +
  geom_point(aes(color=protein_class), size=3.0, alpha=0.85) +
  geom_text_repel(data=examples_df, aes(x=UMAP1, y=UMAP2, label=gene),
    size=8, color="grey10", fontface="bold",
    bg.color="white", bg.r=0.15, force=3, force_pull=0.5,
    box.padding=0.3, point.padding=0.15,
    min.segment.length=Inf, max.overlaps=Inf, seed=42,
    inherit.aes=FALSE) +
  scale_fill_manual(values=pal10, guide="none") +
  scale_color_manual(values=pc_cols, name="Protein class") +
  coord_cartesian(xlim=c(x_range[1]-x_pad, x_range[2]+x_pad),
                  ylim=c(y_range[1]-y_pad, y_range[2]+y_pad), clip="off") +
  guides(color=guide_legend(override.aes=list(size=8, alpha=1), ncol=1)) +
  labs(title=NULL, subtitle=NULL) +
  theme_void(base_size=24) +
  theme(plot.margin=margin(60,60,60,60,"pt"),
        legend.position="right",
        legend.text=element_text(size=26),
        legend.title=element_text(size=28, face="bold"),
        legend.key.size=unit(1.6,"cm"),
        legend.key.spacing.y=unit(0.3,"cm"))

ggsave(out_pdf, p, width=42, height=22)
cat("Saved:", out_pdf, "\n")
