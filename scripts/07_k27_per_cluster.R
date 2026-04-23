suppressMessages({library(tidyverse); library(patchwork)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

umap_df <- read.csv("gene_umap_agent_k10_coordinates.csv", stringsAsFactors=FALSE)
agent <- read.csv("agent_extracted_annotations.csv", stringsAsFactors=FALSE)

cluster_labels <- c(
  "1" = "Nucleic acid binding / Peptidase",
  "2" = "Mixed / Other (uncharacterized)",
  "3" = "Ion channels / transmembrane",
  "4" = "Phosphotransferase / kinase",
  "5" = "Thorax/limb morphogenesis",
  "6" = "Cell projection morphogenesis",
  "7" = "Small molecule / lipid metabolism",
  "8" = "Mixed / Other (scaffolds)",
  "9" = "Transcription regulation",
  "10" = "Cell-cell adhesion / tight junction")

overall_pct <- 100 * sum(umap_df$H3K27me3_bin == 1, na.rm=TRUE) /
  sum(!is.na(umap_df$H3K27me3_bin) & umap_df$cluster > 0)

# Per cluster
per_cl <- umap_df %>%
  dplyr::filter(cluster > 0, !is.na(H3K27me3_bin)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n = n(), n_k27 = sum(H3K27me3_bin == 1), .groups="drop") %>%
  dplyr::mutate(pct_k27 = 100*n_k27/n,
    label = cluster_labels[as.character(cluster)],
    cl_lab = paste0("Cl", cluster, " (n=", n, ")")) %>%
  dplyr::arrange(desc(pct_k27))

per_cl$cl_lab <- factor(per_cl$cl_lab, levels = rev(per_cl$cl_lab))

p_cl <- ggplot(per_cl, aes(x=cl_lab, y=pct_k27, fill=pct_k27)) +
  geom_col(color="grey30", linewidth=0.3, width=0.7) +
  geom_hline(yintercept=overall_pct, linetype="dashed", color="red", linewidth=0.6) +
  annotate("text", x=0.5, y=overall_pct+1.5,
    label=paste0("Overall: ", round(overall_pct,1), "%"),
    color="red", size=3.5, hjust=0, fontface="italic") +
  geom_text(aes(label=paste0(round(pct_k27,0), "%")),
    hjust=-0.15, size=4, fontface="bold") +
  geom_text(aes(label=label, y=2), hjust=0, color="white",
    size=3.2, fontface="italic") +
  scale_fill_gradient(low="#FDBB84", high="#990000", guide="none") +
  coord_flip() +
  scale_y_continuous(expand=expansion(mult=c(0, 0.15)), limits=c(0, 80)) +
  labs(x=NULL, y="% H3K27me3-bound",
    title="H3K27me3 marking by UMAP cluster",
    subtitle="Dashed line = overall rate across diff genes") +
  theme_bw(base_size=12) +
  theme(panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=14),
    plot.subtitle=element_text(face="italic", size=11, color="grey40"))

# Per protein class
df_pc <- umap_df %>%
  dplyr::filter(cluster > 0, !is.na(H3K27me3_bin)) %>%
  dplyr::select(gene, H3K27me3_bin) %>%
  dplyr::left_join(agent %>% dplyr::select(gene, protein_class), by="gene")

per_pc <- df_pc %>%
  dplyr::filter(!is.na(protein_class), protein_class != "NA") %>%
  dplyr::group_by(protein_class) %>%
  dplyr::summarise(n = n(), n_k27 = sum(H3K27me3_bin == 1), .groups="drop") %>%
  dplyr::mutate(pct_k27 = 100*n_k27/n,
    lab = paste0(protein_class, " (n=", n, ")")) %>%
  dplyr::filter(n >= 10) %>%
  dplyr::arrange(desc(pct_k27))

per_pc$lab <- factor(per_pc$lab, levels = rev(per_pc$lab))

p_pc <- ggplot(per_pc, aes(x=lab, y=pct_k27, fill=pct_k27)) +
  geom_col(color="grey30", linewidth=0.3, width=0.7) +
  geom_hline(yintercept=overall_pct, linetype="dashed", color="red", linewidth=0.6) +
  geom_text(aes(label=paste0(round(pct_k27,0), "%")),
    hjust=-0.15, size=4, fontface="bold") +
  scale_fill_gradient(low="#FDBB84", high="#990000", guide="none") +
  coord_flip() +
  scale_y_continuous(expand=expansion(mult=c(0, 0.15)), limits=c(0, 80)) +
  labs(x=NULL, y="% H3K27me3-bound",
    title="H3K27me3 marking by agent-extracted protein class",
    subtitle="Classes with n >= 10 genes; dashed line = overall rate") +
  theme_bw(base_size=12) +
  theme(panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=14),
    plot.subtitle=element_text(face="italic", size=11, color="grey40"))

# Save individual
ggsave("Fig_k27_per_cluster.pdf", p_cl, width=10, height=6)
cat("Saved: Fig_k27_per_cluster.pdf\n")
ggsave("Fig_k27_per_proteinclass.pdf", p_pc, width=10, height=8)
cat("Saved: Fig_k27_per_proteinclass.pdf\n")

# Combined
p_combined <- p_cl / p_pc + plot_layout(heights=c(1, 1.4)) +
  plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(face="bold", size=16))
ggsave("Fig_k27_enrichment_bars.pdf", p_combined, width=10, height=14)
cat("Saved: Fig_k27_enrichment_bars.pdf\n")

write.csv(per_cl, "k27_pct_per_cluster.csv", row.names=FALSE)
write.csv(per_pc, "k27_pct_per_proteinclass.csv", row.names=FALSE)
