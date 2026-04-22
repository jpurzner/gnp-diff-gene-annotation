suppressMessages({
  library(tidyverse)
  library(cluster)
  library(patchwork)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# Load UMAP coordinates (excluding peripheral cluster 0)
umap_df <- read.csv("gene_umap_agent_k10_coordinates.csv", stringsAsFactors=FALSE)
umap_ann <- umap_df %>% filter(cluster > 0)
umap_mat <- as.matrix(umap_ann[, c("UMAP1","UMAP2")])

cat("Computing elbow / silhouette for k = 2 to 20 on", nrow(umap_mat), "genes...\n")

wss <- numeric(); sil <- numeric()
for (k in 2:20) {
  set.seed(42)
  km <- kmeans(umap_mat, centers=k, nstart=50, iter.max=300)
  wss[k-1] <- km$tot.withinss
  sil[k-1] <- mean(silhouette(km$cluster, dist(umap_mat))[,3])
  cat(sprintf("  k=%2d  WSS=%.1f  sil=%.3f\n", k, wss[k-1], sil[k-1]))
}

elbow_df <- data.frame(k=2:20, wss=wss, sil=sil)

p1 <- ggplot(elbow_df, aes(x=k, y=wss)) +
  geom_line(linewidth=0.8) + geom_point(size=2) +
  geom_vline(xintercept=10, linetype="dashed", color="red") +
  annotate("text", x=10.5, y=max(wss)*0.9,
    label="k=10 (used)", color="red", hjust=0, size=4) +
  scale_x_continuous(breaks=2:20) +
  labs(x="Number of clusters (k)",
       y="Total within-cluster sum of squares",
       title="Elbow plot (agent-enriched UMAP)",
       subtitle=paste0(nrow(umap_mat), " genes from combined database + agent features")) +
  theme_bw(base_size=13)

best_k <- elbow_df$k[which.max(elbow_df$sil)]
p2 <- ggplot(elbow_df, aes(x=k, y=sil)) +
  geom_line(linewidth=0.8) + geom_point(size=2) +
  geom_vline(xintercept=best_k, linetype="dashed", color="blue") +
  geom_vline(xintercept=10, linetype="dashed", color="red") +
  annotate("text", x=best_k+0.5, y=max(sil)*0.98,
    label=paste0("k=", best_k, " (best silhouette)"), color="blue", hjust=0, size=4) +
  annotate("text", x=10.5, y=max(sil)*0.88,
    label="k=10 (used)", color="red", hjust=0, size=4) +
  scale_x_continuous(breaks=2:20) +
  labs(x="Number of clusters (k)", y="Mean silhouette width",
       title="Silhouette score") +
  theme_bw(base_size=13)

p_combined <- p1 / p2 +
  plot_annotation(title="Cluster selection — agent-enriched UMAP",
    theme=theme(plot.title=element_text(size=15, face="bold")))

ggsave("Fig_gene_umap_agent_elbow.pdf", p_combined, width=8, height=8)
cat("\nSaved: Fig_gene_umap_agent_elbow.pdf\n")
print(elbow_df)
