suppressMessages(library(tidyverse))
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# Load saved data
cat("Loading saved data...\n")
all_gene_pmids <- readRDS("tiered_pmid_assignments.rds")
all_abstract_text <- readRDS("tiered_abstracts_raw.rds")

master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude))

orth <- read.table("orth_mouse_filt_best.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="\"")
orth_map <- orth %>%
  dplyr::filter(mmusculus_homolog_ensembl_gene %in% diff_genes$ensembl_gene_id) %>%
  dplyr::select(mmusculus_homolog_ensembl_gene, hgnc_symbol) %>%
  dplyr::rename(human_symbol = hgnc_symbol) %>%
  dplyr::filter(!duplicated(mmusculus_homolog_ensembl_gene))

diff_genes <- diff_genes %>%
  dplyr::left_join(orth_map, by = c("ensembl_gene_id" = "mmusculus_homolog_ensembl_gene"))

symbols <- diff_genes$mgi_symbol
human_symbols <- setNames(diff_genes$human_symbol, diff_genes$mgi_symbol)

tiers <- c("Tier1_cerebellum_MB_diff", "Tier2_brain_neuron", "Tier3_brain_cancer",
           "Tier4_cancer", "Tier5_development", "Tier6_general")

cat("Genes:", length(symbols), "\n")
cat("PMIDs in search results:", length(all_gene_pmids), "\n")
cat("Abstracts fetched:", length(all_abstract_text), "\n")

# =============================================================================
# COMPILE PER-GENE
# =============================================================================
cat("\nCompiling per-gene abstract text...\n")

gene_compiled <- list()

for (sym in symbols) {
  pmid_df <- all_gene_pmids[[sym]]
  if (is.null(pmid_df) || nrow(pmid_df) == 0) {
    gene_compiled[[sym]] <- data.frame(
      symbol = sym, n_papers = 0,
      n_tier1 = 0, n_tier2 = 0, n_tier3 = 0, n_tier4 = 0, n_tier5 = 0, n_tier6 = 0,
      combined_text = NA_character_, titles = NA_character_, years = NA_character_,
      stringsAsFactors = FALSE)
    next
  }

  texts <- character(); titles <- character(); years <- character()
  for (j in 1:nrow(pmid_df)) {
    pmid <- pmid_df$pmid[j]
    if (pmid %in% names(all_abstract_text)) {
      art <- all_abstract_text[[pmid]]
      if (!is.null(art$abstract) && nchar(art$abstract) > 10) {
        texts <- c(texts, art$abstract)
        titles <- c(titles, art$title)
        years <- c(years, ifelse(is.na(art$year), "NA", art$year))
      }
    }
  }

  tier_counts <- table(factor(pmid_df$tier, levels = tiers))

  gene_compiled[[sym]] <- data.frame(
    symbol = sym, n_papers = length(texts),
    n_tier1 = as.integer(tier_counts[1]), n_tier2 = as.integer(tier_counts[2]),
    n_tier3 = as.integer(tier_counts[3]), n_tier4 = as.integer(tier_counts[4]),
    n_tier5 = as.integer(tier_counts[5]), n_tier6 = as.integer(tier_counts[6]),
    combined_text = paste(texts, collapse = " ||| "),
    titles = paste(titles, collapse = " ||| "),
    years = paste(years, collapse = ", "),
    stringsAsFactors = FALSE)
}

compiled_df <- do.call(rbind, gene_compiled)
cat("Genes with compiled text:", sum(compiled_df$n_papers > 0), "\n")

# =============================================================================
# EXTRACT KEYWORDS
# =============================================================================
cat("\nExtracting functional keywords...\n")

category_patterns <- list(
  "ion_channel" = "ion channel|voltage.gated|potassium channel|sodium channel|calcium channel|conductance|electrophysiol|gating|pore.forming|channelopathy",
  "transporter" = "transporter|solute carrier|SLC|pump|exchanger|uptake mechanism|efflux|influx|substrate transport",
  "synapse" = "synap|neurotransmit|synaptic vesicle|exocytos|presynap|postsynap|synaptic plast|long.term potent|LTP|LTD|glutamate release|GABA release",
  "receptor_signaling" = "receptor|GPCR|G.protein.coupled|ligand.binding|agonist|antagonist|receptor tyrosine",
  "kinase_phosphorylation" = "kinase|phosphoryl|phosphatase|MAPK|ERK|AKT|PI3K|mTOR|PKC|PKA|CaMK",
  "transcription_regulation" = "transcription factor|transcriptional|DNA.binding|promoter|enhancer|gene expression|chromatin|histone|epigenetic|methylation|acetylation",
  "cell_adhesion_ecm" = "adhesion|cadherin|integrin|junction|extracellular matrix|ECM|laminin|collagen|proteoglycan|cell.cell contact",
  "cytoskeleton_morphology" = "cytoskelet|actin|tubulin|microtubul|migration|motil|morpholog|growth cone|axon guidance|dendrite|spine|neurite",
  "lipid_metabolism" = "lipid|fatty acid|sphingo|phospholipid|cholesterol|ceramide|lipase|phospholipase|biosynthesis of lipid",
  "calcium_signaling" = "calcium signal|calcium.dependent|calmodulin|Ca2\\+|calcium influx|calcium release|calcium homeostasis|calcium store",
  "cerebellar" = "cerebell|granule neuron|granule cell|EGL|IGL|parallel fiber|climbing fiber|cerebellar development",
  "neuronal_differentiation" = "neuronal differentiation|neurite outgrowth|neurogenesis|neural progenitor|neuron maturation|neuronal migration|radial migration",
  "myelination" = "myelin|oligodendrocyte|schwann|demyelinat|white matter|nerve conduction",
  "neurodegenerative" = "neurodegenerat|Alzheimer|Parkinson|Huntington|ALS|ataxia|dementia|neuronal death|neuroprotect",
  "seizure_epilepsy" = "seizure|epilep|convuls|absence seizure|febrile|hyperexcitab|SUDEP",
  "medulloblastoma" = "medulloblastoma|MB subgroup|SHH.MB|WNT.MB|Group 3|Group 4|cerebellar tumor",
  "proliferation_cell_cycle" = "proliferat|cell cycle|mitosis|G1.S|G2.M|cyclin|CDK|cell division|checkpoint",
  "apoptosis_survival" = "apoptosis|cell death|programmed death|caspase|Bcl.2|survival signal|anti.apoptot|pro.apoptot",
  "immune_inflammation" = "immune|inflammat|cytokine|macrophage|T.cell|B.cell|NF.kB|interferon|toll.like",
  "wnt_hedgehog_notch" = "Wnt|hedgehog|Shh|Sonic|Notch|BMP|TGF.beta|Smoothened|Gli|beta.catenin",
  "vesicle_trafficking" = "vesicle|endosom|exosom|trafficking|sorting|secretory|Golgi|endoplasmic reticulum|autophagy"
)

compiled_df$functional_keywords <- NA_character_

for (i in 1:nrow(compiled_df)) {
  txt <- compiled_df$combined_text[i]
  if (is.na(txt) || nchar(txt) < 50) next
  matched <- character()
  for (cat_name in names(category_patterns)) {
    if (grepl(category_patterns[[cat_name]], txt, ignore.case = TRUE))
      matched <- c(matched, cat_name)
  }
  compiled_df$functional_keywords[i] <- paste(matched, collapse = "; ")
}

# =============================================================================
# EXTRACT KEY SENTENCES
# =============================================================================
cat("Extracting key functional sentences...\n")

compiled_df$functional_summary <- NA_character_

for (i in 1:nrow(compiled_df)) {
  sym <- compiled_df$symbol[i]
  txt <- compiled_df$combined_text[i]
  if (is.na(txt) || nchar(txt) < 100) next

  sentences <- unlist(strsplit(txt, "(?<=[.!?])\\s+", perl = TRUE))

  hsym <- human_symbols[sym]
  gene_pattern <- paste0("\\b", sym, "\\b")
  if (!is.na(hsym) && hsym != "" && hsym != sym)
    gene_pattern <- paste0("(\\b", sym, "\\b|\\b", hsym, "\\b)")

  func_words <- "function|role|required|essential|involved|mediat|regulat|critical|necessary|important|contribut|promotes|inhibits|activat|repress|maintain|controls|modulates"

  key_sentences <- sentences[grepl(gene_pattern, sentences, ignore.case = TRUE) &
    grepl(func_words, sentences, ignore.case = TRUE)]

  ko_sentences <- sentences[grepl(gene_pattern, sentences, ignore.case = TRUE) &
    grepl("knockout|knock.out|loss.of.function|null|deficien|deletion|ablat|conditional|cKO", sentences, ignore.case = TRUE)]

  all_key <- unique(c(head(key_sentences, 3), head(ko_sentences, 2)))

  if (length(all_key) > 0) {
    all_key <- sapply(all_key, function(s) {
      if (nchar(s) > 300) paste0(substr(s, 1, 297), "...") else s
    })
    compiled_df$functional_summary[i] <- paste(all_key, collapse = " | ")
  }
}

# =============================================================================
# SAVE
# =============================================================================
cat("\nSaving...\n")

write.csv(compiled_df %>% dplyr::select(-combined_text),
  "diff_gene_tiered_abstract_summaries.csv", row.names = FALSE)
cat("Saved: diff_gene_tiered_abstract_summaries.csv\n")

saveRDS(compiled_df, "diff_gene_tiered_compiled.rds")
cat("Saved: diff_gene_tiered_compiled.rds\n")

# Stats
cat("\n=== FINAL SUMMARY ===\n")
cat("Total genes:", nrow(compiled_df), "\n")
cat("With any papers:", sum(compiled_df$n_papers > 0), "\n")
cat("With Tier 1 (cerebellum/MB/diff):", sum(compiled_df$n_tier1 > 0), "\n")
cat("With Tier 2 (brain/neuron):", sum(compiled_df$n_tier2 > 0), "\n")
cat("With Tier 3 (brain cancer):", sum(compiled_df$n_tier3 > 0), "\n")
cat("With Tier 4 (cancer):", sum(compiled_df$n_tier4 > 0), "\n")
cat("With Tier 5 (development):", sum(compiled_df$n_tier5 > 0), "\n")
cat("With Tier 6 (general):", sum(compiled_df$n_tier6 > 0), "\n")
cat("With functional keywords:", sum(!is.na(compiled_df$functional_keywords) & compiled_df$functional_keywords != ""), "\n")
cat("With functional summary:", sum(!is.na(compiled_df$functional_summary)), "\n")

cat("\n=== KEYWORD FREQUENCY ===\n")
all_kw <- unlist(strsplit(compiled_df$functional_keywords[!is.na(compiled_df$functional_keywords)], "; "))
kw_freq <- sort(table(all_kw), decreasing = TRUE)
for (k in names(kw_freq)) {
  cat(sprintf("  %-30s %3d genes (%2.0f%%)\n", k, kw_freq[k], 100*kw_freq[k]/length(symbols)))
}

cat("\n=== DONE ===\n")
