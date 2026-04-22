suppressMessages({
  library(tidyverse)
  library(rentrez)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# =============================================================================
# 1. LOAD GENES
# =============================================================================
cat("Loading gene data...\n")

master <- read.table("gnp_MB_histone_k_toc_t50.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, quote="\"")
colnames(master) <- gsub("\"","",colnames(master))

diff_genes <- master %>%
  dplyr::filter(max_cat == "P14_P56",
    group_1 == "Dynically expressed in granule neuron lineage",
    group_2 == ">2 log2FC", t50_cat %in% c("Mid","Late"),
    !grepl("_exclude", max_cat_exclude))

# Get human orthologs for searching (papers may use human gene names)
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

cat("Genes to query:", length(symbols), "\n")

# =============================================================================
# 2. DEFINE TIERED SEARCH STRATEGY
# =============================================================================

# Each tier: keywords to AND with the gene name
# We search with both mouse and human symbols
tiers <- list(
  list(name = "Tier1_cerebellum_MB_diff",
    terms = c("cerebellum", "cerebellar", "granule neuron", "granule cell",
              "medulloblastoma", "differentiation")),
  list(name = "Tier2_brain_neuron",
    terms = c("brain", "neuron", "neural", "CNS", "cortex", "hippocampus",
              "synapse", "synaptic", "axon", "dendrite")),
  list(name = "Tier3_brain_cancer",
    terms = c("brain tumor", "glioma", "glioblastoma", "neuroblastoma",
              "brain cancer", "neuro-oncology")),
  list(name = "Tier4_cancer",
    terms = c("cancer", "tumor", "oncogene", "tumor suppressor", "carcinoma",
              "malignant", "metastasis")),
  list(name = "Tier5_development",
    terms = c("development", "embryo", "progenitor", "stem cell", "lineage",
              "morphogenesis", "patterning")),
  list(name = "Tier6_general",
    terms = NULL)  # No filter — just gene name
)

MAX_PAPERS <- 20

# =============================================================================
# 3. TIERED PUBMED SEARCH
# =============================================================================
cat("\nStarting tiered PubMed search...\n")
cat("Strategy: up to", MAX_PAPERS, "papers per gene, prioritized by relevance tier\n")
cat("Estimated time: ~90-120 minutes\n\n")

all_gene_pmids <- list()  # gene -> list of PMIDs with tier info
all_pmids_needed <- character()  # all unique PMIDs to fetch

for (gi in seq_along(symbols)) {
  sym <- symbols[gi]
  hsym <- human_symbols[sym]

  if (gi %% 50 == 0 || gi == 1) {
    cat(sprintf("[%s] Gene %d / %d: %s (human: %s)\n",
      format(Sys.time(), "%H:%M:%S"), gi, length(symbols), sym,
      ifelse(is.na(hsym), "NA", hsym)))
  }

  collected_pmids <- character()
  tier_assignments <- character()

  for (tier_idx in seq_along(tiers)) {
    tier <- tiers[[tier_idx]]
    remaining <- MAX_PAPERS - length(collected_pmids)
    if (remaining <= 0) break

    # Build query: gene name (mouse OR human) AND tier keywords
    gene_part <- paste0("(", sym, "[Title/Abstract]")
    if (!is.na(hsym) && hsym != "" && hsym != sym) {
      gene_part <- paste0(gene_part, " OR ", hsym, "[Title/Abstract])")
    } else {
      gene_part <- paste0(gene_part, ")")
    }

    if (!is.null(tier$terms)) {
      tier_kw <- paste0("(", paste(tier$terms, collapse = " OR "), ")")
      query <- paste(gene_part, "AND", tier_kw)
    } else {
      query <- gene_part
    }

    result <- tryCatch({
      # Sort by relevance (PubMed "best match" which factors in citation impact)
      search <- entrez_search(db = "pubmed", term = query,
        retmax = remaining + 5, sort = "relevance",
        use_history = FALSE)
      search$ids
    }, error = function(e) character())

    # Remove already-collected PMIDs
    new_pmids <- setdiff(result, collected_pmids)
    new_pmids <- head(new_pmids, remaining)

    if (length(new_pmids) > 0) {
      collected_pmids <- c(collected_pmids, new_pmids)
      tier_assignments <- c(tier_assignments, rep(tier$name, length(new_pmids)))
    }

    Sys.sleep(0.35)  # Rate limit
  }

  all_gene_pmids[[sym]] <- data.frame(
    pmid = collected_pmids,
    tier = tier_assignments,
    rank = seq_along(collected_pmids),
    stringsAsFactors = FALSE
  )

  all_pmids_needed <- unique(c(all_pmids_needed, collected_pmids))
}

cat("\nSearch complete.\n")
cat("Total unique PMIDs to fetch:", length(all_pmids_needed), "\n")

# Per-gene stats
pmid_counts <- sapply(all_gene_pmids, nrow)
cat("Papers per gene: median =", median(pmid_counts),
    "mean =", round(mean(pmid_counts), 1),
    "max =", max(pmid_counts), "\n")
cat("Genes with 0 papers:", sum(pmid_counts == 0), "\n")
cat("Genes with 20 papers:", sum(pmid_counts == 20), "\n")

# Tier breakdown
all_tiers <- do.call(rbind, all_gene_pmids)
cat("\nPapers by tier:\n")
print(table(all_tiers$tier))

# Save PMID assignments
saveRDS(all_gene_pmids, "tiered_pmid_assignments.rds")
cat("Saved: tiered_pmid_assignments.rds\n")

# =============================================================================
# 4. BATCH FETCH ABSTRACTS
# =============================================================================
cat("\nFetching abstracts in batches of 200...\n")

batch_size <- 200
n_batches <- ceiling(length(all_pmids_needed) / batch_size)
all_abstract_text <- list()

for (b in 1:n_batches) {
  start_idx <- (b - 1) * batch_size + 1
  end_idx <- min(b * batch_size, length(all_pmids_needed))
  batch_ids <- all_pmids_needed[start_idx:end_idx]

  if (b %% 10 == 0 || b == 1) {
    cat(sprintf("  Batch %d / %d (%d PMIDs)...\n", b, n_batches, length(batch_ids)))
  }

  result <- tryCatch({
    # Fetch as XML for structured parsing
    fetched <- entrez_fetch(db = "pubmed", id = batch_ids, rettype = "xml")

    # Parse XML to extract per-PMID abstracts
    xml_doc <- xml2::read_xml(fetched)
    articles <- xml2::xml_find_all(xml_doc, "//PubmedArticle")

    for (art in articles) {
      pmid <- xml2::xml_text(xml2::xml_find_first(art, ".//PMID"))
      title <- xml2::xml_text(xml2::xml_find_first(art, ".//ArticleTitle"))
      abstract_nodes <- xml2::xml_find_all(art, ".//AbstractText")
      abstract <- paste(xml2::xml_text(abstract_nodes), collapse = " ")

      # Get year
      year <- tryCatch(
        xml2::xml_text(xml2::xml_find_first(art, ".//PubDate/Year")),
        error = function(e) NA)

      all_abstract_text[[pmid]] <- list(
        title = title,
        abstract = abstract,
        year = year
      )
    }
    TRUE
  }, error = function(e) {
    cat(sprintf("  ERROR in batch %d: %s\n", b, e$message))
    FALSE
  })

  Sys.sleep(0.35)
}

cat("Fetched abstracts for", length(all_abstract_text), "/", length(all_pmids_needed), "PMIDs\n")

# Save raw abstracts
saveRDS(all_abstract_text, "tiered_abstracts_raw.rds")
cat("Saved: tiered_abstracts_raw.rds\n")

# =============================================================================
# 5. COMPILE PER-GENE ABSTRACT TEXT
# =============================================================================
cat("\nCompiling per-gene abstract text...\n")

gene_compiled <- list()

for (sym in symbols) {
  pmid_df <- all_gene_pmids[[sym]]
  if (is.null(pmid_df) || nrow(pmid_df) == 0) {
    gene_compiled[[sym]] <- data.frame(
      symbol = sym, n_papers = 0,
      n_tier1 = 0, n_tier2 = 0, n_tier3 = 0, n_tier4 = 0, n_tier5 = 0, n_tier6 = 0,
      combined_text = NA_character_,
      titles = NA_character_,
      years = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  # Gather abstracts in tier order
  texts <- character()
  titles <- character()
  years <- character()
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

  tier_counts <- table(factor(pmid_df$tier, levels = sapply(tiers, "[[", "name")))

  gene_compiled[[sym]] <- data.frame(
    symbol = sym,
    n_papers = length(texts),
    n_tier1 = as.integer(tier_counts[1]),
    n_tier2 = as.integer(tier_counts[2]),
    n_tier3 = as.integer(tier_counts[3]),
    n_tier4 = as.integer(tier_counts[4]),
    n_tier5 = as.integer(tier_counts[5]),
    n_tier6 = as.integer(tier_counts[6]),
    combined_text = paste(texts, collapse = " ||| "),
    titles = paste(titles, collapse = " ||| "),
    years = paste(years, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

compiled_df <- do.call(rbind, gene_compiled)
cat("Genes with compiled text:", sum(compiled_df$n_papers > 0), "\n")
cat("Median papers per gene:", median(compiled_df$n_papers), "\n")

# =============================================================================
# 6. EXTRACT FUNCTIONAL KEYWORDS
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
    if (grepl(category_patterns[[cat_name]], txt, ignore.case = TRUE)) {
      matched <- c(matched, cat_name)
    }
  }
  compiled_df$functional_keywords[i] <- paste(matched, collapse = "; ")
}

cat("Genes with keywords:", sum(!is.na(compiled_df$functional_keywords) & compiled_df$functional_keywords != ""), "\n")

# Keyword frequency
cat("\n=== KEYWORD FREQUENCY ===\n")
all_kw <- unlist(strsplit(compiled_df$functional_keywords[!is.na(compiled_df$functional_keywords)], "; "))
kw_freq <- sort(table(all_kw), decreasing = TRUE)
for (k in names(kw_freq)) {
  cat(sprintf("  %-30s %3d genes (%2.0f%%)\n", k, kw_freq[k], 100*kw_freq[k]/length(symbols)))
}

# =============================================================================
# 7. EXTRACT KEY SENTENCES
# =============================================================================
cat("\nExtracting key functional sentences...\n")

compiled_df$functional_summary <- NA_character_

for (i in 1:nrow(compiled_df)) {
  sym <- compiled_df$symbol[i]
  txt <- compiled_df$combined_text[i]
  if (is.na(txt) || nchar(txt) < 100) next

  # Split into sentences
  sentences <- unlist(strsplit(txt, "(?<=[.!?])\\s+", perl = TRUE))

  # Find sentences about function/role of this gene
  # Use both mouse and human symbol
  hsym <- human_symbols[sym]
  gene_pattern <- paste0("\\b", sym, "\\b")
  if (!is.na(hsym) && hsym != "" && hsym != sym) {
    gene_pattern <- paste0("(\\b", sym, "\\b|\\b", hsym, "\\b)")
  }

  func_words <- "function|role|required|essential|involved|mediat|regulat|critical|necessary|important|contribut|promotes|inhibits|activat|repress|maintain|controls|modulates"

  key_sentences <- sentences[grepl(gene_pattern, sentences, ignore.case = TRUE) &
    grepl(func_words, sentences, ignore.case = TRUE)]

  # Also look for "knockout" or "loss of" sentences
  ko_sentences <- sentences[grepl(gene_pattern, sentences, ignore.case = TRUE) &
    grepl("knockout|knock.out|loss.of.function|null|deficien|deletion|ablat|conditional|cKO", sentences, ignore.case = TRUE)]

  all_key <- unique(c(head(key_sentences, 3), head(ko_sentences, 2)))

  if (length(all_key) > 0) {
    # Truncate each sentence to 300 chars
    all_key <- sapply(all_key, function(s) {
      if (nchar(s) > 300) paste0(substr(s, 1, 297), "...") else s
    })
    compiled_df$functional_summary[i] <- paste(all_key, collapse = " | ")
  }
}

cat("Genes with functional summary:", sum(!is.na(compiled_df$functional_summary)), "\n")

# =============================================================================
# 8. SAVE
# =============================================================================
cat("\nSaving results...\n")

# Save compiled data
write.csv(compiled_df %>% dplyr::select(-combined_text),
  "diff_gene_tiered_abstract_summaries.csv", row.names = FALSE)
cat("Saved: diff_gene_tiered_abstract_summaries.csv\n")

# Save with full text for reference
saveRDS(compiled_df, "diff_gene_tiered_compiled.rds")
cat("Saved: diff_gene_tiered_compiled.rds\n")

# Summary stats
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

cat("\n=== DONE ===\n")
