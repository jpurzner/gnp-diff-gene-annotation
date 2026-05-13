suppressMessages({library(tidyverse); library(rentrez); library(org.Mm.eg.db); library(AnnotationDbi); library(xml2)})
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

cat("=== PROLIFERATION GENE PUBMED FETCH ===\n")
cat("Loading data...\n")

prolif <- read.csv("gene_list_proliferation.csv", stringsAsFactors=FALSE)
symbols <- prolif$mgi_symbol
human_symbols <- setNames(prolif$hgnc_symbol, prolif$mgi_symbol)
cat("Genes to query:", length(symbols), "\n")

# =============================================================================
# TIER STRATEGY — PROLIFERATION-FOCUSED
# =============================================================================
tiers <- list(
  list(name = "Tier1_cellcycle_GNP_MB",
    terms = c("cell cycle", "mitosis", "mitotic", "proliferation", "G1 phase",
              "S phase", "G2 phase", "M phase", "cyclin", "CDK", "CDKN",
              "cell division", "DNA replication", "checkpoint",
              "granule neuron precursor", "GNP", "medulloblastoma",
              "cerebellar progenitor")),
  list(name = "Tier2_brain_progenitor",
    terms = c("brain", "neural progenitor", "neuroepithelium", "neural stem cell",
              "radial glia", "neurogenesis", "cerebellum", "cortex", "rhombic lip",
              "ventricular zone", "subventricular")),
  list(name = "Tier3_brain_cancer",
    terms = c("brain tumor", "glioma", "glioblastoma", "neuroblastoma",
              "brain cancer", "medulloblastoma")),
  list(name = "Tier4_cancer",
    terms = c("cancer", "tumor", "oncogene", "tumor suppressor", "carcinoma",
              "malignant", "metastasis")),
  list(name = "Tier5_development",
    terms = c("development", "embryo", "progenitor", "stem cell", "lineage",
              "morphogenesis", "patterning", "differentiation")),
  list(name = "Tier6_general", terms = NULL))

MAX_PAPERS <- 20

cat("\nStarting tiered PubMed search... (~50 minutes for 438 genes)\n")
all_gene_pmids <- list()
all_pmids_needed <- character()

for (gi in seq_along(symbols)) {
  sym <- symbols[gi]; hsym <- human_symbols[sym]
  if (gi %% 50 == 0 || gi == 1)
    cat(sprintf("[%s] Gene %d / %d: %s\n", format(Sys.time(), "%H:%M:%S"), gi, length(symbols), sym))

  collected_pmids <- character(); tier_assignments <- character()
  for (tier_idx in seq_along(tiers)) {
    tier <- tiers[[tier_idx]]
    remaining <- MAX_PAPERS - length(collected_pmids)
    if (remaining <= 0) break

    gene_part <- paste0("(", sym, "[Title/Abstract]")
    if (!is.na(hsym) && hsym != "" && hsym != sym)
      gene_part <- paste0(gene_part, " OR ", hsym, "[Title/Abstract])")
    else gene_part <- paste0(gene_part, ")")

    if (!is.null(tier$terms)) {
      query <- paste(gene_part, "AND", paste0("(", paste(tier$terms, collapse=" OR "), ")"))
    } else {
      query <- gene_part
    }

    result <- tryCatch({
      search <- entrez_search(db="pubmed", term=query, retmax=remaining + 5,
        sort="relevance", use_history=FALSE)
      search$ids
    }, error = function(e) character())

    new_pmids <- setdiff(result, collected_pmids)
    new_pmids <- head(new_pmids, remaining)
    if (length(new_pmids) > 0) {
      collected_pmids <- c(collected_pmids, new_pmids)
      tier_assignments <- c(tier_assignments, rep(tier$name, length(new_pmids)))
    }
    Sys.sleep(0.35)
  }

  all_gene_pmids[[sym]] <- data.frame(pmid=collected_pmids, tier=tier_assignments,
    rank=seq_along(collected_pmids), stringsAsFactors=FALSE)
  all_pmids_needed <- unique(c(all_pmids_needed, collected_pmids))
}

cat("\nSearch complete. Total unique PMIDs:", length(all_pmids_needed), "\n")
saveRDS(all_gene_pmids, "tiered_pmid_assignments_PROLIF.rds")

# Batch fetch abstracts
cat("\nFetching abstracts in batches of 200...\n")
batch_size <- 200
n_batches <- ceiling(length(all_pmids_needed) / batch_size)
all_abstract_text <- list()

for (b in 1:n_batches) {
  start_idx <- (b-1)*batch_size + 1
  end_idx <- min(b*batch_size, length(all_pmids_needed))
  batch_ids <- all_pmids_needed[start_idx:end_idx]

  if (b %% 10 == 0 || b == 1) cat(sprintf("  Batch %d / %d\n", b, n_batches))

  tryCatch({
    fetched <- entrez_fetch(db="pubmed", id=batch_ids, rettype="xml")
    xml_doc <- xml2::read_xml(fetched)
    articles <- xml2::xml_find_all(xml_doc, "//PubmedArticle")
    for (art in articles) {
      pmid <- xml2::xml_text(xml2::xml_find_first(art, ".//PMID"))
      title <- xml2::xml_text(xml2::xml_find_first(art, ".//ArticleTitle"))
      abstract_nodes <- xml2::xml_find_all(art, ".//AbstractText")
      abstract <- paste(xml2::xml_text(abstract_nodes), collapse=" ")
      year <- tryCatch(xml2::xml_text(xml2::xml_find_first(art, ".//PubDate/Year")),
        error=function(e) NA)
      all_abstract_text[[pmid]] <- list(title=title, abstract=abstract, year=year)
    }
  }, error = function(e) cat("  ERROR batch", b, ":", e$message, "\n"))
  Sys.sleep(0.35)
}

cat("Fetched abstracts for", length(all_abstract_text), "PMIDs\n")
saveRDS(all_abstract_text, "tiered_abstracts_raw_PROLIF.rds")

# Compile per-gene
cat("\nCompiling per-gene text...\n")
tiers_names <- sapply(tiers, "[[", "name")
gene_compiled <- list()
for (sym in symbols) {
  pmid_df <- all_gene_pmids[[sym]]
  if (is.null(pmid_df) || nrow(pmid_df) == 0) {
    gene_compiled[[sym]] <- data.frame(symbol=sym, n_papers=0,
      n_tier1=0, n_tier2=0, n_tier3=0, n_tier4=0, n_tier5=0, n_tier6=0,
      combined_text=NA_character_, titles=NA_character_, years=NA_character_,
      stringsAsFactors=FALSE)
    next
  }
  texts <- character(); titles <- character(); years <- character()
  for (j in 1:nrow(pmid_df)) {
    pmid <- pmid_df$pmid[j]
    if (pmid %in% names(all_abstract_text)) {
      art <- all_abstract_text[[pmid]]
      if (!is.null(art$abstract) && nchar(art$abstract) > 10) {
        texts <- c(texts, art$abstract); titles <- c(titles, art$title)
        years <- c(years, ifelse(is.na(art$year), "NA", art$year))
      }
    }
  }
  tier_counts <- table(factor(pmid_df$tier, levels=tiers_names))
  gene_compiled[[sym]] <- data.frame(symbol=sym, n_papers=length(texts),
    n_tier1=as.integer(tier_counts[1]), n_tier2=as.integer(tier_counts[2]),
    n_tier3=as.integer(tier_counts[3]), n_tier4=as.integer(tier_counts[4]),
    n_tier5=as.integer(tier_counts[5]), n_tier6=as.integer(tier_counts[6]),
    combined_text=paste(texts, collapse=" ||| "),
    titles=paste(titles, collapse=" ||| "),
    years=paste(years, collapse=", "),
    stringsAsFactors=FALSE)
}
compiled_df <- do.call(rbind, gene_compiled)
saveRDS(compiled_df, "tiered_compiled_PROLIF.rds")
cat("Done. ", sum(compiled_df$n_papers > 0), "/", nrow(compiled_df), "genes have abstracts\n")
