suppressMessages(library(tidyverse))
setwd("~/Dropbox/GNP_MB_integrate/Ezh2_2022/")

# Load UniProt
uniprot <- read.table("/tmp/uniprot_mouse_functions.tsv", sep="\t", header=TRUE,
  stringsAsFactors=FALSE, quote="", fill=TRUE, comment.char="")
cat("UniProt entries:", nrow(uniprot), "\n")

# Parse primary gene name
uniprot$primary_gene <- sapply(strsplit(uniprot$Gene.Names, " "), function(x) x[1])
uniprot$primary_gene_upper <- toupper(uniprot$primary_gene)

# Load our diff genes
ann <- read.csv("diff_gene_full_annotations.csv", stringsAsFactors=FALSE)
ann$symbol_upper <- toupper(ann$mgi_symbol)

matched <- ann %>%
  dplyr::left_join(
    uniprot %>% dplyr::select(primary_gene, primary_gene_upper, Protein.names, Function..CC.) %>%
      dplyr::filter(!duplicated(primary_gene_upper)),
    by = c("symbol_upper" = "primary_gene_upper"))

cat("Genes with UniProt function:", sum(!is.na(matched$Function..CC.) & matched$Function..CC. != ""), "/", nrow(matched), "\n")

# Clean function text
matched$uniprot_function <- gsub("\\{ECO:[^\\}]+\\}", "", matched$Function..CC.)
matched$uniprot_function <- gsub("^FUNCTION: ", "", matched$uniprot_function)
matched$uniprot_function <- gsub("\\s+", " ", trimws(matched$uniprot_function))
matched$uniprot_protein <- matched$Protein.names

# Save
write.csv(matched %>% dplyr::select(mgi_symbol, GENENAME, uniprot_protein, uniprot_function),
  "diff_gene_uniprot_functions.csv", row.names=FALSE)
cat("Saved: diff_gene_uniprot_functions.csv\n")

# Show coverage
cat("\nGenes with UniProt function text:", sum(!is.na(matched$uniprot_function) & nchar(matched$uniprot_function) > 10), "\n")
cat("Genes with protein name:", sum(!is.na(matched$uniprot_protein) & matched$uniprot_protein != ""), "\n")
