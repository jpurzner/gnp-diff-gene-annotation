# GNP Differentiation Gene Functional Annotation Pipeline

Agent-based functional annotation of cerebellar granule neuron precursor (GNP) differentiation genes for the EZH2/Polycomb medulloblastoma study.

## Overview

This pipeline combines structured database annotations (GO, KEGG, Reactome, MSigDB, MGI phenotypes, UniProt) with **LLM-based extraction from PubMed abstracts** to produce a tabular, semantically-rich annotation per gene that can be used as features for dimensionality reduction (UMAP) and clustering.

The motivation: traditional GO/pathway enrichment captures surface-level keywords but misses the actual biology — *what the gene does* and *what happens when you perturb it*. Regex-based keyword matching also confuses bioinformatic mentions (e.g., "appears in TCGA screen") with functional studies. Sub-agents that *read* the abstracts can distinguish these and write expert-quality structured annotations.

## Pipeline

```
1. Pull UniProt function descriptions (Swiss-Prot CC FUNCTION field)
            |
            v
2. Tiered PubMed search per gene (cerebellum/MB/diff > brain > cancer > dev > general)
   Up to 20 abstracts per gene, ranked by tier + relevance
            |
            v
3. Group genes by UMAP cluster, prepare per-batch text bundles
   (40 genes per batch × 24 batches = 786 genes)
            |
            v
4. Run 24 sub-agents in parallel — each reads abstracts + UniProt for its batch
   and produces structured TSV with 17 fields per gene
            |
            v
5. Merge results, integrate with database matrix, re-run UMAP
```

## Repository structure

```
.
├── README.md
├── skills/
│   └── gene_extraction_schema.md      # Schema defining 17 extraction fields
├── scripts/
│   ├── 01_parse_uniprot.R             # Parse UniProt function descriptions
│   ├── 02_fetch_tiered_abstracts.R    # Tiered PubMed search via rentrez
│   ├── 02b_process_abstracts.R        # Compile per-gene text bundles
│   ├── 03_prepare_agent_batches.R     # Cluster-grouped batch files
│   ├── 04_merge_agent_results.R       # Merge agent TSVs + Excel build
│   ├── 05_umap_with_agent_features.R  # UMAP with agent + database features
│   ├── umap_semantic_v3.R             # Database-only semantic UMAP
│   ├── umap_k10_dotplots.R            # Per-cluster enrichment dotplots
│   ├── enrichment_umap_overlays.R     # Term overlays on UMAP
│   └── pathway_analysis.R             # GO/KEGG/Reactome enrichment
├── agent_batches/                     # Sample input + output
│   ├── batch_02.txt                   # Sample text bundle for one batch
│   ├── batch_02_results.tsv           # Sample TSV output from agent
│   ├── batch_12.txt
│   └── batch_12_results.tsv
└── data_samples/
    ├── agent_extracted_annotations.csv  # Final merged 786-gene table
    └── differentiation_genes_with_orthologs.csv  # Gene list with human orthologs
```

## The schema (skills/gene_extraction_schema.md)

Each gene gets one row with these 17 columns:

| Field | Description |
|---|---|
| `protein_class` | One of 19 categories (ion channel, GPCR, kinase, TF, ...) |
| `molecular_action` | Brief verb phrase (what it DOES) |
| `key_substrates` | Substrates, ligands, binding partners |
| `subcellular` | Where in the cell |
| `tissues` | Where expressed/studied |
| `cell_types` | Specific cell types |
| `pathways` | Named signaling pathways |
| `ko_phenotype` | Knockout phenotype |
| `ko_severity` | Lethality/viability category |
| `neural_phenotype` | Specific neural phenotype |
| `cerebellar_evidence` | Direct cerebellar data? |
| `cancer_types` | Cancers where functionally studied |
| `cancer_role` | Oncogene / suppressor / biomarker / none |
| `neuro_disease` | Associated neurological diseases |
| `differentiation_role` | Promotes / inhibits / required for / marker |
| `one_line_summary` | Expert single-sentence summary |

## Agent prompt

See `skills/agent_prompt_template.md` for the prompt sent to each sub-agent. Key directives:

- READ and INTERPRET the abstracts — don't keyword match
- Distinguish bioinformatic mentions ("appears in TCGA screen") from functional studies
- Use "NA" for unknown — don't speculate
- Write the one-line summary as if for a biology review paper

## Coverage achieved (786 differentiation genes)

| Field | Coverage |
|---|---|
| Protein class | 100% |
| Molecular action | 95% |
| One-line summary | 96% |
| KO phenotype | 88% |
| Cancer role | 75% |
| Neural phenotype | 62% |
| Cerebellar evidence | 50% |
| Differentiation role | 73% |

## How to use

```bash
# Required R packages
Rscript -e 'install.packages(c("tidyverse","rentrez","openxlsx","uwot","irlba","cluster","ggrepel","patchwork"))'
Rscript -e 'BiocManager::install(c("clusterProfiler","org.Mm.eg.db","msigdbr","enrichplot"))'

# Pipeline (run in order)
Rscript scripts/01_parse_uniprot.R
Rscript scripts/02_fetch_tiered_abstracts.R    # ~90 minutes
Rscript scripts/02b_process_abstracts.R
Rscript scripts/03_prepare_agent_batches.R     # writes 24 .txt batches

# Run sub-agents (one per batch, in parallel — see prompt template)
# Each agent reads agent_batches/batch_NN.txt and writes batch_NN_results.tsv

Rscript scripts/04_merge_agent_results.R       # merge TSVs + build Excel
Rscript scripts/05_umap_with_agent_features.R  # UMAP with agent features
```

## Citation

Part of the EZH2/Polycomb-mediated regulation of GNP differentiation in medulloblastoma study (Purzner et al.).

## License

MIT
