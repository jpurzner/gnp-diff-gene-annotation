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
2. Tiered PubMed search per gene
   (cerebellum/MB/diff > brain > cancer > dev > general)
   Up to 20 abstracts per gene, ranked by tier + relevance
            |
            v
3. Group genes by UMAP cluster, prepare per-batch text bundles
   (~40 genes per batch × 24 batches = 786 genes)
            |
            v
4. Run 24 sub-agents in parallel — each reads abstracts + UniProt
   for its batch and produces structured TSV with 17 fields per gene
            |
            v
4b. Identify under-annotated genes that had abstracts (>= 5 papers
    but < 6 filled fields) and reprocess with a rescue agent
            |
            v
5. Merge results, integrate with database matrix, re-run UMAP
```

## Repository structure

```
.
├── README.md
├── LICENSE
├── skills/
│   ├── gene_extraction_schema.md    # Schema defining 17 extraction fields
│   └── agent_prompt_template.md     # Full prompt sent to each sub-agent
├── scripts/
│   ├── 01_parse_uniprot.R               # Parse UniProt function descriptions
│   ├── 02_fetch_tiered_abstracts.R      # Tiered PubMed search via rentrez
│   ├── 02b_process_abstracts.R          # Compile per-gene text bundles
│   ├── 03_prepare_agent_batches.R       # Cluster-grouped batch files
│   ├── 04_merge_agent_results.R         # Merge agent TSVs (first pass)
│   ├── 04b_fix_and_reprocess.R          # Lenient re-parse + find under-annotated
│   ├── 04c_merge_reprocess.R            # Merge rescue agent results
│   ├── 05_umap_with_agent_features.R    # UMAP with agent + database features
│   ├── umap_semantic_v3.R               # Database-only semantic UMAP
│   ├── umap_k10_dotplots.R              # Per-cluster enrichment dotplots
│   ├── enrichment_umap_overlays.R       # Term overlays on UMAP
│   └── pathway_analysis.R               # GO/KEGG/Reactome enrichment
├── agent_batches/                       # Sample input + output
│   ├── batch_02.txt                     # Sample text bundle (first pass)
│   ├── batch_02_results.tsv
│   ├── batch_12.txt
│   ├── batch_12_results.tsv
│   ├── reprocess_batch_01.txt           # Sample rescue batch
│   └── reprocess_batch_01_results.tsv
├── figures/                             # Generated figures
│   ├── Fig_gene_umap_agent_k10.pdf      # Main UMAP with term labels
│   ├── Fig_gene_umap_agent_k27.pdf      # H3K27me3 overlay
│   └── Fig_gene_umap_agent_proteinclass.pdf  # Protein class overlay
└── data_samples/
    ├── agent_extracted_annotations.csv  # Final merged 786-gene table
    └── differentiation_genes_with_orthologs.csv
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

## Rescue pass for under-annotated genes

The first agent pass correctly identified most RIKEN/Gm-prefixed predicted genes as uncharacterized. But ~9 genes had ≥5 abstracts yet ended up with <6 filled fields — these got a second pass with stronger guidance ("look for aliases, ortholog matches, protein family inference"). This recovered:

- **ST5** → scaffold/adaptor tumor suppressor causing severe neurodevelopmental syndrome
- **A230050P20Rik** → mouse ortholog of SHFL, an interferon-stimulated antiviral factor
- **3110082D06Rik** → mouse ortholog of PTCHD4, a p53-regulated Hedgehog pathway repressor
- **Mansc1** → candidate tumor suppressor at 12p in myeloid malignancies
- **D17H6S56E-3** → likely ortholog of VWA7, putative phospholipase C
- **Chi3l7** → chitinase-3-like family member (inferred function)
- **1700019D03Rik** → hematopoietic regulator (from MGI KO phenotype)
- **D630039A03Rik** → mouse ortholog of human C9orf152

The remaining ~25 predicted genes (Gm9899, Gm5454, Gm7276, etc.) are genuinely uncharacterized with no literature and no ortholog hits — they stay labeled "other" and are flagged as such.

## Coverage achieved (786 differentiation genes)

| Field | Coverage |
|---|---|
| One-line summary | 97% |
| Tissues | 96% |
| Molecular action | 95% |
| Subcellular | 94% |
| Cell types | 94% |
| Pathways | 94% |
| KO phenotype | 86% |
| Protein class | 82% |
| Cancer types | 75% |
| Cancer role | 71% |
| Differentiation role | 69% |
| Neural phenotype | 61% |
| Neuro disease | 56% |
| Cerebellar evidence | 45% |

## Final clustering (k=10 UMAP with agent + database features)

Silhouette score: **0.493**

| Cluster | N | %K27+ | Label | Dominant classes |
|---|---|---|---|---|
| 1 | 92 | 46% | Predicted/uncharacterized + diverse | other (35), scaffold (10) |
| 2 | 114 | **69%** | Ion channels & transporters | ion channel (46), transporter (36) |
| 3 | 60 | **72%** | Signal release / secretion | other (16), GPCR (11) |
| 4 | 72 | 50% | Kinase / phosphorylation signaling | kinase (25), secreted (11) |
| 5 | 54 | 52% | Lipid / small molecule metabolism | **enzyme (49/54 = 91%)** |
| 6 | 64 | **67%** | Cell projection / neuronal morphology | GTPase/reg (18), scaffold (10) |
| 7 | 71 | 58% | Adhesion / junction proteins | other (15), adhesion (11) |
| 8 | 56 | 59% | Morphogenesis / structural | structural (10), TF (7) |
| 9 | 62 | 61% | Transcription regulation | **TF (50/62 = 81%)** |
| 10 | 63 | **67%** | Neuron-associated / metallopeptidase | other (12), enzyme (10) |

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

# Run sub-agents (one per batch, in parallel — see skills/agent_prompt_template.md)
# Each agent reads agent_batches/batch_NN.txt and writes batch_NN_results.tsv

Rscript scripts/04_merge_agent_results.R       # first-pass merge
Rscript scripts/04b_fix_and_reprocess.R        # find under-annotated genes
# Run rescue agent on reprocess_batch_01.txt

Rscript scripts/04c_merge_reprocess.R          # merge rescue results + build Excel
Rscript scripts/05_umap_with_agent_features.R  # UMAP with agent features
```

## Citation

Part of the EZH2/Polycomb-mediated regulation of GNP differentiation in medulloblastoma study (Purzner et al.).

## License

MIT
