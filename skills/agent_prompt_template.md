# Sub-agent Prompt Template

This is the prompt sent to each of the 24 sub-agents for structured gene annotation extraction. Each agent receives one batch file (~30-40 genes), reads the UniProt descriptions and PubMed abstracts, and writes a TSV with 17 structured fields per gene.

## Template

```
You are a molecular biology expert extracting structured annotations from gene literature.

Read the file at /path/to/agent_batches/batch_NN.txt

For EACH gene in the file, read all available text (UniProt function, MGI phenotypes,
abstracts) and extract a TSV table with these columns:

gene, protein_class, molecular_action, key_substrates, subcellular, tissues,
cell_types, pathways, ko_phenotype, ko_severity, neural_phenotype,
cerebellar_evidence, cancer_types, cancer_role, neuro_disease,
differentiation_role, one_line_summary

**protein_class** must be one of: ion channel, transporter, GPCR, receptor, kinase,
phosphatase, transcription factor, enzyme, structural protein, adhesion molecule,
scaffold/adaptor, secreted factor, GTPase/regulator, ubiquitin ligase, protease,
motor protein, channel regulator, ECM protein, other

**Key directives:**

1. READ and INTERPRET the abstracts. Don't just keyword match.

2. Distinguish functional studies from bioinformatic screens:
   - "SLC16A1 inhibits pancreatic β-cell proliferation" (functional) vs
     "SLC16A1 appears in TCGA breast cancer screen" (bioinformatic)
   - For bioinformatic-only mentions, use `cancer_role = biomarker only`,
     not `oncogene` or `tumor suppressor`.

3. Use "NA" for unknown/empty fields. Do not speculate.

4. Focus on the ACTUAL phenotype, not generic terms:
   - Good: "ataxia, impaired cerebellar Purkinje cell arborization"
   - Bad: "abnormal brain morphology"

5. The `one_line_summary` should be written as if for a biology review paper —
   specific, functional, and not generic. Example:
   - Good: "PRKCC is a Purkinje-cell-enriched serine/threonine kinase whose
     gain-of-function mutations cause spinocerebellar ataxia type 14 (SCA14)."
   - Bad: "PRKCC is a kinase involved in signaling."

Write the output as a TSV file to: /path/to/agent_batches/batch_NN_results.tsv
```

## Field-specific guidance

### `molecular_action`
Verb phrase. What the protein does, mechanistically.
- "conducts potassium ions across the plasma membrane"
- "phosphorylates Ser/Thr residues on CREB"
- "activates Rac1 GTPase via nucleotide exchange"

### `subcellular`
Semicolon-separated list. Standard terms only:
plasma membrane, synapse, nucleus, cytoplasm, mitochondria, ER, Golgi, axon, dendrite, growth cone, extracellular

### `pathways`
Named pathways, not generic terms:
Wnt, Hedgehog, Notch, TGF-beta/BMP, MAPK/ERK, PI3K/AKT/mTOR, cAMP/PKA, Calcium, NF-kB, JAK/STAT, Hippo/YAP

### `ko_severity`
One of: embryonic lethal, perinatal lethal, postnatal lethal, subviable, viable, unknown

### `cancer_role`
One of: oncogene, tumor suppressor, biomarker only, none, unknown

### `differentiation_role`
One of: promotes, inhibits, required for, marker of, associated, none known

### `cerebellar_evidence`
One of: granule neuron, cerebellum (general), no, NA

## Input format (per gene in batch file)

```
========== GENE: Sema4f ==========
Description: sema domain, immunoglobulin domain, TM domain...
H3K27me3: 0 | Timing: Mid

--- UniProt Function ---
Probable cell surface receptor that regulates oligodendroglial precursor
cell migration...

--- MGI Knockout Phenotypes ---
decreased bone mineral density | decreased leukocyte cell number | ...

--- Abstracts (20 papers, 7 tier1) ---
Clinical outcomes in neuroblastoma (NB) are closely linked to...
---
The transition of differentiated Schwann cells to support of nerve repair...
---
...
```

## Output format (tab-separated, UTF-8, one row per gene)

```
gene<TAB>protein_class<TAB>molecular_action<TAB>...<TAB>one_line_summary
Sema4f<TAB>receptor<TAB>transmembrane semaphorin binding PlexinB1<TAB>...<TAB>Sema4F is a postsynaptic transmembrane semaphorin at glutamatergic synapses, anchored via PSD-95, that regulates oligodendroglial precursor migration.
```
