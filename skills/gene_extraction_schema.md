# Gene Functional Annotation Extraction Schema

## Purpose
Extract structured, tabular annotations from PubMed abstracts and UniProt
descriptions for GNP differentiation genes. Each gene gets one row with
standardized fields.

## Fields to Extract

### 1. Molecular Function (what the protein DOES)
- **protein_class**: One of: ion channel, transporter, GPCR, receptor,
  kinase, phosphatase, transcription factor, enzyme, structural protein,
  adhesion molecule, scaffold/adaptor, secreted factor, GTPase/regulator,
  ubiquitin ligase, protease, motor protein, other
- **molecular_activity**: Brief verb phrase (e.g., "phosphorylates X",
  "conducts potassium", "binds DNA at E-box motifs")
- **substrates_or_ligands**: Known substrates, ligands, or binding partners
- **pathway_membership**: Named pathways (e.g., Wnt, Hedgehog, Notch,
  MAPK, PI3K/AKT, cAMP, calcium, mTOR)

### 2. Expression & Localization
- **subcellular**: plasma membrane, cytoplasm, nucleus, synapse,
  mitochondria, ER, Golgi, extracellular, axon, dendrite, growth cone
- **tissues_expressed**: Key tissues (brain, cerebellum, heart, muscle,
  liver, kidney, etc.)
- **cell_types**: Specific cell types (granule neuron, Purkinje cell,
  astrocyte, oligodendrocyte, interneuron, etc.)
- **temporal**: developmental timing if mentioned (embryonic, postnatal,
  adult, during differentiation)

### 3. Phenotype (what happens when it's lost)
- **ko_phenotype**: Brief description of knockout/LOF phenotype
- **ko_severity**: viable, subviable, embryonic lethal, perinatal lethal,
  postnatal lethal, conditional phenotype
- **neural_phenotype**: Specific neural phenotype if any (ataxia, seizure,
  learning deficit, axon guidance defect, etc.)
- **cerebellar_phenotype**: Cerebellar-specific phenotype if any

### 4. Disease Association
- **cancer_type**: Cancer types associated (medulloblastoma, glioma,
  breast, colorectal, etc.) or "none known"
- **cancer_role**: oncogene, tumor suppressor, context-dependent, biomarker,
  therapeutic target, none
- **neurological_disease**: Associated neurological diseases (epilepsy,
  autism, intellectual disability, ataxia, neurodegeneration, etc.)
- **other_disease**: Non-neural disease associations

### 5. Relevance to This Study
- **cerebellar_evidence**: Direct evidence in cerebellum (yes/no)
- **gnp_evidence**: Evidence in granule neuron progenitors specifically
- **differentiation_role**: Role in differentiation if known (promotes,
  inhibits, required for, marker of)
- **mb_evidence**: Evidence in medulloblastoma (yes/no)
