# POEM: Prediction of Operons from Metagenomes

POEM is a computational pipeline for predicting operons and identifying conserved core operons from metagenomic sequencing data. It supports paired-end sequencing reads, assembled contigs, and existing gene predictions.

## 1. Installation

Clone the POEM repository:

```bash
git clone <repository-url>
cd MetagenomicOperons
```

Install the dependencies using the provided installation script:

```bash
bash install.sh
```

The installation script creates a Conda environment named `con_env`, installs the required software packages, and downloads the COG2020 reference database.

**Requirements:**
- Conda and Mamba
- Python 3.12 (installed by the script)
- Internet access for installation and database download

Activate the environment:

```bash
conda activate con_env
```

## 2. Running POEM

POEM supports three input options depending on the type of data available and the preprocessing already performed. Users can start with raw paired-end sequencing reads, previously assembled contigs, or existing gene predictions. All three options require an output directory (-o) where the pipeline results will be stored.

### Option 1: Run complete pipeline from paired-end reads

This option runs the complete POEM pipeline starting from raw paired-end metagenomic sequencing reads. POEM performs metagenomic assembly, gene prediction, functional annotation, protein clustering, operon prediction, and core-operon identification.

```bash
./poem -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

**Required parameters:**

| Parameter | Description |
|---|---|
| `-1 <file>` | Forward paired-end FASTQ file |
| `-2 <file>` | Reverse paired-end FASTQ file |
| `-o <directory>` | Output directory for POEM results |

**Optional parameters:**

| Parameter | Description | Default |
|---|---|---|
| `-a <assembler>` | Metagenomic assembly tool | `spades` |
| `--model <name>` | Operon prediction classifier | `all_species` |

POEM uses metaSPAdes as the default assembler. The following assemblers are supported:

- `spades`: metaSPAdes
- `megahit`: MEGAHIT
- `idba`: IDBA-UD

#### Operon prediction classifier selection (`--model`)

POEM provides Random Forest classifiers trained on operon data from nine bacterial species, including individual species-specific classifiers and a pooled classifier trained on the combined dataset.

The following classifiers are supported:

- `all_species`: Pooled classifier trained on all nine bacterial species (default)
- `monocytogenes`: Classifier trained on *Listeria monocytogenes*

The pooled classifier is intended for general metagenomic applications, where the species composition may be unknown.

#### Example with optional parameters

For example, to run POEM using MEGAHIT and the *Listeria monocytogenes* operon prediction classifier:

```bash
./poem \
    -1 reads_R1.fastq.gz \
    -2 reads_R2.fastq.gz \
    -a megahit \
    --model monocytogenes \
    -o results/
```

### Option 2: Run pipeline from assembled contigs

This option is intended for users who have already performed metagenomic assembly using an external assembly tool. POEM skips the assembly stage and begins with gene prediction on the supplied contigs, followed by all subsequent pipeline stages.

```bash
./poem -af contigs.fa -o results/
```

**Required parameters:**

| Parameter | Description |
|---|---|
| `-af <file>` | Assembled contigs in FASTA format |
| `-o <directory>` | Output directory for POEM results |

**Optional parameters:**

| Parameter | Description | Default |
|---|---|---|
| `--model <name>` | Operon prediction classifier | `all_species` |

### Option 3: Run pipeline from existing gene predictions

This option is intended for users who have already performed metagenomic assembly and gene prediction using external tools. POEM skips both assembly and gene prediction, allowing users to retain their existing gene annotations and protein translations.

```bash
./poem \
    -af contigs.fa \
    -gff genes.gff \
    -pf proteins.faa \
    -o results/
```

**Required parameters:**

| Parameter | Description |
|---|---|
| `-af <file>` | Assembled contigs in FASTA format |
| `-gff <file>` | Existing CDS annotations in GFF3 format |
| `-pf <file>` | Protein translations in FASTA format |
| `-o <directory>` | Output directory for POEM results |

**Optional parameters:**

| Parameter | Description | Default |
|---|---|---|
| `--model <name>` | Operon prediction classifier | `all_species` |

The GFF and protein FASTA files must correspond to the same assembled contigs. Both annotation and protein translation files must be provided. POEM standardizes identifiers for downstream processing.

## 3. Pipeline Steps

Depending on the input option, POEM runs the following stages:

1. **Assembly:** Assemble paired-end reads using SPAdes, MEGAHIT, or IDBA-UD.
2. **Gene prediction:** Predict protein-coding genes using Prodigal.
3. **Function prediction:** Assign COG annotations using DIAMOND BLASTP.
4. **Protein clustering:** Cluster predicted proteins using CD-HIT.
5. **Operon prediction:** Predict operons and process their COG associations.
6. **Core-operon identification:** Identify conserved core operons using network-based analysis.
7. **Visualization:** Explore predicted core operons using an interactive Dash interface.

## 4. Output Files

POEM organizes results into the following directories:

```text
results/
├── assembly_output/
│   ├── poem_contigs.fa
│   └── contig_id_mapping.tsv
├── gene_prediction_output/
│   ├── gene.coords.gff
│   ├── gene_coordinates.csv
│   ├── protein.translations.faa
│   └── dna.translations.faa
├── function_prediction_output/
│   └── diamond_cog_reference.tsv
├── cdhit_output/
│   ├── clustered_proteins.fasta
│   ├── clustered_proteins.fasta.clstr
│   ├── cluster_to_genes
│   └── gene_to_clusters
├── operon_prediction_output/
│   ├── gene_pair_names_to_operon.csv
│   ├── gene_to_cogId.pkl
│   ├── cogId_to_gene.pkl
│   └── cog_pairs_df
└── core_operon_output/
    ├── leiden_communities.pkl
    ├── core_operons.csv
    ├── community_edges.csv
    └── core_operon_edges.csv
```

### Assembly output (`assembly_output/`)

| File | Description |
|---|---|
| `poem_contigs.fa` | Assembled contig sequences in FASTA format with standardized identifiers. |
| `contig_id_mapping.tsv` | Mapping between original and standardized contig identifiers. |

### Gene prediction output (`gene_prediction_output/`)

| File | Description |
|---|---|
| `gene.coords.gff` | Predicted or imported gene annotations in GFF format. |
| `gene_coordinates.csv` | Gene identifiers, genomic coordinates, contig identifiers, and strand information used for operon prediction. |
| `protein.translations.faa` | Protein sequences of predicted or imported genes in FASTA format. |
| `dna.translations.faa` | Nucleotide sequences of genes in FASTA format. |

### Function prediction output (`function_prediction_output/`)

| File | Description |
|---|---|
| `diamond_cog_reference.tsv` | DIAMOND BLASTP alignment results against the COG reference protein database. |

### Protein clustering output (`cdhit_output/`)

| File | Description |
|---|---|
| `clustered_proteins.fasta` | Representative protein sequences from CD-HIT clusters. |
| `clustered_proteins.fasta.clstr` | CD-HIT cluster assignments and sequence membership information. |
| `cluster_to_genes` | Mapping of protein clusters to their member genes, stored in Python pickle format. |
| `gene_to_clusters` | Mapping of genes to their assigned protein clusters, stored in Python pickle format. |

### Operon prediction and processing output (`operon_prediction_output/`)

| File | Description |
|---|---|
| `gene_pair_names_to_operon.csv` | Predicted operon membership of adjacent, same-strand gene pairs. |
| `gene_to_cogId.pkl` | Mapping of gene identifiers to their assigned COG identifiers. |
| `cogId_to_gene.pkl` | Mapping of COG identifiers to their associated genes. |
| `cog_pairs_df` | Co-operonic COG pairs obtained after merging predicted gene pairs into operons and filtering redundant gene pairs using CD-HIT cluster assignments. This file is stored in Python pickle format and serves as input for core-operon identification. |

### Core-operon identification output (`core_operon_output/`)

This directory contains the primary results of POEM, including the
identified Leiden communities, predicted core operons, and network
edge files.

| File | Description |
|---|---|
| `leiden_communities.pkl` | Identified Leiden communities containing COG identifiers and weighted edges, stored in Python pickle format. |
| `core_operons.csv` | Final predicted core operons, including their constituent COGs, parent community identifiers, sizes, edge weights, densities, and core-operon scores. Predictions are ranked by their core-operon scores. |
| `community_edges.csv` | Weighted COG associations within Leiden communities, including community identifiers. |
| `core_operon_edges.csv` | Weighted COG associations within predicted core operons, including core-operon and parent community identifiers. |


## 5. Visualization

POEM includes an interactive interface for exploring predicted core operons and their COG associations.

To launch visualization independently:

```bash
./poem -s visualize-core-operon -o results/
```

Open the local URL displayed in the terminal to view the interactive network.

## 6. Example Dataset

Run POEM using the included example dataset:

```bash
./poem --test
```

This runs the pipeline using paths defined in `example_config`.

## 7. Help

Display the available arguments:

```bash
./poem --help
```

POEM also supports individual pipeline stages for testing and debugging:

```bash
./poem -s gene-prediction -o results/
```

Stage-wise execution requires the corresponding upstream input files to be available.