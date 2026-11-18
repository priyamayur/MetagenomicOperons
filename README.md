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

POEM supports three input options depending on the preprocessing already performed.

### Option 1: Paired-end sequencing reads

Run the complete pipeline starting from raw paired-end FASTQ files:

```bash
./poem -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

SPAdes (metaSPAdes) is the default assembler. Alternative assemblers can be selected using `-a`:

```bash
./poem -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -a megahit -o results/
```

Supported assemblers: `spades`, `megahit`, `idba`.

### Option 2: Assembled contigs

Skip assembly and begin with gene prediction:

```bash
./poem -af contigs.fa -o results/
```

POEM normalizes the contig FASTA headers before processing.

### Option 3: Existing gene predictions

Skip assembly and gene prediction by providing assembled contigs, GFF gene annotations, and protein translations:

```bash
./poem \
    -af contigs.fa \
    -gff genes.gff \
    -pf proteins.faa \
    -o results/
```

The GFF and protein FASTA files must correspond to the same assembled contigs. POEM standardizes identifiers for downstream processing.

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
│   └── clustered_proteins.fasta
├── operon_prediction_output/
└── core_operon_output/
    ├── leiden_communities.pkl
    └── path_quasi_hybrid_cores.csv
```

The main core-operon results are stored in `core_operon_output/`.

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