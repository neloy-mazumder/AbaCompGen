# AbaCompGen: *Acinetobacter baumannii* Comparative Genomics Pipeline

A fully reproducible, end-to-end bioinformatics pipeline for comparative genomic analysis of *Acinetobacter baumannii* clinical isolates. This pipeline replicates and extends the analysis published in:

> Kumkar SN, Kamble EE, Chavan NS, Dhotre DP and Pardesi KR (2022)  
> **Diversity of resistant determinants, virulence factors, and mobile genetic elements in *Acinetobacter baumannii* from India: A comprehensive in silico genome analysis.**  
> *Front. Cell. Infect. Microbiol.* 12:997897.  
> doi: [10.3389/fcimb.2022.997897](https://doi.org/10.3389/fcimb.2022.997897)

## Novel Additions Beyond the Paper

| Step | Analysis | Why Added |
|------|----------|-----------|
| 10 | Pan-genome (Roary) | Paper mentions open pan-genome but never quantifies it; provides core/accessory/cloud gene counts and Heaps' law alpha |
| 11 | GTDB-Tk phylogenomics | More robust GTDB-based taxonomic placement than SNP tree alone |
| 11 | Temporal ARG trend analysis | Tracks ARG acquisition and blaOXA-23 prevalence across collection years (2005–2020) |
| 12 | `advanced_analysis.R` — 23 figures | Replaces a basic 7-figure script with a full publication suite: alluvial epidemiology flows, IS–ARG distance, forest plot regression, annotated phylogeny, ridge plots, treemaps, and a complete statistical report |

---

## Table of Contents

- [Repository Name and Structure](#repository-structure)
- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Database Setup](#database-setup)
- [Configuration](#configuration)
- [Usage](#usage)
- [Step-by-Step Description](#step-by-step-description)
- [Outputs](#outputs)
- [Citation](#citation)

---

## Repository Structure

```
acinetobacter-baumannii-genomics/
├── config.sh                       ← USER EDITS THIS (paths + threads)
├── run_pipeline.sh                 ← Master runner (all or individual steps)
├── step00_setup.sh                 ← Directory setup + sample list
├── step01_annotation.sh            ← Genome annotation (Prokka)
├── step02_mlst.sh                  ← MLST — Pasteur + Oxford schemes
├── step03_ARG_prediction.sh        ← ARG prediction (CARD + AMRFinderPlus + ResFinder)
├── step04_VFG_prediction.sh        ← Virulence factor prediction (VFDB)
├── step05_plasmid_genomad.sh       ← Plasmid detection (geNomad + MOB-suite)
├── step06_IS_isescan.sh            ← IS element detection + IS-ARG linkage
├── step07_prophage.sh              ← Prophage detection (geNomad viral mode)
├── step08_SNP_phylogeny.sh         ← SNP phylogeny (snippy + IQ-TREE 2)
├── step09_AbaRI_comM.sh            ← AbaR resistance island detection
├── step10_pangenome_roary.sh       ← Pan-genome analysis [NOVEL]
├── step11_gtdbtk_temporal.sh       ← GTDB-Tk + temporal ARG trends [NOVEL]
├── step12_summary_visualization.sh ← Master table builder + calls advanced_analysis.R
└── advanced_analysis.R             ← 23 publication figures + statistics (called by step12)
```

**Input:** One genome assembly per isolate as `.fna` files in your `assembly/` directory.  
**Output:** Per-sample and merged TSV tables, phylogenetic trees, and publication-quality R figures.

---

## Pipeline Overview

```
 Genome assemblies (.fna)
        │
        ├─ 01 ─ Prokka annotation (GFF, FAA, GBK)
        ├─ 02 ─ MLST (Pasteur + Oxford)
        ├─ 03 ─ ARG prediction (CARD / AMRFinder / ResFinder)
        ├─ 04 ─ Virulence factor genes (VFDB)
        ├─ 05 ─ Plasmid detection (geNomad + MOB-suite)
        ├─ 06 ─ IS elements (ISEScan) + IS-ARG linkage
        ├─ 07 ─ Prophage detection (geNomad viral)
        ├─ 08 ─ SNP phylogeny (snippy → IQ-TREE 2)
        ├─ 09 ─ AbaRI resistance island detection
        ├─ 10 ─ Pan-genome (Roary) ← NEW
        ├─ 11 ─ GTDB-Tk + temporal ARG trends ← NEW
        └─ 12 ─ Master table + R figures
```

---

## Requirements

### Operating System
Linux (Ubuntu 20.04+ recommended). macOS may work but is untested.

### Conda / Mamba
All tools are managed via conda environments. Install [Miniforge](https://github.com/conda-forge/miniforge) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.

```bash
# Recommended: install mamba for faster environment creation
conda install -n base -c conda-forge mamba
```

### Conda Environments and Tool Versions

| Environment | Key Tools | Versions |
|-------------|-----------|----------|
| `annotation_env` | Prokka, Prodigal, Barrnap, tRNAscan-SE | Prokka 1.14.x |
| `clinical_env` | Abricate, AMRFinderPlus, ISEScan, MOB-suite, geNomad, BLAST+, seqkit, entrez-direct | Abricate 1.0.x, AMRFinder 3.12.x, ISEScan 1.7.x, geNomad 1.8.x |
| `snippy2026` | Snippy, samtools, htslib, freebayes, snp-dists, snp-sites | Snippy 4.0.2, samtools 1.23 |
| `phylogeny_env` | IQ-TREE 2, FastTree, MAFFT | IQ-TREE 2.3.x, FastTree 2.1.x |
| `roary_env` | Roary | Roary 3.13.x |
| `assembly_binning_sn_env` | GTDB-Tk | GTDB-Tk 2.6.1 |
| `r_env` | R, tidyverse, ggpubr, pheatmap, viridis, corrplot, ggrepel, ggtree, ggtreeExtra, rstatix, patchwork, ggridges, RColorBrewer, vegan, ggalluvial, treemapify, ggbeeswarm, lme4, broom, effectsize, ape, phytools, ggcorrplot | R 4.4.x |

### Create Environments

```bash
# annotation
mamba create -n annotation_env -c bioconda -c conda-forge prokka

# clinical tools
mamba create -n clinical_env -c bioconda -c conda-forge \
    abricate amrfinderplus isescan mob-suite blast seqkit genomad \
    entrez-direct mmseqs2

# snippy environment
mamba create -n snippy2026 -c bioconda -c conda-forge \
    snippy=4.0.2 samtools=1.23 htslib=1.23 snp-dists snp-sites

# phylogeny
mamba create -n phylogeny_env -c bioconda -c conda-forge \
    iqtree mafft fasttree

# Roary
mamba create -n roary_env -c bioconda -c conda-forge roary

# GTDB-Tk
mamba create -n assembly_binning_sn_env -c bioconda -c conda-forge gtdbtk=2.6.1

# R environment
mamba create -n r_env -c conda-forge r-base=4.4 r-tidyverse r-ggpubr \
    r-pheatmap r-viridis r-corrplot r-ggrepel r-rstatix r-scales \
    r-patchwork r-ggridges r-rcolorbrewer r-vegan r-ggalluvial \
    r-treemapify r-ggbeeswarm r-lme4 r-broom r-effectsize r-ape \
    r-phytools r-ggcorrplot

# Bioconductor packages (ggtree, ggtreeExtra — install inside R after env creation)
conda run -n r_env Rscript -e \
    "BiocManager::install(c('ggtree','ggtreeExtra'), ask=FALSE)"
```

---

## Database Setup

All databases should be placed under a single root directory (set as `DB_BASE` in `config.sh`).

### CARD Database (ARG prediction — Abricate)
```bash
# Abricate downloads CARD automatically:
conda run -n clinical_env abricate --setupdb
# Or manually update:
conda run -n clinical_env abricate --list
```

### AMRFinderPlus Database
```bash
conda run -n clinical_env amrfinder --update
# Default location: $CONDA_PREFIX/share/amrfinderplus/data/latest/
# Copy to: $DB_BASE/amrfinderplus_db/latest/
```

### VFDB (Virulence Factor Database)
```bash
mkdir -p $DB_BASE/vfdb
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz -P $DB_BASE/vfdb/
gunzip $DB_BASE/vfdb/VFDB_setB_nt.fas.gz
```

### geNomad Database (Plasmids + Prophages)
```bash
mkdir -p $DB_BASE/genomad_db
conda run -n clinical_env genomad download-database $DB_BASE/genomad_db/
# Expected path: $DB_BASE/genomad_db/genomad_db/
```

### MOB-suite Database
```bash
conda run -n clinical_env mob_init
# Default: installed to the conda environment
```

### GTDB-Tk Database (R226)
```bash
# ~85 GB — download from GTDB
mkdir -p $DB_BASE/GTDB
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar -xzf gtdbtk_r226_data.tar.gz -C $DB_BASE/GTDB/
```

### Reference Genome (Step 08 — SNP Phylogeny)
The pipeline uses *A. baumannii* ATCC 19606 (GCF_000015425.1) as the SNP reference, matching the paper.

```bash
mkdir -p $PROJECT_DIR/references
# Download from NCBI:
datasets download genome accession GCF_000015425.1 --filename ATCC19606.zip
unzip ATCC19606.zip
cp ncbi_dataset/data/GCF_000015425.1/*.fna $PROJECT_DIR/references/ATCC_19606_genomic.fna
```

Or with conda entrez-direct:
```bash
conda run -n clinical_env efetch \
    -db nuccore -id CP000521.1 -format fasta \
    > $PROJECT_DIR/references/ATCC_19606_genomic.fna
```

---

## Configuration

**Step 1:** Copy and edit the configuration file:
```bash
cp config.sh config.sh.bak   # optional backup
nano config.sh               # or use any text editor
```

**Key variables to set:**

```bash
# Where your project lives (results will be created here)
export PROJECT_DIR="$HOME/acinetobacter_genomics"

# Where your genome assemblies (.fna files) are
export ASSEMBLY_DIR="$PROJECT_DIR/assembly"

# Root directory for all databases
export DB_BASE="$HOME/databases/bacteria_microbial_genomics"

# Path to ATCC 19606 reference genome (needed for step 08)
export REFERENCE_GENOME="$PROJECT_DIR/references/ATCC_19606_genomic.fna"

# Number of CPU threads
export THREADS=8
```

**Step 2:** Place genome assemblies in `$ASSEMBLY_DIR/`:
```
assembly/
├── SAMPLE1.fna
├── SAMPLE2.fna
└── ...
```

**Step 3:** Run setup (creates directory structure and sample list):
```bash
bash step00_setup.sh
```

**Step 4:** Edit `$PROJECT_DIR/metadata.tsv` to fill in sample information (collection year, hospital, location, specimen type). The `sequence_type` column will be filled automatically after step 02.

---

## Usage

```bash
# Setup (always run first)
bash step00_setup.sh

# Run the full pipeline
bash run_pipeline.sh

# Run a single step (e.g., step 3 only)
bash run_pipeline.sh --step 3

# Run from a specific step onwards (e.g., resume from step 6)
bash run_pipeline.sh --from 6

# Dry-run (show what would be executed, without running)
bash run_pipeline.sh --dry-run

# Or run individual step scripts directly
source config.sh
bash step03_ARG_prediction.sh
```

### Quick Start Example

```bash
git clone https://github.com/YOUR_USERNAME/acinetobacter-baumannii-genomics.git
cd acinetobacter-baumannii-genomics

# Edit config.sh with your paths
nano config.sh

# Place assemblies
cp /path/to/your/*.fna $PROJECT_DIR/assembly/

# Run setup
bash step00_setup.sh

# Edit metadata
nano $PROJECT_DIR/metadata.tsv

# Run pipeline
bash run_pipeline.sh
```

---

## Step-by-Step Description

### Step 00 — Setup
Creates all output directories and generates `samples.txt` and `metadata.tsv` template. Also auto-generates `env_config.sh` which is sourced by all downstream steps.

### Step 01 — Genome Annotation (`annotation_env`)
Annotates all genomes with **Prokka 1.14.x** using *Acinetobacter baumannii* genus/species flags. Produces GFF, GBK, FAA, FFN per sample. GFF files are copied to `results/gff/` for use by Roary (Step 10).

**Output:** `results/01_annotation/`, `results/13_summary_tables/annotation_summary.tsv`

### Step 02 — MLST (`clinical_env`)
Runs **mlst (Seemann)** with both the **Pasteur** (abaumannii_2, 7 loci: cpn60, fusA, gltA, pyrG, recA, rplB, rpoB) and **Oxford** (abaumannii, 7 loci: gltA, gyrB, gdhB, recA, cpn60, gpi, rpoD) schemes against PubMLST databases. Matches the paper's Pasteur scheme usage.

**Output:** `results/02_mlst/mlst_pasteur.tsv`, `mlst_oxford.tsv`, `results/13_summary_tables/mlst_combined.tsv`

### Step 03 — ARG Prediction (`clinical_env`)
Three-tool approach for maximum confidence:
- **Abricate + CARD** (≥80% identity, ≥80% coverage) — mirrors paper's RGI/CARD approach
- **AMRFinderPlus** (NCBI, with *A. baumannii* organism flag) — broader database
- **Abricate + ResFinder** (≥90% identity, ≥60% coverage) — acquired genes only

**Output:** `results/03_ARG_RGI/`, `results/04_ARG_amrfinder/`, `results/14_abricate_resfinder/`, `results/13_summary_tables/ARG_counts.tsv`

### Step 04 — VFG Prediction (`clinical_env`)
- **Abricate + VFDB** (≥70% identity, ≥60% coverage) — mirrors paper's VFanalyzer approach
- **BLASTN vs VFDB** (local, same thresholds) — for more sensitive detection
- Generates a category map (adherence, biofilm, iron uptake, immune evasion, etc.)

**Output:** `results/05_VFG_abricate/`, `results/13_summary_tables/VFG_presence_absence.tsv`

### Step 05 — Plasmid Detection (`clinical_env`)
- **geNomad 1.8.x** neural-network-based plasmid/chromosome classification (min-score 0.7)
- **MOB-suite** replicon typing for plasmid classification
- ARG screening on plasmid sequences with Abricate/CARD

**Note:** `mmseqs2` is required by geNomad; the script installs it automatically if missing.

**Output:** `results/06_plasmid_genomad/`, `results/11_mob_suite/`, `results/13_summary_tables/plasmid_summary.tsv`

### Step 06 — IS Element Detection (`clinical_env`)
- **ISEScan** ORF-prediction-based IS element identification
- **IS-ARG linkage analysis**: identifies ARGs within 5 kb of IS elements (key mechanism for ARG mobilisation in *A. baumannii*)
- Includes an automatic patch for the fastcluster/numpy ≥1.24 incompatibility

**Output:** `results/07_IS_isescan/`, `results/13_summary_tables/IS_family_frequency.tsv`, `IS_ARG_linkage.tsv`

### Step 07 — Prophage Detection (`clinical_env`)
- **geNomad** viral-mode classification — reuses the database from Step 05
- Checks prophage regions for embedded ARGs
- Provides an optional PHASTER API submission script for cross-validation with the paper's exact tool

**Output:** `results/08_prophage_genomad/`, `results/13_summary_tables/prophage_summary.tsv`, `prophage_ARG_check.tsv`

### Step 08 — SNP Phylogeny (`snippy2026` + `phylogeny_env`)
- **Snippy 4.0.2** maps each assembly against ATCC 19606 (matching the paper's reference)
- **snippy-core** extracts the core SNP alignment
- **snp-sites** strips constant sites
- **IQ-TREE 2** (GTR+F+ASC model, 1000 UFBoot) builds the ML tree
- **FastTree** provides a rapid pre-visualisation tree
- **snp-dists** computes pairwise SNP distance matrix
- Generates **iTOL** annotation files (ST colours, metadata labels)

**Note:** Includes an automatic patch for the Snippy 4.0.2 samtools version-string comparison bug.

**Output:** `results/09_SNP_phylogeny/iqtree/acinetobacter_snp.treefile`, `core_snps/snp_distance_matrix.tsv`, `iqtree/iTOL_annotations/`

### Step 09 — AbaRI Resistance Island Detection (`clinical_env`)
Implements the flanking-region method to match Kumkar et al. 2022:
1. BLAST *comM* (ATCC 19606 locus A1S_1494) against each genome
2. Extract ±30 kb flanking window around comM
3. BLAST window vs AbaR3, AbGRI1, AbaR4a reference islands
4. BLAST window vs 7 ARG marker sequences (blaOXA-23, tetB, strA, sul2, aadA1, ISAba1, blaOXA-58)
5. Call AbaRI if island hit OR ≥2 markers present

Reference sequences are fetched from NCBI on first run (requires internet, then cached).

**Output:** `results/15_AbaRI_comM/`, `results/13_summary_tables/AbaRI_summary.tsv`

### Step 10 — Pan-genome Analysis (`roary_env`) ⬅ NOVEL
**Roary** (95% identity threshold, core-gene alignment enabled) quantifies:
- Core genome (≥99% of isolates), soft-core (95–99%), shell (15–95%), cloud (<15%)
- Heaps' law alpha (open vs. closed pan-genome character)
- Core genome phylogeny (IQ-TREE 2, GTR+G, 1000 UFBoot)
- Links accessory genes to ARGs from Step 03

**Input:** GFF files from Step 01 (copied to `results/gff/`).

**Output:** `results/10_pangenome_roary/`, `results/13_summary_tables/accessory_ARG_genes.tsv`

### Step 11 — GTDB-Tk + Temporal ARG Analysis (`assembly_binning_sn_env`) ⬅ NOVEL
**Part A — GTDB-Tk 2.6.1** (GTDB R226): Places all genomes in the GTDB taxonomy framework using 120 conserved bacterial marker genes, producing a reference-anchored phylogenomic tree independent of the SNP approach.

**Part B — Temporal ARG trends**: Cross-references ARG data with `metadata.tsv` collection years to calculate:
- Mean acquired ARGs per year (excluding intrinsic genes: blaOXA-51, ADC, adeABC/IJ/K, etc.)
- blaOXA-23 prevalence per year
- Hospital-level ARG burden

**Output:** `results/12_gtdbtk/`, `results/13_summary_tables/temporal_ARG_trends.tsv`

### Step 12 — Summary Tables + Advanced Visualisation (`r_env`)
Merges all step outputs into `master_genome_table.tsv`, then executes `advanced_analysis.R` to produce **23 publication-quality figures** and a complete statistics report.

#### Figures produced

| Figure | Content | Key statistics |
|--------|---------|----------------|
| Fig 01 | ARG burden by sequence type — violin + box + jitter | Kruskal-Wallis + Dunn post-hoc |
| Fig 02 | ARG burden by hospital — barcode range plot | Kruskal-Wallis |
| Fig 03 | ARG burden by specimen type — boxplot | Kruskal-Wallis |
| Fig 04 | Temporal ARG dynamics per hospital — line chart with blaOXA-23 emergence marker | — |
| Fig 05 | ARG distribution shift over years — ridge plot | — |
| Fig 06 | ARG by era (Early/Mid/Recent) + specimen overlay | Kruskal-Wallis + Dunn |
| Fig 07 | Genome size vs ARG count — scatter + regression | Spearman ρ |
| Fig 08 | Plasmid count vs ARG count — scatter + regression | Spearman ρ |
| Fig 09 | IS–ARG physical distance histogram + density | % pairs within 1 kb / 500 bp |
| Fig 10 | IS family prevalence — lollipop (size = mean copy number) | — |
| Fig 11 | IS burden vs ARG count per isolate | Spearman ρ |
| Fig 12 | AMR gene presence/absence heatmap — annotated by ST, hospital, era | Ward D2 clustering |
| Fig 13 | Plasmid mobility class by ST — stacked bar | MOB-suite |
| Fig 14 | Prophage burden per isolate — ordered bar | — |
| Fig 15 | Prophage ARG carriage — donut chart | — |
| Fig 16 | VFG count by ST + specimen overlay — violin | — |
| Fig 17 | ARG vs VFG correlation — scatter + regression | Spearman ρ |
| Fig 18 | Spearman correlation matrix (all genomic features) | p-value masked |
| Fig 19 | Epidemiological flow: Hospital → ST → Specimen — alluvial/Sankey | — |
| Fig 20 | ST distribution over time — bubble chart | — |
| Fig 21 | Isolate distribution by hospital and specimen — treemap | — |
| Fig 22 | Independent ARG predictors — forest plot | Multiple linear regression R² |
| Fig 23 | SNP phylogeny annotated with ST, hospital, year, ARG bar | ggtree + ggtreeExtra |

#### Statistics report (`advanced_statistics.txt`)
Kruskal-Wallis + Dunn tests (ARG ~ ST, era, specimen), Spearman correlations (plasmid × ARG, IS × ARG, VFG × ARG, genome size × ARG, prophage × ARG), conjugative vs mobilizable plasmid size (Wilcoxon), blaOXA-23 χ² by hospital, multivariable linear regression coefficients, Cohen's d effect sizes, and per-ST descriptive summaries.

#### Re-running figures independently
```bash
# After the pipeline has run at least once:
export PROJECT_DIR=/your/project/dir
export RESULTS_DIR=/your/project/dir/results
conda run -n r_env Rscript advanced_analysis.R

# Or with positional arguments:
conda run -n r_env Rscript advanced_analysis.R /your/project /your/project/results
```

#### Extra input file required for Fig 12
`advanced_analysis.R` requires `results/13_summary_tables/amr_matrix.csv` — a wide-format presence/absence matrix (rows = samples, columns = ARG gene names, values = 0/1). Generate it from the Abricate CARD summary:

```bash
# Convert abricate summary to presence/absence matrix
conda run -n clinical_env abricate --summary \
    results/03_ARG_RGI/abricate_card_raw.tsv \
    | awk 'BEGIN{OFS="\t"} NR==1{print} NR>1{
        n=split($1,a,"/"); name=a[n]; sub(/\.fna$/,"",name); $1=name; print}' \
    > results/13_summary_tables/amr_matrix.csv
```

**Output:** `results/13_summary_tables/master_genome_table.tsv`, `results/ADVANCED_FIGURES/*.pdf`, `results/ADVANCED_FIGURES/*.png`, `results/ADVANCED_FIGURES/advanced_statistics.txt`

---

## Outputs

```
$PROJECT_DIR/
├── samples.txt                          ← Auto-generated sample list
├── metadata.tsv                         ← Fill in manually
├── env_config.sh                        ← Auto-generated; do not edit
└── results/
    ├── 00_logs/                         ← Log files for every step
    ├── 01_annotation/                   ← Prokka output per sample
    ├── 02_mlst/                         ← MLST tables (Pasteur + Oxford)
    ├── 03_ARG_RGI/                      ← Abricate CARD results
    ├── 04_ARG_amrfinder/                ← AMRFinderPlus results
    ├── 05_VFG_abricate/                 ← VFDB results + BLAST
    ├── 06_plasmid_genomad/              ← geNomad plasmid output
    ├── 07_IS_isescan/                   ← ISEScan output per sample
    ├── 08_prophage_genomad/             ← geNomad viral output
    ├── 09_SNP_phylogeny/                ← snippy, core SNPs, IQ-TREE
    │   └── iqtree/
    │       └── acinetobacter_snp.treefile  ← used by Fig 23
    ├── 10_pangenome_roary/              ← Roary output + core tree
    ├── 11_mob_suite/                    ← MOB-suite replicon typing
    ├── 12_gtdbtk/                       ← GTDB-Tk classification
    ├── 13_summary_tables/               ← All merged tables + master table
    │   ├── master_genome_table.tsv      ← Master results table (all steps merged)
    │   ├── master_merged_data.csv       ← R-ready merged table (output of step 12)
    │   ├── metadata.tsv                 ← Copy placed here for R convenience
    │   ├── mlst_combined.tsv
    │   ├── ARG_counts.tsv
    │   ├── amr_matrix.csv               ← Required for Fig 12 (see step 12 notes)
    │   ├── VFG_presence_absence.tsv
    │   ├── plasmid_summary.tsv
    │   ├── AbaRI_summary.tsv
    │   ├── IS_family_frequency.tsv
    │   ├── IS_ARG_linkage.tsv
    │   ├── mob_suite_summary.tsv
    │   ├── prophage_summary.tsv
    │   ├── prophage_ARG_check.tsv
    │   └── temporal_ARG_trends.tsv
    ├── 14_abricate_resfinder/           ← ResFinder results
    ├── 15_AbaRI_comM/                   ← AbaRI flanking-region analysis
    └── ADVANCED_FIGURES/                ← 23 figures (PDF + PNG) + statistics
        ├── Fig01_ARG_by_ST.pdf/.png
        ├── Fig02_ARG_by_Hospital.pdf/.png
        ├── Fig03_ARG_by_Specimen.pdf/.png
        ├── Fig04_Temporal_ARG_Dynamics_by_Hospital.pdf/.png
        ├── Fig05_ARG_Ridge_Temporal.pdf/.png
        ├── Fig06_ARG_by_Era.pdf/.png
        ├── Fig07_Genome_vs_ARG.pdf/.png
        ├── Fig08_Plasmid_vs_ARG.pdf/.png
        ├── Fig09_IS_ARG_Distance.pdf/.png
        ├── Fig10_IS_Prevalence.pdf/.png
        ├── Fig11_IS_vs_ARG.pdf/.png
        ├── Fig12_AMR_Heatmap.pdf
        ├── Fig13_Plasmid_Mobility_ST.pdf/.png
        ├── Fig14_Prophage_Burden.pdf/.png
        ├── Fig15_Prophage_ARG_Donut.pdf/.png
        ├── Fig16_VFG_by_ST.pdf/.png
        ├── Fig17_ARG_vs_VFG.pdf/.png
        ├── Fig18_Correlation_Matrix.pdf
        ├── Fig19_Alluvial_Flow.pdf/.png
        ├── Fig20_MLST_Bubble.pdf/.png
        ├── Fig21_Treemap_Hospital_Specimen.pdf/.png
        ├── Fig22_Forest_ARG_Predictors.pdf/.png
        ├── Fig23_SNP_Phylogeny_Annotated.pdf/.png
        ├── Fig23_acinetobacter_midpoint_clean.nwk
        └── advanced_statistics.txt
```

---

## Known Issues and Fixes

The following bugs are patched automatically by the pipeline:

| Tool | Bug | Fix applied in |
|------|-----|----------------|
| ISEScan | `fastcluster.linkage()` passes `copy=None` which is rejected by numpy ≥1.24 | `step06_IS_isescan.sh` |
| Snippy 4.0.2 | Perl string comparison of version numbers causes false rejection of samtools ≥1.10 | `step08_SNP_phylogeny.sh` |

### Extra input file for Fig 12 (`amr_matrix.csv`)

Figure 12 (AMR heatmap) requires a wide-format presence/absence matrix. This is not automatically generated by the pipeline steps — create it with:

```bash
source config.sh
conda run -n clinical_env abricate \
    --summary $RESULTS_DIR/03_ARG_RGI/abricate_card_raw.tsv \
    | awk 'BEGIN{OFS="\t"} NR==1{print} NR>1{
        n=split($1,a,"/"); name=a[n]
        sub(/\.fna$/,"",name); $1=name; print}' \
    > $RESULTS_DIR/13_summary_tables/amr_matrix.csv
```

If this file is absent, `advanced_analysis.R` will skip Fig 12 gracefully with a warning rather than crashing.

---

## Citation

If you use this pipeline, please cite the original paper:

> Kumkar SN, Kamble EE, Chavan NS, Dhotre DP and Pardesi KR (2022)  
> Diversity of resistant determinants, virulence factors, and mobile genetic elements in *Acinetobacter baumannii* from India: A comprehensive in silico genome analysis.  
> *Front. Cell. Infect. Microbiol.* 12:997897. doi: 10.3389/fcimb.2022.997897

And the key tools:

- **Prokka:** Seemann T (2014) *Bioinformatics* 30(14):2068–9
- **mlst:** Seemann T et al. — https://github.com/tseemann/mlst
- **Abricate:** Seemann T — https://github.com/tseemann/abricate
- **AMRFinderPlus:** Feldgarden M et al. (2021) *Sci Rep* 11:16931
- **geNomad:** Camargo AP et al. (2023) *Nat Biotechnol* 41:1alloc
- **ISEScan:** Xie Z & Tang H (2017) *Bioinformatics* 33(18):2850–2852
- **MOB-suite:** Robertson J & Nash JH (2018) *Microb Genomics* 4(8)
- **Snippy:** Seemann T — https://github.com/tseemann/snippy
- **IQ-TREE 2:** Minh BQ et al. (2020) *Mol Biol Evol* 37(5):1530–1534
- **Roary:** Page AJ et al. (2015) *Bioinformatics* 31(22):3691–3
- **GTDB-Tk:** Chaumeil PA et al. (2022) *Bioinformatics* 38(23):5315–5316
- **VFDB:** Chen L et al. (2016) *Nucleic Acids Res* 44:D694–7
- **CARD:** Alcock BP et al. (2023) *Nucleic Acids Res* 51:D690–9

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contributing

Pull requests are welcome. For major changes, please open an issue first.

---

*Pipeline developed for replication and extension of Kumkar et al. 2022. Not affiliated with the original authors.*
