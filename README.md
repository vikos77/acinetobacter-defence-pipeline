# Acinetobacter Defence Systems Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://python.org)
[![R](https://img.shields.io/badge/R-4.4+-blue.svg)](https://r-project.org)

> **A comprehensive, scalable bioinformatics pipeline for analyzing bacterial defence systems, antimicrobial resistance genes, and mobile genetic elements in *Acinetobacter* species.**

## Pipeline Overview

This repository contains a production-ready **Snakemake workflow** that coordinates the analysis of bacterial defence systems across large genomic datasets. The pipeline integrates multiple state-of-the-art bioinformatics tools to provide comprehensive insights into microbial defence mechanisms and their relationship with antimicrobial resistance.

**Research Context**: *Acinetobacter* species are clinically significant pathogens with remarkable adaptability. This pipeline enables systematic investigation of defence system architectures and their potential correlations with resistance acquisition mechanisms.

##  Key Technical Features

### **Scalable Architecture**
- 1000+ genome analysis capability with automatic parallel processing
- Multi-tool integration with dependency management via Conda or Docker
- Modular design allowing individual component execution
- Automated quality control and result consolidation

### **Comprehensive Analysis Suite**
- **Defence System Prediction**: DefenseFinder v2.0.2 + PADLOC v2.0.0
- **CRISPR-Cas Detection**: Integrated CRISPRCasFinder analysis
- **Resistance Gene Screening**: ResFinder with custom databases
- **Mobile Element Analysis**: BLAST-based identification pipelines
- **Statistical Framework**: Correlation analysis with multiple testing correction

### **Production-Ready Implementation**
- Reproducible environments using Conda or Docker containers    
- Error handling with retry logic for network operations
- Resource optimization with configurable CPU/memory allocation
- Automated data consolidation from individual tool outputs

## Pipeline Architecture

```
Input: Genome Accessions → Automated Download → Parallel Analysis → Consolidated Results
    ↓                           ↓                    ↓                    ↓
Config File              NCBI E-utilities      Multi-tool Suite      Statistical Analysis
                                                     ↓
                              ┌─────────────────────────────────────┐
                              │  DefenseFinder  │  PADLOC + CRISPR  │
                              │  ResFinder      │  IME Detection    │
                              │  HMRG Analysis  │  Data Quality QC  │
                              └─────────────────────────────────────┘
```

##  Analysis Scope

| Component | Scale | Technical Implementation |
|-----------|-------|-------------------------|
| **Defence Systems** | 35+ system types | HMM-based prediction + CRISPRCasFinder integration |
| **Resistance Genes** | 2000+ gene database | Sequence similarity with coverage thresholds |
| **Mobile Elements** | ICEberg + BacMet databases | tBLASTn with stringent quality filters |
| **Statistical Analysis** | Correlation matrices | Fisher's exact tests + Spearman correlation |
| **Genome Processing** | 132 complete genomes | Automated download + quality validation |

## Quick Start

### Prerequisites
```bash
# Core requirements
conda >= 4.10
snakemake >= 6.0
R >= 4.4.0
```

### Installation

1. **Clone the repository:**
```bash
git clone https://github.com/vikos77/acinetobacter-defence-pipeline.git
cd acinetobacter-defence-pipeline
```

2. **Choose your deployment method:**

## Method 1: Conda Environment (Traditional)

**Install Snakemake:**
```bash
# Create a new conda environment
conda create -n snakemake_env -c conda-forge -c bioconda snakemake

# Activate the environment
conda activate snakemake_env

# Test the installation
snakemake --help
```

**Configure and run:**
```bash
# Edit the configuration file with your genome accessions
nano config/config.yaml

# Dry run to check workflow
snakemake --dry-run

# Execute with 8 cores
snakemake --cores 8 --use-conda
```

## Method 2: Docker Containers (Recommended for Reproducibility)

**Prerequisites (15 minutes):**
```bash
# Install Docker
sudo apt-get update
sudo apt-get install docker.io docker-compose

# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker
```

**Deployment (30-60 minutes depending on internet):**
```bash
# Build containers (only needs to be done once)
./build_containers.sh
```

**Run Analysis (30 minutes - several hours depending on data):**
```bash
# Configure your analysis
nano config/config.yaml

# Run the pipeline
snakemake --snakefile Snakefile.containerized --cores 8 --configfile config/config.yaml
```

### Configuration

The pipeline is controlled via `config/config.yaml`:

```yaml
samples:
  - "NZ_CP048014.1"  # Example: A. baumannii genome
  - "NZ_CP078045.1"  # Add your genome accessions here

# Tool parameters
resfinder:
  threshold: 0.9    # Identity threshold
  coverage: 0.6     # Minimum coverage

blast:
  evalue: 0.005      # E-value cutoff
  identity: 80       # Minimum identity percentage
```

### **Tool Integration Strategy**
```python
# Example: DefenseFinder with custom model directory
rule run_defensefinder:
    input: genome="genomes/{sample}.fna"
    output: "results/defensefinder/{sample}_systems.tsv"
    conda: "envs/defense.yaml"
    shell: """
        defense-finder run --models-dir resources/models 
                          --out-dir results/defensefinder/{wildcards.sample} 
                          {input.genome}
    """
```


## Repository Structure

```
├── Snakefile                      # Main workflow definition (conda)
├── Snakefile.containerized        # Docker workflow definition
├── containers/                    # Docker container definitions
├── build_containers.sh            # Container build script
├── config/
│   ├── config.yaml               # Pipeline configuration
│   └── samples.txt               # Genome accession list
├── workflow/
│   ├── envs/                     # Conda environment definitions
│   │   ├── defense.yaml          # DefenseFinder + dependencies
│   │   ├── blast.yaml            # BLAST + sequence analysis
│   │   └── r_analysis.yaml       # R statistical environment
│   ├── rules/                    # Modular workflow rules
│   └── scripts/                  # Analysis and visualization scripts
├── results/
│   ├── consolidated/             # Merged tool outputs
│   └── analysis/                 # Statistical results
└── docs/                         # Detailed documentation
└── resources/                    # Databases and Downloaded genomes
```
## Deployment Comparison

| Feature | Conda | Docker |
|---------|-------|--------|
| **Setup Time** | 15 min | 60 min (one-time) |
| **Reproducibility** | Good | Excellent |
| **Portability** | Local | Universal |
| **Maintenance** | Manual updates | Automated |
| **Best For** | Development | Production/Publication |

## **Estimated Runtime**
- **Single genome**: 5-7 minutes (depends on the genome size)

## Scientific Applications

This pipeline enables researchers to:
- Systematically catalog defence system diversity across bacterial species
- Investigate correlations between bacterial immunity and resistance mechanisms 
- Compare tool performance for defence system prediction
- Analyze mobile element association patterns
- Generate publication-ready statistical analyses and visualizations


## Documentation

- **[Installation Guide](docs/config_reference.md)** - Step-by-step analysis workflow


## Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{muthuraman2025acinetobacter,
  title={Adaptive trade-offs between niche-driven defence system selection and horizontal
  gene transfer suggests clinical success in Acinetobacter spp},
  author={Vigneshwaran Muthuraman, Proyash Roy, Paul Dean, Bruno Silvester Lopes,
  Saadlee Shehreen},
  year={2025},
  journal={bioRxiv},
  doi={10.1101/2025.08.01.668115},
  url={https://www.biorxiv.org/content/10.1101/2025.08.01.668115v1}
}
```

## Contact & Support

- **Lead Developer**: Vigneshwaran Muthuraman (vigneshwaran0594@gmail.com)
- **Supervisor**: Dr. Saadlee Shehreen (s.shehreen@tees.ac.uk)
- **Institution**: Teesside University, School of Health & Life Sciences

---

**⭐ If this pipeline is useful for your research, please consider starring the repository!**
