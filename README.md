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
- Multi-tool integration with dependency management via Conda
- Modular design allowing individual component execution
- Automated quality control and result consolidation

### **Comprehensive Analysis Suite**
- **Defence System Prediction**: DefenseFinder v2.0.2 + PADLOC v2.0.0
- **CRISPR-Cas Detection**: Integrated CRISPRCasFinder analysis
- **Resistance Gene Screening**: ResFinder with custom databases
- **Mobile Element Analysis**: BLAST-based identification pipelines
- **Statistical Framework**: Correlation analysis with multiple testing correction

### **Production-Ready Implementation**
- Reproducible environments using Conda
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

2. **Configure your analysis:**
```bash
# Edit the configuration file with your genome accessions
nano config/config.yaml
```

3. **Run the complete pipeline:**
```bash
# Dry run to check workflow
snakemake --dry-run

# Execute with 8 cores
snakemake --cores 8 --use-conda

# For cluster execution
snakemake --cluster "sbatch" --cores 50
```

### Configuration

The pipeline is controlled via `config/config.yaml`:

```yaml
samples:
  - "GCF_000005825.2"  # Example: A. baumannii ATCC 17978
  - "GCF_000746645.1"  # Add your genome accessions here

# Tool parameters
resfinder:
  threshold: 0.9    # Identity threshold
  coverage: 0.6     # Minimum coverage

blast:
  evalue: 0.005      # E-value cutoff
  identity: 80       # Minimum identity percentage
```

## Technical Implementation

### **Workflow Management**
- **Snakemake framework** for reproducible, scalable analysis
- **Conda integration** for automatic environment management
- **Parallel execution** with intelligent resource allocation
- **Checkpointing** for robust long-running analyses

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
├── Snakefile                      # Main workflow definition
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
```

## Computational Requirements

### **Minimum System Requirements**
- **CPU**: 4+ cores (8+ recommended for large datasets)
- **RAM**: 8GB (16GB recommended) 
- **Storage**: 50GB+ (scales with dataset size)
- **OS**: Linux/macOS (tested on Ubuntu 24.04 LTS)

### **Estimated Runtime**
- **Single genome**: 5-7 minutes

## Scientific Applications

This pipeline enables researchers to:
- Systematically catalog defence system diversity across bacterial species
- Investigate correlations between bacterial immunity and resistance mechanisms 
- Compare tool performance for defence system prediction
- Analyze mobile element association patterns
- Generate publication-ready statistical analyses and visualizations


## Documentation

- **[Installation Guide](docs/INSTALLATION.md)** - Detailed setup instructions
- **[Tutorial](docs/TUTORIAL.md)** - Step-by-step analysis workflow
- **[Configuration Reference](docs/CONFIG.md)** - Complete parameter documentation
- **[Troubleshooting](docs/TROUBLESHOOTING.md)** - Common issues and solutions

## 🤝 Contributing

Contributions are welcome! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on:
- Code style conventions
- Testing requirements
- Documentation standards
- Issue reporting procedures

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{muthuraman2025acinetobacter,
  title={Acinetobacter Defense Systems Analysis Pipeline},
  author={Muthuraman, Vigneshwaran and Shehreen, Saadlee},
  year={2025},
  url={https://github.com/vikos77/acinetobacter-defence-pipeline},
  version={1.0.0}
}
```

## Contact & Support

- **Lead Developer**: Vigneshwaran Muthuraman (vigneshwaran0594@gmail.com)
- **Supervisor**: Dr. Saadlee Shehreen (s.shehreen@tees.ac.uk)
- **Institution**: Teesside University, School of Health & Life Sciences

---

**⭐ If this pipeline is useful for your research, please consider starring the repository!**
