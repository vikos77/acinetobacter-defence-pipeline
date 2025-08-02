# Installation Guide

This guide provides detailed instructions for setting up the Acinetobacter Defense Systems Analysis Pipeline on various systems.

## System Requirements

### Hardware Requirements
- **CPU**: Minimum 4 cores, recommended 8+ cores
- **RAM**: Minimum 16GB, recommended 32GB+ for large datasets
- **Storage**: 50GB+ free space (scales with dataset size)
- **Network**: Stable internet connection for database downloads

### Software Requirements
- **Operating System**: Linux (Ubuntu 20.04+), macOS (10.15+), or WSL2 on Windows
- **Python**: 3.8+ (via Conda/Mamba)
- **R**: 4.4.0+
- **Git**: For repository management

## Installation Methods

### Method 1: Quick Setup (Recommended)

```bash
# 1. Clone the repository
git clone https://github.com/vikos77/acinetobacter-defence-pipeline.git
cd acinetobacter-defence-pipeline

# 2. Install Conda/Mamba (if not already installed)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

# 3. Set up bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# 4. Install Snakemake
mamba create -n snakemake snakemake>=6.0
conda activate snakemake

# 5. Test the installation
snakemake --help
```

### Method 2: Container-Based Setup

```bash
# Using Docker (if preferred)
docker pull snakemake/snakemake
docker run -it -v $(pwd):/data snakemake/snakemake
```

## Tool-Specific Installation

### DefenseFinder Setup

```bash
# Create dedicated environment
mamba create -n defensefinder defense-finder>=1.2.0
conda activate defensefinder

# Download and setup models
defense-finder update
defense-finder --help

# Verify installation
defense-finder --version
```

### PADLOC Installation

```bash
# Create PADLOC environment
mamba create -n padloc -c bioconda padloc>=1.0.1
conda activate padloc

# Test PADLOC installation
padloc --help
```

### ResFinder Setup

```bash
# Clone ResFinder database
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git resources/resfinder_db
cd resources/resfinder_db
python3 INSTALL.py
```

### R Environment Setup

```bash
# Install R and required packages
mamba create -n r_analysis r-base>=4.4.0 r-tidyverse r-ggplot2 r-readxl
conda activate r_analysis

# Launch R and install additional packages
R
```

```r
# Inside R console
install.packages(c("patchwork", "ggrepel", "ggupset", "RColorBrewer", 
                   "gridExtra", "viridis", "corrplot", "pheatmap"))
```

## Database Setup

### Required Databases

```bash
# Create database directories
mkdir -p resources/{ime_db,hmrg_db,resfinder_db}

# Download ICEberg database for IME analysis
curl -o resources/ime_db/ICE_proteins.faa \
     https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/IME_aa_all.fas

# Download BacMet database for HMRG analysis
wget -O resources/hmrg_db/BacMet2_EXP_database.fasta \
     https://bacmet.biomedicine.gu.se/downloads/BacMet2_EXP_database.fasta
```

## Configuration

### Configure Pipeline Settings

```bash
# Copy example configuration
cp config/config.yaml.example config/config.yaml

# Edit configuration file
nano config/config.yaml
```

### Example Configuration

```yaml
# Sample genome accessions (replace with your data)
samples:
  - "GCF_000005825.2"  # A. baumannii ATCC 17978
  - "GCF_000746645.1"  # A. baumannii strain

# Tool parameters
resfinder:
  threshold: 90.0      # Identity threshold (%)
  coverage: 60.0       # Coverage threshold (%)

blast:
  evalue: 0.005        # E-value cutoff
  identity: 80         # Minimum identity (%)

hmrg:
  evalue: 0.005
  qcov_hsp_perc: 80    # Query coverage threshold (%)
```

## Validation

### Test Installation

```bash
# Activate Snakemake environment
conda activate snakemake

# Dry run to validate workflow
snakemake --dry-run --cores 1

# Test with small dataset (should complete in ~10 minutes)
snakemake --cores 2 --use-conda results/defensefinder/GCF_000005825.2/
```

### Verify Tool Versions

```bash
# Check all tool versions
conda activate defensefinder && defense-finder --version
conda activate padloc && padloc --version
conda activate r_analysis && R --version
```

## Troubleshooting

### Common Issues

#### Issue 1: Conda Environment Conflicts
```bash
# Solution: Use mamba instead of conda for faster resolution
mamba install -n defensefinder defense-finder
```

#### Issue 2: Database Download Failures
```bash
# Solution: Manual download with retry
for i in {1..3}; do
  wget https://database-url && break
  sleep 5
done
```

#### Issue 3: Permission Errors
```bash
# Solution: Fix permissions
chmod +x workflow/scripts/*.sh
chmod 755 resources/
```

### Performance Optimization

#### For Large Datasets (100+ genomes)
```bash
# Use cluster execution
snakemake --cluster "sbatch --time=4:00:00 --mem=8G" --cores 50

# Or use local execution with resource limits
snakemake --cores 8 --resources mem_mb=32000
```

#### Memory Management
```bash
# Monitor memory usage
htop
# Adjust Snakemake resource allocation if needed
snakemake --cores 4 --resources mem_mb=16000
```

## Advanced Configuration

### Cluster Integration (SLURM)

```bash
# Create cluster configuration
cat > cluster.yaml << 'EOF'
__default__:
  partition: "compute"
  time: "2:00:00"
  mem: "8G"
  cores: 1

run_defensefinder:
  time: "4:00:00"
  mem: "16G"
  cores: 4
EOF

# Execute on cluster
snakemake --cluster-config cluster.yaml --cluster "sbatch" --cores 50
```

### Custom Environment Setup

```bash
# For specific institutional clusters
module load conda/latest
conda create --prefix ./envs/pipeline snakemake
```

## Verification Checklist

- [ ] Conda/Mamba installed and configured
- [ ] All required channels added
- [ ] Snakemake environment created
- [ ] Tool-specific environments created
- [ ] Databases downloaded and processed
- [ ] Configuration file customized
- [ ] Dry run completed successfully
- [ ] Test analysis completed

## Getting Help

If you encounter issues during installation:

1. **Check system requirements** - ensure adequate resources
2. **Verify network connectivity** - database downloads require stable internet
3. **Review log files** - Snakemake provides detailed error messages
4. **Contact support** - See main README for contact information

## Next Steps

After successful installation:
1. Review the [Tutorial](TUTORIAL.md) for usage examples
2. Customize [Configuration](CONFIG.md) for your specific analysis
3. Read [Best Practices](BEST_PRACTICES.md) for optimal performance