# Configuration Reference

This document provides comprehensive details on all configurable parameters in the Acinetobacter Defense Systems Pipeline.

## Configuration File Structure

The pipeline is controlled via `config/config.yaml`. All parameters are organized into logical sections for different analysis components.

## Core Configuration

### Sample Specification

```yaml
samples:
  - "GCF_000005825.2"   # NCBI accession numbers
  - "GCF_000746645.1"   # Add as many as needed
  - "GCF_001022155.1"   # Pipeline scales automatically
```

**Parameter Details:**
- **Type**: List of strings
- **Format**: NCBI GenBank accession numbers (GCF_XXXXXXXXX.X)
- **Validation**: Must be valid, accessible NCBI accessions
- **Scale**: Tested with 1-200 genomes

## Tool-Specific Parameters

### DefenseFinder Configuration

```yaml
defensefinder:
  models_dir: "resources/defensefinder_models"  # Local model directory
  coverage_profile: "DEFAULT"                   # Coverage threshold profile
  custom_models: false                          # Use additional custom models
```

**Parameter Details:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `models_dir` | String | `resources/defensefinder_models` | Path to DefenseFinder HMM models |
| `coverage_profile` | String | `DEFAULT` | Coverage profile: `DEFAULT`, `HIGH_COV`, `LOW_COV` |
| `custom_models` | Boolean | `false` | Enable custom/experimental models |

### ResFinder Configuration

```yaml
resfinder:
  threshold: 90.0        # Identity threshold (%)
  coverage: 60.0         # Coverage threshold (%)
  database_path: "resources/resfinder_db"
  species: "acinetobacter"
  acquired_only: true    # Focus on acquired resistance genes
```

**Parameter Details:**

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `threshold` | Float | 70.0-100.0 | Minimum sequence identity (%) |
| `coverage` | Float | 50.0-100.0 | Minimum query coverage (%) |
| `database_path` | String | - | Path to ResFinder database |
| `species` | String | - | Target species for analysis |
| `acquired_only` | Boolean | - | Exclude intrinsic resistance genes |

**Recommended Values by Use Case:**

| Use Case | Threshold | Coverage | Rationale |
|----------|-----------|----------|-----------|
| **High Sensitivity** | 80.0 | 50.0 | Detect divergent resistance genes |
| **Balanced** | 90.0 | 60.0 | Standard clinical relevance |
| **High Specificity** | 95.0 | 80.0 | Only high-confidence hits |

### BLAST Configuration

```yaml
blast:
  evalue: 0.005          # E-value threshold
  identity: 80           # Minimum identity (%)
  coverage: 70           # Minimum coverage (%)
  max_target_seqs: 1000  # Maximum hits per query
```

**Parameter Details:**

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `evalue` | Float | 1e-10 - 0.1 | Statistical significance threshold |
| `identity` | Integer | 50-100 | Minimum sequence identity (%) |
| `coverage` | Integer | 30-100 | Minimum query coverage (%) |
| `max_target_seqs` | Integer | 100-10000 | Maximum number of hits to report |

### HMRG Analysis Configuration

```yaml
hmrg:
  evalue: 0.005                    # E-value threshold
  qcov_hsp_perc: 80               # Query coverage per HSP (%)
  database: "BacMet2_EXP"         # Database selection
  metal_categories: "all"          # Metal resistance categories
```

**Database Options:**

| Database | Description | Gene Count | Focus |
|----------|-------------|------------|-------|
| `BacMet2_EXP` | Experimentally verified | ~2,700 | High confidence |
| `BacMet2_PRED` | Predicted genes | ~26,000 | Comprehensive |
| `Custom` | User-provided | Variable | Specialized analysis |

### CRISPRCasFinder Configuration

```yaml
crisprcasfinder:
  installation_path: "resources/CRISPRCasFinder"
  combined_analysis: true          # Analyze genomes together
  keep_intermediate: true          # Retain intermediate files
  cas_detection: true              # Enable Cas gene detection
```

## Advanced Configuration

### Resource Management

```yaml
resources:
  default_cores: 1                 # Default CPU cores per job
  default_mem_mb: 4000            # Default memory per job (MB)
  max_concurrent_jobs: 8           # Maximum parallel jobs
  download_retries: 3              # Network retry attempts
```

### Quality Control

```yaml
quality_control:
  min_genome_size: 1000000        # Minimum genome size (bp)
  max_contigs: 1000               # Maximum contigs per genome
  validate_downloads: true         # Verify downloaded genomes
  remove_duplicates: true          # Remove duplicate systems
```

### Output Configuration

```yaml
output:
  consolidate_results: true        # Create merged output files
  generate_summaries: true         # Create summary statistics
  create_visualizations: true      # Generate plots and figures
  compression: "gzip"              # Compress large output files
```

## Environment-Specific Configuration

### Local Workstation

```yaml
# Optimized for single-user workstation
execution:
  cores: 8
  memory_limit: "32GB"
  temp_dir: "/tmp/snakemake"
  
resources:
  defensefinder: {cores: 2, mem_mb: 8000}
  blast: {cores: 4, mem_mb: 4000}
  padloc: {cores: 2, mem_mb: 6000}
```

### HPC Cluster

```yaml
# Optimized for cluster execution
execution:
  cores: 100
  memory_limit: "500GB"
  temp_dir: "/scratch/$USER/snakemake"
  
cluster:
  partition: "compute"
  account: "your_account"
  time_limit: "24:00:00"
```

### Cloud Computing

```yaml
# Optimized for cloud instances
execution:
  cores: 16
  memory_limit: "64GB"
  storage_type: "ssd"
  
cloud:
  instance_type: "compute-optimized"
  auto_scaling: true
  cost_optimization: true
```

## Performance Tuning

### Small Datasets (1-10 genomes)

```yaml
performance:
  cores_per_job: 1
  memory_per_job: 4000
  concurrent_downloads: 2
  blast_threads: 2
```

### Medium Datasets (10-50 genomes)

```yaml
performance:
  cores_per_job: 2
  memory_per_job: 8000
  concurrent_downloads: 4
  blast_threads: 4
```

### Large Datasets (50+ genomes)

```yaml
performance:
  cores_per_job: 4
  memory_per_job: 16000
  concurrent_downloads: 8
  blast_threads: 8
  enable_checkpointing: true
```

## Validation and Testing

### Test Configurations

```yaml
# Quick test configuration
test:
  samples: ["GCF_000005825.2"]     # Single reference genome
  quick_mode: true                  # Skip time-intensive steps
  validate_only: true               # Validate setup without full analysis
```

### Benchmark Configuration

```yaml
# Performance benchmarking
benchmark:
  enable_timing: true               # Record execution times
  enable_memory_tracking: true      # Monitor memory usage
  enable_profiling: true            # Detailed performance profiling
  output_dir: "benchmarks/"        # Benchmark results location
```

## Configuration Validation

### Automatic Validation

The pipeline includes built-in configuration validation:

```bash
# Validate configuration
snakemake --configfile config/config.yaml --dry-run

# Check parameter ranges
python workflow/scripts/validate_config.py config/config.yaml
```

### Manual Validation Checklist

- [ ] All NCBI accessions are valid and accessible
- [ ] Parameter values are within acceptable ranges
- [ ] Required directories exist and are writable
- [ ] Database files are present and readable
- [ ] Resource allocations match system capabilities

## Configuration Templates

### Research-Grade Analysis

```yaml
# High-quality research configuration
samples: [...]  # Your genome list

resfinder: {threshold: 90.0, coverage: 60.0}
blast: {evalue: 0.005, identity: 80}
hmrg: {evalue: 0.005, qcov_hsp_perc: 80}

quality_control:
  validate_downloads: true
  remove_duplicates: true
  min_genome_size: 1000000

output:
  consolidate_results: true
  generate_summaries: true
  create_visualizations: true
```

### High-Throughput Screening

```yaml
# Optimized for large-scale screening
samples: [...]  # Large genome list

# Relaxed parameters for broader detection
resfinder: {threshold: 85.0, coverage: 50.0}
blast: {evalue: 0.01, identity: 75}

# Performance optimization
resources:
  max_concurrent_jobs: 20
  enable_checkpointing: true
  compress_outputs: true
```

## Environment Variables

### System-Level Configuration

```bash
# Set environment variables for optimization
export SNAKEMAKE_CORES=8
export SNAKEMAKE_MEMORY=32000
export TMPDIR=/fast/tmp/directory

# Database locations
export RESFINDER_DB_PATH=/databases/resfinder
export DEFENSEFINDER_MODELS=/models/defensefinder
```

### Tool-Specific Variables

```bash
# DefenseFinder optimization
export DEFENSEFINDER_WORKERS=4
export HMMER_NCPU=2

# BLAST optimization
export BLAST_NUM_THREADS=4
export BLASTDB_LMDB_MAP_SIZE=1000000
```

---

**Related Documentation:**
- [Installation Guide](INSTALLATION.md) - Setup instructions
- [Tutorial](TUTORIAL.md) - Usage examples  
- [Troubleshooting](TROUBLESHOOTING.md) - Common issues