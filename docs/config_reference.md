# Configuration Reference

This document provides comprehensive details on all configurable parameters in the Acinetobacter Defense Systems Pipeline.

## Configuration File Structure

The pipeline is controlled via `config/config.yaml`. All parameters are organized into logical sections for different analysis components.

## Core Configuration

### Sample Specification

```yaml
samples:
  - "NZ_CP019041.1"   # NCBI accession numbers
  - "NZ_CP048849.1"   # Add as many as needed
  - "NZ_CP104347.1"   # Pipeline scales automatically
```

**Parameter Details:**
- **Type**: List of strings
- **Format**: NCBI GenBank accession numbers (NZ_XXXXXXXXX.X)
- **Validation**: Must be valid, accessible NCBI accessions
- **Scale**: Tested with 1-350 genomes

## Tool-Specific Parameters

### ResFinder Configuration

```yaml
resfinder:
  threshold: 0.9        # Identity threshold (%)
  coverage: 0.6         # Coverage threshold (%)
  database_path: "resources/resfinder_db"
  species: "acinetobacter"
  acquired_only: true    # Focus on acquired resistance genes
```

**Parameter Details:**

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `threshold` | Float | 0.7-1.0 | Minimum sequence identity (%) |
| `coverage` | Float | 0.5-1.0 | Minimum query coverage (%) |
| `database_path` | String | - | Path to ResFinder database |
| `species` | String | - | Target species for analysis |
| `acquired_only` | Boolean | - | Exclude intrinsic resistance genes |

**Recommended Values by Use Case:**

| Use Case | Threshold | Coverage | Rationale |
|----------|-----------|----------|-----------|
| **High Sensitivity** | 0.8 | 0.5 | Detect divergent resistance genes |
| **Balanced** | 0.9 | 0.6 | Standard clinical relevance |
| **High Specificity** | 0.95 | 0.8 | Only high-confidence hits |

### BLAST Configuration

```yaml
blast:
  evalue: 1e-6          # E-value threshold
  identity: 80           # Minimum identity (%)
  coverage: 80           # Minimum coverage (%)
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
| `BacMet2_EXP` | Experimentally verified | ~753 | High confidence |
| `BacMet2_PRED` | Predicted genes | ~2,700 | Comprehensive |
| `Custom` | User-provided | Variable | Specialized analysis |

### CRISPRCasFinder Configuration

```yaml
crisprcasfinder:
  installation_path: "resources/CRISPRCasFinder"
  combined_analysis: true          # Analyze genomes together
  keep_intermediate: true          # Retain intermediate files
  cas_detection: true              # Enable Cas gene detection
```

### Automatic Validation

The pipeline includes built-in configuration validation:

```bash
# Validate configuration
snakemake --configfile config/config.yaml --dry-run
```

### Manual Validation Checklist

- [ ] All NCBI accessions are valid and accessible
- [ ] Parameter values are within acceptable ranges
- [ ] Required directories exist and are writable
- [ ] Database files are present and readable
- [ ] Resource allocations match system capabilities

---

**Related Documentation:**
- [Installation Guide](INSTALLATION.md) - Setup instructions
- [Tutorial](TUTORIAL.md) - Usage examples  
- [Troubleshooting](TROUBLESHOOTING.md) - Common issues
