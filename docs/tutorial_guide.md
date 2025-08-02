# Tutorial: Running the Acinetobacter Defense Systems Pipeline

This tutorial walks through the complete pipeline execution, from configuration to results interpretation.

## Quick Start (5 minutes)

For users who want to test the pipeline immediately:

```bash
# 1. Activate environment
conda activate snakemake

# 2. Quick test with example data
snakemake --cores 2 --use-conda \
  results/defensefinder/GCF_000005825.2/GCF_000005825.2_defense_finder_systems.tsv

# 3. View results
head results/defensefinder/GCF_000005825.2/GCF_000005825.2_defense_finder_systems.tsv
```

## Complete Workflow Tutorial

### Step 1: Project Setup

```bash
# Navigate to your project directory
cd acinetobacter-defence-pipeline

# Verify installation
snakemake --dry-run
```

### Step 2: Configure Your Analysis

#### Edit Configuration File

```bash
# Open configuration
nano config/config.yaml
```

#### Example Configuration for Different Use Cases

**Small Test Dataset (3-5 genomes):**
```yaml
samples:
  - "GCF_000005825.2"    # A. baumannii ATCC 17978 (reference)
  - "GCF_000746645.1"    # A. baumannii clinical isolate
  - "GCF_001022155.1"    # A. pittii example

resfinder:
  threshold: 90.0
  coverage: 60.0

blast:
  evalue: 0.005
  identity: 80
```

**Large Dataset (50+ genomes):**
```yaml
samples:
  - "GCF_000005825.2"
  - "GCF_000746645.1"
  # ... add up to 100+ accessions

# Optimized parameters for large datasets
resfinder:
  threshold: 85.0     # Slightly relaxed for broader detection
  coverage: 50.0

blast:
  evalue: 0.001       # More stringent for large datasets
  identity: 85
```

### Step 3: Execute Pipeline Components

#### Run Individual Analysis Steps

```bash
# 1. Download genomes only
snakemake --cores 4 --use-conda \
  resources/genomes/GCF_000005825.2.fna

# 2. Run DefenseFinder analysis
snakemake --cores 4 --use-conda \
  results/defensefinder/GCF_000005825.2/GCF_000005825.2_defense_finder_systems.tsv

# 3. Run PADLOC analysis
snakemake --cores 4 --use-conda \
  results/padloc/GCF_000005825.2/GCF_000005825.2_padloc.csv

# 4. Run resistance gene analysis
snakemake --cores 4 --use-conda \
  results/resfinder/GCF_000005825.2/ResFinder_results_tab.txt
```

#### Run Complete Pipeline

```bash
# Execute all analysis steps
snakemake --cores 8 --use-conda

# For cluster environments
snakemake --cluster "sbatch --time=4:00:00 --mem=8G" --cores 20
```

### Step 4: Monitor Progress

#### Real-time Monitoring

```bash
# Monitor with detailed output
snakemake --cores 8 --use-conda --printshellcmds

# Monitor resource usage
htop    # In separate terminal

# Check Snakemake status
snakemake --summary
```

#### Progress Tracking

```bash
# Check completed rules
snakemake --list-target-rules
snakemake --list-rules

# View workflow graph
snakemake --dag | dot -Tpng > workflow_graph.png
```

## Understanding Pipeline Outputs

### DefenseFinder Results

```bash
# View system-level results
head results/defensefinder/GCF_000005825.2/GCF_000005825.2_defense_finder_systems.tsv
```

**Expected columns:**
- `sys_id`: System identifier
- `type`: Defense system type (e.g., RM, CRISPR-Cas)
- `subtype`: Specific subtype classification
- `genes_count`: Number of genes in system
- `gene_names`: List of component genes

### PADLOC Results

```bash
# View PADLOC predictions
head results/padloc/GCF_000005825.2/GCF_000005825.2_padloc.csv
```

**Key columns:**
- `protein.name`: Protein identifier
- `target.name`: Defense system target
- `system`: Defense system classification
- `confidence`: Prediction confidence score

### Consolidated Results

```bash
# View merged results across all genomes
head results/consolidated/defense_finder_consolidated.tsv
```

## Analysis Examples

### Example 1: Single Genome Analysis

```bash
# Complete analysis for A. baumannii reference genome
snakemake --cores 4 --use-conda \
  results/defensefinder/GCF_000005825.2/GCF_000005825.2_defense_finder_systems.tsv \
  results/padloc/GCF_000005825.2/GCF_000005825.2_padloc.csv \
  results/resfinder/GCF_000005825.2/ResFinder_results_tab.txt

# View defense systems found
echo "Defense systems in A. baumannii ATCC 17978:"
cut -f2 results/defensefinder/GCF_000005825.2/GCF_000005825.2_defense_finder_systems.tsv | tail -n +2 | sort | uniq -c
```

### Example 2: Comparative Analysis

```bash
# Run analysis on multiple genomes for comparison
snakemake --cores 8 --use-conda \
  results/consolidated/defense_finder_consolidated.tsv \
  results/consolidated/padloc_consolidated.tsv

# Generate comparison statistics
Rscript workflow/scripts/compare_tools.R
```

### Example 3: Statistical Analysis

```bash
# Run complete pipeline with statistical analysis
snakemake --cores 8 --use-conda results/analysis/defense_figures

# View generated plots
ls results/analysis/defense_figures/
```

## Advanced Usage

### Custom Parameter Optimization

#### For Sensitivity Analysis
```yaml
# High sensitivity configuration
resfinder:
  threshold: 80.0      # Lower threshold
  coverage: 50.0       # Lower coverage

blast:
  evalue: 0.01         # More permissive
  identity: 70         # Lower identity
```

#### For Specificity Analysis
```yaml
# High specificity configuration
resfinder:
  threshold: 95.0      # Higher threshold
  coverage: 80.0       # Higher coverage

blast:
  evalue: 0.001        # More stringent
  identity: 90         # Higher identity
```

### Cluster Execution Examples

#### SLURM Integration

```bash
# Create cluster profile
mkdir -p ~/.config/snakemake/slurm
cat > ~/.config/snakemake/slurm/config.yaml << 'EOF'
cluster: "sbatch --job-name={rule} --time={resources.time} --mem={resources.mem_mb} --cpus-per-task={threads}"
default-resources: [mem_mb=4000, time="2:00:00"]
max-jobs-per-second: 10
max-status-checks-per-second: 1
EOF

# Execute with cluster profile
snakemake --profile slurm --cores 50
```

#### PBS/Torque Integration

```bash
# PBS submission
snakemake --cluster "qsub -l walltime=4:00:00,mem=8gb" --cores 20
```

### Custom Analysis Scripts

#### Adding New Analysis Rules

```python
# Example: Add custom correlation analysis
rule custom_correlation:
    input:
        defense="results/consolidated/defense_finder_consolidated.tsv",
        resistance="results/consolidated/resfinder_consolidated.tsv"
    output:
        correlation="results/analysis/custom_correlation.png"
    conda:
        "workflow/envs/r_analysis.yaml"
    script:
        "workflow/scripts/custom_analysis.R"
```

## Workflow Customization

### Modifying Analysis Parameters

#### Tool-Specific Customization

```bash
# DefenseFinder: Adjust model sensitivity
# Edit rule in Snakefile:
shell: """
    defense-finder run --models-dir resources/defensefinder_models \
                      --coverage-profile HIGH_COV \
                      --out-dir {params.outdir} \
                      {input.genome}
"""

# PADLOC: Enable additional features
shell: """
    padloc --faa {input.faa} \
           --gff {input.gff} \
           --cpu {threads} \
           --hmm-profile extended
"""
```

### Adding New Genome Sources

```bash
# Modify download rule for different databases
rule download_genome_custom:
    output: "resources/genomes/{accession}.fna"
    shell: """
        # Custom download logic for private databases
        your_custom_download_script.sh {wildcards.accession} {output}
    """
```

## Quality Control

### Validation Steps

```bash
# 1. Check genome integrity
for genome in resources/genomes/*.fna; do
    echo "Checking $genome"
    grep -c "^>" "$genome"  # Count sequences
done

# 2. Validate tool outputs
python workflow/scripts/validate_outputs.py

# 3. Check result completeness
Rscript workflow/scripts/quality_check.R
```

### Performance Monitoring

```bash
# Generate resource usage report
snakemake --report results/pipeline_report.html

# Benchmark analysis
snakemake --benchmark-extended results/benchmarks/
```

## Results Interpretation

### Defense System Classification

Understanding DefenseFinder output types:
- **RM**: Restriction-Modification systems
- **CRISPR-Cas**: Adaptive immunity systems  
- **Abortive infection**: Programmed cell death systems
- **Toxin-antitoxin**: Stress response systems

### Statistical Output

```r
# Load and examine consolidated results
library(tidyverse)

defense_data <- read_tsv("results/consolidated/defense_finder_consolidated.tsv")
resistance_data <- read_tsv("results/consolidated/resfinder_consolidated.tsv")

# Basic statistics
summary(defense_data)
table(defense_data$type)
```

## Best Practices

### Resource Management
- **Start small**: Test with 3-5 genomes before scaling
- **Monitor memory**: Large datasets require substantial RAM
- **Use checkpoints**: Snakemake automatically handles interrupted runs

### Data Organization
- **Keep raw data separate**: Never modify original genome files
- **Version results**: Tag important analysis runs
- **Document parameters**: Always record configuration used

### Reproducibility
- **Use exact versions**: Pin conda package versions
- **Save environments**: Export conda environments for sharing
- **Document changes**: Use git to track configuration modifications

## Troubleshooting Common Issues

### Network-Related Problems
```bash
# Set up retry logic for downloads
# Add to download rule:
max_attempts=3
attempt=1
while [ $attempt -le $max_attempts ]; do
    command && break
    sleep 5
    attempt=$((attempt + 1))
done
```

### Memory Issues
```bash
# Reduce parallel jobs for memory-constrained systems
snakemake --cores 2 --resources mem_mb=8000
```

### Tool-Specific Errors
```bash
# DefenseFinder: Clear model cache
rm -rf ~/.macsyfinder/

# PADLOC: Reset environment
conda env remove -n padloc
mamba create -n padloc padloc
```

---

**Next**: See [Configuration Reference](CONFIG.md) for detailed parameter descriptions.