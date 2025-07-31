#!/bin/bash
# Simple benchmark test for one genome
# Author: Vigneshwaran Muthuraman

echo "=== Single Genome Pipeline Benchmark ==="
echo "Testing individual genome processing time"
echo ""

# Choose one genome to test (change this to one of your actual genome IDs)
TEST_GENOME="CP006963.1"  # Replace with one of your genome accessions

echo "Test genome: $TEST_GENOME"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "Cores available: $(nproc)"
echo ""

# Create a minimal config for just this one genome
cat > config/single_genome_test.yaml << EOF
# Single genome test configuration
samples:
  - "$TEST_GENOME"

# Tool configurations (copy from your main config)
resfinder:
  threshold: "0.9"
  coverage: "0.6"

blast:
  evalue: "1e-6"
  identity: "80"

hmrg:
  evalue: "0.005"
  qcov_hsp_perc: "80"

crisprcasfinder:
  installation_path: "resources/CRISPRCasFinder"
EOF

# Clean any previous results for this genome to ensure fresh timing
echo "Cleaning previous results for $TEST_GENOME..."
rm -rf results/defensefinder/$TEST_GENOME/
rm -rf results/padloc/$TEST_GENOME/
rm -rf results/resfinder/$TEST_GENOME/
rm -rf results/ime_blast/$TEST_GENOME/
rm -rf results/hmrg_blast/$TEST_GENOME/

echo "Starting pipeline benchmark..."
echo "Pipeline start: $(date '+%Y-%m-%d %H:%M:%S')"

# Record start time
START_TIME=$(date +%s)

# Run the pipeline for individual genome tools (excluding CRISPRCasFinder)
time snakemake \
  --configfile config/single_genome_test.yaml \
  --cores 8 \
  --use-conda \
  results/defensefinder/$TEST_GENOME/${TEST_GENOME}_defense_finder_systems.tsv \
  results/padloc/$TEST_GENOME/${TEST_GENOME}_padloc.csv \
  results/resfinder/$TEST_GENOME/ResFinder_results_tab.txt \
  results/ime_blast/$TEST_GENOME/${TEST_GENOME}_ime.blastn \
  results/hmrg_blast/$TEST_GENOME/${TEST_GENOME}_hmrg.tblastn \
  "results/crisprcasfinder_combined/result.json"

# Record end time
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "Pipeline end: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "=== BENCHMARK RESULTS ==="
echo "Genome: $TEST_GENOME"
echo "Total time: $TOTAL_TIME seconds"
echo "Total time: $(($TOTAL_TIME / 60)) minutes $(($TOTAL_TIME % 60)) seconds"

# Calculate genome size for efficiency metrics
if [ -f "resources/genomes/$TEST_GENOME.fna" ]; then
    GENOME_SIZE=$(wc -c < resources/genomes/$TEST_GENOME.fna)
    GENOME_SIZE_MB=$(echo "scale=2; $GENOME_SIZE / 1024 / 1024" | bc -l)
    echo "Genome size: $GENOME_SIZE_MB MB"
    echo "Processing rate: $(echo "scale=2; $GENOME_SIZE_MB / ($TOTAL_TIME / 60)" | bc -l) MB/minute"
fi

# Estimate time for full dataset
echo ""
echo "=== FULL DATASET ESTIMATES ==="
echo "You have 132 genomes in your full dataset"
echo ""
echo "Sequential processing (one at a time):"
SEQUENTIAL_TIME=$((TOTAL_TIME * 132))
echo "  Total time: $((SEQUENTIAL_TIME / 3600)) hours $((SEQUENTIAL_TIME % 3600 / 60)) minutes"

echo ""
echo "Parallel processing estimates:"
echo "  4 cores: $((SEQUENTIAL_TIME / 4 / 3600)) hours $((SEQUENTIAL_TIME / 4 % 3600 / 60)) minutes"
echo "  8 cores: $((SEQUENTIAL_TIME / 8 / 3600)) hours $((SEQUENTIAL_TIME / 8 % 3600 / 60)) minutes"
echo "  16 cores: $((SEQUENTIAL_TIME / 16 / 3600)) hours $((SEQUENTIAL_TIME / 16 % 3600 / 60)) minutes"

echo ""
echo "NOTE: CRISPRCasFinder processes all genomes together (~1.5 hours total)"
echo "So add ~1.5 hours to the above estimates for the complete pipeline"

echo ""
echo "=== RESULTS VERIFICATION ==="
echo "Check that all tools completed successfully:"

# Verify each tool's output
if [ -f "results/defensefinder/$TEST_GENOME/${TEST_GENOME}_defense_finder_systems.tsv" ]; then
    DEFENSE_COUNT=$(wc -l < results/defensefinder/$TEST_GENOME/${TEST_GENOME}_defense_finder_systems.tsv)
    echo "✓ DefenseFinder: $DEFENSE_COUNT lines (including header)"
else
    echo "✗ DefenseFinder: FAILED"
fi

if [ -f "results/padloc/$TEST_GENOME/${TEST_GENOME}_padloc.csv" ]; then
    PADLOC_COUNT=$(wc -l < results/padloc/$TEST_GENOME/${TEST_GENOME}_padloc.csv)
    echo "✓ PADLOC: $PADLOC_COUNT lines (including header)"
else
    echo "✗ PADLOC: FAILED"
fi

if [ -f "results/resfinder/$TEST_GENOME/ResFinder_results_tab.txt" ]; then
    RESFINDER_COUNT=$(wc -l < results/resfinder/$TEST_GENOME/ResFinder_results_tab.txt)
    echo "✓ ResFinder: $RESFINDER_COUNT lines (including header)"
else
    echo "✗ ResFinder: FAILED"
fi

if [ -f "results/ime_blast/$TEST_GENOME/${TEST_GENOME}_ime.blastn" ]; then
    IME_COUNT=$(wc -l < results/ime_blast/$TEST_GENOME/${TEST_GENOME}_ime.blastn)
    echo "✓ IME BLAST: $IME_COUNT hits"
else
    echo "✗ IME BLAST: FAILED"
fi

if [ -f "results/hmrg_blast/$TEST_GENOME/${TEST_GENOME}_hmrg.tblastn" ]; then
    HMRG_COUNT=$(grep -v '^#' results/hmrg_blast/$TEST_GENOME/${TEST_GENOME}_hmrg.tblastn | wc -l)
    echo "✓ HMRG BLAST: $HMRG_COUNT hits"
else
    echo "✗ HMRG BLAST: FAILED"
fi

echo ""
echo "=== SUMMARY ==="
echo "Single genome processing time: $(($TOTAL_TIME / 60)) minutes"
echo "Estimated full dataset time (parallel): $((SEQUENTIAL_TIME / 8 / 60)) - $((SEQUENTIAL_TIME / 4 / 60)) minutes"
echo "Plus CRISPRCasFinder: ~90 minutes"
echo ""
echo "Benchmark complete! 🎉"
