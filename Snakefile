import os
# Main Snakefile
configfile: "config/config.yaml" #Config file with all the required parameters

# Define the final rule that requests all outputs
rule all:
    input:
    	#Genome Downloads
        expand("resources/genomes/{accession}.fna", accession=config["samples"]),
        #Tool outputs
        expand("results/defensefinder/{accession}/{accession}_defense_finder_systems.tsv", accession=config["samples"]),
        expand("results/resfinder/{accession}/ResFinder_results_tab.txt", accession=config["samples"]),
        expand("results/ime_blast/{accession}/{accession}_ime.blastn", accession=config["samples"]),
        expand("results/padloc/{accession}/{accession}_padloc.csv", accession=config["samples"]),
        expand("results/hmrg_blast/{accession}/{accession}_hmrg.tblastn", accession=config["samples"]),
        
        # CRISPRCasFinder - Combined approach
        "results/crisprcasfinder_combined/result.json",  # Combined analysis result
        
        #Consolidated results
        "results/consolidated/defense_finder_consolidated.tsv",
        "results/consolidated/padloc_consolidated.tsv",
        "results/consolidated/resfinder_consolidated.tsv",
        "results/consolidated/ime_consolidated.tsv",
        "results/consolidated/hmrg_consolidated.tsv",
 
        # Analysis results
        "results/analysis/defense_figures"

# Rule to download a genome from NCBI
rule download_genome:
    output:
        "resources/genomes/{accession}.fna"
    params:
        accession = "{accession}"
    shell:
        """
        # Create a directory for the genome if it doesn't exist
        mkdir -p resources/genomes
        
        # Add a random delay between 1-3 seconds to avoid rate limiting
        sleep $((1 + RANDOM % 3))
        
        # Download the genome using NCBI's efetch utility with retry logic
        max_attempts=3
        attempt=1
        
        while [ $attempt -le $max_attempts ]; do
            echo "Attempt $attempt of $max_attempts for {params.accession}"
            
            ~/edirect/efetch -db nucleotide -id {params.accession} -format fasta > {output}
            
            # Check if file has content
            if [ -s {output} ]; then
                echo "Successfully downloaded {params.accession}"
                break
            else
                echo "Download failed, waiting before retry..."
                sleep 5
                attempt=$((attempt + 1))
            fi
        done
        
        # If all attempts failed, exit with error
        if [ ! -s {output} ]; then
            echo "Failed to download {params.accession} after $max_attempts attempts"
            exit 1
        fi
        """
rule setup_defensefinder_models:
    output:
        models_dir = directory("resources/defensefinder_models"),
        casfinder_dir = directory("resources/defensefinder_models/CasFinder"),
        defense_models_dir = directory("resources/defensefinder_models/defense-finder-models"),
        flag = touch("resources/defensefinder_models_ready.flag")
    conda:
        "workflow/envs/defense.yaml"
    shell:
        """
        echo "Setting up DefenseFinder models in project directory..."
        echo "This copies both defense-finder-models AND CasFinder directories"
        
        # Step 1: Update DefenseFinder models in default location
        echo "Downloading/updating DefenseFinder models..."
        defense-finder update
        
        # Step 2: Create project models directory structure
        echo "Creating project models directory structure..."
        mkdir -p {output.models_dir}
        
        # Step 3: Copy defense-finder-models
        echo "Copying defense-finder-models..."
        if [ -d ~/.macsyfinder/models/defense-finder-models ]; then
            cp -r ~/.macsyfinder/models/defense-finder-models {output.models_dir}/
            echo "defense-finder-models copied successfully"
            echo "   Files: $(find {output.models_dir}/defense-finder-models -type f | wc -l)"
        else
            echo "Warning: defense-finder-models not found in ~/.macsyfinder/models/"
            mkdir -p {output.defense_models_dir}
        fi
        
        # Step 4: Copy CasFinder models
        echo "Copying CasFinder models..."
        if [ -d ~/.macsyfinder/models/CasFinder ]; then
            cp -r ~/.macsyfinder/models/CasFinder {output.models_dir}/
            echo "CasFinder models copied successfully"
            echo "   Files: $(find {output.models_dir}/CasFinder -type f | wc -l)"
        else
            echo "Warning: CasFinder not found in ~/.macsyfinder/models/"
            mkdir -p {output.casfinder_dir}
        fi
        
        # Step 5: Verify complete setup
        echo ""
        echo "=== DefenseFinder Models Setup Complete ==="
        echo "Project models directory: {output.models_dir}"
        echo "Directory structure:"
        ls -la {output.models_dir}
        
        echo ""
        echo "Total files copied:"
        echo "  defense-finder-models: $(find {output.models_dir}/defense-finder-models -type f 2>/dev/null | wc -l) files"
        echo "  CasFinder: $(find {output.models_dir}/CasFinder -type f 2>/dev/null | wc -l) files"
        echo "  Total model files (.hmm): $(find {output.models_dir} -name '*.hmm' 2>/dev/null | wc -l)"
        
        echo ""
        echo "DefenseFinder models setup complete!"
        echo "   Models ready for use with: --models-dir {output.models_dir}"
        """
        
rule run_defensefinder_bioconda:
    input:
        genome = "resources/genomes/{accession}.fna",
        models_ready = "resources/defensefinder_models_ready.flag"  # Depends on models setup
    output:
        systems = "results/defensefinder/{accession}/{accession}_defense_finder_systems.tsv",
        genes = "results/defensefinder/{accession}/{accession}_defense_finder_genes.tsv"
    params:
        outdir = "results/defensefinder/{accession}"
    conda:
        "workflow/envs/defense.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        echo "Running DefenseFinder for {wildcards.accession}..."
        echo "Using project models directory: resources/defensefinder_models"
        echo "Input genome: {input.genome}"
        echo "Output directory: {params.outdir}"
        
        # Run DefenseFinder with project models directory (both defense-finder-models and CasFinder available)
        defense-finder run \
            --models-dir resources/defensefinder_models \
            --out-dir {params.outdir} \
            {input.genome}
        
        echo "DefenseFinder analysis complete for {wildcards.accession}"
        """

rule run_padloc_wrapper_style:
    input:
        genome = "resources/genomes/{accession}.fna"
    output:
        csv = "results/padloc/{accession}/{accession}_padloc.csv",
        gff = "results/padloc/{accession}/{accession}_padloc.gff",
        prodigal_faa = "results/padloc/{accession}/{accession}_prodigal.faa",
        prodigal_gff = "results/padloc/{accession}/{accession}_prodigal.gff",
        fixed_gff = "results/padloc/{accession}/{accession}_fixed.gff"
    params:
        accession = "{accession}",
        output_dir = "results/padloc/{accession}"
    threads: 4
    shell:
        """
        # Activate your working environment
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate padloc
        
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Step 1: Run Prodigal to generate protein sequences and GFF
        echo "[$(date +%T)] Running Prodigal for {params.accession}..."
        prodigal -i {input.genome} \
                 -a {output.prodigal_faa} \
                 -f gff \
                 -o {output.prodigal_gff}
        
        if [ $? -ne 0 ]; then
            echo "Error: Prodigal failed. Exiting."
            exit 1
        fi
        
        # Step 2: Extract protein IDs from the FASTA file
        echo "[$(date +%T)] Extracting protein IDs from FASTA..."
        grep ">" {output.prodigal_faa} | sed 's/>//' > {params.output_dir}/protein_ids.txt
        
        # Step 3: Create a properly formatted GFF file
        echo "[$(date +%T)] Creating properly formatted GFF file..."
        awk -v OFS="\\t" '
        BEGIN {{
            # Read protein IDs from file
            while ((getline id < "{params.output_dir}/protein_ids.txt") > 0) {{
                protein_ids[++count] = id;
            }}
            close("{params.output_dir}/protein_ids.txt");
            current_id = 1;
        }}
        /^#/ {{
            print $0;
            next;
        }}
        $3 == "CDS" {{
            # Add ID attribute that exactly matches the FASTA header
            $9 = "ID=" protein_ids[current_id] ";" $9;
            current_id++;
            print $0;
            next;
        }}
        {{
            print $0;  # Print other lines unchanged
        }}' {output.prodigal_gff} > {output.fixed_gff}
        
        # Clean up
        rm {params.output_dir}/protein_ids.txt
        
        # Step 4: Run PADLOC with the fixed files
        echo "[$(date +%T)] Running PADLOC analysis..."
        
        # Change to the output directory
        cd {params.output_dir}
        
        # Run PADLOC without --outdir flag
        padloc --faa {params.accession}_prodigal.faa \
               --gff {params.accession}_fixed.gff \
               --cpu {threads}
        
        if [ $? -eq 0 ]; then
            echo "[$(date +%T)] PADLOC analysis completed successfully!"
            echo "Results saved to: {output.csv}"
            
            # Check if the expected output file was created
            if [ -f "{params.accession}_prodigal_padloc.csv" ]; then
                # Move to the expected output name
                mv {params.accession}_prodigal_padloc.csv {params.accession}_padloc.csv
                echo "Renamed output file to {params.accession}_padloc.csv"
            fi
            
            if [ -f "{params.accession}_prodigal_padloc.gff" ]; then
                # Move to the expected output name  
                mv {params.accession}_prodigal_padloc.gff {params.accession}_padloc.gff
                echo "Renamed GFF output file to {params.accession}_padloc.gff"
            fi
        else
            echo "[$(date +%T)] Error: PADLOC analysis failed."
            exit 1
        fi
        """
        

# Rule to prepare ResFinder database
rule prepare_resfinder_db:
    output:
        db_dir = directory("resources/resfinder_db"),
        flag = touch("resources/resfinder_db_ready.flag")
    conda:
        "workflow/envs/resfinder.yaml"
    shell:
        """
        # Clone the ResFinder database
        mkdir -p resources
        git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git {output.db_dir}
        
        # Run the installation script for the database
        cd {output.db_dir}
        python3 INSTALL.py
        """
        
# Rule to run ResFinder
rule run_resfinder:
    input:
        genome = "resources/genomes/{accession}.fna",
        db_ready = "resources/resfinder_db_ready.flag"
    output:
        results = "results/resfinder/{accession}/ResFinder_results_tab.txt"
    params:
        outdir = os.path.abspath("results/resfinder/{accession}"),
        infile = os.path.abspath("resources/genomes/{accession}.fna"),
        db_path = os.path.abspath("resources/resfinder_db"),
        species = "acinetobacter",
        threshold = config["resfinder"]["threshold"],
        min_cov = config["resfinder"]["coverage"]
    conda:
        "workflow/envs/resfinder.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run ResFinder with explicit database path
        run_resfinder.py \
            -ifa {params.infile} \
            -o {params.outdir} \
            -s "{params.species}" \
            -l {params.min_cov} \
            -t {params.threshold} \
            --acquired \
            -db_res {params.db_path}
        """

#Rule to setup CasFinder
rule setup_crisprcasfinder:
    output:
        flag = touch("resources/crisprcasfinder_ready.flag"),
        repo_dir = directory("resources/CRISPRCasFinder")
    conda:
        "workflow/envs/crisprcasfinder_env.yaml"
    params:
        repo_url = "https://github.com/dcouvin/CRISPRCasFinder.git"
    shell:
        """
        # Create resources directory if it doesn't exist
        mkdir -p resources
        
        # Check if CRISPRCasFinder already exists
        if [ -d "{output.repo_dir}" ]; then
            echo "CRISPRCasFinder already exists at {output.repo_dir}"
            
            # Verify it's a valid installation
            if [ -f "{output.repo_dir}/CRISPRCasFinder.pl" ]; then
                echo "Existing installation verified"
            else
                echo "Invalid installation found. Re-cloning..."
                rm -rf {output.repo_dir}
                git clone {params.repo_url} {output.repo_dir}
            fi
        else
            echo "Cloning CRISPRCasFinder repository..."
            git clone {params.repo_url} {output.repo_dir}
        fi
        
        # Verify the installation
        if [ ! -f "{output.repo_dir}/CRISPRCasFinder.pl" ]; then
            echo "Error: CRISPRCasFinder.pl not found after setup"
            exit 1
        fi
        
        echo "CRISPRCasFinder setup complete at {output.repo_dir}"
        """

# Rule to run CrisprCasFinder

rule run_crisprcasfinder_combined:
    input:
        genomes = expand("resources/genomes/{accession}.fna", accession=config["samples"]),
        setup_ready = "resources/crisprcasfinder_ready.flag",
        repo_dir = "resources/CRISPRCasFinder"
    output:
        combined_results = "results/crisprcasfinder_combined/result.json",
        combined_genomes = "results/crisprcasfinder_combined/combined_genomes.fna"
    params:
        output_dir = "results/crisprcasfinder_combined",
        crispr_dir = config.get("crisprcasfinder", {}).get("installation_path", "resources/CRISPRCasFinder")
    conda:
        "workflow/envs/crisprcasfinder_env.yaml"
    threads: 8
    shell:
        """
        # Get absolute paths before changing directories
        OUTPUT_ABS=$(realpath {params.output_dir})
        COMBINED_GENOMES_ABS=$(realpath {output.combined_genomes})
        RESULT_JSON_ABS=$(realpath {output.combined_results})
        
        # Prepare output directory
        mkdir -p {params.output_dir}
        
        # Create combined genome file
        cat {input.genomes} > {output.combined_genomes}
        
        # Set up working directory in CRISPRCasFinder
        mkdir -p {params.crispr_dir}/temp_analysis
        cp {output.combined_genomes} {params.crispr_dir}/temp_analysis/input.fna
        
        # Execute CRISPRCasFinder analysis
        cd {params.crispr_dir}
        perl CRISPRCasFinder.pl -in temp_analysis/input.fna -cas -keep
        
        # Process results
        RESULTS_DIR=$(find . -maxdepth 1 -type d -name "Result_input_*" | head -1)
        
        if [ -n "$RESULTS_DIR" ] && [ -d "$RESULTS_DIR" ]; then
            # Copy results to output directory
            cp -r "$RESULTS_DIR"/* "$OUTPUT_ABS/" 2>/dev/null || true
            
            # Handle result.json
            if [ -f "$RESULTS_DIR/result.json" ]; then
                cp "$RESULTS_DIR/result.json" "$RESULT_JSON_ABS"
            else
                # Generate summary for completed analysis with results
                echo '{{
                    "status": "completed",
                    "timestamp": "'$(date -Iseconds)'",
                    "input_genomes": '$(echo {input.genomes} | wc -w)',
                    "analysis_method": "CRISPRCasFinder_combined",
                    "results_found": true
                }}' > "$RESULT_JSON_ABS"
            fi
            
            # Clean up temporary files
            rm -rf temp_analysis "$RESULTS_DIR"
        else
            # Handle case where no CRISPRs/Cas genes were found (normal biological result)
            echo "No CRISPR/Cas systems detected - generating summary for completed analysis"
            echo '{{
                "status": "completed",
                "timestamp": "'$(date -Iseconds)'",
                "input_genomes": '$(echo {input.genomes} | wc -w)',
                "analysis_method": "CRISPRCasFinder_combined",
                "results_found": false,
                "note": "No CRISPR arrays or Cas genes detected in input genomes"
            }}' > "$RESULT_JSON_ABS"
            
            # Clean up temporary files
            rm -rf temp_analysis
        fi
        """

            
#IME detection

# Rule to prepare IME protein database
rule prepare_ime_db:
    output:
        db_dir = directory("resources/ime_db"),
        fasta = "resources/ime_db/ime_proteins.fasta",
        flag = touch("resources/ime_db_ready.flag")
    conda:
        "workflow/envs/blast.yaml"
    shell:
        """
        # Create IME database directory
        mkdir -p {output.db_dir}
        
        # Download protein sequences from ICEberg
        # This is the URL for ICEberg protein sequences
        curl -o {output.db_dir}/ICE_proteins_raw.faa https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/IME_aa_all.fas
        
        # Process the file - replace spaces with pipes in headers
        sed '/^>/ s/ /|/g' {output.db_dir}/ICE_proteins_raw.faa > {output.fasta}
        
        # Clean up the raw file
        rm {output.db_dir}/ICE_proteins_raw.faa
        """

# Rule to run tblastn for IMEs
rule run_ime_blast:
    input:
        genome = "resources/genomes/{accession}.fna",
        protein_db = "resources/ime_db/ime_proteins.fasta"
    output:
        blast_results = "results/ime_blast/{accession}/{accession}_ime.blastn"
    params:
        evalue = config["blast"]["evalue"],
        identity = config["blast"]["identity"],
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qtitle"
    conda:
        "workflow/envs/blast.yaml"
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.blast_results})
        
        # Run tblastn (protein query against translated nucleotide database)
        tblastn -query {input.protein_db} \
                -subject {input.genome} \
                -evalue {params.evalue} \
                -outfmt "{params.outfmt}" \
                -out {output.blast_results}
        """


rule prepare_hmrg_db:
    output:
        processed_fasta = "resources/hmrg_db/hmrg_proteins.fasta",
        flag = touch("resources/hmrg_db_ready.flag")
    conda:
        "workflow/envs/blast.yaml"
    shell:
        """
        # Create HMRG database directory if it doesn't exist
        mkdir -p resources/hmrg_db
        
        # Check if BacMet database file exists
        BACMET_FILE="resources/hmrg_db/BacMet2_EXP_database.fasta"
        
        echo "Looking for BacMet file at: $BACMET_FILE"
        echo "Current working directory: $(pwd)"
        
        if [ -f "$BACMET_FILE" ]; then
            echo "Found existing BacMet file: $BACMET_FILE"
            echo "File size: $(du -h $BACMET_FILE | cut -f1)"
            echo "Number of sequences: $(grep -c '^>' $BACMET_FILE)"
        else
            echo "BacMet database file not found. Attempting download..."
            
            # Try to download automatically (with insecure flag due to certificate issue)
            if command -v wget >/dev/null 2>&1; then
                wget --no-check-certificate -O "$BACMET_FILE" 'https://bacmet.biomedicine.gu.se/downloads/BacMet2_EXP_database.fasta' || {{
                    echo "Download failed. Creating placeholder..."
                    echo ">placeholder_protein" > "$BACMET_FILE"
                    echo "MPLACEHOLDERPROTEINSEQUENCE" >> "$BACMET_FILE"
                }}
            else
                echo "wget not available. Creating placeholder..."
                echo ">placeholder_protein" > "$BACMET_FILE"
                echo "MPLACEHOLDERPROTEINSEQUENCE" >> "$BACMET_FILE"
            fi
        fi
        
        # Process the BacMet file - handle CRLF line endings and clean headers
        echo "Processing BacMet database headers..."
        
        # First convert CRLF to LF (handle Windows line endings), then clean headers
        tr -d '\\r' < "$BACMET_FILE" | sed '/^>/ s/ /|/g' > {output.processed_fasta}
        
        echo "HMRG database preparation complete"
        echo "Original file: $BACMET_FILE"
        echo "Processed file: {output.processed_fasta}"
        echo "Number of sequences in processed file: $(grep -c '^>' {output.processed_fasta})"
        
        # Verify the processed file has content
        if [ $(grep -c '^>' {output.processed_fasta}) -eq 0 ]; then
            echo "Warning: Processed file has no sequences!"
            echo "First few lines of original file:"
            head -10 "$BACMET_FILE"
            echo "First few lines of processed file:"
            head -10 {output.processed_fasta}
        fi
        """

# Rule to run HMRG tblastn on individual genomes (like IME approach)
rule run_hmrg_tblastn:
    input:
        genome = "resources/genomes/{accession}.fna",
        protein_db = "resources/hmrg_db/hmrg_proteins.fasta",
        db_ready = "resources/hmrg_db_ready.flag"
    output:
        blast_results = "results/hmrg_blast/{accession}/{accession}_hmrg.tblastn"
    params:
        evalue = config.get("hmrg", {}).get("evalue", "0.005"),
        qcov_hsp_perc = config.get("hmrg", {}).get("qcov_hsp_perc", "80"),
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs qcovhsp"
    conda:
        "workflow/envs/blast.yaml"
    threads: 4
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.blast_results})
        
        # Check if HMRG database has real sequences (not placeholder)
        if grep -q "placeholder_protein" {input.protein_db}; then
            echo "Warning: Using placeholder HMRG database for {wildcards.accession}"
            echo "Please download the real BacMet2_EXP_database.fasta file."
            echo "Creating empty results file..."
            echo "# No HMRG analysis - placeholder database used" > {output.blast_results}
            echo "# Please download BacMet2_EXP_database.fasta from https://bacmet.biomedicine.gu.se/download.html" >> {output.blast_results}
        else
            # Run tblastn with stringent parameters (exactly like your script)
            echo "Running HMRG tblastn for {wildcards.accession}..."
            
            tblastn -query {input.protein_db} \
                   -subject {input.genome} \
                   -evalue {params.evalue} \
                   -qcov_hsp_perc {params.qcov_hsp_perc} \
                   -outfmt "{params.outfmt}" \
                   -out {output.blast_results}
            
            echo "HMRG tblastn completed for {wildcards.accession}"
            echo "Results saved to: {output.blast_results}"
            echo "Number of hits: $(wc -l < {output.blast_results})"
        fi
        """

# Data consolidation rules

# Rule to consolidate DefenseFinder results
rule consolidate_defensefinder:
    input:
        system_files = expand("results/defensefinder/{accession}/{accession}_defense_finder_systems.tsv", 
                              accession=config["samples"])
    output:
        consolidated = "results/consolidated/defense_finder_consolidated.tsv"
    run:
        import pandas as pd
        import os
        
        # Initialize an empty list to store dataframes
        dfs = []
        
        # Process each input file
        for file_path in input.system_files:
            try:
                # Extract accession from file path
                accession = os.path.basename(file_path).split('_defense_finder')[0]
                
                # Read the file
                df = pd.read_csv(file_path, sep='\t')
                
                # Add genome identifier column
                df['Genome_ID'] = accession
                
                # Append to list
                dfs.append(df)
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        # Combine all dataframes
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            
            # Ensure the Genome_ID column is first
            cols = ['Genome_ID'] + [col for col in combined_df.columns if col != 'Genome_ID']
            combined_df = combined_df[cols]
            
            # Write to output file
            combined_df.to_csv(output.consolidated, sep='\t', index=False)
        else:
            # Create empty file with header
            with open(output.consolidated, 'w') as f:
                f.write("Genome_ID\tsys_id\ttype\tsubtype\tgenes_count\tgene_names\tgene_ids\n")

# PADLOC consolidation rule
rule consolidate_padloc:
    input:
        padloc_files = expand("results/padloc/{accession}/{accession}_padloc.csv", 
                             accession=config["samples"])
    output:
        consolidated = "results/consolidated/padloc_consolidated.tsv"
    run:
        import pandas as pd
        import os
        
        # Initialize an empty list to store dataframes
        dfs = []
        
        # Process each input file
        for file_path in input.padloc_files:
            try:
                # Extract accession from file path
                accession = file_path.split('/')[-2]  # Get directory name which is the accession
                
                # Read the file
                df = pd.read_csv(file_path)
                
                # Add genome identifier column
                df['Genome_ID'] = accession
                
                # Append to list
                dfs.append(df)
                print(f"Successfully processed {accession} with {len(df)} rows")
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        # Combine all dataframes
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            
            # Ensure the Genome_ID column is first
            cols = ['Genome_ID'] + [col for col in combined_df.columns if col != 'Genome_ID']
            combined_df = combined_df[cols]
            
            # Write to output file
            combined_df.to_csv(output.consolidated, sep='\t', index=False)
            print(f"Consolidated PADLOC results: {len(combined_df)} total rows from {len(dfs)} genomes")
        else:
            # Create empty file with basic header
            with open(output.consolidated, 'w') as f:
                f.write("Genome_ID\n")
            print("No PADLOC results found - created empty file")

# Rule to consolidate ResFinder results
rule consolidate_resfinder:
    input:
        resfinder_files = expand("results/resfinder/{accession}/ResFinder_results_tab.txt", 
                                accession=config["samples"])
    output:
        consolidated = "results/consolidated/resfinder_consolidated.tsv"
    run:
        import pandas as pd
        import os
        
        # Initialize an empty list to store dataframes
        dfs = []
        
        # Process each input file
        for file_path in input.resfinder_files:
            try:
                # Extract accession from file path
                accession = file_path.split('/')[2]  # Get the accession from path
                
                # Read the file (ResFinder output has a specific format)
                df = pd.read_csv(file_path, sep='\t')
                
                # Add genome identifier column
                df['Genome_ID'] = accession
                
                # Append to list
                dfs.append(df)
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        # Combine all dataframes
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            
            # Ensure the Genome_ID column is first
            cols = ['Genome_ID'] + [col for col in combined_df.columns if col != 'Genome_ID']
            combined_df = combined_df[cols]
            
            # Write to output file
            combined_df.to_csv(output.consolidated, sep='\t', index=False)
        else:
            # Create empty file with appropriate header
            with open(output.consolidated, 'w') as f:
                f.write("Genome_ID\tResistance gene\tIdentity\tAlignment Length\tCoverage\tPosition in reference\tContig\tPosition in contig\tPhenotype\tAccession no.\n")

# Rule to consolidate IME BLAST results
rule consolidate_ime_blast:
    input:
        ime_files = expand("results/ime_blast/{accession}/{accession}_ime.blastn", accession=config["samples"])
        
    output:
        consolidated = "results/consolidated/ime_consolidated.tsv"
    run:
        import pandas as pd
        import os
        
        # Initialize an empty list to store dataframes
        dfs = []
        
        # Process each input file
        for file_path in input.ime_files:
            try:
                # Read the file
                df = pd.read_csv(file_path, sep='\t')
                
                # Append to list
                dfs.append(df)
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        # Combine all dataframes
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            
            # Write to output file
            combined_df.to_csv(output.consolidated, sep='\t', index=False)
        else:
            # Create empty file with header
            with open(output.consolidated, 'w') as f:
                f.write("genome_id\time_id\time_name\time_source\tpercent_identity\talignment_length\tevalue\n")

# Rule to consolidate HMRG results (simple approach like IME)
rule consolidate_hmrg:
    input:
        hmrg_files = expand("results/hmrg_blast/{accession}/{accession}_hmrg.tblastn", 
                           accession=config["samples"])
    output:
        consolidated = "results/consolidated/hmrg_consolidated.tsv"
    run:
        import pandas as pd
        import os
        
        # Initialize an empty list to store dataframes
        dfs = []
        
        # Process each input file
        for file_path in input.hmrg_files:
            try:
                # Extract accession from file path
                accession = file_path.split('/')[2]  # results/hmrg_blast/ACCESSION/...
                
                # Read the file, skip comment lines
                df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
                
                # Only process if file has data
                if not df.empty:
                    # Add column names (matching the outfmt)
                    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                                 'qlen', 'qcovs', 'qcovhsp']
                    
                    # Add genome identifier
                    df['Genome_ID'] = accession
                    
                    # Append to list
                    dfs.append(df)
                    print(f"Successfully processed {accession} with {len(df)} HMRG hits")
                else:
                    print(f"No HMRG hits found for {accession}")
                    
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        # Combine all dataframes
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            
            # Ensure the Genome_ID column is first
            cols = ['Genome_ID'] + [col for col in combined_df.columns if col != 'Genome_ID']
            combined_df = combined_df[cols]
            
            # Write to output file
            combined_df.to_csv(output.consolidated, sep='\t', index=False)
            print(f"Consolidated HMRG results: {len(combined_df)} total hits from {len(dfs)} genomes with hits")
        else:
            # Create empty file with header (matching the BLAST output format)
            with open(output.consolidated, 'w') as f:
                f.write("Genome_ID\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tqcovs\tqcovhsp\n")
            print("No HMRG results found - created empty file with headers")


rule defense_distribution_analysis:
    input:
        defense = "results/consolidated/defense_finder_consolidated.tsv"
    output:
        figures_dir = directory("results/analysis/defense_figures")
    conda:
        "workflow/envs/r_analysis.yaml"
    script:
        "workflow/scripts/defense_distribution.R"
        


