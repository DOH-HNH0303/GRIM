process RESOLVE_PHOENIX_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::biopython=1.81 bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0':
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"


    input:
    tuple val(meta), path(phoenix_run_dir), path(ont_genome), path(gfa_file)

    output:
    tuple val(meta), path("gamma_ar_file"), path("amrfinder_report"), path("phoenix_assembly"), path(ont_genome), path(gfa_file), emit: resolved_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_id = meta.id
    """
    echo "DEBUG: Working directory contents:"
    ls -la .

    echo "DEBUG: Phoenix run directory variable: ${phoenix_run_dir}"
    echo "DEBUG: Sample ID: ${sample_id}"

    # Determine if we're working with a traditional directory structure or S3 staged files
    # When S3 directories are staged, Nextflow flattens the directory structure
    
    # First, check if we have the traditional nested structure
    if [ -d "${phoenix_run_dir}" ]; then
        echo "Found traditional Phoenix directory structure"
        search_dir="${phoenix_run_dir}"
        summary_file="${phoenix_run_dir}/Phoenix_Summary.tsv"
        sample_dir="${phoenix_run_dir}/${sample_id}"
    else
        echo "Phoenix directory not found as nested structure, checking for flattened S3 staging"
        # S3 staging flattens directories - files are staged at the current working directory
        search_dir="."
        
        # Look for Phoenix summary file - try both lowercase and capitalized versions
        if [ -f "Phoenix_Summary.tsv" ]; then
            summary_file="Phoenix_Summary.tsv"
        elif [ -f "Phoenix_Summary.tsv" ]; then
            summary_file="Phoenix_Summary.tsv"
        else
            echo "ERROR: Could not find phoenix_summary.tsv or Phoenix_Summary.tsv"
            echo "Available files in working directory:"
            ls -la .
            exit 1
        fi
        
        # In flattened structure, sample directory is at root level
        sample_dir="${sample_id}"
    fi

    echo "Using search directory: \$search_dir"
    echo "Using summary file: \$summary_file"
    echo "Using sample directory: \$sample_dir"

    # Verify Phoenix summary file exists
    if [ ! -f "\$summary_file" ]; then
        echo "ERROR: Phoenix summary file not found: \$summary_file"
        echo "Available files in search directory:"
        ls -la "\$search_dir"
        exit 1
    fi

    # Check if sample directory exists
    if [ ! -d "\$sample_dir" ]; then
        echo "ERROR: Could not find sample directory: \$sample_dir"
        echo "Available directories in search area:"
        ls -la "\$search_dir"
        exit 1
    fi

    # Find GAMMA file (search in sample subdirectories, 2 levels deep)
    gamma_file=\$(find "\$sample_dir"/*/ -name "${sample_id}_ResGANNCBI_*_srst2.gamma" 2>/dev/null | head -1)
    if [ -z "\$gamma_file" ]; then
        echo "ERROR: Could not find GAMMA file for sample ${sample_id}"
        echo "Searched in: \$sample_dir/*/"
        echo "Expected pattern: ${sample_id}_ResGANNCBI_*_srst2.gamma"
        echo "Available GAMMA files in \$sample_dir/:"
        find "\$sample_dir/" -type f -name "*.gamma" 2>/dev/null || echo "No .gamma files found"
        echo "Sample directory structure:"
        ls -la "\$sample_dir/"
        echo "Subdirectories in sample directory:"
        find "\$sample_dir/" -type d -mindepth 1 -maxdepth 1 2>/dev/null | while read subdir; do
            echo "Contents of \$subdir:"
            ls -la "\$subdir/" 2>/dev/null || echo "  Cannot access \$subdir"
        done
        exit 1
    fi
    ln -s "\$gamma_file" gamma_ar_file

    # Find AMRFinder report (search in sample subdirectories, 1 level deep)
    amrfinder_file=\$(find "\$sample_dir"/*/ -name "${sample_id}_all_genes.tsv" 2>/dev/null | head -1)
    if [ -z "\$amrfinder_file" ]; then
        echo "ERROR: Could not find AMRFinder report for sample ${sample_id}"
        echo "Searched in: \$sample_dir/*/"
        echo "Expected pattern: ${sample_id}_all_genes.tsv"
        echo "Available TSV files in \$sample_dir/:"
        find "\$sample_dir/" -type f -name "*.tsv" 2>/dev/null || echo "No .tsv files found"
        echo "Sample directory structure:"
        ls -la "\$sample_dir/"
        echo "Subdirectories in sample directory:"
        find "\$sample_dir/" -type d -mindepth 1 -maxdepth 1 2>/dev/null | while read subdir; do
            echo "Contents of \$subdir:"
            ls -la "\$subdir/" 2>/dev/null || echo "  Cannot access \$subdir"
        done
        exit 1
    fi
    ln -s "\$amrfinder_file" amrfinder_report

    # Find Phoenix assembly (search in sample subdirectories)
    assembly_file=\$(find "\$sample_dir"/*/ -name "${sample_id}.scaffolds.fa.gz" 2>/dev/null | head -1)
    if [ -z "\$assembly_file" ]; then
        echo "ERROR: Could not find Phoenix assembly for sample ${sample_id}"
        echo "Searched in: \$sample_dir/*/"
        echo "Expected patterns: ${sample_id}.scaffolds.fa.gz"
        echo "Available FASTA files in \$sample_dir/:"
        find "\$sample_dir/" -type f -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fa.gz" 2>/dev/null || echo "No FASTA files found"
        echo "Sample directory structure:"
        ls -la "\$sample_dir/"
        echo "Subdirectories in sample directory:"
        find "\$sample_dir/" -type d -mindepth 1 -maxdepth 1 2>/dev/null | while read subdir; do
            echo "Contents of \$subdir:"
            ls -la "\$subdir/" 2>/dev/null || echo "  Cannot access \$subdir"
        done
        exit 1
    fi
    ln -s "\$assembly_file" phoenix_assembly

    echo "Successfully resolved files for sample ${sample_id}:"
    echo "  Phoenix summary: \$summary_file"
    echo "  GAMMA file: \$gamma_file"
    echo "  AMRFinder report: \$amrfinder_file" 
    echo "  Phoenix assembly: \$assembly_file"
    echo "  ONT genome: ${ont_genome}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    """
    touch gamma_ar_file
    touch amrfinder_report  
    touch phoenix_assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //g' | sed 's/ .*//g')
    END_VERSIONS
    """
}