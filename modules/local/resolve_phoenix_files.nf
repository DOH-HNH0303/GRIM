process RESOLVE_PHOENIX_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::biopython=1.81 bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0':
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"


    input:
    tuple val(meta), path(phoenix_run_dir), path(ont_genome)

    output:
    tuple val(meta), path("gamma_ar_file"), path("amrfinder_report"), path("phoenix_assembly"), path(ont_genome), emit: resolved_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_id = meta.id
    """
    # Check if phoenix run directory exists and is accessible
    if [ ! -d "${phoenix_run_dir}" ]; then
        echo "ERROR: Phoenix run directory does not exist or is not accessible: ${phoenix_run_dir}"
        echo "Please verify the phoenix_outdir path in your samplesheet points to the correct Phoenix run directory."
        exit 1
    fi

    # Check if phoenix_summary.tsv exists in the run directory
    if [ ! -f "${phoenix_run_dir}/phoenix_summary.tsv" ]; then
        echo "ERROR: Could not find phoenix_summary.tsv in run directory ${phoenix_run_dir}/"
        echo "Expected structure: ${phoenix_run_dir}/phoenix_summary.tsv"
        echo "Available files in run directory:"
        ls -la "${phoenix_run_dir}/"
        exit 1
    fi

    # Check if sample directory exists in the run directory
    if [ ! -d "${phoenix_run_dir}/${sample_id}" ]; then
        echo "ERROR: Could not find sample directory ${phoenix_run_dir}/${sample_id}/"
        echo "Available directories in ${phoenix_run_dir}:"
        ls -la "${phoenix_run_dir}/"
        exit 1
    fi

    # Find GAMMA file (search in sample subdirectories, 2 levels deep)
    gamma_file=\$(find ${phoenix_run_dir}/${sample_id}/*/ -name "${sample_id}_ResGANNCBI_*_srst2.gamma" 2>/dev/null | head -1)
    if [ -z "\$gamma_file" ]; then
        echo "ERROR: Could not find GAMMA file for sample ${sample_id}"
        echo "Searched in: ${phoenix_run_dir}/${sample_id}/*/"
        echo "Expected pattern: ${sample_id}_ResGANNCBI_*_srst2.gamma"
        echo "Available GAMMA files in ${phoenix_run_dir}/${sample_id}/:"
        find "${phoenix_run_dir}/${sample_id}/" -type f -name "*.gamma" 2>/dev/null || echo "No .gamma files found"
        echo "Sample directory structure:"
        ls -la "${phoenix_run_dir}/${sample_id}/"
        echo "Subdirectories in sample directory:"
        find "${phoenix_run_dir}/${sample_id}/" -type d -mindepth 1 -maxdepth 1 2>/dev/null | while read subdir; do
            echo "Contents of \$subdir:"
            ls -la "\$subdir/" 2>/dev/null || echo "  Cannot access \$subdir"
        done
        exit 1
    fi
    ln -s "\$gamma_file" gamma_ar_file

    # Find AMRFinder report (search in sample subdirectories, 1 level deep)
    amrfinder_file=\$(find ${phoenix_run_dir}/${sample_id}/*/ -name "${sample_id}_all_genes.tsv" 2>/dev/null | head -1)
    if [ -z "\$amrfinder_file" ]; then
        echo "ERROR: Could not find AMRFinder report for sample ${sample_id}"
        echo "Searched in: ${phoenix_run_dir}/${sample_id}/*/"
        echo "Expected pattern: ${sample_id}_all_genes.tsv"
        echo "Available TSV files in ${phoenix_run_dir}/${sample_id}/:"
        find "${phoenix_run_dir}/${sample_id}/" -type f -name "*.tsv" 2>/dev/null || echo "No .tsv files found"
        echo "Sample directory structure:"
        ls -la "${phoenix_run_dir}/${sample_id}/"
        echo "Subdirectories in sample directory:"
        find "${phoenix_run_dir}/${sample_id}/" -type d -mindepth 1 -maxdepth 1 2>/dev/null | while read subdir; do
            echo "Contents of \$subdir:"
            ls -la "\$subdir/" 2>/dev/null || echo "  Cannot access \$subdir"
        done
        exit 1
    fi
    ln -s "\$amrfinder_file" amrfinder_report

    # Find Phoenix assembly (search in sample subdirectories)
    assembly_file=\$(find ${phoenix_run_dir}/${sample_id}/*/ -name "${sample_id}.scaffolds.fa.gz" 2>/dev/null | head -1)
    if [ -z "\$assembly_file" ]; then
        echo "ERROR: Could not find Phoenix assembly for sample ${sample_id}"
        echo "Searched in: ${phoenix_run_dir}/${sample_id}/*/"
        echo "Expected patterns: ${sample_id}.scaffolds.fa.gz"
        echo "Available FASTA files in ${phoenix_run_dir}/${sample_id}/:"
        find "${phoenix_run_dir}/${sample_id}/" -type f -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fa.gz" 2>/dev/null || echo "No FASTA files found"
        echo "Sample directory structure:"
        ls -la "${phoenix_run_dir}/${sample_id}/"
        echo "Subdirectories in sample directory:"
        find "${phoenix_run_dir}/${sample_id}/" -type d -mindepth 1 -maxdepth 1 2>/dev/null | while read subdir; do
            echo "Contents of \$subdir:"
            ls -la "\$subdir/" 2>/dev/null || echo "  Cannot access \$subdir"
        done
        exit 1
    fi
    ln -s "\$assembly_file" phoenix_assembly

    echo "Successfully resolved files for sample ${sample_id}:"
    echo "  Phoenix summary: ${phoenix_run_dir}/phoenix_summary.tsv"
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