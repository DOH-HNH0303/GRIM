process RESOLVE_PHOENIX_FILES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(phoenix_outdir), path(ont_genome)

    output:
    tuple val(meta), path("gamma_ar_file"), path("amrfinder_report"), path("phoenix_assembly"), path(ont_genome), emit: resolved_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_id = meta.id
    """
    # Find GAMMA file (with wildcard pattern)
    gamma_file=\$(find ${phoenix_outdir}/${sample_id}/ -name "${sample_id}_ResGANNCBI_*_srst2.gamma" | head -1)
    if [ -z "\$gamma_file" ]; then
        echo "ERROR: Could not find GAMMA file for sample ${sample_id} in ${phoenix_outdir}/${sample_id}/"
        exit 1
    fi
    ln -s "\$gamma_file" gamma_ar_file

    # Find AMRFinder report
    amrfinder_file="${phoenix_outdir}/${sample_id}/${sample_id}_amrfinder_report.tsv"
    if [ ! -f "\$amrfinder_file" ]; then
        echo "ERROR: Could not find AMRFinder report for sample ${sample_id} at \$amrfinder_file"
        exit 1
    fi
    ln -s "\$amrfinder_file" amrfinder_report

    # Find Phoenix assembly
    assembly_file="${phoenix_outdir}/${sample_id}/${sample_id}_scaffolds.fasta"
    if [ ! -f "\$assembly_file" ]; then
        echo "ERROR: Could not find Phoenix assembly for sample ${sample_id} at \$assembly_file"
        exit 1
    fi
    ln -s "\$assembly_file" phoenix_assembly

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