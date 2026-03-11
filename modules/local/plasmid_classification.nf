process PLASMID_CLASSIFICATION {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::mob_suite=3.1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_0':
        'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ont_genome)

    output:
    tuple val(meta), path("${prefix}_contig_classification.tsv"), emit: classification
    tuple val(meta), path("${prefix}_plasmid_replicons.tsv"), emit: replicons
    path "mob_output", emit: mob_results, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plasmid_classification.py \\
        --sample_id ${meta.id} \\
        --ont_genome ${ont_genome} \\
        --output_classification ${prefix}_contig_classification.tsv \\
        --output_replicons ${prefix}_plasmid_replicons.tsv \\
        --threads ${task.cpus} \\
        --outdir mob_output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        mob_suite: \$(mob_recon --version 2>&1 | grep -oP 'MOB-suite \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_contig_classification.tsv
    touch ${prefix}_plasmid_replicons.tsv
    mkdir -p mob_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        mob_suite: \$(mob_recon --version 2>&1 | grep -oP 'MOB-suite \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
