process BANDAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bandage:0.8.1--hc9558a2_2':
        'quay.io/biocontainers/bandage:0.8.1--hc9558a2_2' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.png")             , emit: png      , optional: true
    tuple val(meta), path("*.svg")             , emit: svg      , optional: true
    tuple val(meta), path("*_bandage_info.tsv"), emit: info_tsv
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args_image = task.ext.args_image ?: ''
    def args_info  = task.ext.args_info  ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    # Generate assembly graph visualization
    Bandage image \\
        ${gfa} \\
        ${prefix}.png \\
        ${args_image}

    # Generate SVG version
    Bandage image \\
        ${gfa} \\
        ${prefix}.svg \\
        ${args_image}

    # Extract assembly graph information as TSV
    Bandage info \\
        ${gfa} \\
        ${args_info} \\
        --tsv > ${prefix}_bandage_info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bandage: \$(echo \$(Bandage --version 2>&1) | sed 's/^Version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}.svg
    touch ${prefix}_bandage_info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bandage: \$(echo \$(Bandage --version 2>&1) | sed 's/^Version: //')
    END_VERSIONS
    """
}
