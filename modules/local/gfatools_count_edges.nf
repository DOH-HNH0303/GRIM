process GFATOOLS_COUNT_EDGES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfatools:0.5--h9a82719_3' :
        'biocontainers/gfatools:0.5--h9a82719_3' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("${prefix}_edge_counts.tsv"), emit: edge_counts
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Count occurrences of each edge_N in the GFA file
    # Extract all edge names and count their occurrences
    
    echo -e "contig\\tedge_count" > ${prefix}_edge_counts.tsv
    
    gfatools view ${gfa} | \\
        grep -oE 'edge_[0-9]+' | \\
        sort | \\
        uniq -c | \\
        awk '{print \$2"\\t"\$1}' | \\
        sort -k1,1V >> ${prefix}_edge_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfatools: \$(gfatools version 2>&1 | sed '1!d; s/.*: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_edge_counts.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/GNU Awk //' | sed 's/,.*//')
    END_VERSIONS
    """
}
