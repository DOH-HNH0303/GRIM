process COUNT_GFA_EDGES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*_edge_counts.tsv"), emit: edge_counts
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 $projectDir/bin/count_gfa_edges.py \\
        ${gfa} \\
        ${prefix}_edge_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "contig\\tedge_count" > ${prefix}_edge_counts.tsv
    echo -e "contig_1\\t4" >> ${prefix}_edge_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
