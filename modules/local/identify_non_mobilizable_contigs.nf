process IDENTIFY_NON_MOBILIZABLE_CONTIGS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(mobtyper_tsv), path(platon_tsv), path(edge_counts, stageAs: 'edge_counts.tsv')

    output:
    tuple val(meta), path("${prefix}_non_mobilizable_contigs.txt"), emit: non_mobilizable_list
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Only add --bandage argument if edge_counts file exists and has content
    // Note: Check is done in bash since file is staged as 'edge_counts.tsv'
    """
    # Build command with optional --bandage argument
    CMD="python $projectDir/bin/identify_contigs.py \\
        --mobtyper ${mobtyper_tsv} \\
        --platon ${platon_tsv}"
    
    # Add --bandage argument if edge_counts.tsv exists and has content
    if [ -f "edge_counts.tsv" ] && [ -s "edge_counts.tsv" ]; then
        CMD="\$CMD --bandage edge_counts.tsv"
    fi
    
    CMD="\$CMD --output ${prefix}_non_mobilizable_contigs.txt"
    
    # Execute the command
    eval \$CMD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_non_mobilizable_contigs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
