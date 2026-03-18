process PLATON {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/platon:1.6--pyhdfd78af_0':
        'biocontainers/platon:1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}.plasmid.fasta")    , emit: plasmids, optional: true
    tuple val(meta), path("${prefix}.chromosome.fasta") , emit: chromosomes, optional: true
    tuple val(meta), path("${prefix}.json")             , emit: json, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def db_arg = params.platon_db ? "--db ${params.platon_db}" : ''
    """
    # Check if FASTA is empty or has no sequences
    if [ ! -s ${fasta} ]; then
        echo "Warning: Empty FASTA file, creating empty Platon results" >&2
        touch ${prefix}.tsv
        echo "ID\\t# Chromosome\\tLength\\tCoverage\\t# ORFs\\tRDS\\tCircular\\tInc Type(s)\\t# Replication\\t# Mobilization\\t# ORFs on Conjugation\\t# Conjugation\\tAnnotation File" > ${prefix}.tsv
    else
        # Decompress if needed
        if [ "$is_compressed" == "true" ]; then
            gzip -c -d ${fasta} > ${fasta_name}
        fi

        # Determine which file to use
        FASTA_INPUT="${is_compressed ? fasta_name : fasta}"

        # Platon will use database from params.platon_db or download automatically
        platon \\
            --output . \\
            ${db_arg} \\
            --prefix ${prefix} \\
            --threads ${task.cpus} \\
            ${args} \\
            \${FASTA_INPUT} || {
                echo "Platon failed, creating empty results" >&2
                touch ${prefix}.tsv
                echo "ID\\t# Chromosome\\tLength\\tCoverage\\t# ORFs\\tRDS\\tCircular\\tInc Type(s)\\t# Replication\\t# Mobilization\\t# ORFs on Conjugation\\t# Conjugation\\tAnnotation File" > ${prefix}.tsv
            }
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platon: \$(platon --version 2>&1 | sed 's/platon //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.plasmid.fasta
    touch ${prefix}.chromosome.fasta
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platon: \$(platon --version 2>&1 | sed 's/platon //')
    END_VERSIONS
    """
}