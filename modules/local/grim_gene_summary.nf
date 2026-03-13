process GRIM_GENE_SUMMARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'public.ecr.aws/o8h2f0o1/illumina_amr_locator:1.0.0' }"

    input:
    tuple val(meta), path(gene_mappings), path(mobtyper_results)

    output:
    tuple val(meta), path("${prefix}_gene_summary.csv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    grim_gene_summary.py \\
        --sample_id ${meta.id} \\
        --gene_mappings ${gene_mappings} \\
        --mobtyper_results ${mobtyper_results} \\
        --output_summary ${prefix}_gene_summary.csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_gene_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
