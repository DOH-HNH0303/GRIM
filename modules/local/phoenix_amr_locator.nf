process PHOENIX_AMR_LOCATOR {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1a6f9b2c2e9b8c7d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7/data':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(gamma_ar_file), path(amrfinder_report), path(assembly_fasta)

    output:
    tuple val(meta), path("${prefix}_gene_locations.tsv"), emit: locations
    tuple val(meta), path("${prefix}_detailed_amr.tsv"), emit: detailed_amr
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    phoenix_amr_locator.py \\
        --sample_id ${meta.id} \\
        --gamma_ar ${gamma_ar_file} \\
        --amrfinder_report ${amrfinder_report} \\
        --assembly ${assembly_fasta} \\
        --output_locations ${prefix}_gene_locations.tsv \\
        --output_detailed ${prefix}_detailed_amr.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_gene_locations.tsv
    touch ${prefix}_detailed_amr.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}