process PHOENIX_AMR_LOCATOR {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::biopython=1.81 bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0':
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(gamma_ar_file), path(amrfinder_report), path(phoenix_assembly_fasta), path(ont_complete_genome)

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
        --phoenix_assembly ${phoenix_assembly_fasta} \\
        --ont_genome ${ont_complete_genome} \\
        --output_locations ${prefix}_gene_locations.tsv \\
        --output_detailed ${prefix}_detailed_amr.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        blast: \$(blastn -version | head -n1 | sed 's/blastn: //')
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
        blast: \$(blastn -version | head -n1 | sed 's/blastn: //')
    END_VERSIONS
    """
}