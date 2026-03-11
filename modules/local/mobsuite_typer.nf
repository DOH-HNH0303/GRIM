process MOBSUITE_TYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_0':
        'biocontainers/mob_suite:3.1.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*mobtyper_results.tsv"), emit: mobtyper_results, optional: false
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mob_typer \\
        --infile $fasta_name \\
        $args \\
        --num_threads $task.cpus \\
        --out_file ${prefix}_mobtyper_results.tsv \\
        --sample_id ${prefix} \\
        --multi 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_typer--version 2>&1) | sed 's/^.*mob_typer/; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results

    touch results/chromosome.fasta
    touch results/contig_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_typer --version 2>&1) | sed 's/^.*mob_typer //; s/ .*\$//')
    END_VERSIONS
    """
}