process RESOLVE_PHOENIX_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::biopython=1.81 bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0':
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"


    input:
    tuple val(meta), path(phoenix_summary), path(gamma_file), path(amrfinder_file), path(assembly_file), path(ont_genome), path(gfa_file, stageAs: 'input.gfa')

    output:
    tuple val(meta), path("gamma_ar_file"), path("amrfinder_report"), path("phoenix_assembly"), path("ont_genome.fasta"), path("assembly.gfa"), emit: resolved_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_id = meta.id
    """
    echo "Processing sample: ${sample_id}"
    echo "Files staged by Nextflow:"
    ls -lah

    # Verify all required input files exist
    if [ ! -f "${gamma_file}" ]; then
        echo "ERROR: GAMMA file not found: ${gamma_file}"
        exit 1
    fi

    if [ ! -f "${amrfinder_file}" ]; then
        echo "ERROR: AMRFinder file not found: ${amrfinder_file}"
        exit 1
    fi

    if [ ! -f "${assembly_file}" ]; then
        echo "ERROR: Assembly file not found: ${assembly_file}"
        exit 1
    fi

    if [ ! -f "${ont_genome}" ]; then
        echo "ERROR: ONT genome file not found: ${ont_genome}"
        exit 1
    fi

    # Create symlinks with standardized names
    ln -s "${gamma_file}" gamma_ar_file
    ln -s "${amrfinder_file}" amrfinder_report
    ln -s "${assembly_file}" phoenix_assembly
    ln -s "${ont_genome}" ont_genome.fasta

    # Handle GFA file - check if input.gfa exists and has content
    if [ -f "input.gfa" ] && [ -s "input.gfa" ]; then
        echo "GFA file provided with content"
        ln -s "input.gfa" assembly.gfa
    elif [ -f "input.gfa" ]; then
        echo "GFA file staged but empty - creating placeholder"
        touch assembly.gfa
    else
        echo "No GFA file provided, creating empty placeholder"
        touch assembly.gfa
    fi

    echo "Successfully resolved files for sample ${sample_id}:"
    echo "  GAMMA file: ${gamma_file}"
    echo "  AMRFinder report: ${amrfinder_file}"
    echo "  Phoenix assembly: ${assembly_file}"
    echo "  ONT genome: ${ont_genome}"
    if [ -f "assembly.gfa" ] && [ -s "assembly.gfa" ]; then
        GFA_SIZE=\$(stat -f%z input.gfa 2>/dev/null || stat -c%s input.gfa)
        echo "  GFA file: input.gfa (\${GFA_SIZE} bytes)"
    else
        echo "  GFA file: none (placeholder created)"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    """
    touch gamma_ar_file
    touch amrfinder_report  
    touch phoenix_assembly
    touch ont_genome.fasta
    touch assembly.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //g' | sed 's/ .*//g')
    END_VERSIONS
    """
}