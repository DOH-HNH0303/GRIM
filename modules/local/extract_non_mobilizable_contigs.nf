process EXTRACT_NON_MOBILIZABLE_CONTIGS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(assembly_fasta), path(contig_list)

    output:
    tuple val(meta), path("${prefix}_non_mobilizable.fasta"), emit: non_mobilizable_fasta, optional: true
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > extract_contigs.py <<'EOF'
from Bio import SeqIO
import sys

# Read contig IDs to extract
contig_ids = set()
try:
    with open('${contig_list}', 'r') as f:
        for line in f:
            contig_id = line.strip()
            if contig_id:
                contig_ids.add(contig_id)
except Exception as e:
    print(f"Warning: Could not read contig list: {e}", file=sys.stderr)

# Extract sequences
extracted_count = 0
if contig_ids:
    with open('${prefix}_non_mobilizable.fasta', 'w') as out:
        for record in SeqIO.parse('${assembly_fasta}', 'fasta'):
            if record.id in contig_ids:
                SeqIO.write(record, out, 'fasta')
                extracted_count += 1
else:
    # Create empty file if no contigs to extract
    open('${prefix}_non_mobilizable.fasta', 'w').close()

print(f"Extracted {extracted_count} non-mobilizable contigs", file=sys.stderr)
EOF

    python3 extract_contigs.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_non_mobilizable.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "1.78")
    END_VERSIONS
    """
}
