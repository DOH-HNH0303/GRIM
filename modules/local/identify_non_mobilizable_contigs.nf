process IDENTIFY_NON_MOBILIZABLE_CONTIGS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(mobtyper_tsv), path(platon_tsv)

    output:
    tuple val(meta), path("${prefix}_non_mobilizable_contigs.txt"), emit: non_mobilizable_list
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > identify_contigs.py <<'EOF'
import csv
import sys

# Read MOB-typer results
mobtyper_contigs = set()
try:
    with open('${mobtyper_tsv}', 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            if row.get('num_contigs'):  # MOB-typer found mobile elements
                # Get contig IDs from MOB-typer output
                # This assumes contig names are in the first column
                contig_id = row.get('sample_id', '')
                if contig_id:
                    mobtyper_contigs.add(contig_id)
except Exception as e:
    print(f"Warning: Could not parse MOB-typer TSV: {e}", file=sys.stderr)

# Read Platon results
platon_contigs = set()
try:
    with open('${platon_tsv}', 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            # Platon marks plasmids in the 'ID' column
            contig_id = row.get('ID', '')
            if contig_id and contig_id != 'ID':  # Skip header if repeated
                platon_contigs.add(contig_id)
except Exception as e:
    print(f"Warning: Could not parse Platon TSV: {e}", file=sys.stderr)

# Find contigs identified by Platon but not by MOB-typer
non_mobilizable = platon_contigs - mobtyper_contigs

# Write output
with open('${prefix}_non_mobilizable_contigs.txt', 'w') as out:
    for contig in sorted(non_mobilizable):
        out.write(f"{contig}\\n")

print(f"Found {len(platon_contigs)} contigs in Platon results", file=sys.stderr)
print(f"Found {len(mobtyper_contigs)} contigs in MOB-typer results", file=sys.stderr)
print(f"Identified {len(non_mobilizable)} non-mobilizable contigs", file=sys.stderr)
EOF

    python3 identify_contigs.py

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
