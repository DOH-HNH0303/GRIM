#!/usr/bin/env python3
"""
Count edges for each contig in a GFA assembly graph.
Outputs a TSV file with contig name and edge count.
"""
import sys
from collections import Counter
from pathlib import Path


def count_gfa_edges(gfa_file, output_file):
    """
    Parse GFA file and count edges for each contig.
    
    Args:
        gfa_file: Path to input GFA file
        output_file: Path to output TSV file
    """
    edges = Counter()
    contigs = set()

    with open(gfa_file) as f:
        for line in f:
            if line.startswith("S"):
                # Segment line: S <name> <sequence> [tags]
                fields = line.rstrip().split("\t")
                contigs.add(fields[1])

            elif line.startswith("L"):
                # Link line: L <from> <from_orient> <to> <to_orient> <overlap>
                fields = line.rstrip().split("\t")
                c1 = fields[1]
                c2 = fields[3]
                edges[c1] += 1
                edges[c2] += 1

    # Write output TSV with header
    with open(output_file, 'w') as out:
        out.write("contig\tedge_count\n")
        # Print all contigs, including those with zero edges
        for contig in sorted(contigs):
            out.write(f"{contig}\t{edges[contig]}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.gfa> <output.tsv>", file=sys.stderr)
        sys.exit(1)
    
    gfa_file = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    
    if not gfa_file.exists():
        print(f"Error: GFA file not found: {gfa_file}", file=sys.stderr)
        sys.exit(1)
    
    count_gfa_edges(gfa_file, output_file)
    print(f"Edge counts written to {output_file}")
