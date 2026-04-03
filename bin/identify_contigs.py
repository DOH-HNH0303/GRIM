#!/usr/bin/env python3
"""
Identify non-mobilizable contigs by comparing MOB-typer and Platon results.

This script identifies contigs that are marked as plasmids by Platon but not
identified as mobile by MOB-typer, suggesting they are non-mobilizable plasmids.
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Identify non-mobilizable contigs from MOB-typer, Platon, and Bandage results"
    )
    parser.add_argument(
        "--mobtyper",
        type=Path,
        required=True,
        help="Path to MOB-typer TSV output file"
    )
    parser.add_argument(
        "--platon",
        type=Path,
        required=True,
        help="Path to Platon TSV output file"
    )
    parser.add_argument(
        "--bandage",
        type=Path,
        required=False,
        help="Path to Bandage info TSV output file (optional)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to output file for non-mobilizable contigs"
    )
    return parser.parse_args()


def read_mobtyper_contigs(mobtyper_file):
    """
    Read MOB-typer results and extract contig IDs.
    
    Args:
        mobtyper_file: Path to MOB-typer TSV file
        
    Returns:
        Set of contig IDs identified by MOB-typer
    """
    mobtyper_contigs = set()
    try:
        with open(mobtyper_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('num_contigs'):  # MOB-typer found mobile elements
                    # Get contig IDs from MOB-typer output
                    # This assumes contig names are in the first column
                    contig_id = row.get('sample_id', '')
                    if contig_id:
                        mobtyper_contigs.add(contig_id)
    except Exception as e:
        print(f"Warning: Could not parse MOB-typer TSV: {e}", file=sys.stderr)
    
    return mobtyper_contigs


def read_platon_contigs(platon_file):
    """
    Read Platon results and extract contig IDs.
    
    Args:
        platon_file: Path to Platon TSV file
        
    Returns:
        Set of contig IDs identified by Platon as plasmids
    """
    platon_contigs = set()
    try:
        with open(platon_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Platon marks plasmids in the 'ID' column
                contig_id = row.get('ID', '')
                if contig_id and contig_id != 'ID':  # Skip header if repeated
                    platon_contigs.add(contig_id)
    except Exception as e:
        print(f"Warning: Could not parse Platon TSV: {e}", file=sys.stderr)
    
    return platon_contigs


def read_bandage_info(bandage_file):
    """
    Read Bandage info results and extract contig information.
    
    Args:
        bandage_file: Path to Bandage info TSV file
        
    Returns:
        Dictionary mapping contig IDs to their properties (length, depth, etc.)
    """
    bandage_info = {}
    if not bandage_file:
        return bandage_info
    
    try:
        with open(bandage_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Bandage info TSV typically has columns like:
                # Node, Length, Depth, etc.
                contig_id = row.get('Node', row.get('Name', ''))
                if contig_id:
                    bandage_info[contig_id] = {
                        'length': row.get('Length', ''),
                        'depth': row.get('Depth', ''),
                        'coverage': row.get('Coverage', '')
                    }
    except Exception as e:
        print(f"Warning: Could not parse Bandage TSV: {e}", file=sys.stderr)
    
    return bandage_info


def write_non_mobilizable_contigs(contigs, output_file):
    """
    Write non-mobilizable contigs to output file.
    
    Args:
        contigs: Set of contig IDs
        output_file: Path to output file
    """
    with open(output_file, 'w') as out:
        for contig in sorted(contigs):
            out.write(f"{contig}\n")


def main():
    """Main function."""
    args = parse_args()
    
    # Read input files
    mobtyper_contigs = read_mobtyper_contigs(args.mobtyper)
    platon_contigs = read_platon_contigs(args.platon)
    bandage_info = read_bandage_info(args.bandage) if args.bandage else {}
    
    # Find contigs identified by Platon but not by MOB-typer
    non_mobilizable = platon_contigs - mobtyper_contigs
    
    # Write output
    write_non_mobilizable_contigs(non_mobilizable, args.output)
    
    # Print summary statistics
    print(f"Found {len(platon_contigs)} contigs in Platon results", file=sys.stderr)
    print(f"Found {len(mobtyper_contigs)} contigs in MOB-typer results", file=sys.stderr)
    print(f"Identified {len(non_mobilizable)} non-mobilizable contigs", file=sys.stderr)
    
    if bandage_info:
        print(f"Bandage info available for {len(bandage_info)} contigs", file=sys.stderr)
        # Count how many non-mobilizable contigs have Bandage info
        contigs_with_info = sum(1 for c in non_mobilizable if c in bandage_info)
        print(f"Bandage info available for {contigs_with_info}/{len(non_mobilizable)} non-mobilizable contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
