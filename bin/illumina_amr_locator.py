#!/usr/bin/env python3
"""
Phoenix AMR Locator (Refactored)
Maps AMR genes from Phoenix Illumina assembly to ONT complete genome using BLAST.
Classification logic has been moved to separate plasmid_classification.py module.
"""

import argparse
import pandas as pd
import sys
import os
import gzip
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML
import tempfile

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Map Phoenix AMR genes onto ONT complete genomes (gene mapping only)'
    )
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--gamma_ar', required=True, help='GAMMA AR output file (.gamma)')
    parser.add_argument('--amrfinder_report', required=False, help='AMRFinder report file')
    parser.add_argument('--illumina_assembly', required=True, help='illumina assembly FASTA file')
    parser.add_argument('--ont_genome', required=True, help='ONT complete genome FASTA file')
    parser.add_argument('--output_mappings', required=True, help='Output gene mappings TSV')
    parser.add_argument('--min_identity', type=float, default=95.0, 
                       help='Minimum BLAST identity percentage (default: 95.0)')
    parser.add_argument('--min_coverage', type=float, default=90.0,
                       help='Minimum BLAST coverage percentage (default: 90.0)')
    return parser.parse_args()

def parse_gamma_ar_file(gamma_file):
    """
    Parse GAMMA AR file to extract gene information with coordinates
    
    Args:
        gamma_file: Path to GAMMA .gamma file
        
    Returns:
        List of gene dictionaries with illumina assembly coordinates
    """
    genes = []
    
    try:
        with open(gamma_file, 'r') as f:
            header = next(f).strip().split('\t')
            
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                
                # Parse the gene identifier (format: database__version__gene__accession__category)
                gene_id_parts = fields[0].split('__')
                if len(gene_id_parts) >= 5:
                    database = gene_id_parts[0]
                    version = gene_id_parts[1] 
                    gene_name = gene_id_parts[2]
                    accession = gene_id_parts[3]
                    category = gene_id_parts[4]
                else:
                    gene_name = fields[0]
                    category = 'Unknown'
                    database = 'Unknown'
                    accession = 'Unknown'
                
                # Extract relevant information
                illumina_contig = fields[1] if len(fields) > 1 else 'Unknown'
                start_pos = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
                end_pos = int(fields[3]) if len(fields) > 3 and fields[3].isdigit() else 0
                
                # Quality metrics from GAMMA
                percent_identity = float(fields[9]) * 100 if len(fields) > 9 and fields[9] != '' else 0
                percent_length = float(fields[11]) * 100 if len(fields) > 11 and fields[11] != '' else 0
                
                # Apply Phoenix's filtering criteria (90% length, 98% identity)
                if percent_length >= 90 and percent_identity >= 98:
                    genes.append({
                        'gene_name': gene_name,
                        'gene_id': fields[0],
                        'illumina_contig': illumina_contig,
                        'phoenix_start': start_pos,
                        'phoenix_end': end_pos,
                        'category': category,
                        'database': database,
                        'accession': accession,
                        'gamma_identity': round(percent_identity, 2),
                        'gamma_coverage': round(percent_length, 2),
                        'is_beta_lactam': 'LACTAM' in category.upper(),
                        'source': 'GAMMA'
                    })
    
    except Exception as e:
        print(f"Error parsing GAMMA file {gamma_file}: {e}", file=sys.stderr)
        return []
    
    return genes

def parse_amrfinder_report(amrfinder_file):
    """
    Parse AMRFinder report for additional AMR information
    
    Args:
        amrfinder_file: Path to AMRFinder report
        
    Returns:
        List of AMR entry dictionaries with Phoenix assembly coordinates
    """
    amr_data = []
    
    if not amrfinder_file or not Path(amrfinder_file).exists():
        return amr_data
    
    try:
        with open(amrfinder_file, 'r') as f:
            header = next(f).strip().split('\t')
            
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue
                
                # AMRFinder format: Protein identifier, Contig id, Start, Stop, Strand, Gene symbol, ...
                illumina_contig = fields[1] if len(fields) > 1 else 'Unknown'
                start_pos = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
                end_pos = int(fields[3]) if len(fields) > 3 and fields[3].isdigit() else 0
                gene_symbol = fields[5] if len(fields) > 5 else 'Unknown'
                
                # Check if it's a point mutation
                is_point_mutation = 'POINT' in line
                
                # Extract class if available
                gene_class = fields[10] if len(fields) > 10 else 'Unknown'
                
                amr_data.append({
                    'gene_name': gene_symbol,
                    'gene_id': f"AMRFinder_{gene_symbol}",
                    'illumina_contig': illumina_contig,
                    'phoenix_start': start_pos,
                    'phoenix_end': end_pos,
                    'category': gene_class,
                    'is_point_mutation': is_point_mutation,
                    'is_beta_lactam': 'BETA-LACTAM' in gene_class.upper() if gene_class != 'Unknown' else False,
                    'source': 'AMRFinder'
                })
    
    except Exception as e:
        print(f"Error parsing AMRFinder file {amrfinder_file}: {e}", file=sys.stderr)
        return []
    
    return amr_data


def open_maybe_gzip(path):
    # Read first two bytes to detect gzip magic number
    with open(path, "rb") as f:
        start = f.read(2)

    if start == b"\x1f\x8b":  # gzip magic number
        return gzip.open(path, "rt")
    else:
        return open(path, "r")



def get_contig_info(assembly_file):
    """
    Get contig information from assembly FASTA
    
    Args:
        assembly_file: Path to FASTA file
        
    Returns:
        Dictionary mapping contig names to their info (length, sequence)
    """
    contig_info = {}

    try:
        # Detect gzip by file extension

        with open_maybe_gzip(assembly_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contig_info[record.id] = {
                    'length': len(record.seq),
                    'description': record.description,
                    'sequence': str(record.seq)
                }


    except Exception as e:
        print(f"Error reading assembly file {assembly_file}: {e}", file=sys.stderr)
        return {}
    
    return contig_info

def extract_gene_sequence(illumina_contigs, contig_name, start_pos, end_pos):
    """
    Extract gene sequence from Phoenix assembly
    
    Args:
        illumina_contigs: Dictionary of contig information
        contig_name: Name of contig containing gene
        start_pos: Gene start position (1-based)
        end_pos: Gene end position (1-based, inclusive)
        
    Returns:
        Gene sequence string or None if extraction fails
    """
    if contig_name not in illumina_contigs:
        print(f"{contig_name} not in contigs")
        return None
    
    contig_seq = illumina_contigs[contig_name]['sequence']
    
    # Ensure coordinates are within bounds
    start_pos = max(0, start_pos - 1)  # Convert to 0-based indexing
    end_pos = min(len(contig_seq), end_pos)
    
    if start_pos >= end_pos:
        return None
    
    return contig_seq[start_pos:end_pos]

def run_blast_search(query_seq, ont_genome_file, temp_dir, gene_name="query"):
    """
    Run BLAST search to find gene location in ONT genome
    
    Args:
        query_seq: Gene sequence to search for
        ont_genome_file: Path to ONT genome FASTA
        temp_dir: Temporary directory for BLAST files
        gene_name: Name of gene (for logging)
        
    Returns:
        Dictionary with BLAST results or None if no hit found
    """
    
    # Create temporary query file
    query_file = os.path.join(temp_dir, f"query_{gene_name}.fasta")
    with open(query_file, 'w') as f:
        f.write(f">{gene_name}\n{query_seq}\n")
    
    # Create BLAST database (only once per temp_dir)
    db_file = os.path.join(temp_dir, "ont_db")
    if not os.path.exists(f"{db_file}.nhr"):
        makeblastdb_cmd = f"makeblastdb -in {ont_genome_file} -dbtype nucl -out {db_file} -parse_seqids"
        result = subprocess.run(makeblastdb_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error creating BLAST database: {result.stderr}", file=sys.stderr)
            return None
    
    # Run BLAST search
    blast_output = os.path.join(temp_dir, f"blast_{gene_name}.xml")
    blastn_cmd = [
        "blastn",
        "-query", query_file,
        "-db", db_file,
        "-out", blast_output,
        "-outfmt", "5",
        "-evalue", "1e-10",
        "-max_target_seqs", "1"
    ]
    
    result = subprocess.run(blastn_cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running BLAST for {gene_name}: {result.stderr}", file=sys.stderr)
        return None
    
    # Parse BLAST results
    try:
        with open(blast_output, 'r') as f:
            blast_records = NCBIXML.parse(f)
            for blast_record in blast_records:
                if blast_record.alignments:
                    # Get best hit
                    alignment = blast_record.alignments[0]
                    hsp = alignment.hsps[0]
                    
                    return {
                        'ont_contig': alignment.title.split()[0].replace('>', ''),
                        'ont_start': hsp.sbjct_start,
                        'ont_end': hsp.sbjct_end,
                        'blast_identity': round(hsp.identities / hsp.align_length * 100, 2),
                        'blast_coverage': round(hsp.align_length / blast_record.query_length * 100, 2),
                        'blast_evalue': hsp.expect,
                        'blast_bitscore': hsp.bits
                    }
    except Exception as e:
        print(f"Error parsing BLAST results for {gene_name}: {e}", file=sys.stderr)
        return None
    
    return None

def map_amr_genes_to_ont(all_genes, illumina_contigs, ont_genome_file, sample_id, min_identity=95.0, min_coverage=90.0):
    """
    Map AMR genes from Phoenix assembly to ONT complete genome using BLAST
    
    Args:
        all_genes: List of gene dictionaries (from GAMMA and AMRFinder)
        illumina_contigs: Dictionary of Phoenix assembly contigs
        ont_genome_file: Path to ONT genome FASTA
        min_identity: Minimum BLAST identity percentage
        min_coverage: Minimum BLAST coverage percentage
        
    Returns:
        List of mapped gene dictionaries with ONT coordinates
    """
    mapped_genes = []
    with tempfile.TemporaryDirectory() as temp_dir:
        #illumina_contigs = {"_".join(k.split("_")[1:-2]): v for k, v in illumina_contigs.items()}
        illumina_contigs = {"_".join(k.replace(f"{sample_id}_", "")
                                    .replace("NODE_", "")
                                    .replace("EDGE_", "")
                                    .split("_")[:-2]): v for k, v in illumina_contigs.items()}

        for idx, gene in enumerate(all_genes, 1):
            gene_name = gene['gene_name']
            gene["illumina_contig"] = gene["illumina_contig"].replace(f"{sample_id}_", "")
            print("here2")
            print(gene["illumina_contig"])

            print(f"[{idx}/{len(all_genes)}] Processing gene: {gene_name}")
            
            # Extract gene sequence from Phoenix assembly
            gene_seq = extract_gene_sequence(
                illumina_contigs, 
                gene['illumina_contig'], 
                gene['phoenix_start'], 
                gene['phoenix_end']
            )
            if gene_name != "blaOXY-1-1_NG_049841.1":
                #print(illumina_contigs)
                print("here")
                print(gene['illumina_contig'], gene['phoenix_start'], gene['phoenix_end'])
                
            
            if not gene_seq:
                print(f"  ⚠ Could not extract sequence for gene {gene_name}", file=sys.stderr)
                # Create entry for unmapped gene
                mapped_gene = gene.copy()
                mapped_gene.update({
                    'ont_contig': None,
                    'ont_start': None,
                    'ont_end': None,
                    'blast_identity': None,
                    'blast_coverage': None,
                    'blast_evalue': None,
                    'blast_bitscore': None,
                    'mapping_status': 'not_mapped',
                    'mapping_failure_reason': 'sequence_extraction_failed'
                })
                mapped_genes.append(mapped_gene)
                continue
            
            # BLAST against ONT genome
            blast_result = run_blast_search(gene_seq, ont_genome_file, temp_dir, gene_name)
            
            if blast_result and blast_result['blast_identity'] >= min_identity and blast_result['blast_coverage'] >= min_coverage:
                # Successfully mapped
                print(f"  ✓ Mapped to {blast_result['ont_contig']} "
                      f"({blast_result['blast_identity']:.1f}% identity, "
                      f"{blast_result['blast_coverage']:.1f}% coverage)")
                
                mapped_gene = gene.copy()
                mapped_gene.update({
                    'ont_contig': blast_result['ont_contig'],
                    'ont_start': blast_result['ont_start'],
                    'ont_end': blast_result['ont_end'],
                    'blast_identity': blast_result['blast_identity'],
                    'blast_coverage': blast_result['blast_coverage'],
                    'blast_evalue': blast_result['blast_evalue'],
                    'blast_bitscore': blast_result['blast_bitscore'],
                    'mapping_status': 'mapped',
                    'mapping_failure_reason': None
                })
                mapped_genes.append(mapped_gene)
                
            elif blast_result:
                # Hit found but below thresholds
                print(f"  ✗ Hit found but below thresholds "
                      f"({blast_result['blast_identity']:.1f}% identity, "
                      f"{blast_result['blast_coverage']:.1f}% coverage)")
                
                mapped_gene = gene.copy()
                mapped_gene.update({
                    'ont_contig': blast_result['ont_contig'],
                    'ont_start': blast_result['ont_start'],
                    'ont_end': blast_result['ont_end'],
                    'blast_identity': blast_result['blast_identity'],
                    'blast_coverage': blast_result['blast_coverage'],
                    'blast_evalue': blast_result['blast_evalue'],
                    'blast_bitscore': blast_result['blast_bitscore'],
                    'mapping_status': 'below_threshold',
                    'mapping_failure_reason': f'identity={blast_result["blast_identity"]:.1f}%_coverage={blast_result["blast_coverage"]:.1f}%'
                })
                mapped_genes.append(mapped_gene)
                
            else:
                # No BLAST hit found
                print(f"  ✗ No BLAST hit found in ONT genome")
                
                mapped_gene = gene.copy()
                mapped_gene.update({
                    'ont_contig': None,
                    'ont_start': None,
                    'ont_end': None,
                    'blast_identity': None,
                    'blast_coverage': None,
                    'blast_evalue': None,
                    'blast_bitscore': None,
                    'mapping_status': 'not_mapped',
                    'mapping_failure_reason': 'no_blast_hit'
                })
                mapped_genes.append(mapped_gene)
    
    return mapped_genes

def main():
    args = parse_args()
    
    print(f"=" * 80)
    print(f"Phoenix AMR Locator (Refactored)")
    print(f"Sample: {args.sample_id}")
    print(f"=" * 80)
    
    # Parse input files
    print("\n[1/5] Parsing GAMMA AR file...")
    gamma_genes = parse_gamma_ar_file(args.gamma_ar)
    print(f"  Found {len(gamma_genes)} genes from GAMMA")
    
    print("\n[2/5] Parsing AMRFinder report...")
    amrfinder_data = parse_amrfinder_report(args.amrfinder_report) if args.amrfinder_report else []
    print(f"  Found {len(amrfinder_data)} entries from AMRFinder")
    
    # Deduplicate: remove AMRFinder genes that are already in GAMMA
    amrfinder_unique = []
    for amr_gene in amrfinder_data:
        # Check if gene already exists in GAMMA results
        is_duplicate = any(
            gamma_gene['illumina_contig'] == amr_gene['illumina_contig'] and 
            gamma_gene['gene_name'].lower() == amr_gene['gene_name'].lower()
            for gamma_gene in gamma_genes
        )
        if not is_duplicate:
            amrfinder_unique.append(amr_gene)
    
    print(f"  Unique AMRFinder genes (not in GAMMA): {len(amrfinder_unique)}")
    
    # Combine all genes
    all_genes = gamma_genes + amrfinder_unique
    print(f"  Total unique genes to process: {len(all_genes)}")
    
    print("\n[3/5] Loading assembly files...")
    illumina_contigs = get_contig_info(args.illumina_assembly)
    print(f"  Phoenix assembly contigs: {len(illumina_contigs)}")
    
    ont_contigs = get_contig_info(args.ont_genome)
    print(f"  ONT genome contigs: {len(ont_contigs)}")
    
    # Map AMR genes to ONT genome
    print(f"\n[4/5] Mapping AMR genes to ONT genome (min_identity={args.min_identity}%, min_coverage={args.min_coverage}%)...")
    mapped_genes = map_amr_genes_to_ont(
        all_genes, 
        illumina_contigs, 
        args.ont_genome,
        args.sample_id,
        args.min_identity,
        args.min_coverage
    )
    
    # Save results
    print(f"\n[5/5] Saving results...")
    results_df = pd.DataFrame(mapped_genes)
    
    # Reorder columns for clarity
    column_order = [
        'gene_name', 'gene_id', 'source', 'category', 'is_beta_lactam',
        'illumina_contig', 'phoenix_start', 'phoenix_end',
        'ont_contig', 'ont_start', 'ont_end',
        'blast_identity', 'blast_coverage', 'blast_evalue', 'blast_bitscore',
        'mapping_status', 'mapping_failure_reason'
    ]
    
    # Add any additional columns not in the predefined order
    for col in results_df.columns:
        if col not in column_order:
            column_order.append(col)
    
    # Reorder columns (only include columns that exist)
    column_order = [col for col in column_order if col in results_df.columns]
    results_df = results_df[column_order]
    
    results_df.to_csv(args.output_mappings, sep='\t', index=False)
    print(f"  ✓ Saved mappings to: {args.output_mappings}")
    
    # Print summary
    print(f"\n" + "=" * 80)
    print(f"Summary for sample {args.sample_id}:")
    print(f"=" * 80)
    print(f"  Total AMR genes processed: {len(mapped_genes)}")
    
    if mapped_genes:
        status_counts = results_df['mapping_status'].value_counts()
        for status, count in status_counts.items():
            print(f"  {status}: {count}")
        
        # Count mapped genes only
        mapped_df = results_df[results_df['mapping_status'] == 'mapped']
        if not mapped_df.empty:
            print(f"\nMapped genes breakdown:")
            
            # By source
            source_counts = mapped_df['source'].value_counts()
            for source, count in source_counts.items():
                print(f"  From {source}: {count}")
            
            # Beta-lactam count
            beta_lactam_count = mapped_df['is_beta_lactam'].sum() if 'is_beta_lactam' in mapped_df.columns else 0
            print(f"  Beta-lactam genes: {beta_lactam_count}")
            
            # By category (top 5)
            if 'category' in mapped_df.columns:
                category_counts = mapped_df['category'].value_counts().head(5)
                print(f"\nTop 5 categories:")
                for category, count in category_counts.items():
                    print(f"  {category}: {count}")
    
    print(f"\n✓ Pipeline complete!")

if __name__ == '__main__':
    main()
