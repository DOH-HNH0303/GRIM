#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import tempfile
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Map Phoenix AMR genes onto ONT complete genomes')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--gamma_ar', required=True, help='GAMMA AR output file (.gamma)')
    parser.add_argument('--amrfinder_report', required=False, help='AMRFinder report file')
    parser.add_argument('--phoenix_assembly', required=True, help='Phoenix assembly FASTA file')
    parser.add_argument('--ont_genome', required=True, help='ONT complete genome FASTA file')
    parser.add_argument('--output_locations', required=True, help='Output gene locations TSV')
    parser.add_argument('--output_detailed', required=True, help='Output detailed AMR TSV')
    return parser.parse_args()

def classify_ont_contig(contig_name, contig_length):
    """Classify ONT contig as chromosome or plasmid based on naming and size"""
    contig_name_lower = contig_name.lower()
    
    # Common plasmid indicators in contig names
    plasmid_indicators = ['plasmid', 'plas', 'p1', 'p2', 'p3', 'p4', 'p5', 'unnamed']
    chromosome_indicators = ['chromosome', 'chr', 'chrom', 'complete', 'genome']
    
    # Check name-based classification first
    for indicator in plasmid_indicators:
        if indicator in contig_name_lower:
            return 'plasmid', contig_name
    
    for indicator in chromosome_indicators:
        if indicator in contig_name_lower:
            return 'chromosome', contig_name
    
    # Size-based heuristic for complete genomes (adjust thresholds as needed)
    if contig_length > 2000000:  # > 2 Mb likely chromosome
        return 'chromosome', contig_name
    elif contig_length < 1000000:  # < 1 Mb likely plasmid
        return 'plasmid', contig_name
    else:
        return 'unknown', contig_name

def parse_gamma_ar_file(gamma_file):
    """Parse GAMMA AR file to extract gene information with coordinates"""
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
                phoenix_contig = fields[1] if len(fields) > 1 else 'Unknown'
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
                        'phoenix_contig': phoenix_contig,
                        'phoenix_start_pos': start_pos,
                        'phoenix_end_pos': end_pos,
                        'category': category,
                        'database': database,
                        'accession': accession,
                        'percent_identity': round(percent_identity, 2),
                        'percent_length': round(percent_length, 2),
                        'is_beta_lactam': 'LACTAM' in category.upper()
                    })
    
    except Exception as e:
        print(f"Error parsing GAMMA file {gamma_file}: {e}", file=sys.stderr)
        return []
    
    return genes

def parse_amrfinder_report(amrfinder_file):
    """Parse AMRFinder report for additional AMR information"""
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
                phoenix_contig = fields[1] if len(fields) > 1 else 'Unknown'
                start_pos = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
                end_pos = int(fields[3]) if len(fields) > 3 and fields[3].isdigit() else 0
                gene_symbol = fields[5] if len(fields) > 5 else 'Unknown'
                
                # Check if it's a point mutation
                is_point_mutation = 'POINT' in line
                
                amr_data.append({
                    'gene_symbol': gene_symbol,
                    'phoenix_contig': phoenix_contig,
                    'phoenix_start_pos': start_pos,
                    'phoenix_end_pos': end_pos,
                    'is_point_mutation': is_point_mutation,
                    'full_line': line.strip()
                })
    
    except Exception as e:
        print(f"Error parsing AMRFinder file {amrfinder_file}: {e}", file=sys.stderr)
        return []
    
    return amr_data

def get_contig_info(assembly_file):
    """Get contig information from assembly FASTA"""
    contig_info = {}
    
    try:
        for record in SeqIO.parse(assembly_file, "fasta"):
            contig_info[record.id] = {
                'length': len(record.seq),
                'description': record.description,
                'sequence': str(record.seq)
            }
    except Exception as e:
        print(f"Error reading assembly file {assembly_file}: {e}", file=sys.stderr)
        return {}
    
    return contig_info

def extract_gene_sequence(phoenix_contigs, contig_name, start_pos, end_pos):
    """Extract gene sequence from Phoenix assembly"""
    if contig_name not in phoenix_contigs:
        return None
    
    contig_seq = phoenix_contigs[contig_name]['sequence']
    
    # Ensure coordinates are within bounds
    start_pos = max(0, start_pos - 1)  # Convert to 0-based indexing
    end_pos = min(len(contig_seq), end_pos)
    
    if start_pos >= end_pos:
        return None
    
    return contig_seq[start_pos:end_pos]

def run_blast_search(query_seq, ont_genome_file, temp_dir):
    """Run BLAST search to find gene location in ONT genome"""
    
    # Create temporary query file
    query_file = os.path.join(temp_dir, "query.fasta")
    with open(query_file, 'w') as f:
        f.write(f">query\n{query_seq}\n")
    
    # Create BLAST database
    db_file = os.path.join(temp_dir, "ont_db")
    makeblastdb_cmd = f"makeblastdb -in {ont_genome_file} -dbtype nucl -out {db_file}"
    subprocess.run(makeblastdb_cmd, shell=True, capture_output=True)
    
    # Run BLAST search
    blast_output = os.path.join(temp_dir, "blast_results.xml")
    blastn_cmd = f"blastn -query {query_file} -db {db_file} -out {blast_output} -outfmt 5 -evalue 1e-10"
    result = subprocess.run(blastn_cmd, shell=True, capture_output=True)
    
    if result.returncode != 0:
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
                        'identity': hsp.identities / hsp.align_length * 100,
                        'coverage': hsp.align_length / blast_record.query_length * 100,
                        'evalue': hsp.expect
                    }
    except Exception as e:
        print(f"Error parsing BLAST results: {e}", file=sys.stderr)
        return None
    
    return None

def map_amr_genes_to_ont(gamma_genes, amrfinder_data, phoenix_contigs, ont_contigs, ont_genome_file):
    """Map AMR genes from Phoenix assembly to ONT complete genome"""
    mapped_genes = []
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Process GAMMA genes
        for gene in gamma_genes:
            print(f"Processing gene: {gene['gene_name']}")
            
            # Extract gene sequence from Phoenix assembly
            gene_seq = extract_gene_sequence(
                phoenix_contigs, 
                gene['phoenix_contig'], 
                gene['phoenix_start_pos'], 
                gene['phoenix_end_pos']
            )
            
            if not gene_seq:
                print(f"Could not extract sequence for gene {gene['gene_name']}")
                continue
            
            # BLAST against ONT genome
            blast_result = run_blast_search(gene_seq, ont_genome_file, temp_dir)
            
            if blast_result and blast_result['identity'] >= 95 and blast_result['coverage'] >= 90:
                # Get ONT contig info
                ont_contig_name = blast_result['ont_contig']
                ont_contig_length = ont_contigs.get(ont_contig_name, {}).get('length', 'unknown')
                
                # Classify location on ONT genome
                location_type, location_name = classify_ont_contig(ont_contig_name, ont_contig_length if ont_contig_length != 'unknown' else 0)
                
                mapped_gene = gene.copy()
                mapped_gene.update({
                    'ont_contig': ont_contig_name,
                    'ont_start_pos': blast_result['ont_start'],
                    'ont_end_pos': blast_result['ont_end'],
                    'ont_contig_length': ont_contig_length,
                    'location_type': location_type,
                    'blast_identity': round(blast_result['identity'], 2),
                    'blast_coverage': round(blast_result['coverage'], 2),
                    'blast_evalue': blast_result['evalue'],
                    'mapping_status': 'mapped'
                })
                mapped_genes.append(mapped_gene)
            else:
                # Gene not found in ONT genome
                mapped_gene = gene.copy()
                mapped_gene.update({
                    'ont_contig': 'Not found',
                    'ont_start_pos': 'N/A',
                    'ont_end_pos': 'N/A',
                    'ont_contig_length': 'N/A',
                    'location_type': 'not_found',
                    'blast_identity': blast_result['identity'] if blast_result else 0,
                    'blast_coverage': blast_result['coverage'] if blast_result else 0,
                    'blast_evalue': blast_result['evalue'] if blast_result else 'N/A',
                    'mapping_status': 'not_mapped'
                })
                mapped_genes.append(mapped_gene)
        
        # Process AMRFinder data similarly
        for amr_entry in amrfinder_data:
            # Check if this gene is already in GAMMA results
            gene_in_gamma = any(g['phoenix_contig'] == amr_entry['phoenix_contig'] and 
                               g['gene_name'].lower() == amr_entry['gene_symbol'].lower() 
                               for g in gamma_genes)
            
            if not gene_in_gamma:
                print(f"Processing AMRFinder gene: {amr_entry['gene_symbol']}")
                
                # Extract gene sequence from Phoenix assembly
                gene_seq = extract_gene_sequence(
                    phoenix_contigs, 
                    amr_entry['phoenix_contig'], 
                    amr_entry['phoenix_start_pos'], 
                    amr_entry['phoenix_end_pos']
                )
                
                if not gene_seq:
                    continue
                
                # BLAST against ONT genome
                blast_result = run_blast_search(gene_seq, ont_genome_file, temp_dir)
                
                if blast_result and blast_result['identity'] >= 95 and blast_result['coverage'] >= 90:
                    # Get ONT contig info
                    ont_contig_name = blast_result['ont_contig']
                    ont_contig_length = ont_contigs.get(ont_contig_name, {}).get('length', 'unknown')
                    
                    # Classify location on ONT genome
                    location_type, location_name = classify_ont_contig(ont_contig_name, ont_contig_length if ont_contig_length != 'unknown' else 0)
                    
                    mapped_gene = {
                        'gene_name': amr_entry['gene_symbol'],
                        'gene_id': f"AMRFinder_{amr_entry['gene_symbol']}",
                        'phoenix_contig': amr_entry['phoenix_contig'],
                        'phoenix_start_pos': amr_entry['phoenix_start_pos'],
                        'phoenix_end_pos': amr_entry['phoenix_end_pos'],
                        'ont_contig': ont_contig_name,
                        'ont_start_pos': blast_result['ont_start'],
                        'ont_end_pos': blast_result['ont_end'],
                        'ont_contig_length': ont_contig_length,
                        'location_type': location_type,
                        'blast_identity': round(blast_result['identity'], 2),
                        'blast_coverage': round(blast_result['coverage'], 2),
                        'blast_evalue': blast_result['evalue'],
                        'mapping_status': 'mapped',
                        'category': 'AMRFinder',
                        'database': 'AMRFinder',
                        'accession': 'N/A',
                        'percent_identity': 'N/A',
                        'percent_length': 'N/A',
                        'is_beta_lactam': 'Unknown',
                        'is_point_mutation': amr_entry['is_point_mutation']
                    }
                    mapped_genes.append(mapped_gene)
    
    return mapped_genes

def main():
    args = parse_args()
    
    print(f"Processing sample: {args.sample_id}")
    print(f"Mapping Phoenix AMR results to ONT complete genome")
    
    # Parse input files
    gamma_genes = parse_gamma_ar_file(args.gamma_ar)
    amrfinder_data = parse_amrfinder_report(args.amrfinder_report) if args.amrfinder_report else []
    phoenix_contigs = get_contig_info(args.phoenix_assembly)
    ont_contigs = get_contig_info(args.ont_genome)
    
    print(f"Found {len(gamma_genes)} genes from GAMMA")
    print(f"Found {len(amrfinder_data)} entries from AMRFinder")
    print(f"Found {len(phoenix_contigs)} contigs in Phoenix assembly")
    print(f"Found {len(ont_contigs)} contigs in ONT genome")
    
    # Map AMR genes to ONT genome
    mapped_genes = map_amr_genes_to_ont(gamma_genes, amrfinder_data, phoenix_contigs, ont_contigs, args.ont_genome)
    
    # Prepare results
    results = []
    detailed_results = []
    
    for gene in mapped_genes:
        # Create result entry
        result = {
            'sample_id': args.sample_id,
            'gene_name': gene['gene_name'],
            'gene_id': gene['gene_id'],
            'location_type': gene['location_type'],
            'ont_contig_name': gene['ont_contig'],
            'ont_contig_length': gene['ont_contig_length'],
            'ont_start_position': gene['ont_start_pos'],
            'ont_end_position': gene['ont_end_pos'],
            'phoenix_contig_name': gene['phoenix_contig'],
            'phoenix_start_position': gene['phoenix_start_pos'],
            'phoenix_end_position': gene['phoenix_end_pos'],
            'blast_identity': gene.get('blast_identity', 'N/A'),
            'blast_coverage': gene.get('blast_coverage', 'N/A'),
            'mapping_status': gene['mapping_status'],
            'is_beta_lactam': gene.get('is_beta_lactam', 'Unknown'),
            'source': 'GAMMA' if 'AMRFinder' not in gene['gene_id'] else 'AMRFinder'
        }
        
        results.append(result)
        
        # Add to detailed results
        detailed_result = result.copy()
        detailed_result.update({
            'category': gene.get('category', 'Unknown'),
            'database': gene.get('database', 'Unknown'),
            'accession': gene.get('accession', 'Unknown'),
            'phoenix_percent_identity': gene.get('percent_identity', 'N/A'),
            'phoenix_percent_length': gene.get('percent_length', 'N/A'),
            'blast_evalue': gene.get('blast_evalue', 'N/A')
        })
        detailed_results.append(detailed_result)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output_locations, sep='\t', index=False)
    
    detailed_df = pd.DataFrame(detailed_results)
    detailed_df.to_csv(args.output_detailed, sep='\t', index=False)
    
    # Print summary
    print(f"\nResults for sample {args.sample_id}:")
    print(f"  Total AMR genes processed: {len(mapped_genes)}")
    
    if results:
        mapped_count = sum(1 for r in results if r['mapping_status'] == 'mapped')
        print(f"  Successfully mapped to ONT genome: {mapped_count}")
        print(f"  Not found in ONT genome: {len(results) - mapped_count}")
        
        location_counts = results_df[results_df['mapping_status'] == 'mapped']['location_type'].value_counts()
        for location, count in location_counts.items():
            print(f"  {location}: {count}")
        
        beta_lactam_count = sum(1 for r in results if r.get('is_beta_lactam', False) and r['mapping_status'] == 'mapped')
        print(f"  Beta-lactam genes mapped: {beta_lactam_count}")
        
        source_counts = results_df[results_df['mapping_status'] == 'mapped']['source'].value_counts()
        for source, count in source_counts.items():
            print(f"  From {source}: {count}")

if __name__ == '__main__':
    main()