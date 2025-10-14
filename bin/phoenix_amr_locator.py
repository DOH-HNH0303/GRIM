#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
from pathlib import Path
from Bio import SeqIO
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Locate AMR genes using existing Phoenix output files')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--gamma_ar', required=True, help='GAMMA AR output file (.gamma)')
    parser.add_argument('--amrfinder_report', required=False, help='AMRFinder report file')
    parser.add_argument('--assembly', required=True, help='Assembly FASTA file')
    parser.add_argument('--output_locations', required=True, help='Output gene locations TSV')
    parser.add_argument('--output_detailed', required=True, help='Output detailed AMR TSV')
    return parser.parse_args()

def classify_contig(contig_name, contig_length):
    """Classify contig as chromosome or plasmid based on naming and size"""
    contig_name_lower = contig_name.lower()
    
    # Common plasmid indicators in contig names
    plasmid_indicators = ['plasmid', 'plas', 'p1', 'p2', 'p3', 'p4', 'p5']
    chromosome_indicators = ['chromosome', 'chr', 'chrom']
    
    # Check name-based classification first
    for indicator in plasmid_indicators:
        if indicator in contig_name_lower:
            return 'plasmid', contig_name
    
    for indicator in chromosome_indicators:
        if indicator in contig_name_lower:
            return 'chromosome', contig_name
    
    # Size-based heuristic (adjust thresholds as needed)
    if contig_length > 1000000:  # > 1 Mb likely chromosome
        return 'chromosome', contig_name
    elif contig_length < 500000:  # < 500 kb likely plasmid
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
                contig = fields[1] if len(fields) > 1 else 'Unknown'
                start_pos = fields[2] if len(fields) > 2 else 'Unknown'
                end_pos = fields[3] if len(fields) > 3 else 'Unknown'
                
                # Quality metrics from GAMMA
                percent_identity = float(fields[9]) * 100 if len(fields) > 9 and fields[9] != '' else 0
                percent_length = float(fields[11]) * 100 if len(fields) > 11 and fields[11] != '' else 0
                
                # Apply Phoenix's filtering criteria (90% length, 98% identity)
                if percent_length >= 90 and percent_identity >= 98:
                    genes.append({
                        'gene_name': gene_name,
                        'gene_id': fields[0],
                        'contig': contig,
                        'start_pos': start_pos,
                        'end_pos': end_pos,
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
                contig = fields[1] if len(fields) > 1 else 'Unknown'
                start_pos = fields[2] if len(fields) > 2 else 'Unknown'
                end_pos = fields[3] if len(fields) > 3 else 'Unknown'
                gene_symbol = fields[5] if len(fields) > 5 else 'Unknown'
                
                # Check if it's a point mutation
                is_point_mutation = 'POINT' in line
                
                amr_data.append({
                    'gene_symbol': gene_symbol,
                    'contig': contig,
                    'start_pos': start_pos,
                    'end_pos': end_pos,
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
                'description': record.description
            }
    except Exception as e:
        print(f"Error reading assembly file {assembly_file}: {e}", file=sys.stderr)
        return {}
    
    return contig_info

def main():
    args = parse_args()
    
    print(f"Processing sample: {args.sample_id}")
    
    # Parse input files
    gamma_genes = parse_gamma_ar_file(args.gamma_ar)
    amrfinder_data = parse_amrfinder_report(args.amrfinder_report) if args.amrfinder_report else []
    contig_info = get_contig_info(args.assembly)
    
    print(f"Found {len(gamma_genes)} genes from GAMMA")
    print(f"Found {len(amrfinder_data)} entries from AMRFinder")
    print(f"Found {len(contig_info)} contigs in assembly")
    
    # Process gene locations
    results = []
    detailed_results = []
    
    # Process GAMMA genes (primary source)
    for gene in gamma_genes:
        contig_name = gene['contig']
        contig_length = contig_info.get(contig_name, {}).get('length', 'unknown')
        
        # Classify location
        location_type, location_name = classify_contig(contig_name, contig_length if contig_length != 'unknown' else 0)
        
        # Create result entry
        result = {
            'sample_id': args.sample_id,
            'gene_name': gene['gene_name'],
            'gene_id': gene['gene_id'],
            'location_type': location_type,
            'contig_name': contig_name,
            'contig_length': contig_length,
            'start_position': gene['start_pos'],
            'end_position': gene['end_pos'],
            'percent_identity': gene['percent_identity'],
            'percent_length': gene['percent_length'],
            'is_beta_lactam': gene['is_beta_lactam'],
            'source': 'GAMMA'
        }
        
        results.append(result)
        
        # Add to detailed results
        detailed_result = result.copy()
        detailed_result.update({
            'category': gene['category'],
            'database': gene['database'],
            'accession': gene['accession']
        })
        detailed_results.append(detailed_result)
    
    # Process AMRFinder data for additional context
    for amr_entry in amrfinder_data:
        contig_name = amr_entry['contig']
        contig_length = contig_info.get(contig_name, {}).get('length', 'unknown')
        location_type, location_name = classify_contig(contig_name, contig_length if contig_length != 'unknown' else 0)
        
        # Check if this gene is already in GAMMA results
        gene_in_gamma = any(g['contig'] == contig_name and 
                           g['gene_name'].lower() == amr_entry['gene_symbol'].lower() 
                           for g in gamma_genes)
        
        if not gene_in_gamma:
            # Add AMRFinder-only genes
            result = {
                'sample_id': args.sample_id,
                'gene_name': amr_entry['gene_symbol'],
                'gene_id': f"AMRFinder_{amr_entry['gene_symbol']}",
                'location_type': location_type,
                'contig_name': contig_name,
                'contig_length': contig_length,
                'start_position': amr_entry['start_pos'],
                'end_position': amr_entry['end_pos'],
                'percent_identity': 'N/A',
                'percent_length': 'N/A',
                'is_beta_lactam': 'Unknown',
                'source': 'AMRFinder'
            }
            
            results.append(result)
            
            detailed_result = result.copy()
            detailed_result.update({
                'category': 'AMRFinder',
                'database': 'AMRFinder',
                'accession': 'N/A',
                'is_point_mutation': amr_entry['is_point_mutation']
            })
            detailed_results.append(detailed_result)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output_locations, sep='\t', index=False)
    
    detailed_df = pd.DataFrame(detailed_results)
    detailed_df.to_csv(args.output_detailed, sep='\t', index=False)
    
    # Print summary
    print(f"Results for sample {args.sample_id}:")
    print(f"  Total AMR genes found: {len(results)}")
    
    if results:
        location_counts = results_df['location_type'].value_counts()
        for location, count in location_counts.items():
            print(f"  {location}: {count}")
        
        beta_lactam_count = sum(1 for r in results if r.get('is_beta_lactam', False))
        print(f"  Beta-lactam genes: {beta_lactam_count}")
        
        source_counts = results_df['source'].value_counts()
        for source, count in source_counts.items():
            print(f"  From {source}: {count}")

if __name__ == '__main__':
    main()