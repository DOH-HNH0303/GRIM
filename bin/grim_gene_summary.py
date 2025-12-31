#!/usr/bin/env python3
"""
GRIM Gene Summary
Generates per-gene CSV with location and plasmid information.
Combines gene mappings, contig classifications, and replicon data.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate per-gene summary CSV for GRIM pipeline'
    )
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--gene_mappings', required=True, 
                       help='Gene mappings TSV from phoenix_amr_locator')
    parser.add_argument('--contig_classification', required=True,
                       help='Contig classification TSV from plasmid_classification')
    parser.add_argument('--plasmid_replicons', required=True,
                       help='Plasmid replicons TSV from plasmid_classification')
    parser.add_argument('--output_summary', required=True,
                       help='Output gene summary CSV')
    return parser.parse_args()

def load_gene_mappings(mappings_file):
    """Load gene mappings from TSV"""
    try:
        df = pd.read_csv(mappings_file, sep='\t')
        print(f"  Loaded {len(df)} gene mappings")
        return df
    except Exception as e:
        print(f"Error loading gene mappings: {e}", file=sys.stderr)
        sys.exit(1)

def load_contig_classification(classification_file):
    """Load contig classifications from TSV"""
    try:
        df = pd.read_csv(classification_file, sep='\t')
        print(f"  Loaded {len(df)} contig classifications")
        return df
    except Exception as e:
        print(f"Error loading contig classifications: {e}", file=sys.stderr)
        sys.exit(1)

def load_plasmid_replicons(replicons_file):
    """Load plasmid replicons from TSV"""
    try:
        if Path(replicons_file).stat().st_size == 0:
            print(f"  No plasmid replicons found (empty file)")
            return pd.DataFrame()
        
        df = pd.read_csv(replicons_file, sep='\t')
        print(f"  Loaded {len(df)} plasmid replicons")
        return df
    except Exception as e:
        print(f"Error loading plasmid replicons: {e}", file=sys.stderr)
        # Return empty DataFrame instead of exiting
        return pd.DataFrame()

def generate_gene_summary(sample_id, gene_mappings_df, classification_df, replicons_df):
    """
    Generate per-gene summary with location and plasmid information
    
    Args:
        sample_id: Sample identifier
        gene_mappings_df: DataFrame with gene mappings
        classification_df: DataFrame with contig classifications
        replicons_df: DataFrame with plasmid replicon information
        
    Returns:
        DataFrame with per-gene summary
    """
    print(f"\n[Processing] Generating gene summary...")
    
    summary_rows = []
    
    for _, gene in gene_mappings_df.iterrows():
        # Base information
        summary_row = {
            'sample_id': sample_id,
            'gene_name': gene['gene_name'],
            'gene_category': gene.get('category', 'Unknown'),
            'source': gene.get('source', 'Unknown'),
            'mapping_status': gene.get('mapping_status', 'unknown')
        }
        
        # Check if gene was successfully mapped
        if gene.get('mapping_status') == 'mapped' and pd.notna(gene.get('ont_contig')):
            ont_contig = gene['ont_contig']
            
            # Get contig classification
            contig_info = classification_df[
                classification_df['original_contig_name'] == ont_contig
            ]
            
            if not contig_info.empty:
                contig_type = contig_info.iloc[0]['classification']
                summary_row['ont_contig_name'] = ont_contig
                summary_row['ont_contig_type'] = contig_type
                
                # Get plasmid information if contig is a plasmid
                if contig_type == 'plasmid':
                    plasmid_info = replicons_df[
                        replicons_df['original_contig_name'] == ont_contig
                    ]
                    
                    if not plasmid_info.empty:
                        summary_row['plasmid_name'] = plasmid_info.iloc[0].get('plasmid_name', 'N/A')
                        summary_row['replicon_type'] = plasmid_info.iloc[0].get('replicon_type', 'N/A')
                        summary_row['relaxase_type'] = plasmid_info.iloc[0].get('relaxase_type', 'N/A')
                        summary_row['predicted_mobility'] = plasmid_info.iloc[0].get('predicted_mobility', 'N/A')
                    else:
                        # Plasmid but no replicon info (possible for untyped plasmids)
                        summary_row['plasmid_name'] = 'Unknown_plasmid'
                        summary_row['replicon_type'] = 'N/A'
                        summary_row['relaxase_type'] = 'N/A'
                        summary_row['predicted_mobility'] = 'N/A'
                else:
                    # Chromosomal gene
                    summary_row['plasmid_name'] = 'N/A'
                    summary_row['replicon_type'] = 'N/A'
                    summary_row['relaxase_type'] = 'N/A'
                    summary_row['predicted_mobility'] = 'N/A'
            else:
                # Contig not found in classification (shouldn't happen)
                summary_row['ont_contig_name'] = ont_contig
                summary_row['ont_contig_type'] = 'unknown'
                summary_row['plasmid_name'] = 'N/A'
                summary_row['replicon_type'] = 'N/A'
                summary_row['relaxase_type'] = 'N/A'
                summary_row['predicted_mobility'] = 'N/A'
            
            # Add position information
            ont_start = gene.get('ont_start')
            ont_end = gene.get('ont_end')
            
            if pd.notna(ont_start) and pd.notna(ont_end):
                summary_row['ont_position'] = f"{int(ont_start)}-{int(ont_end)}"
            else:
                summary_row['ont_position'] = 'N/A'
            
            # Add Phoenix information
            summary_row['phoenix_contig'] = gene.get('phoenix_contig', 'N/A')
            summary_row['blast_identity'] = gene.get('blast_identity', 'N/A')
            summary_row['blast_coverage'] = gene.get('blast_coverage', 'N/A')
            
        else:
            # Gene not mapped
            summary_row['ont_contig_name'] = 'Not mapped'
            summary_row['ont_contig_type'] = 'N/A'
            summary_row['plasmid_name'] = 'N/A'
            summary_row['replicon_type'] = 'N/A'
            summary_row['relaxase_type'] = 'N/A'
            summary_row['predicted_mobility'] = 'N/A'
            summary_row['ont_position'] = 'N/A'
            summary_row['phoenix_contig'] = gene.get('phoenix_contig', 'N/A')
            summary_row['blast_identity'] = gene.get('blast_identity', 'N/A')
            summary_row['blast_coverage'] = gene.get('blast_coverage', 'N/A')
        
        summary_rows.append(summary_row)
    
    # Create DataFrame
    summary_df = pd.DataFrame(summary_rows)
    
    # Sort by category, then gene name
    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            by=['gene_category', 'gene_name'],
            ascending=[True, True]
        ).reset_index(drop=True)
    
    return summary_df

def main():
    args = parse_args()
    
    print(f"=" * 80)
    print(f"GRIM Gene Summary Generator")
    print(f"Sample: {args.sample_id}")
    print(f"=" * 80)
    
    # Load input files
    print("\n[1/4] Loading input files...")
    gene_mappings_df = load_gene_mappings(args.gene_mappings)
    classification_df = load_contig_classification(args.contig_classification)
    replicons_df = load_plasmid_replicons(args.plasmid_replicons)
    
    # Generate summary
    print("\n[2/4] Generating gene summary...")
    summary_df = generate_gene_summary(
        args.sample_id,
        gene_mappings_df,
        classification_df,
        replicons_df
    )
    
    # Define column order
    column_order = [
        'sample_id',
        'gene_name',
        'gene_category',
        'source',
        'mapping_status',
        'ont_contig_name',
        'ont_contig_type',
        'plasmid_name',
        'replicon_type',
        'relaxase_type',
        'predicted_mobility',
        'ont_position',
        'phoenix_contig',
        'blast_identity',
        'blast_coverage'
    ]
    
    # Reorder columns
    summary_df = summary_df[[col for col in column_order if col in summary_df.columns]]
    
    # Save results
    print("\n[3/4] Saving results...")
    output_file = args.output_summary
    summary_df.to_csv(output_file, index=False)
    print(f"  ✓ Saved gene summary to: {output_file}")
    
    # Print summary
    print(f"\n[4/4] Summary statistics:")
    print(f"=" * 80)
    print(f"  Total genes: {len(summary_df)}")
    
    if not summary_df.empty:
        # Mapping status
        status_counts = summary_df['mapping_status'].value_counts()
        print(f"\n  Mapping status:")
        for status, count in status_counts.items():
            print(f"    - {status}: {count}")
        
        # Mapped genes breakdown
        mapped_df = summary_df[summary_df['mapping_status'] == 'mapped']
        if not mapped_df.empty:
            print(f"\n  Mapped genes breakdown:")
            
            # By location type
            type_counts = mapped_df['ont_contig_type'].value_counts()
            print(f"    Location types:")
            for loc_type, count in type_counts.items():
                print(f"      - {loc_type}: {count}")
            
            # Plasmid genes
            plasmid_genes = mapped_df[mapped_df['ont_contig_type'] == 'plasmid']
            if not plasmid_genes.empty:
                print(f"\n    Plasmid-associated genes: {len(plasmid_genes)}")
                
                # By mobility
                mobility_counts = plasmid_genes['predicted_mobility'].value_counts()
                print(f"      Predicted mobility:")
                for mobility, count in mobility_counts.items():
                    if mobility != 'N/A':
                        print(f"        - {mobility}: {count}")
                
                # Top replicon types
                replicon_counts = plasmid_genes['replicon_type'].value_counts()
                print(f"      Top replicon types:")
                for replicon, count in list(replicon_counts.items())[:5]:
                    if replicon != 'N/A':
                        print(f"        - {replicon}: {count}")
            
            # Chromosomal genes
            chr_genes = mapped_df[mapped_df['ont_contig_type'] == 'chromosome']
            if not chr_genes.empty:
                print(f"\n    Chromosomal genes: {len(chr_genes)}")
            
            # By category
            category_counts = mapped_df['gene_category'].value_counts()
            print(f"\n    Top gene categories:")
            for category, count in list(category_counts.items())[:5]:
                print(f"      - {category}: {count}")
    
    print(f"\n✓ Pipeline complete!")

if __name__ == '__main__':
    main()
