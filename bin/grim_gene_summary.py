#!/usr/bin/env python3
"""
GRIM Gene Summary - Optimized Version
Generates per-gene CSV with location and plasmid information.
Combines gene mappings and MOB-typer results using vectorized pandas operations.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path



def id_cleanup(col):
    """Cleans up columns with values """
    col = col.split(" ")[0].replace(">", "")
    return col


def ge_type(col):
    # (Conjugative, Mobilizable, Non-mobilizable)
    if col == "conjugative":
        return "plasmid"
    elif col == "mobilizable":
        return "plasmid"
    elif col == "non-mobilizable":
        # Check to see if GE is chromosome-like
        pass
    else:
        return "ambiguous"


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate per-gene summary CSV for GRIM pipeline'
    )
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--gene_mappings', required=True, 
                       help='Gene mappings TSV from illumina_amr_locator')
    parser.add_argument('--mobtyper_results', required=True,
                       help='MOB-typer results TSV from mobsuite_typer')
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


def load_mobtyper_results(mobtyper_file):
    """Load MOB-typer results from TSV"""
    try:
        if Path(mobtyper_file).stat().st_size == 0:
            print(f"  No MOB-typer results found (empty file)")
            return pd.DataFrame()
        
        df = pd.read_csv(mobtyper_file, sep='\t')
        df['sample_id'] = df['sample_id'].apply(lambda x: id_cleanup(x))

        print(f"  Loaded {len(df)} MOB-typer results")
        return df
    except Exception as e:
        print(f"Error loading MOB-typer results: {e}", file=sys.stderr)
        # Return empty DataFrame instead of exiting
        return pd.DataFrame()


def generate_gene_summary(sample_id, gene_mappings_df, mobtyper_df):
    print(f"\n[Processing] Generating gene summary...")

    summary_df = gene_mappings_df.copy()
    plasmid_df = mobtyper_df.copy()

    # Add sample_id column to summary
    summary_df['sample_id'] = sample_id

    # Rename MOB-typer sample_id to match gene mappings
    plasmid_df.rename(columns={'sample_id': 'ont_contig'}, inplace=True)

    mobtyper_to_keep = [
        'ont_contig',
        'predicted_mobility',
        'mash_nearest_neighbor',
        'mash_neighbor_distance',
        'mash_neighbor_identification',
        'primary_cluster_id',
        'secondary_cluster_id',
        'predicted_host_range_overall_rank',
        'predicted_host_range_overall_name'
    ]

    # Merge on ont_contig
    summary_df = summary_df.merge(
        plasmid_df[mobtyper_to_keep],
        on='ont_contig',
        how='left'
    )


    print(summary_df)


    

    return summary_df



def main():
    args = parse_args()
    
    print(f"=" * 80)
    print(f"GRIM Gene Summary Generator (Optimized)")
    print(f"Sample: {args.sample_id}")
    print(f"=" * 80)
    
    # Load input files
    print("\n[1/4] Loading input files...")
    gene_mappings_df = load_gene_mappings(args.gene_mappings)
    mobtyper_df = load_mobtyper_results(args.mobtyper_results)
    
    # Generate summary
    print("\n[2/4] Generating gene summary...")
    summary_df = generate_gene_summary(
        args.sample_id,
        gene_mappings_df,
        mobtyper_df
    )
    
    # Define column order
    column_order = [
        'sample_id',
        'gene_name',
        'gene_category',
        'source',
        'ont_contig_name',
        'ont_contig_type',
        'contig_size',
        'gc_content',
        'replicon_type',
        'relaxase_type',
        'mpf_type',
        'predicted_mobility',
        'primary_cluster_id',
        'secondary_cluster_id',
        'mash_nearest_neighbor',
        'mash_neighbor_distance',
        'predicted_host_range_overall_name',
        'ont_position',
        'blast_identity',
        'blast_coverage'
    ]
    
    # Reorder columns
    summary_df.rename(columns={'blast_identity': 'illumina_mapping_blast_identity', 
                               'blast_coverage': 'illumina_mapping_blast_coverage'}, inplace=True)
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
    
    ##################################################
    #### Add in mapped info if wanted, otherwise, redundant
    ##################################################
    # if not summary_df.empty:
    #     # Mapping status
    #     status_counts = summary_df['mapping_status'].value_counts()
    #     print(f"\n  Mapping status:")
    #     for status, count in status_counts.items():
    #         print(f"    - {status}: {count}")
        
    #     # Mapped genes breakdown
    #     mapped_df = summary_df[summary_df['mapping_status'] == 'mapped']
    #     if not mapped_df.empty:
    #         print(f"\n  Mapped genes breakdown:")
            
    #         # By location type
    #         type_counts = mapped_df['ont_contig_type'].value_counts()
    #         print(f"    Location types:")
    #         for loc_type, count in type_counts.items():
    #             print(f"      - {loc_type}: {count}")
            
    #         # Plasmid genes
    #         plasmid_genes = mapped_df[mapped_df['ont_contig_type'] == 'plasmid']
    #         if not plasmid_genes.empty:
    #             print(f"\n    Plasmid-associated genes: {len(plasmid_genes)}")
                
    #             # By mobility
    #             mobility_counts = plasmid_genes['predicted_mobility'].value_counts()
    #             print(f"      Predicted mobility:")
    #             for mobility, count in mobility_counts.items():
    #                 if mobility != 'N/A':
    #                     print(f"        - {mobility}: {count}")
                
    #             # Top replicon types
    #             replicon_counts = plasmid_genes['replicon_type'].value_counts()
    #             print(f"      Top replicon types:")
    #             for replicon, count in list(replicon_counts.items())[:5]:
    #                 if replicon != 'N/A':
    #                     print(f"        - {replicon}: {count}")
            
    #         # Chromosomal genes
    #         chr_genes = mapped_df[mapped_df['ont_contig_type'] == 'chromosome']
    #         if not chr_genes.empty:
    #             print(f"\n    Chromosomal genes: {len(chr_genes)}")
            
    #         # By category
    #         category_counts = mapped_df['gene_category'].value_counts()
    #         print(f"\n    Top gene categories:")
    #         for category, count in list(category_counts.items())[:5]:
    #             print(f"      - {category}: {count}")
    
    print(f"\n✓ Pipeline complete!")


if __name__ == '__main__':
    main()
