#!/usr/bin/env python3
"""
Plasmid Classification using MOB-suite
Classifies ONT contigs (edge_*) as chromosome/plasmid using marker-based detection
and identifies replicon types, Inc groups, and plasmid names.
Uses: mob_recon with --run_typer flag
"""

import argparse
import pandas as pd
import sys
import os
import subprocess
from pathlib import Path
from Bio import SeqIO
import glob
import shutil

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Classify ONT contigs and identify plasmid replicons using MOB-suite'
    )
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--ont_genome', required=True, help='ONT complete genome FASTA file')
    parser.add_argument('--output_classification', required=True, 
                       help='Output contig classification TSV')
    parser.add_argument('--output_replicons', required=True, 
                       help='Output plasmid replicons TSV')
    parser.add_argument('--threads', type=int, default=4, 
                       help='Number of threads for MOB-suite (default: 4)')
    parser.add_argument('--outdir', default='mob_output',
                       help='Output directory for MOB-suite results (default: mob_output)')
    return parser.parse_args()

def check_mob_suite_installed():
    """Check if MOB-suite is installed and accessible"""
    try:
        result = subprocess.run(['mob_recon', '--version'], 
                              capture_output=True, text=True, check=False)
        if result.returncode == 0:
            print(f"  ✓ MOB-suite found: {result.stdout.strip() if result.stdout else result.stderr.strip()}")
            return True
        else:
            print(f"  ✗ MOB-suite not found or error: {result.stderr}", file=sys.stderr)
            return False
    except FileNotFoundError:
        print("  ✗ MOB-suite (mob_recon) not found in PATH", file=sys.stderr)
        return False

def run_mob_recon(ont_genome_file, outdir, threads=4):
    """
    Run MOB-suite mob_recon with --run_typer to classify contigs and identify replicons
    
    Args:
        ont_genome_file: Path to ONT genome FASTA
        outdir: Output directory for MOB-suite results
        threads: Number of threads
        
    Returns:
        Boolean indicating success
    """
    print(f"\n[MOB-suite] Running mob_recon with --run_typer...")
    print(f"  Input: {ont_genome_file}")
    print(f"  Output directory: {outdir}")
    print(f"  Threads: {threads}")
    
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # Build mob_recon command
    cmd = [
        'mob_recon',
        '--infile', ont_genome_file,
        '--outdir', outdir,
        '--num_threads', str(threads),
        '--run_typer',
        '--force'
    ]
    
    print(f"  Command: {' '.join(cmd)}")
    
    # Run mob_recon
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        if result.returncode == 0:
            print(f"  ✓ mob_recon completed successfully")
            if result.stdout:
                print(f"  {result.stdout}")
            return True
        else:
            print(f"  ✗ mob_recon failed with return code {result.returncode}", file=sys.stderr)
            print(f"  STDOUT: {result.stdout}", file=sys.stderr)
            print(f"  STDERR: {result.stderr}", file=sys.stderr)
            return False
            
    except Exception as e:
        print(f"  ✗ Error running mob_recon: {e}", file=sys.stderr)
        return False

def read_original_contigs(ont_genome_file):
    """
    Read original contig names and lengths from ONT genome FASTA
    
    Args:
        ont_genome_file: Path to ONT genome FASTA
        
    Returns:
        Dictionary mapping contig names to their lengths
    """
    original_contigs = {}
    
    try:
        for record in SeqIO.parse(ont_genome_file, "fasta"):
            # Handle contig names with spaces - take first word only
            contig_name = record.id.split()[0]
            original_contigs[contig_name] = len(record.seq)
    except Exception as e:
        print(f"Error reading ONT genome file: {e}", file=sys.stderr)
        return {}
    
    return original_contigs

def parse_contig_report(contig_report_path):
    """
    Parse mob_recon contig_report.txt to get classification
    
    Args:
        contig_report_path: Path to contig_report.txt
        
    Returns:
        Dictionary mapping contig names to classification info
    """
    contig_classifications = {}
    
    if not os.path.exists(contig_report_path):
        print(f"  ⚠ Warning: contig_report.txt not found at {contig_report_path}", file=sys.stderr)
        return contig_classifications
    
    try:
        # Read contig report
        df = pd.read_csv(contig_report_path, sep='\t')
        
        for _, row in df.iterrows():
            contig_id = row.get('contig_id', row.get('primary_cluster_id', 'unknown'))
            
            contig_classifications[contig_id] = {
                'molecule_type': row.get('molecule_type', 'unknown'),
                'rep_type': row.get('rep_type(s)', '-'),
                'relaxase_type': row.get('relaxase_type(s)', '-'),
                'size': row.get('size', 0),
                'gc_content': row.get('gc_content', 0),
                'circular': row.get('circular', 'unknown')
            }
    
    except Exception as e:
        print(f"  ⚠ Error parsing contig_report.txt: {e}", file=sys.stderr)
        return contig_classifications
    
    return contig_classifications

def parse_mobtyper_results(mobtyper_path):
    """
    Parse mobtyper_results.txt for replicon typing details
    
    Args:
        mobtyper_path: Path to mobtyper_results.txt
        
    Returns:
        Dictionary mapping contig/file names to replicon info
    """
    replicon_info = {}
    
    if not os.path.exists(mobtyper_path):
        print(f"  ⚠ Warning: mobtyper_results.txt not found at {mobtyper_path}", file=sys.stderr)
        return replicon_info
    
    try:
        df = pd.read_csv(mobtyper_path, sep='\t')
        
        for _, row in df.iterrows():
            # mob_typer uses different column names - try multiple
            file_id = row.get('file_id', row.get('sample_id', 'unknown'))
            
            replicon_info[file_id] = {
                'num_contigs': row.get('num_contigs', 1),
                'size': row.get('size', row.get('total_length', 0)),
                'gc': row.get('gc', 0),
                'md5': row.get('md5', 'unknown'),
                'rep_type': row.get('rep_type(s)', '-'),
                'rep_type_accession': row.get('rep_type_accession(s)', '-'),
                'relaxase_type': row.get('relaxase_type(s)', '-'),
                'relaxase_type_accession': row.get('relaxase_type_accession(s)', '-'),
                'mpf_type': row.get('mpf_type', '-'),
                'mpf_type_accession': row.get('mpf_type_accession(s)', '-'),
                'orit_type': row.get('orit_type(s)', '-'),
                'orit_accession': row.get('orit_accession(s)', '-'),
                'predicted_mobility': row.get('predicted_mobility', 'unknown'),
                'mash_nearest_neighbor': row.get('mash_nearest_neighbor', 'N/A'),
                'mash_neighbor_distance': row.get('mash_neighbor_distance', 'N/A'),
                'mash_neighbor_identification': row.get('mash_neighbor_identification', 'N/A'),
                'primary_cluster_id': row.get('primary_cluster_id', 'unknown'),
                'secondary_cluster_id': row.get('secondary_cluster_id', '-')
            }
    
    except Exception as e:
        print(f"  ⚠ Error parsing mobtyper_results.txt: {e}", file=sys.stderr)
        return replicon_info
    
    return replicon_info

def identify_plasmid_contigs_from_fasta(mob_outdir):
    """
    Identify which original contigs are plasmids by parsing mob_recon FASTA outputs
    
    Args:
        mob_outdir: MOB-suite output directory
        
    Returns:
        Dictionary mapping original contig names to their plasmid/chromosome classification
    """
    contig_to_type = {}
    
    # Parse chromosome.fasta (if exists)
    chr_fasta = os.path.join(mob_outdir, 'chromosome.fasta')
    if os.path.exists(chr_fasta):
        try:
            for record in SeqIO.parse(chr_fasta, 'fasta'):
                contig_name = record.id.split()[0]
                contig_to_type[contig_name] = 'chromosome'
        except Exception as e:
            print(f"  ⚠ Error parsing chromosome.fasta: {e}", file=sys.stderr)
    
    # Parse plasmid_*.fasta files
    plasmid_files = glob.glob(os.path.join(mob_outdir, 'plasmid_*.fasta'))
    for plasmid_file in plasmid_files:
        plasmid_id = os.path.basename(plasmid_file).replace('.fasta', '')
        try:
            for record in SeqIO.parse(plasmid_file, 'fasta'):
                contig_name = record.id.split()[0]
                contig_to_type[contig_name] = 'plasmid'
        except Exception as e:
            print(f"  ⚠ Error parsing {plasmid_file}: {e}", file=sys.stderr)
    
    return contig_to_type

def count_replicons(rep_type_str):
    """Count number of replicons in comma-separated string"""
    if not rep_type_str or rep_type_str == '-' or rep_type_str == 'unknown':
        return 0
    return len([r.strip() for r in rep_type_str.split(',') if r.strip()])

def generate_plasmid_name(primary_cluster, replicons, relaxase):
    """
    Generate plasmid name from MOB-suite typing results
    
    Priority:
    1. Use primary replicon + relaxase if available
    2. Use primary_cluster_id if no replicon
    3. Use "Unknown_plasmid" if nothing found
    
    Examples:
        IncFIB + MOBF → "IncFIB(MOBF)"
        IncI1 + MOBC → "IncI1(MOBC)"
        Multiple: IncFIB,IncFII + MOBF → "IncFIB/IncFII(MOBF)"
        No replicon + cluster AA087 → "MOB_AA087"
    """
    if replicons and replicons != '-' and replicons != 'unknown':
        # Parse replicons (may be comma-separated)
        rep_list = [r.strip() for r in replicons.split(',') if r.strip()]
        
        if relaxase and relaxase != '-' and relaxase != 'unknown':
            relax_list = [r.strip() for r in relaxase.split(',') if r.strip()]
            primary_relax = relax_list[0] if relax_list else ''
            
            if len(rep_list) == 1:
                return f"{rep_list[0]}({primary_relax})"
            else:
                rep_str = '/'.join(rep_list[:2])  # Limit to 2 for readability
                return f"{rep_str}({primary_relax})"
        else:
            return '/'.join(rep_list[:2])
    
    elif primary_cluster and primary_cluster != 'unknown' and primary_cluster != '-':
        return f"MOB_{primary_cluster}"
    
    else:
        return "Unknown_plasmid"

def calculate_confidence(classification, num_replicons):
    """
    Calculate confidence score for classification
    
    High confidence: Multiple replicons found + clear plasmid markers
    Medium confidence: Single replicon or chromosomal classification
    Low confidence: No clear markers (ambiguous)
    """
    if classification == 'plasmid' and num_replicons > 1:
        return 'high'
    elif classification == 'plasmid' and num_replicons == 1:
        return 'medium'
    elif classification == 'chromosome':
        return 'medium'
    else:
        return 'low'

def create_classification_and_replicon_tables(sample_id, original_contigs, mob_outdir):
    """
    Create classification and replicon TSV files from MOB-suite outputs
    
    Args:
        sample_id: Sample identifier
        original_contigs: Dictionary of original contig names and lengths
        mob_outdir: MOB-suite output directory
        
    Returns:
        Tuple of (classification_df, replicon_df)
    """
    print(f"\n[Parsing] Extracting results from MOB-suite outputs...")
    
    # Identify plasmid/chromosome contigs
    contig_to_type = identify_plasmid_contigs_from_fasta(mob_outdir)
    print(f"  Found {sum(1 for v in contig_to_type.values() if v == 'chromosome')} chromosome(s)")
    print(f"  Found {sum(1 for v in contig_to_type.values() if v == 'plasmid')} plasmid(s)")
    
    # Parse contig report
    contig_report_path = os.path.join(mob_outdir, 'contig_report.txt')
    contig_info = parse_contig_report(contig_report_path)
    
    # Parse mobtyper results
    mobtyper_path = os.path.join(mob_outdir, 'mobtyper_results.txt')
    replicon_info = parse_mobtyper_results(mobtyper_path)
    
    # Build classification table
    classification_rows = []
    replicon_rows = []
    
    for contig_name, contig_length in original_contigs.items():
        # Get classification
        classification = contig_to_type.get(contig_name, 'ambiguous')
        
        # Get additional info from contig report
        info = contig_info.get(contig_name, {})
        rep_type = info.get('rep_type', '-')
        relaxase_type = info.get('relaxase_type', '-')
        num_replicons = count_replicons(rep_type)
        
        # Determine circularity
        circular = 'yes' if classification == 'plasmid' else info.get('circular', 'unknown')
        
        # Calculate confidence
        confidence = calculate_confidence(classification, num_replicons)
        
        # Add to classification table
        classification_rows.append({
            'sample_id': sample_id,
            'original_contig_name': contig_name,
            'contig_length': contig_length,
            'classification': classification,
            'circular': circular,
            'confidence_score': confidence,
            'mob_recon_cluster_id': 'N/A',  # MOB-suite doesn't directly provide this per contig
            'num_replicons_found': num_replicons
        })
        
        # If plasmid, add to replicon table
        if classification == 'plasmid':
            # Try to find matching replicon info
            # mob_typer output might use plasmid_*.fasta as file_id
            plasmid_file_key = None
            for key in replicon_info.keys():
                if contig_name in key or key in contig_name:
                    plasmid_file_key = key
                    break
            
            # If not found by name matching, use first plasmid entry (for single-contig plasmids)
            if not plasmid_file_key and replicon_info:
                # Get plasmid file this contig belongs to
                for plasmid_file in glob.glob(os.path.join(mob_outdir, 'plasmid_*.fasta')):
                    try:
                        for record in SeqIO.parse(plasmid_file, 'fasta'):
                            if record.id.split()[0] == contig_name:
                                plasmid_file_key = os.path.basename(plasmid_file)
                                break
                    except:
                        continue
            
            rep_data = replicon_info.get(plasmid_file_key, {}) if plasmid_file_key else {}
            
            # Generate plasmid name
            plasmid_name = generate_plasmid_name(
                rep_data.get('primary_cluster_id', 'unknown'),
                rep_data.get('rep_type', rep_type),
                rep_data.get('relaxase_type', relaxase_type)
            )
            
            replicon_rows.append({
                'sample_id': sample_id,
                'original_contig_name': contig_name,
                'replicon_type': rep_data.get('rep_type', rep_type),
                'relaxase_type': rep_data.get('relaxase_type', relaxase_type),
                'mob_cluster_id': rep_data.get('primary_cluster_id', 'unknown'),
                'predicted_mobility': rep_data.get('predicted_mobility', 'unknown'),
                'mash_nearest_neighbor': rep_data.get('mash_nearest_neighbor', 'N/A'),
                'plasmid_name': plasmid_name,
                'num_replicons': count_replicons(rep_data.get('rep_type', rep_type))
            })
    
    # Convert to DataFrames
    classification_df = pd.DataFrame(classification_rows)
    replicon_df = pd.DataFrame(replicon_rows) if replicon_rows else pd.DataFrame(columns=[
        'sample_id', 'original_contig_name', 'replicon_type', 'relaxase_type',
        'mob_cluster_id', 'predicted_mobility', 'mash_nearest_neighbor',
        'plasmid_name', 'num_replicons'
    ])
    
    return classification_df, replicon_df

def main():
    args = parse_args()
    
    print(f"=" * 80)
    print(f"Plasmid Classification using MOB-suite")
    print(f"Sample: {args.sample_id}")
    print(f"=" * 80)
    
    # Check MOB-suite installation
    print("\n[1/5] Checking MOB-suite installation...")
    if not check_mob_suite_installed():
        print("\nERROR: MOB-suite is not installed or not in PATH")
        print("Please install MOB-suite: conda install -c conda-forge -c bioconda mob_suite")
        sys.exit(1)
    
    # Read original contigs
    print("\n[2/5] Reading original ONT genome...")
    original_contigs = read_original_contigs(args.ont_genome)
    print(f"  Found {len(original_contigs)} contigs in ONT genome")
    for contig_name, length in list(original_contigs.items())[:5]:
        print(f"    - {contig_name}: {length:,} bp")
    if len(original_contigs) > 5:
        print(f"    ... and {len(original_contigs) - 5} more")
    
    # Run mob_recon
    print("\n[3/5] Running MOB-suite mob_recon...")
    success = run_mob_recon(args.ont_genome, args.outdir, args.threads)
    
    if not success:
        print("\nERROR: mob_recon failed")
        sys.exit(1)
    
    # Parse results
    print("\n[4/5] Parsing MOB-suite results...")
    classification_df, replicon_df = create_classification_and_replicon_tables(
        args.sample_id,
        original_contigs,
        args.outdir
    )
    
    # Save results
    print("\n[5/5] Saving results...")
    classification_df.to_csv(args.output_classification, sep='\t', index=False)
    print(f"  ✓ Saved classification to: {args.output_classification}")
    
    replicon_df.to_csv(args.output_replicons, sep='\t', index=False)
    print(f"  ✓ Saved replicons to: {args.output_replicons}")
    
    # Print summary
    print(f"\n" + "=" * 80)
    print(f"Summary for sample {args.sample_id}:")
    print(f"=" * 80)
    print(f"  Total contigs: {len(classification_df)}")
    
    if not classification_df.empty:
        class_counts = classification_df['classification'].value_counts()
        for classification, count in class_counts.items():
            print(f"  {classification}: {count}")
        
        # Plasmid details
        plasmid_df = classification_df[classification_df['classification'] == 'plasmid']
        if not plasmid_df.empty:
            print(f"\nPlasmid details:")
            print(f"  Total plasmids: {len(plasmid_df)}")
            print(f"  Circular: {plasmid_df[plasmid_df['circular'] == 'yes'].shape[0]}")
            
            if not replicon_df.empty:
                print(f"  With replicons: {len(replicon_df)}")
                
                # Replicon types
                rep_types = replicon_df['replicon_type'].value_counts()
                print(f"\n  Top replicon types:")
                for rep_type, count in rep_types.head(5).items():
                    print(f"    - {rep_type}: {count}")
                
                # Mobility
                mobility_counts = replicon_df['predicted_mobility'].value_counts()
                print(f"\n  Predicted mobility:")
                for mobility, count in mobility_counts.items():
                    print(f"    - {mobility}: {count}")
    
    print(f"\n✓ Pipeline complete!")

if __name__ == '__main__':
    main()
