# GRIM Pipeline Samplesheet Examples

The GRIM pipeline supports two input formats for maximum flexibility.

## Format 1: Individual File Paths

Use this format when you want to specify exact paths to each Phoenix output file:

```csv
sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome
ECOLI_001,/data/phoenix/ECOLI_001/ECOLI_001_ResGANNCBI_20230417_srst2.gamma,/data/phoenix/ECOLI_001/ECOLI_001_amrfinder_report.tsv,/data/phoenix/ECOLI_001/ECOLI_001_scaffolds.fasta,/data/ont/ECOLI_001/ECOLI_001_complete_genome.fasta
KPNEU_002,/data/phoenix/KPNEU_002/KPNEU_002_ResGANNCBI_20230417_srst2.gamma,/data/phoenix/KPNEU_002/KPNEU_002_amrfinder_report.tsv,/data/phoenix/KPNEU_002/KPNEU_002_scaffolds.fasta,/data/ont/KPNEU_002/KPNEU_002_complete_genome.fasta
STAPH_003,/data/phoenix/STAPH_003/STAPH_003_ResGANNCBI_20230417_srst2.gamma,/data/phoenix/STAPH_003/STAPH_003_amrfinder_report.tsv,/data/phoenix/STAPH_003/STAPH_003_scaffolds.fasta,/data/ont/STAPH_003/STAPH_003_complete_genome.fasta
```

## Format 2: Phoenix Output Directory (Recommended)

Use this format when you have Phoenix output directories and want GRIM to automatically find the required files:

```csv
sample,phoenix_outdir,ont_complete_genome
ECOLI_001,/data/phoenix/output/ECOLI_001,/data/ont/ECOLI_001/ECOLI_001_complete_genome.fasta
KPNEU_002,/data/phoenix/output/KPNEU_002,/data/ont/KPNEU_002/KPNEU_002_complete_genome.fasta
STAPH_003,/data/phoenix/output/STAPH_003,/data/ont/STAPH_003/STAPH_003_complete_genome.fasta
```

## Expected Phoenix Directory Structure

When using Format 2, GRIM expects the following Phoenix directory structure:

```
phoenix_outdir/
└── SAMPLE_NAME/
    ├── SAMPLE_NAME_ResGANNCBI_YYYYMMDD_srst2.gamma
    ├── SAMPLE_NAME_amrfinder_report.tsv
    └── SAMPLE_NAME_scaffolds.fasta
```

For example:
```
/data/phoenix/output/ECOLI_001/
├── ECOLI_001_ResGANNCBI_20230417_srst2.gamma
├── ECOLI_001_amrfinder_report.tsv
└── ECOLI_001_scaffolds.fasta
```

## File Requirements

- **GAMMA files**: Must have `.gamma` extension
- **AMRFinder reports**: Must have `.tsv` or `.txt` extension
- **Assembly files**: Must have `.fasta`, `.fa`, or `.fna` extension
- **ONT genomes**: Must have `.fasta`, `.fa`, or `.fna` extension (optionally gzipped with `.gz`)
- **Sample names**: Must not contain spaces

## Usage Examples

### Using Format 1:
```bash
nextflow run nf-core/grim \
    -profile docker \
    --input samplesheet_individual_files.csv \
    --outdir results
```

### Using Format 2:
```bash
nextflow run nf-core/grim \
    -profile docker \
    --input samplesheet_phoenix_outdir.csv \
    --outdir results
```