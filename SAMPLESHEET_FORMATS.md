# GRIM Pipeline Samplesheet Formats

The GRIM pipeline now supports two flexible input formats to accommodate different use cases.

## Format 1: Individual File Paths (Original)

**Use when:** You want precise control over which files are used, or when Phoenix files are in non-standard locations.

**Columns:** `sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome`

**Example:**
```csv
sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome
ECOLI_001,/data/phoenix/ECOLI_001/ECOLI_001_ResGANNCBI_20230417_srst2.gamma,/data/phoenix/ECOLI_001/ECOLI_001_amrfinder_report.tsv,/data/phoenix/ECOLI_001/ECOLI_001_scaffolds.fasta,/data/ont/ECOLI_001/ECOLI_001_complete_genome.fasta
```

## Format 2: Phoenix Output Directory (New, Recommended)

**Use when:** You have standard Phoenix output directories and want simplified input preparation.

**Columns:** `sample,phoenix_outdir,ont_complete_genome`

**Example:**
```csv
sample,phoenix_outdir,ont_complete_genome
ECOLI_001,/data/phoenix/output/ECOLI_001,/data/ont/ECOLI_001/ECOLI_001_complete_genome.fasta
```

## How Format 2 Works

When you use Format 2, GRIM automatically looks for these files in the Phoenix output directory:

1. **GAMMA file**: `{phoenix_outdir}/{sample}/{sample}_ResGANNCBI_*_srst2.gamma`
2. **AMRFinder report**: `{phoenix_outdir}/{sample}/{sample}_amrfinder_report.tsv`
3. **Phoenix assembly**: `{phoenix_outdir}/{sample}/{sample}_scaffolds.fasta`

The `phoenix_outdir` should be the same directory that was specified as `--outdir` when running the Phoenix pipeline.

## Expected Phoenix Directory Structure

```
phoenix_outdir/
└── SAMPLE_NAME/
    ├── SAMPLE_NAME_ResGANNCBI_YYYYMMDD_srst2.gamma  # Date varies
    ├── SAMPLE_NAME_amrfinder_report.tsv
    └── SAMPLE_NAME_scaffolds.fasta
```

## Benefits of Format 2

- **Simpler**: Only 3 columns instead of 5
- **Less error-prone**: Automatic file discovery reduces typos
- **Flexible**: Handles varying GAMMA file dates automatically
- **Consistent**: Works with standard Phoenix output structure

## Pipeline Detection

The pipeline automatically detects which format you're using based on the number of columns:
- **5 columns**: Individual file paths format
- **3 columns**: Phoenix output directory format

## Error Handling

If using Format 2 and required files are not found, the pipeline will:
1. Report exactly which files are missing
2. Show the expected file paths
3. Exit with a clear error message

This helps you quickly identify and fix any issues with file locations or naming.