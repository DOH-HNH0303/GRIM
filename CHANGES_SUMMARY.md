# GRIM Pipeline Changes Summary

## New Feature: Dual Samplesheet Format Support

The GRIM pipeline now supports two input formats for maximum flexibility and ease of use.

### Changes Made

#### 1. Updated Input Schema (`assets/schema_input.json`)
- Modified to support both input formats using JSON Schema `oneOf`
- Format 1: Individual file paths (5 columns)
- Format 2: Phoenix output directory (3 columns)
- Automatic validation for both formats

#### 2. New Process Module (`modules/local/resolve_phoenix_files.nf`)
- Automatically discovers Phoenix files from output directory
- Handles wildcard patterns for GAMMA files (date variations)
- Provides clear error messages if files are not found
- Creates symbolic links to maintain file integrity

#### 3. Enhanced Workflow Logic (`workflows/grim.nf`)
- Added `RESOLVE_PHOENIX_FILES` import
- Implemented channel branching to separate input formats
- Automatic format detection based on column count
- Seamless integration of both formats into single processing channel

#### 4. Updated Documentation
- Enhanced README with detailed samplesheet format explanations
- Created comprehensive examples for both formats
- Added tips and recommendations for format selection
- Clear column descriptions and requirements

#### 5. Example Files Created
- `assets/samplesheet_phoenix_outdir.csv`: Phoenix directory format example
- `example_individual_files.csv`: Individual files format example
- `example_phoenix_outdir.csv`: Phoenix directory format example
- `SAMPLESHEET_FORMATS.md`: Detailed format documentation
- `example_samplesheets.md`: Comprehensive usage examples

### Input Format Details

#### Format 1: Individual File Paths
```csv
sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome
SAMPLE_1,/path/to/gamma.gamma,/path/to/amrfinder.tsv,/path/to/assembly.fasta,/path/to/ont.fasta
```

#### Format 2: Phoenix Output Directory (Recommended)
```csv
sample,phoenix_outdir,ont_complete_genome
SAMPLE_1,/path/to/phoenix/output/SAMPLE_1,/path/to/ont.fasta
```

### Benefits

1. **Simplified Input**: Format 2 reduces columns from 5 to 3
2. **Error Reduction**: Automatic file discovery prevents path typos
3. **Flexibility**: Handles varying GAMMA file dates automatically
4. **Backward Compatibility**: Original format still fully supported
5. **Clear Error Messages**: Helpful feedback when files are missing

### File Discovery Logic

For Format 2, the pipeline automatically finds:
- GAMMA file: `{phoenix_outdir}/{sample}/{sample}_ResGANNCBI_*_srst2.gamma`
- AMRFinder report: `{phoenix_outdir}/{sample}/{sample}_amrfinder_report.tsv`
- Phoenix assembly: `{phoenix_outdir}/{sample}/{sample}_scaffolds.fasta`

### Usage Examples

Both formats work with the same command:

```bash
# Using Format 1 (individual files)
nextflow run nf-core/grim -profile docker --input samplesheet_individual.csv --outdir results

# Using Format 2 (phoenix outdir)
nextflow run nf-core/grim -profile docker --input samplesheet_phoenix_outdir.csv --outdir results
```

### Technical Implementation

- **Automatic Detection**: Pipeline detects format based on column count
- **Channel Branching**: Separates formats into different processing paths
- **File Resolution**: New process handles Phoenix directory structure
- **Channel Mixing**: Combines both formats into unified processing stream
- **Error Handling**: Clear messages for missing files or invalid formats

This enhancement makes the GRIM pipeline much more user-friendly while maintaining full backward compatibility with existing samplesheets.