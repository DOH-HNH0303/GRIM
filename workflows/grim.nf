/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ILLUMINA_AMR_LOCATOR } from '../modules/local/illumina_amr_locator'
include { MOBSUITE_TYPER } from '../modules/local/mobsuite_typer'
include { PLATON } from '../modules/local/platon'
include { BANDAGE } from '../modules/local/bandage'
include { GFATOOLS_COUNT_EDGES } from '../modules/local/gfatools_count_edges'
include { IDENTIFY_NON_MOBILIZABLE_CONTIGS } from '../modules/local/identify_non_mobilizable_contigs'
include { EXTRACT_NON_MOBILIZABLE_CONTIGS } from '../modules/local/extract_non_mobilizable_contigs'
include { GRIM_GENE_SUMMARY } from '../modules/local/grim_gene_summary'
include { RESOLVE_PHOENIX_FILES } from '../modules/local/resolve_phoenix_files'
include { MULTIQC             } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc} from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_GRIM_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW - PHOENIX STYLE APPROACH (Direct S3 URL Handling)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GRIM {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Parse samplesheet and separate into two formats
    // Format 1: sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome[,hybrid_assembly_gfa]
    //           - hybrid_assembly_gfa (optional): GFA file from hybrid assembly (Illumina + ONT) for Bandage visualization
    // Format 2: sample,phoenix_outdir,ont_complete_genome (no GFA - Bandage will not run)
    //
    ch_samplesheet
        .branch { row ->
            individual_files: row[0] == 'individual_files'
                // Format: [meta, gamma, amrfinder, assembly, ont, gfa (optional)]
                return [row[1], row[2], row[3], row[4], row[5], row[6]]
            phoenix_outdir: row[0] == 'phoenix_outdir'
                // Format: [meta, phoenix_dir, ont, gfa (optional)]
                return [row[1], row[2], row[3], row[4]]
            invalid: true
                error "Invalid samplesheet format. Could not determine format from row: ${row}"
        }
        .set { ch_input_formats }

    //
    // For phoenix_outdir format, resolve file paths
    //
    RESOLVE_PHOENIX_FILES (
        ch_input_formats.phoenix_outdir
    )
    ch_versions = ch_versions.mix(RESOLVE_PHOENIX_FILES.out.versions)

    //
    // Combine both input formats into a single channel
    // Both formats now have 6 elements: [meta, gamma, amr, assembly, ont, gfa]
    // GFA may be null for samples without GFA files
    //
    ch_phoenix_files = ch_input_formats.individual_files
        .mix(RESOLVE_PHOENIX_FILES.out.resolved_files)

    //
    // MODULE: Process each sample using existing Phoenix AMR indexing files
    // This leverages the pre-computed GAMMA and AMRFinder results
    // No need to re-parse Phoenix summary files or re-run BLAST!
    // Drop the GFA element (6th element) before passing to ILLUMINA_AMR_LOCATOR
    //
    ch_illumina_amr_input = ch_phoenix_files
        .map { meta, gamma, amr, assembly, ont, gfa ->
            tuple(meta, gamma, amr, assembly, ont)
        }
    
    ILLUMINA_AMR_LOCATOR (
        ch_illumina_amr_input
    )
    ch_versions = ch_versions.mix(ILLUMINA_AMR_LOCATOR.out.versions.first())

    //
    // MODULE: Classify plasmids and chromosomes using MOB-suite recon
    // Extract ONT genome from phoenix_files channel
    //
    ch_ont_genomes = ch_phoenix_files.map { meta, _gamma, _amrfinder, _phoenix_asm, ont, _gfa ->
        tuple(meta, ont)
    }
    
    //
    // Extract GFA files for Bandage (only samples with non-empty GFA files)
    //
    ch_gfa_files = ch_phoenix_files
        .map { meta, _gamma, _amrfinder, _phoenix_asm, _ont, gfa ->
            tuple(meta, gfa)
        }
        .filter { meta, gfa -> 
            // Filter out null and empty files (created by touch in RESOLVE_PHOENIX_FILES)
            gfa != null && gfa.size() > 0
        }
    
    //
    // MODULE: Type plasmids and MGEs using MOB-suite typer
    //
    MOBSUITE_TYPER (
        ch_ont_genomes
    )
    ch_versions = ch_versions.mix(MOBSUITE_TYPER.out.versions.first())

    //
    // MODULE: Run Platon on ONT genomes to classify plasmids and chromosomes
    //
    // Prepare platon database channel - collect all files from the database directory
    ch_platon_db = params.platon_db 
        ? Channel.fromPath("${params.platon_db}/**", checkIfExists: false).collect()
        : Channel.value(file('NO_DB_FILE'))
    
    PLATON (
        ch_ont_genomes,
        ch_platon_db
    )
    ch_versions = ch_versions.mix(PLATON.out.versions.first())

    //
    // MODULE: Run Bandage on GFA files (optional - only for samples with GFA files)
    //
    BANDAGE (
        ch_gfa_files
    )
    ch_versions = ch_versions.mix(BANDAGE.out.versions.first())


    //
    // MODULE: Count number of edges in assembly with gfatools
    //
       GFATOOLS_COUNT_EDGES(
        ch_gfa_files  // Channel with [meta, gfa_file]
    )

    //
    // MODULE: Identify non-mobilizable contigs from MOB-typer and Platon results
    //
    ch_mobtyper_for_identification = MOBSUITE_TYPER.out.mobtyper_results
    .join(PLATON.out.tsv, remainder: true)
    .join(GFATOOLS_COUNT_EDGES.out.edge_counts, remainder: true)
    


    IDENTIFY_NON_MOBILIZABLE_CONTIGS (
        ch_mobtyper_for_identification
    )

    //
    // MODULE: Extract non-mobilizable contigs from ONT genome FASTA
    //
    ch_extract_input = ch_ont_genomes
        .join(IDENTIFY_NON_MOBILIZABLE_CONTIGS.out.non_mobilizable_list)
    
    EXTRACT_NON_MOBILIZABLE_CONTIGS (
        ch_extract_input
    )

    //
    // MODULE: Generate per-gene summary CSV with location and plasmid information
    // Combine outputs from ILLUMINA_AMR_LOCATOR, MOBSUITE_TYPER, and PLATON
    //
    ch_summary_input = ILLUMINA_AMR_LOCATOR.out.mappings
        .join(MOBSUITE_TYPER.out.mobtyper_results)
        .join(PLATON.out.tsv)
    
    GRIM_GENE_SUMMARY (
        ch_summary_input
    )
    ch_versions = ch_versions.mix(GRIM_GENE_SUMMARY.out.versions.first())

    //
    // Collect results for MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA_AMR_LOCATOR.out.mappings.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(GRIM_GENE_SUMMARY.out.summary.collect{it[1]})

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'grim_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    gene_mappings         = ILLUMINA_AMR_LOCATOR.out.mappings        // channel: [ meta, gene_mappings.tsv ]
    unmapped_genes        = ILLUMINA_AMR_LOCATOR.out.unmapped        // channel: [ meta, unmapped_genes.txt ]
    mobtyper_results      = MOBSUITE_TYPER.out.mobtyper_results      // channel: [ meta, mobtyper_results.tsv ]
    gene_summary          = GRIM_GENE_SUMMARY.out.summary           // channel: [ meta, gene_summary.csv ]
    multiqc_report        = MULTIQC.out.report.toList()             // channel: /path/to/multiqc_report.html
    versions              = ch_versions                              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/