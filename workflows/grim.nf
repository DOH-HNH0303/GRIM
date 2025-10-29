/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PHOENIX_AMR_LOCATOR } from '../modules/local/phoenix_amr_locator'
include { RESOLVE_PHOENIX_FILES } from '../modules/local/resolve_phoenix_files'
include { MULTIQC             } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc} from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_grim_pipeline'

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
    // Format 1: sample,gamma_ar_file,amrfinder_report,phoenix_assembly_fasta,ont_complete_genome
    // Format 2: sample,phoenix_outdir,ont_complete_genome
    //
    ch_samplesheet//.brach_samplesheet
        .branch { row ->
            individual_files: row[0] == 'individual_files'
                return [row[1], row[2], row[3], row[4], row[5]]
            phoenix_outdir: row[0] == 'phoenix_outdir'
                return [row[1], row[2], row[3]]  // ← Correct indices!
                // Format 2: [meta, phoenix_outdir, ont_complete_genome]
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
    //
    ch_phoenix_files = ch_input_formats.individual_files
        .mix(RESOLVE_PHOENIX_FILES.out.resolved_files)

    //
    // MODULE: Process each sample using existing Phoenix AMR indexing files
    // This leverages the pre-computed GAMMA and AMRFinder results
    // No need to re-parse Phoenix summary files or re-run BLAST!
    //
    PHOENIX_AMR_LOCATOR (
        ch_phoenix_files
    )
    ch_versions = ch_versions.mix(PHOENIX_AMR_LOCATOR.out.versions.first())

    //
    // Collect results for MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(PHOENIX_AMR_LOCATOR.out.locations.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(PHOENIX_AMR_LOCATOR.out.detailed_amr.collect{it[1]})

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
    gene_locations = PHOENIX_AMR_LOCATOR.out.locations    // channel: [ meta, gene_locations.tsv ]
    detailed_amr   = PHOENIX_AMR_LOCATOR.out.detailed_amr // channel: [ meta, detailed_amr.tsv ]
    multiqc_report = MULTIQC.out.report.toList()         // channel: /path/to/multiqc_report.html
    versions       = ch_versions                          // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/