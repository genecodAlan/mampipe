#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { busco_precheck; check_busco_dup; read_dup_score; purge_duplicates; scaffold_assembly; busco_final } from './Post_Assembly.nf'

// Default parameters for local execution
params.reads = 'data/*.fastq'
params.ref = 'refSeq/*.fna'
params.outdir = 'WGS_results'
params.busco_threshold = 6.0
params.lineage = 'mammalia_odb10'

// Thread parameters
params.threads = 8
params.threads_assembly = 4
params.threads_polish = 4

// Medaka model parameter (adjust based on your sequencing chemistry)
params.medaka_model = 'r1041_e82_400bps_fast_g632'

// For Google Cloud Batch, override with:
// params.reads = 'gs://YOUR_BUCKET/data/*.fastq'
// params.ref = 'gs://YOUR_BUCKET/refSeq/*.fna'
// params.outdir = 'gs://YOUR_BUCKET/results'

workflow {
    // Validate inputs exist
    if (!file(params.reads).exists()) {
        error "ERROR: Reads not found at ${params.reads}"
    }
    if (!file(params.ref).exists()) {
        error "ERROR: Reference genome not found at ${params.ref}"
    }
    
    reads = Channel.fromPath(params.reads)
    
    // Quality control and filtering
    pre_chopped_results = plot_reads(reads)
    chopped_reads = chopper(reads)
    chopped_plot = chopped_plot_reads(chopped_reads)
    
    // Assembly and polishing
    assembled_reads = assemble_long_reads(chopped_reads)
    polished_assembly = polish_assembly(chopped_reads, assembled_reads)
    reference_genome = Channel.fromPath(params.ref)
    
    // Run BUSCO on initial assembly
    busco_results_ch = busco_precheck(polished_assembly)
    
    // Extract duplication score
    dup_score_file = check_busco_dup(busco_results_ch)
    
    // Read the actual score value from the file
    dup_score_value = read_dup_score(dup_score_file)
    
    // Create combined channel with assembly, reads, and score for decision making
    decision_ch = polished_assembly
        .combine(chopped_reads)
        .combine(dup_score_value)
        .map { assembly, reads, score -> tuple(assembly, reads, score) }
    
    // Split into two paths based on threshold
    decision_ch
        .branch {
            high_dup: it[2].toDouble() > params.busco_threshold
            low_dup: it[2].toDouble() <= params.busco_threshold
        }
        .set { branched_ch }
    
    // Path 1: High duplicates - run purge_duplicates
    purged_assembly_ch = purge_duplicates(branched_ch.high_dup)
    
    // Path 2: Low duplicates - use original assembly
    original_assembly_ch = branched_ch.low_dup
        .map { assembly, reads, score -> assembly }
    
    // Combine both possible final assemblies
    final_assembly_ch = purged_assembly_ch.mix(original_assembly_ch)
    
    // Scaffold final assembly
    scaffold_assembly(final_assembly_ch, reference_genome)
}

process chopper{
    label 'low_mem'
    conda './envs/Chopper.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/chopped", mode: 'copy'

    input:
    path sample

    output:
    path "${sample.baseName}_chopped.fastq"

    script:
    """
    chopper -q 10 -l 500 -i $sample > ${sample.baseName}_chopped.fastq"
    """
}

process plot_reads{
    label 'low_mem'
    conda './envs/Nanoplot.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/read_plots", mode: 'copy'

    input:
    path sample
    
    output:
    path "read_length_distribution"
    
    script:
    """
    NanoPlot --fastq $sample --plots dot -o read_length_distribution    
    """
}

process chopped_plot_reads{
    label 'low_mem'
    conda './envs/Nanoplot.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/read_plots", mode: 'copy'

    input:
    path sample
    
    output:
    path "chopped_read_length_distribution"
    
    script:
    """
    NanoPlot --fastq $sample --plots dot -o chopped_read_length_distribution    
    """
}

process assemble_long_reads {
    label 'high_mem'
    conda './envs/Flye.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path sample

    output:
    path "assembly"

    script:
    """
    flye --nano-raw $sample --out-dir assembly --threads ${params.threads_assembly}
    """
}

process polish_assembly {
    label 'medium_mem'
    conda './envs/Medaka.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/polished", mode: 'copy'

    input:
    path sample
    path assembled_fasta

    output:
    path "polished_assembly"

    script:
    """
    medaka_consensus -i $sample -d $assembled_fasta -o polished_assembly -t ${params.threads_polish} -m ${params.medaka_model}
    """
}
