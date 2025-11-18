#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { busco_precheck; check_busco_dup; read_dup_score; purge_duplicates; scaffold_assembly; busco_final } from './Post_Assembly.nf'

// Default parameters for local execution
params.reads = 'data/*.fastq'
params.ref = 'refSeq/*.fna'
params.outdir = 'WGS_results'
params.busco_threshold = 6.0
params.lineage = 'mammalia_odb10'

// For Google Cloud Batch, override with:
// params.reads = 'gs://YOUR_BUCKET/data/*.fastq'
// params.ref = 'gs://YOUR_BUCKET/refSeq/*.fna'
// params.outdir = 'gs://YOUR_BUCKET/results'

workflow {
    reads = Channel.fromPath(params.reads)
    
    //List of all methods = Plot, Chopper, Filter, Assemble, Polish, Purge, Scaffold, Assess Quality, Annotate Genome, Variant Calling
    //List of all tools = Nanoplot, Chopper, flye, medaka, Purge Dups, juicer, (Quast,BUSCO), Braker, (gatk4,bcftools,bwa)
    // Step 1: Run fastqc on raw reads
    

    pre_chopped_results = plot_reads(reads)
    chopped_reads = chopper(reads)
    chopped_plot = chopped_plot_reads(chopped_reads)
    assembled_reads = assemble_long_reads(chopped_reads)
    polished_assembly = polish_assembly(chopped_reads, assembled_reads)
    reference_genome = Channel.fromPath(params.ref)

    //polished_assembly = Channel.fromPath("results/polished/polished_assembly/consensus.fasta")
    //chopped_reads = Channel.fromPath("results/chopped/*_chopped.fastq")
    
    // Run BUSCO on initial assembly
    busco_results_ch = busco_precheck(polished_assembly)
    
    // Extract duplication score
    dup_score_file = check_busco_dup(busco_results_ch)  // emits path to score file
    
    // Read the actual score value from the file
    dup_score_value = read_dup_score(dup_score_file)  // emits val dup_score
    
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

/*
    final_score = read_dup_score(dup_score)

    // Run purge_dups only if dup % exceeds threshold
    purged_assembly_ch = purge_duplicates(polished_assembly, chopped_reads, final_score)

    // Choose final assembly based on duplication score
    final_assembly_ch = purged_assembly_ch.ifEmpty { polished_assembly }
    // Proceed with scaffolding
    scaffolded_assembly_ch = scaffold_assembly(final_assembly_ch, reference_genome)
    //get the only fasta file from the scaffolded assembly channel, dont use first since the fasta file is not guaranteed to be the first
    scaffolded = scaffolded_assembly_ch.flatMap { dir -> file("$dir").listFiles().findAll { it.name.endsWith('.fasta') }}
    busco_final_results = busco_final(scaffolded)


}*/

/*
process fastqc_trimmed {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'

    input:
    path sample

    output:
    path "*_fastqc.{html,zip}"

    

    script:
    sample_id = sample.getBaseName()
    """
    fastqc $sample
    """ 
}

process fastqc_filt {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_filt", mode: 'copy'

    input:
    path sample

    output:
    path "*_fastqc.{html,zip}"

    script:
    sample_id = sample.getBaseName()
    """
    fastqc $sample
    """ 
}

process fastqc_raw {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path sample

    output:
    path "*_fastqc.{html,zip}"

    script:
    sample_id = sample.getBaseName()
    """
    fastqc $sample
    """
}

process multiqc_raw {
    publishDir "${params.outdir}/multiqc_raw", mode: 'copy'

    input:
    path fastqc_results

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

process multiqc_trimmed {
    publishDir "${params.outdir}/multiqc_trimmed", mode: 'copy'

    input:
    path fastqc_results

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

process trimmomatic {
    tag { file(sample).baseName }
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path sample

    output:
    path "${sample.baseName}_trimmed.fastq"

    script:
    """
    trimmomatic SE -phred33 $sample ${sample.baseName}_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30
    """
}
*/
process chopper{
    conda './envs/Chopper.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/chopped", mode: 'copy'

    input:
    path sample

    output:
    path "${sample.baseName}_chopped.fastq"

    script:
    """
    chopper -q 10 -l 500 -i $sample > ${sample.baseName}_chopped.fastq
    """
}

process plot_reads{
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

/*
process filt_plot_reads{
    tag { file(sample).baseName }
    publishDir "${params.outdir}/read_plots", mode: 'copy'

    input:
    path sample
    output:
    path "filt_read_length_distribution"
    script:
    """
    NanoPlot --fastq $sample --plots dot -o filt_read_length_distribution    
    """
    


}


process nfilt_reads{
    tag { file(sample).baseName }
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'

    input:
    path sample

    output:
    path "${sample.baseName}_filtered.fastq"

    script:
    """
    NanoFilt -q 8 -l 200 --headcrop 30 --tailcrop 30 $sample > ${sample.baseName}_filtered.fastq
    """
}


process filter_reads {
    tag { file(reads).baseName }
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path reads

    output:
    path "filtered_reads.fastq"
    script:
    """
    filtlong --min_length 1000 --keep_percent 90 $reads > filtered_reads.fastq
    """
}

*/ 

process assemble_long_reads {
    conda './envs/Flye.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path sample

    output:
    path "assembly"

    script:
    """
    flye --nano-raw $sample --out-dir assembly --threads 4
    """
}

process polish_assembly {
    tag { file(sample).baseName }
    publishDir "${params.outdir}/polished", mode: 'copy'

    input:
    path sample
    path assembled_fasta

    output:
    path "polished_assembly"

    script:
    """
    medaka_consensus -i $sample -d $assembled_fasta -o polished_assembly -t 4 -m r1041_e82_400bps_fast_g632
    """

}


// BWA Index
process bwa_index {
    label 'Index_Reference_Genome_Map_BWA'
    conda './envs/vcf_env.yml'
    tag "$ref.baseName"
    
    input:
        path ref
    
    output:
        path ref, emit: indexed_ref
        path "${ref}.*", emit: index_files
    
    script:
    """
    bwa index $ref
    """
}




process assess_quality{ 
    conda './envs/Quast.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/quality_assessment", mode: 'copy'  
    input:
    path sample
    path reference_genome
    output:
    path "quast_results"
    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    quast $sample -r $reference_genome -o quast_results --threads 4
    """
}

process annotate_genome {
    conda './envs/Braker.yml'
    tag { file(assembly).baseName }
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path assembly

    output:
    path "annotated_genome.gff3"


    script:
    """
    braker.p1 --genome=$assembly --species=mammal --softmasking --cores=4
    """
}


