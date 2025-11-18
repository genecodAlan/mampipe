#!/usr/bin/env nextflow

/*
 * Optional processes that can be included in custom workflows
 * These are not part of the main pipeline but available for advanced users
 */

nextflow.enable.dsl=2

// Quality assessment with QUAST
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
    quast $sample -r $reference_genome -o quast_results --threads ${params.threads_assembly}
    """
}

// Genome annotation with BRAKER
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
    braker.p1 --genome=$assembly --species=mammal --softmasking --cores=${params.threads_assembly}
    """
}

// BWA Index (used by vcf.nf)
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
