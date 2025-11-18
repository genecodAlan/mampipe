#!/usr/bin/env nextflow

/*
# Conda environment setup:
conda create -n vfc_env bwa samtools gatk bcftools  -c bioconda -c conda-forge
*/

/* Lets build and effective VCF pipeline that follows the GATK standards of the broad institute.
    This pipeline will include the following steps:
    1. Index the reference genome using BWA and samtools
    2. Align the reads to the indexed reference genome using BWA mem
    3. Sort the BAM files using samtools
    4. Mark duplicates using GATK
    5. Call variants using GATK HaplotypeCaller
    6. Filter the variants using bcftools
    7. Output the filtered VCF file
*/

// We are skipping basequality recalibration as it is optional and not always necessary for nanopore data. It can always be added later if needed.
// Base quality recalibration is only possible when knowing areas of variation which you could bootstrap as well, but for now we will skip it.
// We will also skip the joint genotyping step as it is not necessary for single sample variant

nextflow.enable.dsl=2

// Default parameters
params.ref = 'refSeq/*.fna'
params.outdir = 'vfc_results'
params.chopped_reads = 'results/chopped/*.fastq'
params.threads = 4

workflow {
    // Validate inputs exist
    if (!file(params.ref).exists()) {
        error "ERROR: Reference genome not found at ${params.ref}"
    }
    if (!file(params.chopped_reads).exists()) {
        error "ERROR: Chopped reads not found at ${params.chopped_reads}"
    }
    
    //Establish Path to reference Sequence and Raw trimmed reads
    reference = Channel.fromPath(params.ref).first()
    reads  = Channel.fromPath(params.chopped_reads, checkIfExists: true)

    //Set up the reference index and dictionary for GATK
    bwa_index_out = bwa_index(reference)
    faidx_out     = faidx_index(reference)
    dict_out      = mk_dict(reference)

    //Align reads to reference 
    sam_files = bwaAlign(reads, reference, bwa_index_out.index_files)

    //Sort the bam file and correct for duplicates (allows for efficient indexing)
    sorted_bam = sortBam(sam_files)
    dedup_bam  = markDuplicates(sorted_bam)
    //Index the final bam file for GATK
    bamIndexed = bamIndex(dedup_bam)
    //Call variants 
    raw_vcf  = callVariants(dedup_bam, bamIndexed, reference, faidx_out.fai, dict_out.dict)

    // Filter
    filtered_vcf = filterVCF(raw_vcf)
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

//Alignment
process bwaAlign {
    label "Map_to_Reference_Genome_BWAmem"
    conda './envs/vcf_env.yml'
    tag "$reads.baseName"
    publishDir "${params.outdir}/bwa", mode: 'copy'
    
    input:
        path reads
        path ref
        path index_files
    
    output:
        path "${reads.baseName}.sam"
    
    script:
    """
    bwa mem -t ${params.threads} -R "@RG\\tID:${reads.baseName}\\tPL:ONT\\tSM:${reads.baseName}" $ref $reads > ${reads.baseName}.sam
    """
}


//Sort BAM
process sortBam {
    label "Sort_BAM_Files_Samtools"
    conda './envs/vcf_env.yml'
    tag "$sam.baseName"
    publishDir "${params.outdir}/sorted_bam", mode: 'copy'
    input:
        path sam
    output:
        path "${sam.baseName}.sorted.bam"
    script:
    """
    samtools sort -o ${sam.baseName}.sorted.bam $sam
    """
}

//Index final BAM
process bamIndex {
    label "Index_BAM_Files_Samtools"
    conda './envs/vcf_env.yml'
    tag "$bam.baseName"
    publishDir "${params.outdir}/indexed_bam", mode: 'copy'

    input:
        path bam

    output:
        path("${bam}.bai")

    script:
    """
    samtools index $bam
    """
}


// Mark Duplicates
process markDuplicates {
    label "Mark_Duplicates_GATK"
    conda './envs/gatk4.yml'
    tag "$bam.baseName"
    publishDir "${params.outdir}/dedup_bam", mode: 'copy'

    input:
        path bam
    output:
        path "${bam.baseName}.dedup.bam"
    script:
    """
    gatk MarkDuplicates -I $bam -O ${bam.baseName}.dedup.bam -M ${bam.baseName}.metrics.txt
    """
}

//Get fai index for ref 
process faidx_index {
    label "Index_Reference_Genome_Samtools"
    conda './envs/vcf_env.yml'
    tag "$ref.baseName"

    input:
        path ref

    output:
        path "${ref}.fai", emit: fai

    script:
    """
    samtools faidx $ref
    """
}

// Create dict for reference 
process mk_dict {
    label "Create_RefSeq_Dictionary_GATK"
    conda './envs/gatk4.yml'
    tag "$ref.baseName"

    input:
        path ref

    output:
        path "${ref.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R $ref -O ${ref.baseName}.dict
    """
}



// GATK4 variant Calling
process callVariants {
    label "Call_Variants_GATK_HaplotypeCaller"
    conda './envs/gatk4.yml'
    tag "$bam.baseName"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        path bam
        path bai 
        path ref 
        path fai
        path dict

    output:
        path "${bam.baseName}.raw.vcf"

    script:
    """
    gatk HaplotypeCaller -R $ref -I $bam -O ${bam.baseName}.raw.vcf
    """
}




//Filter Variants
process filterVCF {
    label "Filter_Variants_bcftools"
    conda './envs/vcf_env.yml'
    tag "$vcf.baseName"
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'
    
    input:
        path vcf
    output:
        path "${vcf.baseName}.filtered.vcf"
    script:
    """
    bcftools filter -e 'INFO/DP<10 || QUAL<20' $vcf > ${vcf.baseName}.filtered.vcf
    """
}