

nextflow.enable.dsl=2

// Default parameters for local execution
params.reads = 'data/*.fastq'
params.ref = 'refSeq/*.fna'
params.outdir = 'DUP_BUSCO_RESULTS'
params.busco_threshold = 6.0
params.lineage = 'mammalia_odb10'
params.assembly_input = 'results/polished/polished_assembly/consensus.fasta'

// Thread parameters
params.threads = 8

// For Google Cloud Batch, override with:
// params.reads = 'gs://YOUR_BUCKET/data/*.fastq'
// params.ref = 'gs://YOUR_BUCKET/refSeq/*.fna'
// params.outdir = 'gs://YOUR_BUCKET/results'
// params.assembly_input = 'gs://YOUR_BUCKET/results/polished/polished_assembly/consensus.fasta'




workflow {
    // Validate inputs exist
    if (!file(params.reads).exists()) {
        error "ERROR: Reads not found at ${params.reads}"
    }
    if (!file(params.ref).exists()) {
        error "ERROR: Reference genome not found at ${params.ref}"
    }
    if (!file(params.assembly_input).exists()) {
        error "ERROR: Assembly not found at ${params.assembly_input}"
    }
    
    // Define input channels
    reads = Channel.fromPath(params.reads).first()
    reference_genome = Channel.fromPath(params.ref)
    assembly = Channel.fromPath(params.assembly_input)


// Run BUSCO on initial assembly
busco_results_ch = busco_precheck(assembly)
// Extract duplication score
dup_score = check_busco_dup(busco_results_ch)  // emits val dup_score

final_score = read_dup_score(dup_score)



// Run purge_dups only if dup % exceeds threshold
purged_assembly_ch = purge_duplicates(assembly, reads, final_score)

// Choose final assembly based on duplication score
final_assembly_ch = purged_assembly_ch.ifEmpty { assembly }
// Proceed with scaffolding
scaffolded_assembly_ch = scaffold_assembly(final_assembly_ch, reference_genome)

// Extract the scaffolded fasta from the output directory
scaffolded = scaffolded_assembly_ch.flatMap { dir -> 
    file("${dir}/ragtag.scaffold.fasta")
}
busco_final_results = busco_final(scaffolded)

}

process busco_precheck {
    label 'medium_mem'
    conda './envs/BUSCO.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/busco_precheck", mode: 'copy'
    input:
    path sample
    output:
    path "busco_results/short_summary.specific.${params.lineage}.busco_results.txt" 
    script: 
    """
    busco -i $sample -l ${params.lineage} -o busco_results --mode genome
    """
}


process check_busco_dup {
    conda './envs/BUSCO.yml'
    publishDir "${params.outdir}/busco_dup_check", mode: 'copy'
    input:
    path summary_txt
    output:
    path "dup_score.txt"

    script:
    """
    grep -o 'D:[0-9.]*%' ${summary_txt} | cut -d ':' -f2 | tr -d '%' >> dup_score.txt || touch dup_score.txt
    """
    
}



process read_dup_score {
    input:
        path dup_score_txt

    output:
        stdout

    script:
    """
    cat ${dup_score_txt}
    """
}

process test_result {
    input:
    val result

    when:
    result.toDouble() < params.busco_threshold

    output:
    stdout

    script:
    //print out the score and say data recieved to console
    """
    echo "Duplication score received: $result"
    echo "Data received successfully."
    """
}

process purge_duplicates {
    label 'medium_mem'
    conda './envs/Purge.yml'
   
    input:
        tuple path(assembly), path(reads), val(dup_score)
    
    output:
        path "purged_output/purged.fa"
   
    script:
    """
    minimap2 -x map-ont -t ${params.threads} $assembly $reads | gzip -c > aln.paf.gz
    pbcstat aln.paf.gz
    calcuts PB.stat > cutoffs  
    split_fa $assembly > split.fasta
    minimap2 -x asm5 -DP -t ${params.threads} split.fasta split.fasta | gzip -c > self.paf.gz
    purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed
    get_seqs dups.bed $assembly > purged.fa
    mkdir -p purged_output && mv purged.fa purged_output/
    """
}

process scaffold_assembly {
    label 'medium_mem'
    conda './envs/RagTag.yml'
    tag { file(assembly).baseName }
    publishDir "${params.outdir}/scaffolded", mode: 'copy'

    input:
    path assembly
    path ref

    output:
    path "ragtag_output"

    script:
    """
    ragtag.py scaffold $ref $assembly
    """
}   

process busco_final {
    label 'medium_mem'
    conda './envs/BUSCO.yml'
    tag { file(sample).baseName }
    publishDir "${params.outdir}/FINAL_BUSCO", mode: 'copy'
    input:
    //only take the assmebly ending in fasta
    path sample 
    output:
    path "FINAL_BUSCO_QC" 
    script: 
    """
    busco -i $sample -l ${params.lineage} -o FINAL_BUSCO_QC --mode genome 
    """
}