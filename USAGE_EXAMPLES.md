# Usage Examples - MamPipe

## Basic Usage

### Example 1: Default Run (Simplest)
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna'
```
Uses all default parameters:
- threads: 8
- threads_assembly: 4
- threads_polish: 4
- medaka_model: r1041_e82_400bps_fast_g632
- busco_threshold: 6.0
- lineage: mammalia_odb10

---

## Advanced Usage

### Example 2: Custom Thread Counts
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --threads 16 \
  --threads_assembly 8 \
  --threads_polish 8
```
**Use case**: You have a high-performance server with many cores

### Example 3: Different Sequencing Chemistry
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --medaka_model 'r941_min_high_g360'
```
**Use case**: You used R9.4.1 flowcells instead of R10.4.1

Common Medaka models:
- `r1041_e82_400bps_fast_g632` - R10.4.1, fast basecalling
- `r1041_e82_400bps_hac_g632` - R10.4.1, high accuracy
- `r1041_e82_400bps_sup_g615` - R10.4.1, super accuracy
- `r941_min_high_g360` - R9.4.1, high accuracy
- `r941_prom_high_g360` - R9.4.1 PromethION

### Example 4: Custom BUSCO Lineage and Threshold
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --lineage 'vertebrata_odb10' \
  --busco_threshold 8.0
```
**Use case**: Non-mammalian vertebrate with higher duplication tolerance

Common lineages:
- `mammalia_odb10` - Mammals
- `vertebrata_odb10` - Vertebrates
- `aves_odb10` - Birds
- `actinopterygii_odb10` - Ray-finned fishes
- `metazoa_odb10` - Animals

### Example 5: Low-Resource System
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --threads 4 \
  --threads_assembly 2 \
  --threads_polish 2
```
**Use case**: Running on a laptop or small VM

---

## Post-Assembly Processing

### Example 6: Process Existing Assembly
```bash
nextflow run Post_Assembly.nf \
  --assembly_input 'results/polished/polished_assembly/consensus.fasta' \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --lineage 'mammalia_odb10'
```
**Use case**: You already have an assembly and want to run BUSCO, purging, and scaffolding

### Example 7: Post-Assembly with Custom Paths
```bash
nextflow run Post_Assembly.nf \
  --assembly_input 'my_assembly/genome.fasta' \
  --reads 'my_reads/sample.fastq' \
  --ref 'my_ref/reference.fna' \
  --outdir 'post_assembly_results' \
  --threads 12
```
**Use case**: Custom directory structure

---

## Variant Calling

### Example 8: Call Variants from Assembly
```bash
nextflow run vcf.nf \
  --ref 'refSeq/*.fna' \
  --chopped_reads 'results/chopped/*.fastq' \
  --outdir 'variants' \
  --threads 8
```
**Use case**: After running main.nf, call variants against reference

### Example 9: Variant Calling with Custom Paths
```bash
nextflow run vcf.nf \
  --ref 'my_ref/reference.fna' \
  --chopped_reads 'my_filtered_reads/*.fastq' \
  --outdir 'my_variants' \
  --threads 4
```

---

## Google Cloud Batch

### Example 10: Run on Google Cloud
```bash
# First, configure nextflow.config for your project
nextflow run batch_main.nf \
  -c nextflow.config \
  --reads 'gs://my-bucket/data/*.fastq' \
  --ref 'gs://my-bucket/refSeq/*.fna' \
  --outdir 'gs://my-bucket/results' \
  --threads 16 \
  --threads_assembly 8
```
**Use case**: Large genome, need cloud compute

### Example 11: Cloud with Custom Resources
Edit `nextflow.config`:
```groovy
process {
    withLabel: 'high_mem' {
        memory = '64 GB'
        cpus = 16
        machineType = 'n1-highmem-16'
    }
}
```
Then run:
```bash
nextflow run batch_main.nf -c nextflow.config
```

---

## Combining Parameters

### Example 12: Full Custom Run
```bash
nextflow run main.nf \
  --reads 'data/sample_*.fastq' \
  --ref 'refSeq/mouse_genome.fna' \
  --outdir 'mouse_assembly' \
  --lineage 'mammalia_odb10' \
  --busco_threshold 5.0 \
  --threads 24 \
  --threads_assembly 12 \
  --threads_polish 8 \
  --medaka_model 'r1041_e82_400bps_sup_g615'
```
**Use case**: High-quality mouse genome assembly with super-accurate basecalling

### Example 13: Quick Test Run
```bash
nextflow run main.nf \
  --reads 'test_data/small_sample.fastq' \
  --ref 'test_data/small_ref.fna' \
  --outdir 'test_results' \
  --threads 2 \
  --threads_assembly 1 \
  --threads_polish 1
```
**Use case**: Testing pipeline with small dataset

---

## Resume Failed Runs

### Example 14: Resume After Failure
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  -resume
```
**Use case**: Pipeline failed partway through, resume from last successful step

---

## Generating Reports

### Example 15: Run with Execution Report
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  -with-report report.html \
  -with-timeline timeline.html \
  -with-dag flowchart.html
```
**Use case**: Generate visual reports of pipeline execution

---

## Using Optional Processes

### Example 16: Include Quality Assessment
```groovy
// In a custom workflow file
include { assess_quality } from './optional_processes.nf'

workflow {
    // ... your main workflow ...
    
    // Add quality assessment
    quality_results = assess_quality(final_assembly, reference_genome)
}
```

### Example 17: Include Genome Annotation
```groovy
include { annotate_genome } from './optional_processes.nf'

workflow {
    // ... your main workflow ...
    
    // Add annotation
    annotations = annotate_genome(final_assembly)
}
```

---

## Troubleshooting Examples

### Example 18: Check Input Files
```bash
# The pipeline will now validate inputs and give clear errors
nextflow run main.nf \
  --reads 'wrong_path/*.fastq' \
  --ref 'refSeq/*.fna'

# Output: ERROR: Reads not found at wrong_path/*.fastq
```

### Example 19: Dry Run (Check Configuration)
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  -preview
```

---

## Performance Tuning

### Example 20: Optimize for Speed
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --threads 32 \
  --threads_assembly 16 \
  --threads_polish 16 \
  --medaka_model 'r1041_e82_400bps_fast_g632'  # Fast model
```

### Example 21: Optimize for Quality
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --threads 16 \
  --threads_assembly 8 \
  --threads_polish 8 \
  --medaka_model 'r1041_e82_400bps_sup_g615'  # Super-accurate model
  --busco_threshold 3.0  # More aggressive purging
```

---

## Real-World Scenarios

### Scenario 1: Mouse Genome Assembly
```bash
nextflow run main.nf \
  --reads 'data/mouse_nanopore.fastq' \
  --ref 'refSeq/GRCm39.fna' \
  --outdir 'mouse_results' \
  --lineage 'mammalia_odb10' \
  --threads 16 \
  --medaka_model 'r1041_e82_400bps_hac_g632'
```

### Scenario 2: Human Chromosome Assembly
```bash
nextflow run main.nf \
  --reads 'data/chr22_reads.fastq' \
  --ref 'refSeq/chr22_reference.fna' \
  --outdir 'chr22_assembly' \
  --lineage 'mammalia_odb10' \
  --busco_threshold 8.0 \
  --threads 32 \
  --threads_assembly 16
```

### Scenario 3: De Novo Assembly (No Reference)
```bash
# Skip scaffolding by modifying workflow or use partial pipeline
nextflow run main.nf \
  --reads 'data/novel_species.fastq' \
  --ref 'refSeq/closest_relative.fna' \
  --outdir 'denovo_results' \
  --lineage 'vertebrata_odb10' \
  --threads 24
```

---

## Tips

1. **Start small**: Test with a subset of your data first
2. **Monitor resources**: Use `-with-report` to see resource usage
3. **Use resume**: Always use `-resume` if rerunning after changes
4. **Check BUSCO**: Review BUSCO results to decide if purging is needed
5. **Adjust thresholds**: If assembly quality is poor, try different `--busco_threshold` values
6. **Match chemistry**: Always use the correct `--medaka_model` for your data
7. **Scale threads**: More threads isn't always better - match to your data size

---

## Getting Help

If you encounter issues:
1. Check the error message (now more descriptive!)
2. Verify input files exist
3. Check `.nextflow.log` for details
4. Review BUSCO results for quality issues
5. Try with fewer threads if memory errors occur
6. Open an issue on GitHub with your command and error
