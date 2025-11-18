# Example Usage

## Example 1: Basic Local Run

### Directory Structure
```
mampipe/
├── data/
│   └── sample_reads.fastq
├── refSeq/
│   └── reference_genome.fna
└── main.nf
```

### Command
```bash
nextflow run main.nf \
  --reads 'data/sample_reads.fastq' \
  --ref 'refSeq/reference_genome.fna' \
  --outdir 'results' \
  --busco_threshold 6.0 \
  --lineage 'mammalia_odb10'
```

## Example 2: Multiple Samples

### Directory Structure
```
mampipe/
├── data/
│   ├── sample1.fastq
│   ├── sample2.fastq
│   └── sample3.fastq
├── refSeq/
│   └── reference.fna
└── main.nf
```

### Command
```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --outdir 'results' \
  --busco_threshold 6.0 \
  --lineage 'mammalia_odb10'
```

## Example 3: Google Cloud Batch

### Setup
```bash
# 1. Upload data to Google Cloud Storage
gsutil -m cp data/*.fastq gs://my-bucket/data/
gsutil -m cp refSeq/*.fna gs://my-bucket/refSeq/

# 2. Build and push Docker image
docker build -t gcr.io/my-project/mampipe:latest .
docker push gcr.io/my-project/mampipe:latest

# 3. Configure nextflow.config
cat > nextflow.config << 'EOF'
process {
    executor = 'google-batch'
    container = 'gcr.io/my-project/mampipe:latest'
    workDir = 'gs://my-bucket/work'
}

google {
    project = 'my-project'
    region = 'us-central1'
    batch {
        spot = true
        bootDiskSize = '50 GB'
    }
}

conda {
  enabled = true
  autoActivate = true
  cacheable = true
  cacheDir = '/tmp/conda-cache'
}
EOF
```

### Run
```bash
# Copy and configure batch template
cp batch_main.nf.template batch_main.nf

# Edit batch_main.nf to set:
# params.reads = 'gs://my-bucket/data/*.fastq'
# params.ref = 'gs://my-bucket/refSeq/*.fna'
# params.outdir = 'gs://my-bucket/results'

# Run pipeline
nextflow run batch_main.nf -c nextflow.config
```

## Example 4: Post-Assembly Processing Only

If you already have a polished assembly and want to run BUSCO, purging, and scaffolding:

```bash
# Place your assembly at: results/polished/polished_assembly/consensus.fasta

nextflow run Post_Assembly.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --outdir 'DUP_BUSCO_RESULTS' \
  --busco_threshold 6.0 \
  --lineage 'mammalia_odb10'
```

## Example 5: Variant Calling

After assembly, call variants:

```bash
nextflow run vcf.nf \
  --ref 'refSeq/*.fna' \
  --outdir 'vfc_results'
```

This assumes:
- Chopped reads are in `results/chopped/*.fastq`
- Reference genome is in `refSeq/*.fna`

## Example 6: Custom BUSCO Threshold

Adjust the duplication threshold for purging:

```bash
# Lower threshold (more aggressive purging)
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --busco_threshold 3.0 \
  --lineage 'mammalia_odb10'

# Higher threshold (less purging)
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --busco_threshold 10.0 \
  --lineage 'mammalia_odb10'
```

## Example 7: Different Lineages

### For a mouse genome
```bash
nextflow run main.nf \
  --reads 'data/mouse_reads.fastq' \
  --ref 'refSeq/mouse_ref.fna' \
  --lineage 'mammalia_odb10'
```

### For a bird genome
```bash
nextflow run main.nf \
  --reads 'data/bird_reads.fastq' \
  --ref 'refSeq/bird_ref.fna' \
  --lineage 'aves_odb10'
```

### For a fish genome
```bash
nextflow run main.nf \
  --reads 'data/fish_reads.fastq' \
  --ref 'refSeq/fish_ref.fna' \
  --lineage 'actinopterygii_odb10'
```

## Expected Output

After a successful run, you'll find:

```
results/
├── read_plots/
│   ├── read_length_distribution/
│   │   ├── NanoPlot-report.html
│   │   └── *.png
│   └── chopped_read_length_distribution/
│       ├── NanoPlot-report.html
│       └── *.png
├── chopped/
│   └── sample_chopped.fastq
├── assembly/
│   ├── assembly.fasta
│   ├── assembly_info.txt
│   └── assembly_graph.gfa
├── polished/
│   └── polished_assembly/
│       └── consensus.fasta
├── busco_precheck/
│   └── short_summary.specific.mammalia_odb10.busco_results.txt
├── purged_output/           # Only if duplication > threshold
│   └── purged.fa
├── scaffolded/
│   └── ragtag_output/
│       ├── ragtag.scaffold.fasta
│       └── ragtag.scaffold.agp
└── FINAL_BUSCO/
    └── FINAL_BUSCO_QC/
        └── short_summary.specific.mammalia_odb10.FINAL_BUSCO_QC.txt
```

## Monitoring Progress

### Local Execution
```bash
# Watch the log
tail -f .nextflow.log

# View the execution report
nextflow run main.nf -with-report report.html
```

### Google Cloud Batch
```bash
# Monitor in Google Cloud Console
# Navigate to: Batch > Jobs

# Or use gcloud CLI
gcloud batch jobs list --location=us-central1

# View logs
gcloud logging read "resource.type=batch_job" --limit 50
```

## Resuming Failed Runs

Nextflow supports resuming from the last successful step:

```bash
nextflow run main.nf -resume
```

This is especially useful for long-running pipelines or if a step fails.
