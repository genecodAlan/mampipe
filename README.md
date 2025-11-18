# MamPipe - Mammalian Whole Genome Sequencing Pipeline

A comprehensive Nextflow pipeline for mammalian whole genome sequencing (WGS) analysis using Oxford Nanopore long-read data. This pipeline can run locally or on Google Cloud Batch for scalable compute.

## Pipeline Overview

MamPipe performs end-to-end genome assembly and quality assessment:

1. **Quality Control** - NanoPlot visualization of read quality
2. **Read Filtering** - Chopper for quality and length filtering
3. **Assembly** - Flye for de novo genome assembly
4. **Polishing** - Medaka for consensus polishing
5. **Quality Assessment** - BUSCO for completeness evaluation
6. **Duplicate Purging** - Purge_dups (conditional, based on BUSCO duplication score)
7. **Scaffolding** - RagTag for reference-guided scaffolding
8. **Final QC** - BUSCO assessment of final assembly

### Optional Modules

- **Variant Calling** (`vcf.nf`) - BWA alignment, GATK variant calling, and filtering

## Requirements

### Local Execution
- Nextflow (≥22.04)
- Conda or Mamba
- 16+ GB RAM recommended
- Multi-core CPU (4+ cores recommended)

### Google Cloud Batch Execution
- Google Cloud Project with Batch API enabled
- Google Cloud Storage bucket
- Docker image (provided Dockerfile)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/genecodAlan/mampipe.git
cd mampipe
```

2. Install Nextflow:
```bash
curl -s https://get.nextflow.io | bash
```

3. Ensure Conda/Mamba is installed for environment management

## Usage

### Local Execution

```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --outdir 'results' \
  --busco_threshold 6.0 \
  --lineage 'mammalia_odb10'
```

### Google Cloud Batch Execution

1. Build and push the Docker image:
```bash
cd Github_Upload
docker build -t gcr.io/YOUR_PROJECT/mampipe:latest .
docker push gcr.io/YOUR_PROJECT/mampipe:latest
```

2. Configure your Google Cloud settings in `nextflow.config`:
```groovy
process {
    executor = 'google-batch'
    container = 'gcr.io/YOUR_PROJECT/mampipe:latest'
    workDir = 'gs://YOUR_BUCKET/work'
}

google {
    project = 'YOUR_PROJECT_ID'
    region = 'us-central1'
    batch {
        spot = true
        bootDiskSize = '50 GB'
    }
}
```

3. Run the pipeline:
```bash
nextflow run batch_main.nf -c nextflow.config
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `data/*.fastq` | Path to input FASTQ files |
| `--ref` | `refSeq/*.fna` | Path to reference genome |
| `--outdir` | `WGS_results` | Output directory |
| `--busco_threshold` | `6.0` | BUSCO duplication threshold (%) |
| `--lineage` | `bacteria_odb10` | BUSCO lineage dataset |

## Pipeline Workflows

### Main Assembly Pipeline (`main.nf` / `batch_main.nf`)
Complete genome assembly from raw reads to scaffolded assembly

### Post-Assembly Processing (`Post_Assembly.nf` / `batch_Post_assembly.nf`)
BUSCO assessment, conditional purging, and scaffolding for existing assemblies

### Variant Calling (`vcf.nf`)
SNP/indel calling from aligned reads

## Output Structure

```
results/
├── read_plots/              # NanoPlot QC reports
├── chopped/                 # Filtered reads
├── assembly/                # Flye assembly output
├── polished/                # Medaka polished assembly
├── busco_precheck/          # Initial BUSCO results
├── purged_output/           # Purge_dups output (if triggered)
├── scaffolded/              # RagTag scaffolded assembly
└── FINAL_BUSCO/             # Final quality assessment
```

## Tools & Versions

- **NanoPlot** - Read quality visualization
- **Chopper** - Read filtering
- **Flye** - Genome assembly
- **Medaka** - Assembly polishing
- **BUSCO** - Completeness assessment
- **Purge_dups** - Haplotig purging
- **RagTag** - Reference-guided scaffolding
- **minimap2** - Read alignment
- **GATK4** - Variant calling (optional)

## Citation

If you use this pipeline, please cite the individual tools used in your analysis.

## License

MIT License

## Contact

For issues and questions, please open an issue on GitHub.
