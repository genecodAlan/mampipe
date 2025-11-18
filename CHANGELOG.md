# Changelog

All notable changes to MamPipe will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Parameterized thread counts**: `--threads`, `--threads_assembly`, `--threads_polish` for flexible resource allocation
- **Medaka model parameter**: `--medaka_model` to support different sequencing chemistries
- **Assembly input parameter**: `--assembly_input` for Post_Assembly.nf to specify custom assembly paths
- **Resource management**: Process labels (`low_mem`, `medium_mem`, `high_mem`, `very_high_mem`) for better scheduling
- **Input validation**: Checks if input files exist before starting pipeline with clear error messages
- **Retry strategy**: Automatic retry on transient failures (max 2 retries)
- **Nanoplot.yml dependencies**: Fixed empty environment file with proper dependencies
- **optional_processes.nf**: Separated unused processes (assess_quality, annotate_genome, bwa_index) for advanced users

### Changed
- **Post_Assembly.nf**: Removed hardcoded paths, now uses parameters for flexibility
- **Scaffolded output**: Now uses process output directly instead of hardcoded path
- **Thread counts**: All hardcoded thread values replaced with parameters
- **Default lineage**: Changed from `bacteria_odb10` to `mammalia_odb10` for mammalian genomes
- **vcf.nf paths**: Fixed inconsistent parameter names (`ref_Seq` → `refSeq/*.fna`)
- **vcf.nf reads**: Parameterized chopped reads path with `--chopped_reads`

### Fixed
- **Missing conda environment**: Added `conda './envs/Medaka.yml'` to polish_assembly process
- **Nanoplot.yml**: Added missing dependencies (was empty)
- **Path consistency**: Unified reference genome paths across all workflows

### Removed
- **Commented code blocks**: Removed 150+ lines of unused commented processes (fastqc, multiqc, trimmomatic, etc.)
- **Unused processes from main.nf**: Moved assess_quality, annotate_genome, bwa_index to optional_processes.nf
- **Redundant comments**: Cleaned up old workflow logic comments

## [1.0.0] - 2025-01-XX

### Added
- Initial release of MamPipe
- Complete WGS pipeline for mammalian genomes
- Support for Oxford Nanopore long-read data
- Google Cloud Batch integration
- Conda environment management
- Docker support

### Pipeline Features
- Read quality control with NanoPlot
- Read filtering with Chopper
- De novo assembly with Flye
- Consensus polishing with Medaka
- Quality assessment with BUSCO
- Conditional duplicate purging with Purge_dups
- Reference-guided scaffolding with RagTag
- Optional variant calling with GATK

---

## Migration Guide

### From Original Version to Refactored Version

If you were using the original version, here's what changed:

#### 1. Thread Parameters (Optional - Backward Compatible)
You can now customize thread counts:
```bash
# Old way (still works):
nextflow run main.nf --reads data/*.fastq --ref refSeq/*.fna

# New way (with custom threads):
nextflow run main.nf \
  --reads data/*.fastq \
  --ref refSeq/*.fna \
  --threads 16 \
  --threads_assembly 8 \
  --threads_polish 8
```

#### 2. Medaka Model (Optional - Backward Compatible)
Default model is still `r1041_e82_400bps_fast_g632`, but you can change it:
```bash
nextflow run main.nf \
  --reads data/*.fastq \
  --ref refSeq/*.fna \
  --medaka_model r941_min_high_g360
```

#### 3. Post_Assembly.nf Changes (Breaking Change)
If you use Post_Assembly.nf, you now need to specify the assembly path:
```bash
# Old way (hardcoded path):
nextflow run Post_Assembly.nf

# New way (specify path):
nextflow run Post_Assembly.nf \
  --assembly_input results/polished/polished_assembly/consensus.fasta \
  --reads data/*.fastq \
  --ref refSeq/*.fna
```

#### 4. vcf.nf Changes (Breaking Change)
Parameter names changed for consistency:
```bash
# Old way:
# params.ref = 'ref_Seq'

# New way:
nextflow run vcf.nf \
  --ref refSeq/*.fna \
  --chopped_reads results/chopped/*.fastq
```

#### 5. Resource Management
The pipeline now has better resource management. If you were setting resources manually in your config, the new labels will override them. You can adjust in `nextflow.config`:
```groovy
process {
    withLabel: 'high_mem' {
        memory = '64 GB'  // Increase if needed
        cpus = 16
    }
}
```

---

## Upgrade Instructions

1. **Pull latest changes**:
   ```bash
   git pull origin main
   ```

2. **Update your run commands** (if using Post_Assembly.nf or vcf.nf):
   - Add `--assembly_input` parameter to Post_Assembly.nf
   - Update `--ref` parameter in vcf.nf if you were using `ref_Seq`

3. **Test with your data**:
   ```bash
   nextflow run main.nf --reads your_data/*.fastq --ref your_ref/*.fna
   ```

4. **Optional**: Customize thread counts and Medaka model for your setup

---

## Support

For questions about changes or migration issues, please open an issue on GitHub.
