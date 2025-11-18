# Setup Guide for MamPipe

## Local Setup

### 1. Install Dependencies

#### Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

#### Install Conda/Mamba
If you don't have Conda installed:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Install Mamba for faster environment resolution:
```bash
conda install -n base -c conda-forge mamba
```

### 2. Prepare Your Data

Create the following directory structure:
```
mampipe/
├── data/              # Place your FASTQ files here
│   └── sample.fastq
├── refSeq/            # Place your reference genome here
│   └── reference.fna
└── envs/              # Conda environment files (included)
```

### 3. Run the Pipeline

```bash
nextflow run main.nf \
  --reads 'data/*.fastq' \
  --ref 'refSeq/*.fna' \
  --outdir 'results' \
  --busco_threshold 6.0 \
  --lineage 'mammalia_odb10'
```

## Google Cloud Batch Setup

### 1. Prerequisites

- Google Cloud Project with billing enabled
- Google Cloud SDK installed
- Docker installed locally

### 2. Enable Required APIs

```bash
gcloud services enable batch.googleapis.com
gcloud services enable compute.googleapis.com
gcloud services enable storage.googleapis.com
```

### 3. Create a Google Cloud Storage Bucket

```bash
gsutil mb -l us-central1 gs://YOUR_BUCKET_NAME
```

### 4. Upload Your Data

```bash
# Upload reads
gsutil -m cp data/*.fastq gs://YOUR_BUCKET_NAME/data/

# Upload reference genome
gsutil -m cp refSeq/*.fna gs://YOUR_BUCKET_NAME/refSeq/
```

### 5. Build and Push Docker Image

```bash
# Authenticate with Google Container Registry
gcloud auth configure-docker

# Build the Docker image
cd Github_Upload
docker build -t gcr.io/YOUR_PROJECT_ID/mampipe:latest .

# Push to Google Container Registry
docker push gcr.io/YOUR_PROJECT_ID/mampipe:latest
```

### 6. Configure Nextflow

Copy the template configuration:
```bash
cp nextflow.config.template nextflow.config
```

Edit `nextflow.config` and replace:
- `YOUR_PROJECT_ID` with your Google Cloud project ID
- `YOUR_BUCKET_NAME` with your bucket name

Example configuration:
```groovy
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
```

### 7. Run on Google Cloud Batch

```bash
# Copy the template
cp batch_main.nf.template batch_main.nf

# Edit batch_main.nf to set your bucket paths
# Then run:
nextflow run batch_main.nf -c nextflow.config
```

## BUSCO Lineage Datasets

Choose the appropriate lineage for your organism:

- `mammalia_odb10` - Mammals
- `vertebrata_odb10` - Vertebrates
- `metazoa_odb10` - Animals
- `eukaryota_odb10` - Eukaryotes
- `bacteria_odb10` - Bacteria

Full list: https://busco-data.ezlab.org/v5/data/lineages/

## Troubleshooting

### Conda Environment Issues

If conda environments fail to create:
```bash
# Clear conda cache
conda clean --all

# Manually create an environment
conda env create -f envs/BUSCO.yml
```

### Google Cloud Authentication

```bash
# Login to Google Cloud
gcloud auth login

# Set your project
gcloud config set project YOUR_PROJECT_ID

# Configure application default credentials
gcloud auth application-default login
```

### Memory Issues

If you encounter out-of-memory errors, increase the machine type in your config:
```groovy
process {
    machineType = 'n1-highmem-8'  // 8 CPUs, 52 GB RAM
}
```

### Spot Instance Preemption

If using spot instances and jobs are preempted:
```groovy
google {
    batch {
        spot = false  // Use regular instances
    }
}
```

## Cost Optimization

1. **Use Spot Instances**: Save up to 80% on compute costs
2. **Choose Appropriate Regions**: Some regions are cheaper
3. **Clean Up Work Directory**: Delete intermediate files after successful runs
   ```bash
   gsutil -m rm -r gs://YOUR_BUCKET/work/
   ```
4. **Use Lifecycle Policies**: Auto-delete old files in your bucket

## Support

For issues, please open a GitHub issue with:
- Your Nextflow version (`nextflow -version`)
- Error messages
- Configuration file (with sensitive info removed)
