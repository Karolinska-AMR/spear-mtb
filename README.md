## <img src="https://raw.githubusercontent.com/xrazmo/spear-mtb/main/assets/report/img/logo.png" width="50" height="50" > SPEAR-MTB: en*S*emble *P*r*E*diction of *A*ntibiotic *R*esistance in *M*ycobacterium *T*u*B*erculosis

This repository has employed the following pipelines in order to arrive at a consensus on drug-resistant _Mycobacterium tuberculosis_ predictions:

- **[CRyPTIC](https://github.com/iqbal-lab-org)**: The pipeline involves decontaminating reads, performing quality controls and variant calling using various tools and the the H37Rv reference genome. We have incorporated the pipeline into the nextflow DSL2 and employed precompiled reference genome databases (e.g., H37Rv, nontuberculous Mycobacterium, human, etc.) to facilitate analyzing input genomes. Moreover, the [gnomonicus](https://github.com/oxfordmmm/gnomonicus) repository is used to predict durg reistance based on the detected mutations and catalogues available. The spear-MTB utilizes _WHO_ and _CRyPTIC_ catalogues, which are presented in [GARC1](https://fowlerlab.org/2018/11/25/goarc-a-general-ontology-for-antimicrobial-resistance-catalogues/) format and maintained [here](https://github.com/oxfordmmm/tuberculosis_amr_catalogues).

- **[TB Profiler](https://github.com/jodyphelan/TBProfiler)**:
  The pipeline includes aligning reads to the H37Rv reference genome, utilizing a pairwise aligner and then identifying variations with the use of bcftools. The catalogue of mutations is in [hgvs nomenclature](http://varnomen.hgvs.org/bg-material/simple/) and is maintained [here](https://github.com/jodyphelan/tbdb).

> &#x26A0; **Warning:**
> The spear-MTB pipeline is tested using Illumina paired-end reads and on Unix-based system having Singularity.

## Prerequisites

- **Nextflow**
- **Conda**
- **Singularity**
-**git** 

## **Installation**

After the successful installations of prerequisites, proceed with running the following commands:

```
git clone https://github.com/Karolinska-AMR/spear-mtb.git
cd ./spear-mtb | chmod +x ./setup.sh | ./setup.sh
```

These commands perform the following steps:

- Downloading the repository from Github.
- Creating a conda environment named _spear-mtb_.
- Downloading and extracting the precompiled assets, which encompasses indexed genomes, the reference H37rv, and a KRAKEN2 database containing genomes belong to _Mycobacteriaceae_ family.


## **Usage**

### - Nextflow config file:

Before running the pipeline, please review and modify the configuration file (nextflow.config) to suit your specific requirements and preferences. Note that the current configuration file includes predefined profiles and resource allocation.

### Basic Usage

```bash
./spear-mtb.sh --input_dir <path> --output_dir <path>
```

### Full Command Syntax

```bash
./spear-mtb.sh --input_dir <path> --output_dir <path> [OPTIONS]
```

## Parameters

### Required Parameters

| Parameter | Short | Description |
|-----------|-------|-------------|
| `--input_dir` | `-i` | Input directory path containing the data to process |
| `--output_dir` | `-o` | Output directory path where results will be saved |

### Optional Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--working_dir` | `-w` | Working directory for temporary files | `./[source_dir]/.tmp` |
| `--config_file` | `-c` | Nextflow configuration file | `./nextflow.config` |
| `--assets_dir` | `-a` | Assets directory path | `./assets` |
| `--profile` | `-p` | Nextflow profile to use | `slurm` |
| `--archive_dir` | `-A` | Directory for archiving input files | None |
| `--resume` | `-r` | Resume previous pipeline run | `false` |
| `--skip_cryptic` | `-s` | Skip cryptic workflow steps | `false` |
| `--ticket` | `-t` | Ticket ID for tracking | Auto-generated |
| `--singularity_cache` | `-S` | Singularity cache directory | None |
| `--help` | `-h` | Show help message | - |

## Examples

### Basic execution:
```bash
./spear-mtb.sh --input_dir /data/input --output_dir /data/output
```

### With custom profile and resume:
```bash
./spear-mtb.sh -i /data/input -o /data/output --profile slurm --resume
```

### Full parameter example:
```bash
./spear-mtb.sh \
  --input_dir /data/input \
  --output_dir /data/output \
  --working_dir /tmp/work \
  --config_file /path/to/custom.config \
  --profile slurm \
  --archive_dir /data/archive \
  --singularity_cache /cache/singularity \
  --ticket ABC12345 \
  --resume
```

## Output Structure

The pipeline generates a structured output directory with the following organization:

```
{output-directory}/
├── logs/                                    # Execution logs and traces
│   ├── {ticket}_trace.txt                  # Nextflow execution trace
│   ├── {ticket}_timeline.html              # Timeline visualization
│   └── {ticket}_report.html                # Execution report
└── results/                                 # Analysis results
    ├── {ticket}.json                       # Summary results in JSON format
    ├── {ticket}_spear-mtb_report.html      # Comprehensive HTML report
    └── intermediate/                        # Per-sample analysis results
        └── {fastq_file_name}/              # Individual sample directory
            ├── cryptic/*                    # CRyPTIC pipeline results
            └── tbprofiler/*                 # TB-Profiler pipeline results
              
```

### Key Output Files

- **`{ticket}.json`**: Machine-readable summary of all results
- **`{ticket}_spear-mtb_report.html`**: Interactive HTML report with visualizations
- **`logs/{ticket}_*`**: Nextflow execution logs, traces, and performance metrics
- **`intermediate/{sample}/cryptic/`**: CRyPTIC pipeline outputs for drug resistance analysis
- **`intermediate/{sample}/tbprofiler/`**: TB-Profiler results for genotype-phenotype predictions

> **Note**: `{ticket}` represents your unique run identifier, and `{fastq_file_name}` corresponds to each input sample processed.
