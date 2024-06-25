# RNAseq-Workflow-Snakemake-Singularity-ETHZ-LeoMed-SLURM

This repository contains an RNAseq workflow implemented using Snakemake. The workflow utilizes Singularity containers, which are verified by Nexus, to ensure reproducibility and ease of deployment across different computing environments.

## Overview

The workflow includes the following steps:
1. Quality control with FastQC
2. Trimming of adapter sequences
3. Alignment using HISAT2
4. Conversion, sorting, and indexing of BAM files with SAMtools
5. Gene counting with HTSeq
6. Creation of a count table

## Requirements

- Snakemake
- Singularity
- Slurm (for running on a cluster)

## Directory Structure

- `input/`: Directory for input FASTQ files.
- `output/`: Directory for output files including FastQC reports, trimmed reads, BAM files, and count tables.
- `rules/`: Directory containing Snakemake rule definitions (if applicable).
- `envs/`: Directory containing environment definition files for Singularity (if applicable).

## Running the Workflow

### On the Login Node

To perform a dry run on the login node to see what commands will be executed:

\`\`\`bash
snakemake -n -p
\`\`\`

To execute the workflow on the login node (not recommended for heavy computation):

\`\`\`bash
snakemake --use-singularity
\`\`\`

### On a Slurm Cluster

To submit jobs to the Slurm cluster, use the following command:

\`\`\`bash
snakemake --cluster "sbatch -N 1 -n 1 -t 60" --jobs 1 --use-singularity
\`\`\`

This command will submit the jobs to the cluster, requesting 1 node and 1 task per job, with a runtime limit of 60 minutes per job.

## Singularity Containers

The workflow uses Singularity containers to encapsulate the software dependencies. The containers are verified by Nexus to ensure they are secure and reproducible.

Ensure Singularity is installed and properly configured on your system. You can pull the required container images using:

\`\`\`bash
singularity pull [container_image_url]
\`\`\`

## Configuration

You can configure the samples and other parameters directly in the \`Snakefile\`:

\`\`\`python
SAMPLES = ["A", "B"]
TRIMFILE = "/path/to/adapters.fa"
\`\`\`

Adjust the paths and sample names as necessary for your specific dataset.

## Contributing

Contributions are welcome. Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License.