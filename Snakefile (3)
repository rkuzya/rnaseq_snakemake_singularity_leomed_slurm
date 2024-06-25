### import packages

import glob 

import os

import logging


SAMPLES = ["A","B"]

TRIMFILE="/cluster/home/rkuzyakiv/snakemake_test5/adapters.fa"

#localrules: all, download_dataset, run_fastqc, trim, hisat_map, samtools_convert, samtools_sort, samtools_index, count_genes_htseq, create_count_table, upload_dataset 
rule all:
    input:
        expand("input/{sample}.fastq.gz", sample=SAMPLES) +
        expand("output/{sample}_fastqc.html", sample=SAMPLES) +
        expand("output/{sample}_trimmed.fastq.gz", sample=SAMPLES) +
        expand("output/{sample}.sam", sample=SAMPLES) +
        expand("output/{sample}_converted.bam", sample=SAMPLES) +
        expand("output/{sample}_sorted.bam", sample=SAMPLES) +
        expand("output/{sample}_sorted.bam.bai", sample=SAMPLES) +
        expand("output/{sample}_htseq_counts.txt", sample=SAMPLES) +
        expand("output/gene_counts.csv") +
        expand("output/uploaded_to_openbis.txt")


### Download FASTQ file from openBIS with pyBIS         
rule download_dataset:
    output:
        "input/{sample}.fastq.gz"

#    shell: "python3 /cluster/home/rkuzyakiv/snakemake_test5/openbis_pybis_dataset_download.py"
    shell: "python3 /cluster/home/rkuzyakiv/snakemake_test5/openbis_pybis_MULTIPLE_datasets_download.py"

### Running FastQC
rule run_fastqc:
    input: "input/{sample}.fastq.gz"

    output: "output/{sample}_fastqc.html"

    singularity: "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/fastqc:0.11.9--0"
    shell: "fastqc {input} --outdir=output/"


### Running TRIMMOMATIC
rule trim:
    input: "input/{sample}.fastq.gz"
    output: "output/{sample}_trimmed.fastq.gz"
    
    singularity: "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/trimmomatic:0.39--1"
    shell: "trimmomatic SE -phred33 {input} {output} ILLUMINACLIP:"+TRIMFILE+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"


### Running HISAT2 aligner
rule hisat_map:
    input:
        "output/{sample}_trimmed.fastq.gz"
    output:
        "output/{sample}.sam"
        
    singularity: "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/hisat2:2.2.1--h87f3376_5"
    shell:
        "hisat2 -x /cluster/home/rkuzyakiv/snakemake_test5/hisat/grch38/genome -U {input} -S {output}"


### Converting SAM to BAM
rule samtools_convert:
    input:
        "output/{sample}.sam"
    output:
        "output/{sample}_converted.bam"

    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/samtools:1.9--h10a08f8_12"
    shell:
        "samtools view -@ 24 -bS {input} > {output}"


### Sorting BAM file 
rule samtools_sort:
    input:
        "output/{sample}_converted.bam"
    output:
        "output/{sample}_sorted.bam"
    
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/samtools:1.9--h10a08f8_12"
    shell:
        "samtools sort -@ 24 -T ./  -O bam {input} > {output}"


### Indexing BAM file
rule samtools_index:
    input:
        "output/{sample}_sorted.bam"
    output:
        "output/{sample}_sorted.bam.bai"
    
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/samtools:1.9--h10a08f8_12"
    shell:
        "samtools index {input}"


### Counting reads using HTSeq
rule count_genes_htseq:
    input:
        bam="output/{sample}_sorted.bam"
    output:
        counts="output/{sample}_htseq_counts.txt"
    
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/htseq:2.0.5--py39hd5189a5_1"
    params:
        gtf="/cluster/home/rkuzyakiv/snakemake_test5/hg38_annotattion/Homo_sapiens.GRCh38.111.gtf",
        strandedness="reverse", # Adjust as necessary: 'yes', 'no', or 'reverse'
        id_attribute="gene_id", # or use 'gene_id' depending on your preference
        feature_type="exon", # Typically 'exon' for gene-level counts
        mode="union" # 'union', 'intersection-strict', or 'intersection-nonempty'
    
    shell:
        "htseq-count -f bam -r pos -s {params.strandedness} -i {params.id_attribute} -t {params.feature_type} --additional-attr=gene_name -m {params.mode} {input.bam} {params.gtf} > {output.counts}"


### Create read counts table
rule create_count_table:
    input:
        expand("output/{sample}_htseq_counts.txt", sample = SAMPLES)
    output:
        "output/gene_counts.csv"
    run:
#        import pandas
#        import glob

#        filez = glob.glob('output/*.txt')
#        t1 = pandas.read_table(filez[1], header=1)
#        tout = t1.iloc[:,1]
#        for f in filez:
#           t1=pandas.read_table(f, header=1)
#           tout=pandas.concat([tout, t1.iloc[:,2]], axis=1)
#           print(f)

#        tout.to_csv('output/counts.csv')
        
        import pandas as pd
        import glob

        filez = glob.glob('output/*_counts.txt')

        # Read the first file to initialize the 'tout' DataFrame
        t1 = pd.read_table(filez[0], header=1)
        tout = t1.iloc[:, 1]

        # Iterate over all files in 'filez' and concatenate desired columns
        for f in filez:
            t1 = pd.read_table(f, header=1)
            tout = pd.concat([tout, t1.iloc[:, 2]], axis=1)
            print(f)

        # Write the concatenated DataFrame to CSV
        tout.to_csv('output/gene_counts.csv', index=False)

### Upload the FASTQC output files back to openBIS with pyBIS
rule upload_dataset:
    input:
        expand("output/{sample}_fastqc.html", sample = SAMPLES)
    output:
        "output/uploaded_to_openbis.txt"  # Output file to indicate successful upload
        
    shell: 
        """
        python3 /cluster/home/rkuzyakiv/snakemake_test5/openbis_pybis_dataset_upload.py {input} &&
        touch {output}
        """




