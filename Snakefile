import pandas as pd
configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"])
bioreps = sorted(samples["biorep"].unique())
lanes = sorted(samples["lane"].unique())
lr = ["R1", "R2"]
strains = list(set(samples["strain"]))

# Load output directories from config
OUTPUT_DIRS = config["output_dirs"]

def get_raw_file(wildcards):
    row = samples.loc[
        (samples['strain'] == wildcards.strain) &
        (samples['biorep'] == int(wildcards.biorep)) &
        (samples['lane'] == wildcards.lane)
    ]
    if wildcards.lr == "R1":
        return row["fq1"].iloc[0]
    elif wildcards.lr == "R2":
        return row["fq2"].iloc[0]
    else:
        raise ValueError("Invalid read direction (lr): {}".format(wildcards.lr))

rule all:
    input:
        #expand(f"{OUTPUT_DIRS['qc']}/untrimmed_multiqc/{{strain}}_multiqc_report.html", strain=strains),
        #expand(f"{OUTPUT_DIRS['qc']}/trimmed_multiqc/{{strain}}_multiqc_report.html", strain=strains),
        expand(f"{OUTPUT_DIRS['htseq']}/{{strain}}/{{strain}}_{{biorep}}.count", strain=strains, biorep=bioreps),
        expand(f"{OUTPUT_DIRS['qc']}/fastp_reports/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_fastp.html", strain=strains, biorep=bioreps, lane=lanes)

rule fastp_trim:
    input:
        fq1=lambda wildcards: get_raw_file(wildcards)._replace(lr="R1"),
        fq2=lambda wildcards: get_raw_file(wildcards)._replace(lr="R2")
    output:
        trimmed_fq1=f"{OUTPUT_DIRS['trimmed']}/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_R1_trimmed.fq.gz",
        trimmed_fq2=f"{OUTPUT_DIRS['trimmed']}/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_R2_trimmed.fq.gz",
        html_report=f"{OUTPUT_DIRS['qc']}/fastp_reports/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_fastp.html",
        json_report=f"{OUTPUT_DIRS['qc']}/fastp_reports/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_fastp.json"
    conda:
        "workflow/envs/fastp.yaml"
    threads: 4
    shell:
        "mkdir -p {OUTPUT_DIRS['qc']}/fastp_reports/{{wildcards.strain}} && "
        "fastp --in1 {input.fq1} --in2 {input.fq2} "
        "--out1 {output.trimmed_fq1} --out2 {output.trimmed_fq2} "
        "--html {output.html_report} --json {output.json_report} "
        "--thread {threads}"

rule untrimmed_fastqc:
    input:
        get_raw_file
    output:
        f"{OUTPUT_DIRS['qc']}/untrimmed_fastqc/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_{{lr}}_fastqc.html",
        f"{OUTPUT_DIRS['qc']}/untrimmed_fastqc/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_{{lr}}_fastqc.zip"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        "mkdir -p {OUTPUT_DIRS['qc']}/untrimmed_fastqc/{{wildcards.strain}} && "
        "fastqc -t 1 {input} -o {OUTPUT_DIRS['qc']}/untrimmed_fastqc/{{wildcards.strain}}"

rule trim_galore:
    input:
        get_raw_file
    output:
        f"{OUTPUT_DIRS['trimmed']}/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_{{lr}}_trimmed.fq.gz"
    conda:
        "workflow/envs/rna_seq_processing.yaml"
    threads: 4
    shell:
        "trim_galore -j 4 --paired {input} -o {OUTPUT_DIRS['trimmed']}/{{wildcards.strain}} && "
        "mv {OUTPUT_DIRS['trimmed']}/{{wildcards.strain}}/{{wildcards.strain}}_{{wildcards.biorep}}_{{wildcards.lane}}_R1_val_1.fq.gz "
        "{output} && "
        "mv {OUTPUT_DIRS['trimmed']}/{{wildcards.strain}}/{{wildcards.strain}}_{{wildcards.biorep}}_{{wildcards.lane}}_R2_val_2.fq.gz "
        "{output.replace('R1', 'R2')}"

rule samtools_commands:
    input:
        f"{OUTPUT_DIRS['mapped']}/{{strain}}/{{strain}}_{{biorep}}_Aligned.out.sam"
    output:
        f"{OUTPUT_DIRS['sorted']}/{{strain}}/{{strain}}_{{biorep}}_mapped_sorted_reads.bam",
        f"{OUTPUT_DIRS['sorted']}/{{strain}}/{{strain}}_{{biorep}}_mapped_sorted_reads.bam.bai"
    conda:
        "workflow/envs/rna_seq_processing.yaml"
    threads: 4
    resources:
        mem_gb="30G",
        threads="4"
    shell:
        "samtools view -b {input} | samtools sort - -o {output[0]} && samtools index {output[0]} {output[1]}"

rule htseq_quantification:
    input:
        f"{OUTPUT_DIRS['sorted']}/{{strain}}/{{strain}}_{{biorep}}_mapped_sorted_reads.bam"
    output:
        f"{OUTPUT_DIRS['htseq']}/{{strain}}/{{strain}}_{{biorep}}.count"
    conda:
        "workflow/envs/rna_seq_processing.yaml"
    resources:
        mem_gb="20G"
    shell:
        "htseq-count -r pos -s reverse -t gene -i locus_tag -f bam {input} resources/genomes/h_vol/genomic.gff > {output}"

rule map_trimmed_reads_star:
    input:
        read_group_1=expand(f"{OUTPUT_DIRS['trimmed']}/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_R1_trimmed.fq.gz", lane=lanes, allow_missing=True),
        read_group_2=expand(f"{OUTPUT_DIRS['trimmed']}/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_R2_trimmed.fq.gz", lane=lanes, allow_missing=True)
    output:
        temp(f"{OUTPUT_DIRS['mapped']}/{{strain}}/{{strain}}_{{biorep}}_Aligned.out.sam")
    params:
        read_group_1_string = lambda wildcards, input: ",".join(input.read_group_1),
        read_group_2_string = lambda wildcards, input: ",".join(input.read_group_2)
    threads: 12
    resources:
        mem_gb="20G",
        threads="12"
    conda:
        "workflow/envs/rna_seq_processing.yaml"
    shell:
        "STAR --runThreadN 12 --genomeDir resources/genomes/h_vol/indexes --outFileNamePrefix {OUTPUT_DIRS['mapped']}/{{wildcards.strain}}/{{wildcards.strain}}_{{wildcards.biorep}}_ --readFilesCommand gunzip -c --readFilesIn '{params.read_group_1_string}' '{params.read_group_2_string}'"

rule trimmed_multiqc:
    input:
        expand(f"{OUTPUT_DIRS['qc']}/trimmed_fastqc/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_{{lr}}_trimmed_fastqc.html", biorep=bioreps, lane=lanes, lr=lr, allow_missing=True)
    output:
        f"{OUTPUT_DIRS['qc']}/trimmed_multiqc/{{strain}}_multiqc_report.html"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        "multiqc {OUTPUT_DIRS['qc']}/trimmed_fastqc/{{wildcards.strain}}/* -n {output}"

rule untrimmed_multiqc:
    input:
        expand(f"{OUTPUT_DIRS['qc']}/untrimmed_fastqc/{{strain}}/{{strain}}_{{biorep}}_{{lane}}_{{lr}}_fastqc.html", biorep=bioreps, lane=lanes, lr=lr, allow_missing=True)
    output:
        f"{OUTPUT_DIRS['qc']}/untrimmed_multiqc/{{strain}}_multiqc_report.html"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        "multiqc {OUTPUT_DIRS['qc']}/untrimmed_fastqc/{{wildcards.strain}}/* -n {output}"