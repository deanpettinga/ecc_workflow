import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version

### set minimum snakemake version ###
min_version("5.14.0")

configfile: "bin/config.yaml"

units = pd.read_table(config["units"]).set_index("sample", drop=False)

rule all:
    input:
                # ref_index:
                expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa"]),
                # align
                expand("analysis/align/{units.sample}.coordSorted.bam", units=units.itertuples()),
                expand("analysis/align/{units.sample}.coordSorted.bam.bai", units=units.itertuples()),
                # circleMap_Repeats
                expand("analysis/circleMap_Repeats/{units.sample}.circleMap_Repeats.bed", units=units.itertuples()),

rule ref_index:
    input:
                config["reference_genome"]
    output:
                # bwa index
                expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa"]),
                # samtools index
                expand("{ref}.fai", ref=config["reference_genome"]),
    log:
                bwa = expand("logs/ref_index/{ref}.bwa_index.log", ref=config["reference_genome"]),
                samtools = expand("logs/ref_index/{ref}.samtools_faidx.log", ref=config["reference_genome"]),
    conda:
                "envs/circle-map.yaml",
    resources:
                threads = 1,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                # Circle-Map indexing
                samtools faidx {input} 2> {log.samtools}

                # aligner indexing
                bwa index {input} 2> {log.bwa}
                """

rule trim_galore:
    input:
                R1 = "raw_data/{sample}-R1.fastq.gz",
                R2 = "raw_data/{sample}-R2.fastq.gz",
    output:
                "analysis/trim_galore/{sample}-R1_val_1.fq.gz",
                "analysis/trim_galore/{sample}-R1_val_1_fastqc.html",
                "analysis/trim_galore/{sample}-R1_val_1_fastqc.zip",
                "analysis/trim_galore/{sample}-R1.fastq.gz_trimming_report.txt",
                "analysis/trim_galore/{sample}-R2_val_2.fq.gz",
                "analysis/trim_galore/{sample}-R2_val_2_fastqc.html",
                "analysis/trim_galore/{sample}-R2_val_2_fastqc.zip",
                "analysis/trim_galore/{sample}-R2.fastq.gz_trimming_report.txt"
    params:
                outdir="analysis/trim_galore/"
    log:
                stdout="logs/trim_galore/{sample}.o",
                stderr="logs/trim_galore/{sample}.e"
    benchmark:
                "benchmarks/trim_galore/{sample}.txt"
    conda:
                "envs/trim-galore.yaml"
    resources:
                threads = 4,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                trim_galore \
                --paired \
                {input} \
                --output_dir {params.outdir} \
                --cores {threads} \
                -q 20 \
                --fastqc \
                2> {log.stderr} 1> {log.stdout}
                """

rule align:
    input:
                ref = config["reference_genome"],
                ref_index = expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa","fai"]),
                R1 = "analysis/trim_galore/{sample}-R1_val_1.fq.gz",
                R2 = "analysis/trim_galore/{sample}-R2_val_2.fq.gz",
    params:
                tmp = "tmp",

    output:
                sam = temp("analysis/align/{sample}.sam"),
                bam_coordSorted = "analysis/align/{sample}.coordSorted.bam",
                bai_coordSorted = "analysis/align/{sample}.coordSorted.bam.bai",
    log:
                bwa = "logs/align/{sample}.bwa.log",
                coordSort = "logs/align/{sample}.coordSort.log",
                coordSort_index = "logs/align/{sample}.coordSort_index.log",

    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                # align with bwa mem (to SAM)
                bwa mem -M -q -t {resources.threads} {input.ref} {input.R1} {input.R2} 2> {log.bwa} 1> {output.sam}
                # coordSort the BAM
                samtools sort -T {params.tmp} -O BAM -o {output.bam_coordSorted} {output.sam} 2> {log.coordSort}
                # coordSort index
                samtools index -b -@ {resources.threads} {output.bam_coordSorted} 2> {log.coordSort_index}
                """

rule circleMap_Repeats:
    input:
                bam_coordSorted = "analysis/align/{sample}.coordSorted.bam",
                bai_coordSorted = "analysis/align/{sample}.coordSorted.bam.bai",
    output:
                "analysis/circleMap_Repeats/{sample}.circleMap_Repeats.bed"
    log:
                "logs/circleMap_Repeats/{sample}.circleMap_Repeats.log",

    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
                "Circle-Map Repeats -i {input.bam_coordSorted} -o {output} 2> {log}"

# rule circleMap_realign:
#     input:
#                 circleMap_coordSorted = "analysis/circle-map/{sample}.circleMapReadExtractor.coordSorted.bam",
#                 bam_nameSorted = "analysis/align/{sample}.nameSorted.bam",
#                 bam_coordSorted = "analysis/align/{sample}.coordSorted.bam",
#                 bai_coordSorted = "analysis/align/{sample}.coordSorted.bam.bai",
#                 ref = config["reference_genome"],
#     output:
#                 "analysis/circle-map/{sample}.circleMap.bed",
#     log:
#                 "logs/circleMap_realign/{sample}.realign.log",
#     conda:
#                 "envs/circle-map.yaml"
#     resources:
#                 threads = 8,
#                 nodes =   1,
#                 mem_gb =  64,
#     shell:
#         """
#         Circle-Map Realign \
#         -i {input.circleMap_coordSorted} \
#         -qbam {input.bam_nameSorted} \
#         -sbam {input.bam_coordSorted} \
#         -fasta {input.ref} \
#         -o {output} |
#         2> {log}
#         """
