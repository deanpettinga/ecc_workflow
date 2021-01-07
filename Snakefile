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
                # # ref_index:
                # expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa"]),
                # # trim_galore
                # expand("analysis/trim_galore/{units.sample}-R1_val_1.fq.gz", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R1_val_1_fastqc.html", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R1_val_1_fastqc.zip", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R1.fastq.gz_trimming_report.txt", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R2_val_2.fq.gz", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R2_val_2_fastqc.html", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R2_val_2_fastqc.zip", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R2.fastq.gz_trimming_report.txt", units=units.itertuples()),
                # # align
                # expand("analysis/align/{units.sample}.coordSorted.bam", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.bai", units=units.itertuples()),
                # # align_stats
                # expand("analysis/align/{units.sample}.coordSorted.bam.stats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.idxstats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.flagstat", units=units.itertuples()),
                # # circleMap_Repeats
                # expand("analysis/circleMap_Repeats/{units.sample}.circleMap_Repeats.bed", units=units.itertuples()),
                # getfasta
                expand("analysis/circleMap_Repeats/{units.sample}.circleMap_Repeats.fa", units=units.itertuples()),
                # multiqc
                # "analysis/multiqc/multiqc_report.html",
                # "analysis/multiqc/multiqc_report_data/multiqc.log",
                # "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
                # QC
                # align magnaporthe
                expand("analysis/align/{units.sample}.coordSorted.magna.bam", units=units.itertuples()),
                expand("analysis/align/{units.sample}.coordSorted.magna.bam.bai", units=units.itertuples()),
                expand("analysis/align/{units.sample}.coordSorted.magna.bam.stats", units=units.itertuples()),
                expand("analysis/align/{units.sample}.coordSorted.magna.bam.idxstats", units=units.itertuples()),
                expand("analysis/align/{units.sample}.coordSorted.magna.bam.flagstat", units=units.itertuples()),
                # plotPCA
                "analysis/deeptools/multiBamSummary.pca.png",

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
                "envs/trim_galore.yaml"
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

rule align_stats:
    input:
                bam_coordSorted = "analysis/align/{sample}.coordSorted.bam",
    output:
                stats = "analysis/align/{sample}.coordSorted.bam.stats",
                idxstats = "analysis/align/{sample}.coordSorted.bam.idxstats",
                flagstat = "analysis/align/{sample}.coordSorted.bam.flagstat",
    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 1,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                samtools stats {input.bam_coordSorted} > {output.stats}
                samtools idxstats {input.bam_coordSorted} > {output.idxstats}
                samtools flagstat {input.bam_coordSorted} > {output.flagstat}
                """

# rule circleMap_Repeats:
#     input:
#                 bam_coordSorted = "analysis/align/{sample}.coordSorted.bam",
#                 bai_coordSorted = "analysis/align/{sample}.coordSorted.bam.bai",
#     output:
#                 "analysis/circleMap_Repeats/{sample}.circleMap_Repeats.bed"
#     log:
#                 "logs/circleMap_Repeats/{sample}.circleMap_Repeats.log",
#     conda:
#                 "envs/circle-map.yaml"
#     resources:
#                 threads = 8,
#                 nodes =   1,
#                 mem_gb =  64,
#     shell:
#                 "Circle-Map Repeats -i {input.bam_coordSorted} -o {output} 2> {log}"

rule getfasta:
    input:
                bed = "analysis/circleMap_Repeats/{sample}.circleMap_Repeats.bed",
                ref = config["reference_genome"],
    output:
                fa = "analysis/circleMap_Repeats/{sample}.circleMap_Repeats.fa",
    log:
                "logs/getfasta/{sample}.getfasta.log"
    conda:
                "envs/bedtools.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                bedtools getfasta \
                -fi {input.ref} \
                -bed {input.bed} \
                -fo {output.fa} \
                2> {log}
                """

rule multiqc:
        input:
                    # trim_galore
                    expand("analysis/trim_galore/{units.sample}-R1_val_1_fastqc.html", units=units.itertuples()),
                    expand("analysis/trim_galore/{units.sample}-R2_val_2_fastqc.html", units=units.itertuples()),
                    expand("analysis/trim_galore/{units.sample}-R1.fastq.gz_trimming_report.txt", units=units.itertuples()),
                    expand("analysis/trim_galore/{units.sample}-R2.fastq.gz_trimming_report.txt", units=units.itertuples()),
                    # align_stats
                    expand("analysis/align/{units.sample}.coordSorted.bam.stats", units=units.itertuples()),
                    expand("analysis/align/{units.sample}.coordSorted.bam.idxstats", units=units.itertuples()),
                    expand("analysis/align/{units.sample}.coordSorted.bam.flagstat", units=units.itertuples()),
        params:
                    "analysis/align/",
                    "analysis/trim_galore/",
        output:
                    "analysis/multiqc/multiqc_report.html",
                    "analysis/multiqc/multiqc_report_data/multiqc.log",
                    "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
                    "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
                    "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
                    "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
        log:
                    "logs/multiqc/multiqc.log",
        conda:
                    "envs/multiqc.yaml"
        resources:
                    threads = 1,
                    nodes =   1,
                    mem_gb =  64,
        shell:
                    """
                    multiqc -f {params} \
                    -o analysis/multiqc \
                    -n multiqc_report.html \
                    2> {log}
                    """


# QC for suspected contamination

# align to magnaporthae to see if it is the source of contamination
rule align_magna:
    input:
                ref = 'refs/Magnaporthe_oryzae.MG8.dna_sm.toplevel.fa',
                ref_index = expand("refs/Magnaporthe_oryzae.MG8.dna_sm.toplevel.fa.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa","fai"]),
                R1 = "analysis/trim_galore/{sample}-R1_val_1.fq.gz",
                R2 = "analysis/trim_galore/{sample}-R2_val_2.fq.gz",
    params:
                tmp = "tmp",
    output:
                sam = temp("analysis/align/{sample}.sam.magna"),
                bam_coordSorted = "analysis/align/{sample}.coordSorted.magna.bam",
                bai_coordSorted = "analysis/align/{sample}.coordSorted.magna.bam.bai",
    log:
                bwa = "logs/align/{sample}.bwa.magna.log",
                coordSort = "logs/align/{sample}.coordSort.magna.log",
                coordSort_index = "logs/align/{sample}.coordSort_index.magna.log",
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

rule align_magna_stats:
    input:
                bam_coordSorted = "analysis/align/{sample}.coordSorted.magna.bam",
    output:
                stats = "analysis/align/{sample}.coordSorted.magna.bam.stats",
                idxstats = "analysis/align/{sample}.coordSorted.magna.bam.idxstats",
                flagstat = "analysis/align/{sample}.coordSorted.magna.bam.flagstat",
    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 1,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                samtools stats {input.bam_coordSorted} > {output.stats}
                samtools idxstats {input.bam_coordSorted} > {output.idxstats}
                samtools flagstat {input.bam_coordSorted} > {output.flagstat}
                """


# make a PCA from the alignments
rule plotPCA:
    input:
                bam = "analysis/align/{sample}.coordSorted.bam",
    output:
                multiBamSummary = "analysis/deeptools/multiBamSummary.npz",
                pca = "analysis/deeptools/multiBamSummary.pca.png",
    conda:
                "envs/deeptools.yaml"
    resources:
                threads = 1,
                nodes =   1,
                mem_gb =  64,
    shell:
                """
                # compute read coverage for full genome
                multiBamSummary bins -p {resources.threads} \
                --bamfiles {input.bam} \
                --smartLabels \
                --outFileName {output}

                plotPCA -in {output.multiBamSummary} \
                -o {output.pca} \
                -T "PCA of read counts"
                """
