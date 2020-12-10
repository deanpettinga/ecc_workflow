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
                expand("analysis/bwa/{units.sample}.nameSorted.bam", units=units.itertuples()),
                expand("analysis/bwa/{units.sample}.nameSorted.bam.bai", units=units.itertuples()),
                expand("analysis/bwa/{units.sample}.coordSorted.bam", units=units.itertuples()),
                expand("analysis/bwa/{units.sample}.coordSorted.bam.bai", units=units.itertuples()),
                # circleMap_readExtractor
                expand("analysis/circle-map/{units.sample}.circleMapReadExtractor.bam", units=units.itertuples()),
                expand("analysis/circle-map/{units.sample}.circleMapReadExtractor.coordSorted.bam", units=units.itertuples()),
                expand("analysis/circle-map/{units.sample}.circleMapReadExtractor.coordSorted.bam.bai", units=units.itertuples()),
                # circleMap_realign
                expand("analysis/circle-map/{units.sample}.circleMap.bed", units=units.itertuples()),

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

rule align:
    input:
                ref = config["reference_genome"],
                ref_index = expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa","fai"]),
                R1 = "raw_data/{sample}-R1.fastq.gz",
                R2 = "raw_data/{sample}-R2.fastq.gz",
    params:
                tmp = "tmp"
    output:
                bam_nameSorted = "analysis/bwa/{sample}.nameSorted.bam",
                bai_nameSorted = "analysis/bwa/{sample}.nameSorted.bam.bai",
                bam_coordSorted = "analysis/bwa/{sample}.coordSorted.bam",
                bai_coordSorted = "analysis/bwa/{sample}.coordSorted.bam.bai",
    log:
                bwa = "logs/bwa/bwa.{sample}.log",
                samtools_sort = "logs/bwa/samtools-sort.{sample}.log",
                samtools_index = "logs/bwa/samtools-index.{sample}.log",
                coordSort_samtools_sort = "logs/bwa/coordSort_samtools-sort.{sample}.log",
                coordSort_samtools_index = "logs/bwa/coordSort_samtools-index.{sample}.log",
    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
        """
        # align with bwa mem
        bwa mem -M -q -t {resources.threads} {input.ref} {input.R1} {input.R2} 2> {log.bwa} |\
        # pipe sam output to samtools view into bam format
        samtools view -b -@ {resources.threads} - |\
        # pipe into nameSort bam
        samtools sort -n -T {params.tmp} -o {output.bam_nameSorted} - 2> {log.samtools_sort}
        # index the nameSort
        samtools index -b -@ {resources.threads} {output.bam_nameSorted} 2> {log.samtools_index}

        # coordSort bam
        samtools sort -T {params.tmp} -o {output.bam_coordSorted} {output.bam_nameSorted} 2> {log.coordSort_samtools_sort}
        # coordSort bam index
        samtools index -b -@ {resources.threads} {output.bam_coordSorted} 2> {log.coordSort_samtools_index}
        """

rule circleMap_readExtractor:
    input:
                bam_nameSorted = "analysis/bwa/{sample}.nameSorted.bam",
    output:
                circleMap = "analysis/circle-map/{sample}.circleMapReadExtractor.bam",
                circleMap_coordSorted = "analysis/circle-map/{sample}.circleMapReadExtractor.coordSorted.bam",
                circleMap_coordSorted_bai = "analysis/circle-map/{sample}.circleMapReadExtractor.coordSorted.bam.bai",
    log:
                circleMap = "logs/circle-map/circleMap.{sample}.log",
                sort = "logs/circle-map/circleMap.{sample}.sort.log",
    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
        """
        # extract circle reads
        Circle-Map ReadExtractor -i {input.bam_nameSorted} -o {output.circleMap} 2> {log.circleMap}
        # sort bam
        samtools sort -@ {resources.threads} -o {output.circleMap_coordSorted} {output.circleMap} 2> {log.sort}
        # index coordSorted
        samtools index {output.circleMap_coordSorted}
        """

rule circleMap_realign:
    input:
                circleMap_coordSorted = "analysis/circle-map/{sample}.circleMapReadExtractor.coordSorted.bam",
                ref = config["reference_genome"],
    output:
                "analysis/circle-map/{sample}.circleMap.bed",
    log:
                "logs/circle-map/circleMap.{sample}.realign.log",
    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
        """
        Circle-Map Realign \
        -i sort_circular_read_candidates.bam \
        -qbam qname_unknown_circle.bam \
        -sbam sorted_unknown_circle.bam \
        -fasta {input.ref} \
        -o {output} |
        2> {log}
        """
