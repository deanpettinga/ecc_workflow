import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version

### set minimum snakemake version ###
min_version("5.14.0")

configfile: "config.yaml"

units = pd.read_table(config["units"]).set_index("sample", drop=False)

print(units)
rule all:
    input:
                # ref_index:
                expand("refs/hg38.fa.{suffix}", suffix=["amb","ann","bwt","pac","sa"]),
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


# rule download:
#     output:
#                 R1 = "raw_data/unknown_circle_reads_1.fastq",
#                 R2 = "raw_data/unknown_circle_reads_2.fastq",
#                 hg38 =  "refs/hg38.fa",
#     shell:
#         """
#         wget https://raw.githubusercontent.com/iprada/Circle-Map/master/tutorial/unknown_circle_reads_1.fastq {output.R1}
#         wget https://raw.githubusercontent.com/iprada/Circle-Map/master/tutorial/unknown_circle_reads_2.fastq {output.R2}
#         wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz {output.hg38}.gz
#         gunzip -d {output.hg38}.gz
#         """

rule ref_index:
    input:
                config["reference_genome"]
    output:
                # bwa index
                expand("refs/hg38.fa.{suffix}", suffix=["amb","ann","bwt","pac","sa"]),
                # samtools index
                "refs/hg38.fa.fai",
    log:
                bwa = "logs/ref_index.bwa.log",
                samtools = "logs/ref_index.bwa.log",
    conda:
                "envs/circle-map.yaml",
    resources:
                threads = 1,
                nodes =   1,
                mem_gb =  64,
    shell:
        """
        # aligner indexing
        bwa index {input} 2> {log.bwa}

        # Circle-Map indexing
        samtools faidx {input} 2> {log.samtools}
        """

rule align:
    input:
                R1 = "raw_data/{sample}-R1.fastq",
                R2 = "raw_data/{sample}-R2.fastq",
    params:
                ref = "refs/hg38.fa",
                tmp = "tmp"
    output:
                bam_nameSorted = "analysis/bwa/{sample}.nameSorted.bam",
                bai_nameSorted = "analysis/bwa/{sample}.nameSorted.bam.bai",
                bam_coordSorted = "analysis/bwa/{sample}.coordSorted.bam",
                bai_coordSorted = "analysis/bwa/{sample}.coordSorted.bam.bai",
    log:
                bwa = "logs/bwa/bwa.{sample}.log",
                samtools_view = "logs/bwa/samtools-view.{sample}.log",
                samtools_sort = "logs/bwa/samtools-sort.{sample}.log",
                samtools_index = "logs/bwa/samtools-index.{sample}.log",
    conda:
                "envs/circle-map.yaml"
    resources:
                threads = 8,
                nodes =   1,
                mem_gb =  64,
    shell:
        """
        # align with bwa mem
        bwa mem -q -t {resources.threads} {params.ref} {input.R1} {input.R2} 2> {log.bwa} |\
        # pipe sam output to samtools view into bam format
        samtools view -@ {resources.threads} -b - 2> {log.samtools_view} |\
        # nameSort bam
        samtools sort -n -T {params.tmp} -o {output.bam_nameSorted} - 2> {log.samtools_sort}
        # nameSort bam index
        samtools index -b -@ {resources.threads} {output.bam_nameSorted} 2> {log.samtools_index}

        # coordSort bam
        samtools sort -n -T {params.tmp} -o {output.bam_coordSorted} {output.bam_nameSorted}
        # coordSort bam index
        samtools index -b -@ {resources.threads} {output.bam_coordSorted} 2> {log.samtools_index}
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
    params:
                    ref = "refs/hg38.fa",
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
        -fasta {params.ref} \
        -o {output} |
        2> {log}
        """
