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
                # expand("analysis/align/{units.sample}.nameSorted.bam", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.bai", units=units.itertuples()),
                # # # align_stats
                # expand("analysis/align/{units.sample}.coordSorted.bam.stats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.idxstats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.flagstat", units=units.itertuples()),
                # # circleMap_Repeats
                # expand("analysis/circleMap_Repeats/{units.sample}.circleMap_Repeats.bed", units=units.itertuples()),
                # circleMap_ReadExtractor
                # expand("analysis/circleMap_Realign/{units.sample}.circular_read_candidates.unsorted.bam", units=units.itertuples()),
                # expand("analysis/circleMap_Realign/{units.sample}.circular_read_candidates.coordSorted.bam", units=units.itertuples()),
                # expand("analysis/circleMap_Realign/{units.sample}.circular_read_candidates.coordSorted.bam.bai", units=units.itertuples()),
                # circleMap_Realign
                # expand("analysis/circleMap_Realign/{units.sample}.circles.bed", units=units.itertuples()),
                # getfasta
                # expand("analysis/R/circles.collapsed.{condition}.fa", condition=["IF","RC"]),
                # multiqc
                # "analysis/multiqc/multiqc_report.html",
                # "analysis/multiqc/multiqc_report_data/multiqc.log",
                # "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
                # QC
                # # align magnaporthe
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.bai", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.stats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.idxstats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.flagstat", units=units.itertuples()),
                # plotPCA
                # "analysis/deeptools/multiBamSummary.pca.png",
                # ecc_caller_createMapfile
                # "analysis/ecc_caller/mapfile",
                # ecc_caller_align
                # expand("analysis/ecc_caller/{units.sample}.filtered.sorted.bam", units=units.itertuples()),
                # call_ecc_regions
                # expand("analysis/ecc_caller/{units.sample}.confirmedsplitreads.bed", units=units.itertuples()),
                # assign_confidence
                # expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.details.tsv", units=units.itertuples()),
                # filter_by_conf
                expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.details.conf.tsv", units=units.itertuples()),
                expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.conf.bed", units=units.itertuples()),
                # merge_tech_reps
                expand("analysis/ecc_caller/{treatment}.merged.bed",treatment=["IF","RC"]),
                # prepare_epi_datasets
                expand("analysis/epi_marks/{mark}.bw", mark=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                # get_fasta_from_bed
                expand("analysis/meme/{condition}.fa", condition=["RC","IF"]),
                # plotHeatmap
                expand("analysis/deeptools/{mark}.epi_marks.mat.gz", mark=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                expand("analysis/deeptools/{mark}.epi_marks.heatmap.png", mark=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                # get_centromere_bed
                "refs/Oryza_sativa.IRGSP-1.0.dna.centromere_clones.bed",
                # get_TE_annotations
                "refs/Oryza_sativa.IRGSP-1.0.dna.TEs.bed",
                "refs/Oryza_sativa.IRGSP-1.0.dna.TEs.bw",
                # get_wgbs_data
                # expand("analysis/epi_marks/{SRR}_{read}.fastq", SRR=["SRR8427224"], read=[1,2]),
                # # biscuit_align
                # expand("analysis/epi_marks/{SRR}.biscuit.bam", SRR=["SRR8427224"]),
                # expand("analysis/epi_marks/{SRR}.biscuit.bam.bai", SRR=["SRR8427224"]),
                # # biscuit_pileup
                # expand("analysis/epi_marks/{SRR}.biscuit.pileup.vcf.gz", SRR=["SRR8427224"]),
                # biscuit_vcf2bed
                expand("analysis/epi_marks/{SRR}.biscuit.pileup.bed", SRR="SRR8427224"),
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
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "{sample}.ref_index",
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
                "envs/trim_galore.yaml",
    resources:
                threads =   4,
                nodes =     1,
                mem_gb =    64,
                name =      "{sample}.trimgalore",
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


# pierre's ecc_caller rules follow below.
# reference: https://github.com/pierrj/ecc_caller
rule gunzip_reads:
    input:
                R1 = "raw_data/{sample}-R1.fastq.gz",
                R2 = "raw_data/{sample}-R2.fastq.gz",
    output:
                R1 = "raw_data/{sample}-R1.fastq",
                R2 = "raw_data/{sample}-R2.fastq",
    log:
                R1 = "logs/gunzip_reads/{sample}-R1.log",
                R2 = "logs/gunzip_reads/{sample}-R2.log",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "{sample}.gunzip_reads",
    shell:
                """
                gunzip -c {input.R1} 1> {output.R1} 2> {log.R1}
                gunzip -c {input.R2} 1> {output.R2} 2> {log.R2}
                """

rule ecc_caller_createMapfile:
    # ecc_caller step 1: prepare a 'mapfile'
    input:
                ref = config["reference_genome"],
    output:
                "analysis/ecc_caller/mapfile",
    log:
                "logs/ecc_caller/createMapfile.log",
    conda:
                "envs/ecc_caller.yaml",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "ecc_caller_createMapfile",
    shell:
                "grep '>' {input} | grep chromosome | awk '{{print substr($1,2)}}' 1> {output} 2> {log}"

rule ecc_caller_align:
    # ecc_caller step 2: align reads to reference
    input:
                ref = config["reference_genome"],
                bwa_index = expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa"]),
                samtools_index = expand("{ref}.fai", ref=config["reference_genome"]),
                mapfile = "analysis/ecc_caller/mapfile",
                R1 = "raw_data/{sample}-R1.fastq",
                R2 = "raw_data/{sample}-R2.fastq",
    params:
                outname = "analysis/ecc_caller/{sample}",
    output:
                "analysis/ecc_caller/{sample}.filtered.sorted.bam",
    log:
                "logs/ecc_caller/{sample}.ecc_caller_align.log",
    conda:
                "envs/ecc_caller.yaml",
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "{sample}.ecc_caller_align",
    shell:
                """
                export ECC_CALLER_PYTHON_SCRIPTS=envs/ecc_caller/python_scripts

                envs/ecc_caller/generate_bam_file.sh \
                -g {input.ref} \
                -1 {input.R1} \
                -2 {input.R2} \
                -s {params.outname} \
                -t {resources.threads} \
                -m {input.mapfile} \
                &> {log}
                """

rule call_ecc_regions:
    # ecc_caller step 3: call eccDNAs from alignment
    input:
                bam = "analysis/ecc_caller/{sample}.filtered.sorted.bam",
                mapfile = "analysis/ecc_caller/mapfile",
    params:
                sample = "analysis/ecc_caller/{sample}",
    output:
                "analysis/ecc_caller/{sample}.confirmedsplitreads.bed",
    log:
                "logs/ecc_caller/{sample}.call_ecc_regions.log",
    conda:
                "envs/ecc_caller.yaml",
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "{sample}.call_ecc_regions",
    shell:
                """
                export ECC_CALLER_PYTHON_SCRIPTS=envs/ecc_caller/python_scripts

                envs/ecc_caller/call_ecc_regions.sh \
                -m {input.mapfile} \
                -s {params.sample} \
                -t {resources.threads} \
                -b {input.bam} \
                &> {log}
                """

rule assign_confidence:
    # ecc_caller step 3: assign confidence to eccDNA calls
    input:
                mapfile = "analysis/ecc_caller/mapfile",
                bed = "analysis/ecc_caller/{sample}.confirmedsplitreads.bed",
                bam = "analysis/ecc_caller/{sample}.filtered.sorted.bam",
    params:
                sample = "analysis/ecc_caller/{sample}",
    output:
                "analysis/ecc_caller/{sample}.ecccaller_output.renamed.details.tsv",
                "analysis/ecc_caller/{sample}.ecccaller_output.renamed.bed",
    log:
                "logs/ecc_caller/{sample}.assign_confidence.log"
    conda:
                "envs/ecc_caller.yaml"
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "{sample}.assign_confidence",
    shell:
                """
                export ECC_CALLER_PYTHON_SCRIPTS=envs/ecc_caller/python_scripts

                envs/ecc_caller/assign_confidence.sh \
                -m {input.mapfile} \
                -s {params.sample} \
                -t {resources.threads} \
                -b {input.bam} \
                -r {input.bed} \
                &> {log}
                """

rule filter_by_conf:
    # filters the output (bed and tsv) eccDNAs by confidence assigned by ecc_caller.
    # using the color codes provided in the input.tsv awk pulls hconf/conf and removes lowq calls
    input:
                tsv = "analysis/ecc_caller/{sample}.ecccaller_output.renamed.details.tsv",
                bed = "analysis/ecc_caller/{sample}.ecccaller_output.renamed.bed",
    output:
                filtered_tsv = "analysis/ecc_caller/{sample}.ecccaller_output.renamed.details.conf.tsv",
                filtered_bed = "analysis/ecc_caller/{sample}.ecccaller_output.renamed.conf.bed",
    log:
                "logs/ecc_caller/{sample}.filter_by_conf.log"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "filter_by_conf.{sample}",
    shell:
                """
                # filter tsv based on the string denoting quality (exclude lowq)
                awk -v OFS="\t" '$6 ~ /hconf/ {{print}}' {input.tsv} 1> {output.filtered_tsv} 2> {log}

                # then filter the bed based on the RGB values (exclude red).
                # also remove redundant start and end cols (7-8) and append stderr to log.
                awk -v OFS="\t" '$9 ~ /0,255,0/ {{print}}' {input.bed} 1> {output.filtered_bed} 2>> {log}
                """

rule merge_tech_reps:
    input:
                A = "analysis/ecc_caller/{treatment}_{bio_rep,[0-9]}A.ecccaller_output.renamed.conf.bed",
                B = "analysis/ecc_caller/{treatment}_{bio_rep,[0-9]}B.ecccaller_output.renamed.conf.bed",
                C = "analysis/ecc_caller/{treatment}_{bio_rep,[0-9]}C.ecccaller_output.renamed.conf.bed",
    params:
                pct_overlap = "0.7",
    output:
                "analysis/ecc_caller/{treatment}_{bio_rep}.cat.bed",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "merge_tech_reps",
    conda:
                "envs/bedtools.yaml",
    shell:
                """
                cat {input.A} {input.B} {input.C} > {output}
                """

rule merge_bio_reps:
    input:
                rep1 = "analysis/ecc_caller/{treatment}_1.cat.bed",
                rep2 = "analysis/ecc_caller/{treatment}_2.cat.bed",
                rep3 = "analysis/ecc_caller/{treatment}_3.cat.bed",
    params:
                pct_overlap = "0.7",
    output:
                temp1_2 = temp("analysis/ecc_caller/tmp.{treatment}.12.bed"),
                temp2_3 = temp("analysis/ecc_caller/tmp.{treatment}.23.bed"),
                temp1_3 = temp("analysis/ecc_caller/tmp.{treatment}.13.bed"),
                merged = "analysis/ecc_caller/{treatment}.merged.bed",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "merge_bio_reps",
    conda:
                "envs/bedtools.yaml",
    shell:
                """
                # reciprocal overlap
                bedtools intersect -a {input.rep1} -b {input.rep2} -f {params.pct_overlap} -r -wa -wb > {output.temp1_2}
                bedtools intersect -a {input.rep2} -b {input.rep3} -f {params.pct_overlap} -r -wa -wb > {output.temp2_3}
                bedtools intersect -a {input.rep1} -b {input.rep3} -f {params.pct_overlap} -r -wa -wb > {output.temp1_3}

                # now merge them.
                cat {output.temp1_2} {output.temp2_3} {output.temp1_3} |
                # and select only the first 4 columns
                awk -v OFS='\t' '{{print $1, $2, $3, $4}}' | sort | uniq > {output.merged}
                """

# chip datasets
rule prepare_chip_datasets:
    # downloading data from GEO: GSE79033
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79033
    input:
                config["reference_genome"]+".fai",
    params:
                outdir = "analysis/epi_marks/"
    output:
                bed = expand("analysis/epi_marks/{entry}_rice_leaves.macs14_peaks.bed", entry=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                wig = expand("analysis/epi_marks/{entry}_rice_leaves.macs14_treat_afterfiting_all.wig", entry=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                # bw = expand("analysis/epi_marks/{entry}_rice_leaves.macs14_treat_afterfiting_all.wig.bw", entry=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
    log:
                "logs/ecc_caller/prepare_chip_datasets.log"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "prepare_chip_datasets",
    conda:
                "envs/wigtobigwig.yaml"
    shell:
                """
                mkdir -p {params.outdir} 2> {log}
                wget -O {params.outdir}epi_marks.tar https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79033/suppl/GSE79033_RAW.tar  2>> {log}
                tar -xf {params.outdir}epi_marks.tar -C {params.outdir} 2>> {log}
                gunzip {params.outdir}*.gz 2>> {log}

                # remove '^Chr' prefix from chromosome coordinates so beds match genome
                for FILE in {output.bed}
                do
                    awk '{{gsub(/^Chr/, ""); print }}' $FILE > $FILE.new
                    mv $FILE.new $FILE 2>> {log}
                done
                """

rule wig_to_bw:
    input:
                chromsizes = "analysis/epi_marks/{mark}.chromsizes",
                wig = "analysis/epi_marks/{mark}_rice_leaves.macs14_treat_afterfiting_all.wig",
    output:
                chrRemoved = "analysis/epi_marks/{mark}_rice_leaves.macs14_treat_afterfiting_all.wig.chrRemoved",
                bw = "analysis/epi_marks/{mark}.bw",
    log:
                "logs/wig_to_bw/{mark}.wig_to_bw.log",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "wig_to_bw.{mark}",
    conda:
                "envs/wigtobigwig.yaml"
    shell:
                """
                # remove 'Chr prefix from wig file so it is compatible with my reference'
                awk '{{gsub(/Chr/,""); print}}' {input.wig} |
                # and subtract 1 from position to make it 0-based
                awk -v OFS='\t' '{{if ($1=="variableStep" || $1=="track")  print; else if($1 != 1) print $1-1, $2}}' > {output.chrRemoved}
                wigToBigWig {output.chrRemoved} {input.chromsizes} {output.bw} 2>> {log}
                """

rule ecc_bed_to_fasta:
    input:
                bed = "analysis/ecc_caller/{condition}.merged.bed"
    params:
                ref = config["reference_genome"],
    output:
                fa = "analysis/meme/{condition}.fa",
    log:
                "logs/getfasta/{condition}.getfasta.log"
    conda:
                "envs/bedtools.yaml"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "ecc_bed_to_fasta.{condition}",
    shell:
                """
                bedtools getfasta \
                -fi {params.ref} \
                -bed {input.bed} \
                -fo {output.fa} \
                2> {log}
                """

rule plotHeatmap:
    # use merged IF/RC regions to compare and contrast with epigenetic marks in heatmap
    input:
                beds = expand("analysis/ecc_caller/{condition}.merged.bed",condition=["IF","RC"]),
                mark = "analysis/epi_marks/{mark}.bw",
    output:
                matrix = "analysis/deeptools/{mark}.epi_marks.mat.gz",
                heatmap = "analysis/deeptools/{mark}.epi_marks.heatmap.png",
    log:
                "logs/plotHeatmap/{mark}.plotHeatmap.log"
    conda:
                "envs/deeptools.yaml"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "plotHeatmap.{mark}",
    shell:
                """
                computeMatrix scale-regions \
                -S {input.mark} \
                -R {input.beds} \
                --beforeRegionStartLength 3000 \
                --regionBodyLength 5000 \
                --afterRegionStartLength 3000 \
                --skipZeros \
                -out {output.matrix} \
                > {log}

                plotHeatmap \
                -m {output.matrix} \
                -out {output.heatmap} \
                2>> {log}
                """

rule get_centromere_annotation:
    # centromeres are annotated in nipponbare here: http://rice.plantbiology.msu.edu/annotation_pseudo_centromeres.shtml
    input:
                bacs = "refs/Oryza_sativa.IRGSP-1.0.dna.centromere_clones.tsv",
    output:
                bac_tiling = "refs/rice_r7_all_tiling_path.gff3",
                centromere_bed = "refs/Oryza_sativa.IRGSP-1.0.dna.centromere_clones.bed",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "get_centromere_annotation",
    shell:
                """
                wget -O {output.bac_tiling} http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/rice_r7_all_tiling_path.gff3
                for BAC in `awk 'NR > 1 {{print $3}}' {input.bacs}`;
                    do
                        grep ${{BAC}} {output.bac_tiling} |\
                        awk -v OFS='\t' '{{ split($1, a, /Chr/); print a[2], $4, $5, $7, $9}}' \
                        >> {output.centromere_bed}
                    done
                """

rule get_TE_annotation:
    input:
                ref = config["reference_genome"],
    output:
                brief_info = "refs/all.locus_brief_info.7.0",
                TE_bed = "refs/Oryza_sativa.IRGSP-1.0.dna.TEs.bed",
                TE_bedgraph = "refs/Oryza_sativa.IRGSP-1.0.dna.TEs.bedgraph",
                chrom_sizes = config["reference_genome"]+".chromsizes",
                TE_bw = "refs/Oryza_sativa.IRGSP-1.0.dna.TEs.bw",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "get_TE_annotation",
    conda:
                "envs/bedgraphtobigwig.yaml",
    shell:
                """
                # make bed
                wget -O {output.brief_info} http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.locus_brief_info.7.0
                awk -v FS='\t' -v OFS='\t' '{{gsub(/Chr/,""); if ($7 == "Y") print $1, $4, $5, $2, "-", $6, $10}}' {output.brief_info} > {output.TE_bed}

                # bed to bedgraph
                awk '{{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,1}}' {output.TE_bed} |\
                sort -k1,1 -k2,2n - > {output.TE_bedgraph}

                # make chrom.sizes
                samtools faidx {input.ref}
                cut -f1,2 {input.ref}.fai > {output.chrom_sizes}

                # bedgraph to bw
                bedGraphToBigWig {output.TE_bedgraph} {output.chrom_sizes} {output.TE_bw}
                """

### Methylation Analysis ###

rule get_wgbs_reads:
    # download reads from SRA
    input:
    params:
                sra = "{SRR}",
                outdir = "analysis/epi_marks/",
    output:
                sra = temp("analysis/epi_marks/{SRR}/{SRR}.sra"),
                R1 = "analysis/epi_marks/{SRR}_1.fastq",
                R2 = "analysis/epi_marks/{SRR}_2.fastq",
    log:
                "logs/get_wgbs_reads.{SRR}.log"
    conda:
                "envs/wgbs.yaml",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "get_wgbs_reads.{SRR}",
    shell:
                """
                prefetch -O {params.outdir} {params.sra} > {log}
                fastq-dump --outdir {params.outdir} --split-files {output.sra}
                """

rule biscuit_index:
    input:
                ref = config["reference_genome"],
    output:
                expand(config["reference_genome"]+"{suffix}", suffix=[".bis.amb",".bis.ann",".bis.pac",".dau.bwt",".dau.sa",".par.bwt",".par.sa"]),
    log:
                "logs/biscuit_index.log"
    conda:
                "envs/wgbs.yaml",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "biscuit_prep",
    shell:
                """
                biscuit index {input.ref} 2> {log}
                """

rule biscuit_align:
    input:
                R1 = "analysis/epi_marks/{SRR}_1.fastq",
                R2 = "analysis/epi_marks/{SRR}_2.fastq",
                ref = config["reference_genome"],
                ref_index = expand(config["reference_genome"]+"{suffix}", suffix=[".bis.amb",".bis.ann",".bis.pac",".dau.bwt",".dau.sa",".par.bwt",".par.sa"]),
    params:
                tempdir = "analysis/epi_marks/{SRR}.temp_dir",
                sample = "{SRR}"
    output:
                # split and discordant SAMs and heavily clipped reads
                # clipped = "analysis/epi_marks/{SRR}.clipped.fastq",
                # disc = "analysis/epi_marks/{SRR}.disc.sam",
                # split = "analysis/epi_marks/{SRR}.split.sam",
                # duplicate marked, sorted, indexed bam
                bam = "analysis/epi_marks/{SRR}.biscuit.bam",
                bai = "analysis/epi_marks/{SRR}.biscuit.bam.bai",
    log:
                "logs/biscuit_alig.{SRR}.log",
    conda:
                "envs/wgbs.yaml",
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "biscuit_map.{SRR}",
    shell:
                """
                biscuit align -b 1 -t {resources.threads} {input.ref} {input.R1} {input.R2} 2> {log} |\
                samblaster 2>> {log} |\
                samtools sort -@ {resources.threads} -o {output.bam} -O BAM - 2>> {log}

                samtools index {output.bam} 2>> {log}
                """

# rule biscuit_qc:
#     input:
#                 bam = "analysis/epi_marks/{SRR}.biscuit.bam",
#                 ref = config["reference_genome"],
#     params:
#                 sample = "{SRR}",
#                 path_to_assets = ""
#     output:
#
#     log:
#                 "logs/biscuit_map.log"
#     conda:
#                 "envs/wgbs.yaml",
#     resources:
#                 threads =   1,
#                 nodes =     1,
#                 mem_gb =    64,
#                 name =      "biscuit_map.{SRR}",
#     shell:
#                 """
#                 """

rule biscuit_pileup:
    input:
                bam = "analysis/epi_marks/{SRR}.biscuit.bam",
                bai = "analysis/epi_marks/{SRR}.biscuit.bam.bai",
                ref = config["reference_genome"],
    output:
                vcf = temp("analysis/epi_marks/{SRR}.biscuit.pileup.vcf"),
                vcf_gz = "analysis/epi_marks/{SRR}.biscuit.pileup.vcf.gz",
    log:
                "logs/biscuit_pileup.{SRR}.log"
    conda:
                "envs/wgbs.yaml",
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "biscuit_pileup.{SRR}",
    shell:
                """
                biscuit pileup \
                -v 1 \
                -q {resources.threads} \
                -o {output.vcf} \
                {input.ref} {input.bam} 2> {log}

                bgzip {output.vcf}
                tabix -p vcf {output.vcf_gz}
                """

rule biscuit_vcf2bed:
    input:
                vcf_gz = "analysis/epi_marks/{SRR}.biscuit.pileup.vcf.gz",
    params:
                t = "cg",
    output:
                bed = "analysis/epi_marks/{SRR}.biscuit.pileup.bed"
    log:
                "logs/biscuit_pileup.{SRR}.log"
    conda:
                "envs/wgbs.yaml",
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "biscuit_bed.{SRR}",
    shell:
                """
                biscuit vcf2bed -t {params.t} {input.vcf_gz} > {output.bed}
                """

#
# rule plotPCA:
#     input:
#                 bam = expand("analysis/align/{units.sample}.coordSorted.bam", units=units.itertuples()),
#     output:
#                 multiBamSummary = "analysis/deeptools/multiBamSummary.npz",
#                 pca = "analysis/deeptools/multiBamSummary.pca.png",
#     log:
#                 "logs/plotPCA.log"
#     conda:
#                 "envs/deeptools.yaml"
#     resources:
#                 threads =   8,
#                 nodes =     1,
#                 mem_gb =    64,
#                 name =      "plotPCA",
#     shell:
#                 """
#                 # make output dir if it doesnt exist
#                 mkdir -p analysis/deeptools
#                 # compute read coverage for full genome
#                 multiBamSummary bins -p {resources.threads} \
#                 --bamfiles {input.bam} \
#                 --smartLabels \
#                 --outFileName {output.multiBamSummary} \
#                 2> {log}
#
#                 plotPCA -in {output.multiBamSummary} \
#                 -o {output.pca} \
#                 -T "PCA of read counts" \
#                 2>> {log}
#                 """
#
# rule eccDNA_analysis_R:
#     input:
#                 expand("analysis/circleMap_Repeats/{units.sample}.circleMap_Repeats.bed", units=units.itertuples()),
#     params:
#                 Rmd = "analysis/R/eccDNA-analysis.Rmd"
#     output:
#                 "analysis/R/eccDNA-analysis.html",
#                 expand("analysis/R/circles.collapsed.{condition}.bed", condition=["RC","IF"]),
#     log:
#                 "logs/R/eccDNA-analysis-R.log"
#     conda:
#                 "envs/r.yaml"
#     resources:
#                 threads = 1,
#                 nodes =   1,
#                 mem_gb =  64,
#     shell:
#                 """
#                 Rscript -e "rmarkdown::render('{params.Rmd}',output_format='html_document')"
#                 """
#
#
# rule multiqc:
#         input:
#                     # trim_galore
#                     expand("analysis/trim_galore/{units.sample}-R1_val_1_fastqc.html", units=units.itertuples()),
#                     expand("analysis/trim_galore/{units.sample}-R2_val_2_fastqc.html", units=units.itertuples()),
#                     expand("analysis/trim_galore/{units.sample}-R1.fastq.gz_trimming_report.txt", units=units.itertuples()),
#                     expand("analysis/trim_galore/{units.sample}-R2.fastq.gz_trimming_report.txt", units=units.itertuples()),
#                     # align_stats
#                     expand("analysis/align/{units.sample}.coordSorted.bam.stats", units=units.itertuples()),
#                     expand("analysis/align/{units.sample}.coordSorted.bam.idxstats", units=units.itertuples()),
#                     expand("analysis/align/{units.sample}.coordSorted.bam.flagstat", units=units.itertuples()),
#         params:
#                     "analysis/align/",
#                     "analysis/trim_galore/",
#         output:
#                     "analysis/multiqc/multiqc_report.html",
#                     "analysis/multiqc/multiqc_report_data/multiqc.log",
#                     "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
#                     "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
#                     "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
#                     "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
#         log:
#                     "logs/multiqc/multiqc.log",
#         conda:
#                     "envs/multiqc.yaml"
#         resources:
#                     threads =   1,
#                     nodes =     1,
#                     mem_gb =    64,
#                     name =      "multiqc"
#         shell:
#                     """
#                     multiqc -f {params} \
#                     -o analysis/multiqc \
#                     -n multiqc_report.html \
#                     2> {log}
#                     """

# # QC for suspected contamination
# # align to magnaporthae to see if it is the source of contamination
# rule align_magna:
#     input:
#                 ref = 'refs/Magnaporthe_oryzae.MG8.dna_sm.toplevel.fa',
#                 ref_index = expand("refs/Magnaporthe_oryzae.MG8.dna_sm.toplevel.fa.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa","fai"]),
#                 R1 = "analysis/trim_galore/{sample}-R1_val_1.fq.gz",
#                 R2 = "analysis/trim_galore/{sample}-R2_val_2.fq.gz",
#     params:
#                 tmp = "tmp",
#     output:
#                 sam = temp("analysis/align/{sample}.sam.magna"),
#                 bam_coordSorted = "analysis/align/{sample}.coordSorted.magna.bam",
#                 bai_coordSorted = "analysis/align/{sample}.coordSorted.magna.bam.bai",
#     log:
#                 bwa = "logs/align/{sample}.bwa.magna.log",
#                 coordSort = "logs/align/{sample}.coordSort.magna.log",
#                 coordSort_index = "logs/align/{sample}.coordSort_index.magna.log",
#     conda:
#                 "envs/circle-map.yaml"
#     resources:
#                 threads =   8,
#                 nodes =     1,
#                 mem_gb =    64,
#                 name =      "{sample}.align_magna",
#     shell:
#                 """
#                 # align with bwa mem (to SAM)
#                 bwa mem -M -q -t {resources.threads} {input.ref} {input.R1} {input.R2} 2> {log.bwa} 1> {output.sam}
#                 # coordSort the BAM
#                 samtools sort -T {params.tmp} -O BAM -o {output.bam_coordSorted} {output.sam} 2> {log.coordSort}
#                 # coordSort index
#                 samtools index -b -@ {resources.threads} {output.bam_coordSorted} 2> {log.coordSort_index}
#                 """
#
# rule align_magna_stats:
#     input:
#                 bam_coordSorted = "analysis/align/{sample}.coordSorted.magna.bam",
#     output:
#                 stats = "analysis/align/{sample}.coordSorted.magna.bam.stats",
#                 idxstats = "analysis/align/{sample}.coordSorted.magna.bam.idxstats",
#                 flagstat = "analysis/align/{sample}.coordSorted.magna.bam.flagstat",
#     conda:
#                 "envs/circle-map.yaml"
#     resources:
#                 threads =   1,
#                 nodes =     1,
#                 mem_gb =    64,
#                 name =      "align_magna_stats",
#     shell:
#                 """
#                 samtools stats {input.bam_coordSorted} > {output.stats}
#                 samtools idxstats {input.bam_coordSorted} > {output.idxstats}
#                 samtools flagstat {input.bam_coordSorted} > {output.flagstat}
#                 """
