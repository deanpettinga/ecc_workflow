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
                # config["reference_genome"]+".fai",
                # # trim_galore
                # expand("analysis/trim_galore/{units.sample}-R1_val_1_fastqc.html", units=units.itertuples()),
                # expand("analysis/trim_galore/{units.sample}-R2_val_2_fastqc.zip", units=units.itertuples()),
                # # align_stats
                # expand("analysis/align/{units.sample}.coordSorted.bam.stats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.idxstats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.bam.flagstat", units=units.itertuples()),
                # # multiqc
                # "analysis/multiqc/multiqc_report.html",
                # "analysis/multiqc/multiqc_report_data/multiqc.log",
                # "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
                # "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
                # ## QC --------------------------------------------------------
                # # align magnaporthe
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.bai", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.stats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.idxstats", units=units.itertuples()),
                # expand("analysis/align/{units.sample}.coordSorted.magna.bam.flagstat", units=units.itertuples()),
                # # plotPCA
                # "analysis/deeptools/multiBamSummary.pca.png",
                # ## ECC_CALLER ------------------------------------------------
                # # ecc_caller_createMapfile
                # "analysis/ecc_caller/mapfile",
                # # ecc_caller_align
                # expand("analysis/ecc_caller/{units.sample}.filtered.sorted.bam", units=units.itertuples()),
                # # call_ecc_regions
                # expand("analysis/ecc_caller/{units.sample}.confirmedsplitreads.bed", units=units.itertuples()),
                # # assign_confidence
                # expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.details.tsv", units=units.itertuples()),
                # # filter_by_conf
                # expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.details.conf.tsv", units=units.itertuples()),
                # expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.conf.bed", units=units.itertuples()),
                # # merge_tech_reps
                # expand("analysis/ecc_caller/{treatment}.merged.bed",treatment=["IF","RC"]),
                # #
                # expand("analysis/ecc_caller/{condition}.merged.shuf.bed",(treatment=["IF","RC"]),
                ## HEATMAPS ----------------------------------------------------
                # # prepare_epi_datasets
                # expand("analysis/epi_marks/{mark}.bw", mark=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                # # get_centromere_bw
                # "analysis/deeptools/centromeres.heatmap.bw",
                # # get_TE_bw
                # "analysis/deeptools/TEs.heatmap.bw",
                # # get_GC_bw
                # "analysis/deeptools/50bps.GC.bw",
                # # plotHeatmap
                expand("analysis/deeptools/{feature}.heatmap.png", feature=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3","TEs","50bps.GC"]),
                # # plotProfile
                expand("analysis/deeptools/both_vs_background.{feature}.profile.png", feature=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3","TEs","50bps.GC","centromere"]),
                ### MOTIF ANALYSIS ---------------------------------------------
                # homer
                expand("analysis/homer/{condition}/knownResults.html", condition=["IF","RC"]),
                # ARS
                expand("analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random{number}", condition=["IF","RC"], number=[1,2,3,4,5,6,7,8,9,10]),

rule download_ref:
    input:
    output:
                gz = temp(config["reference_genome"]+".gz"),
                ref = config["reference_genome"],
    log:
                "logs/download_ref.log",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "download_ref",
    shell:
                """
                wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/oryza_sativa/dna//Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz {output.gz} -o {log}
                gunzip {output.gz}
                """
rule ref_index:
    input:
                ancient(config["reference_genome"]),
    output:
                # bwa index
                expand("{ref}.{suffix}", ref=config["reference_genome"], suffix=["amb","ann","bwt","pac","sa"]),
                # samtools index
                config["reference_genome"]+".fai"
    log:
                bwa = expand("logs/ref_index/{ref}.bwa_index.log", ref=config["reference_genome"]),
                samtools = expand("logs/ref_index/{ref}.samtools_faidx.log", ref=config["reference_genome"]),
    conda:
                "envs/circle-map.yaml",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "ref_index",
    shell:
                """
                # make .fai
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
                "grep '>' {input} | grep chromosome | awk '{{print substr($1,$2)}}' 1> {output} 2> {log}"

rule ecc_caller_align:
    # ecc_caller step 2: align reads to reference
    input:
                ref = config["reference_genome"],
                ref_index = config["reference_genome"]+".fai",
                bwa_index = expand(config["reference_genome"]+".{suffix}", suffix=["amb","ann","bwt","pac","sa"]),
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

rule ecc_caller_callEccRegions:
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
                # and select only the first 3 columns
                awk -v OFS='\t' '{{print $1, $2, $3}}' | sort | uniq > {output.merged}
                """

rule samtools_flagstat:
    input:
                bam = "analysis/ecc_caller/{sample}.filtered.sorted.bam",
    output:
                flagstat = "analysis/ecc_caller/{sample}.filtered.sorted.bam.flagstat",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "samtools_flagstat",
    conda:
                "envs/samtools.yaml",
    shell:
                "samtools flagstat {input.bam} > {output.flagstat}"

rule multiqc:
    input:
                expand("analysis/ecc_caller/{units.sample}.filtered.sorted.bam.flagstat", units=units.itertuples()),
    params:
                "analysis/ecc_caller/",
    output:
                "analysis/multiqc/multiqc_report.html",
    log:
                "logs/multiqc.log"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "multiqc",
    conda:
                "envs/multiqc.yaml",
    shell:
                """
                multiqc -f {params} \
                -o analysis/multiqc \
                -n multiqc_report.html \
                2> {log}
                """

# chip datasets-----------------------------------------------------------------
rule prepare_chip_datasets:
    # downloading data from GEO: GSE79033
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79033
    input:
                config["reference_genome"]+".fai",
    params:
                outdir = "analysis/epi_marks/"
    output:
                bed = expand("analysis/epi_marks/{mark}_rice_leaves.macs14_peaks.bed", mark=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
                wig = expand("analysis/epi_marks/{mark}_rice_leaves.macs14_treat_afterfiting_all.wig", mark=["GSM2084216_H3K9me1","GSM2084217_H3K4ac","GSM2084218_H3K27me3","GSM2084219_H3K27ac","GSM2084220_H3K9ac","GSM2084221_H3K9me3"]),
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
                bw = "analysis/deeptools/{mark}.bw",
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

rule plotHeatmap:
    # use merged IF/RC regions to compare and contrast with epigenetic marks in heatmap
    input:
                beds = expand("analysis/ecc_caller/{condition}.merged.bed",condition=["IF","RC"]),
                mark = "analysis/deeptools/{mark}.bw",
    output:
                matrix = "analysis/deeptools/{mark}.mat.gz",
                heatmap = "analysis/deeptools/{mark}.heatmap.png",
    log:
                "logs/plotHeatmap/{mark}.plotHeatmap.log"
    conda:
                "envs/deeptools.yaml"
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "plotHeatmap",
    shell:
                """
                computeMatrix scale-regions \
                -p {resources.threads} \
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

rule plotProfile_background:
    # use merged IF/RC regions to compare and contrast with epigenetic marks in heatmap
    input:
                IF = "analysis/ecc_caller/IF.merged.bed",
                IF_bkg = "analysis/ecc_caller/IF.merged.shuf.bed",
                RC = "analysis/ecc_caller/RC.merged.bed",
                RC_bkg = "analysis/ecc_caller/RC.merged.shuf.bed",
    output:
                matrix = "analysis/deeptools/both_vs_background.{mark}.mat.gz",
                profile = "analysis/deeptools/both_vs_background.{mark}.profile.png",
    log:
                "logs/plotProfile_background/both_vs_background.{mark}.plotProfile_background.log"
    conda:
                "envs/deeptools.yaml"
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "plotProfile_background",
    shell:
                """
                computeMatrix scale-regions \
                -p {resources.threads} \
                -S {input.mark} \
                -R {input.IF} {input.IF_bkg} {input.RC} {input.RC_bkg} \
                --beforeRegionStartLength 3000 \
                --regionBodyLength 5000 \
                --afterRegionStartLength 3000 \
                --skipZeros \
                -out {output.matrix} \
                > {log}

                plotProfile \
                -m {output.matrix} \
                -out {output.profile} \
                –regionsLabel IF IF_bkg RC RC_bkg \
                2>> {log}
                """

rule get_centromere_bw:
    # centromeres are annotated in nipponbare here: http://rice.plantbiology.msu.edu/annotation_pseudo_centromeres.shtml
    input:
                ref = config["reference_genome"],
                ref_fai = config["reference_genome"]+".fai",
                bacs = config["reference_genome"]+".centromere_clones.tsv",
    output:
                bac_tiling = "refs/rice_r7_all_tiling_path.gff3",
                centromere_bed = config["reference_genome"]+".centromeres.bed",
                centromere_bdg = temp(config["reference_genome"]+".centromeres.bdg"),
                centromere_bdg_sorted = temp(config["reference_genome"]+".centromeres.sorted.bdg"),
                centromere_bdg_sorted_merged = temp(config["reference_genome"]+".centromeres.sorted.merged.bdg"),
                chrom_sizes = temp(config["reference_genome"]+".chromsizes"),
                centromere_bw = config["reference_genome"]+"centromeres.bw",
                centromere_bw_moved = "analysis/deeptools/centromeres.bw",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "get_centromere_annotation",
    conda:
                "envs/bedgraphtobigwig.yaml",
    shell:
                """
                wget -O {output.bac_tiling} http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/rice_r7_all_tiling_path.gff3
                for BAC in `awk 'NR > 1 {{print $3}}' {input.bacs}`;
                    do
                        grep ${{BAC}} {output.bac_tiling} |\
                        awk -v OFS='\t' '{{ split($1, a, /Chr/); print a[2], $4, $5, $7, $9}}' \
                        >> {output.centromere_bed}
                    done

                # convert BED to BDG and merge
                awk -v OFS='\t' '{{print $1,$2,$3,1}}' {output.centromere_bed} > {output.centromere_bdg}
                echo 'extracted bdg'
                sort -k1,1 -k2,2n {output.centromere_bdg} > {output.centromere_bdg_sorted}
                bedtools merge -i {output.centromere_bdg_sorted} -c 1 -o count > {output.centromere_bdg_sorted_merged}
                echo 'bedgraph converted.'

                # make chrom.sizes
                cut -f1,2 {input.ref_fai} | awk -v FS='\t' -v OFS='\t' '{{gsub(/chr/,""); print}}' > {output.chrom_sizes}
                echo 'chrom.sizes made'
                cat {output.chrom_sizes}

                # bedgraph to bw
                bedGraphToBigWig {output.centromere_bdg_sorted_merged} {output.chrom_sizes} {output.centromere_bw}
                cp {output.centromere_bw} {output.centromere_bw_moved}
                """

rule get_TE_bw:
    input:
                ref = config["reference_genome"],
                ref_fai = config["reference_genome"]+".fai",
    output:
                brief_info = temp("refs/all.locus_brief_info.7.0"),
                TE_bed = config["reference_genome"]+".TEs.bed",
                TE_bdg = temp(config["reference_genome"]+".TEs.bdg"),
                TE_bdg_sorted = temp(config["reference_genome"]+".TEs.sorted.bdg"),
                TE_bdg_sorted_merged = temp(config["reference_genome"]+".TEs.sorted.merged.bdg"),
                chrom_sizes = temp(config["reference_genome"]+".chromsizes"),
                TE_bw = config["reference_genome"]+".TEs.bw",
                TE_bw_moved = "analysis/deeptools/TEs.bw",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "get_TE_annotation",
    conda:
                "envs/bedgraphtobigwig.yaml",
    shell:
                """
                wget -O {output.brief_info} http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.locus_brief_info.7.0
                awk -v FS='\t' -v OFS='\t' '{{gsub(/Chr/,""); if ($7 == "Y") print $1, $4, $5, $2, "-", $6, $10}}' {output.brief_info} > {output.TE_bed}
                echo 'refs/all.locus_brief_info.7.0 downloaded.'

                # bed to bedgraph and merge
                awk -v OFS='\t' '{{print $1,$2,$3,1}}' {output.TE_bed} > {output.TE_bdg}
                echo 'extracted bdg'
                sort -k1,1 -k2,2n {output.TE_bdg} > {output.TE_bdg_sorted}
                bedtools merge -i {output.TE_bdg_sorted} -c 1 -o count > {output.TE_bdg_sorted_merged}
                echo 'bedgraph converted.'

                # make chrom.sizes
                cut -f1,2 {input.ref_fai} | awk -v FS='\t' -v OFS='\t' '{{gsub(/chr/,""); print}}' > {output.chrom_sizes}
                echo 'chrom.sizes made'
                cat {output.chrom_sizes}

                # bedgraph to bw
                bedGraphToBigWig {output.TE_bdg_sorted_merged} {output.chrom_sizes} {output.TE_bw}
                cp {output.TE_bw} {output.TE_bw_moved}
                """

rule get_GC_bw:
    input:
                ref = config["reference_genome"],
                ref_fai = config["reference_genome"]+".fai",
    params:
                width = "{width}",
    output:
                chrom_sizes = temp(config["reference_genome"]+".{width}.chromsizes"),
                width_bed = temp(config["reference_genome"]+".{width}bps.bed"),
                GC_bdg = temp(config["reference_genome"]+".{width}bps.GC.bdg"),
                GC_bdg_sorted = config["reference_genome"]+".{width}bps.GC.sorted.bdg",
                GC_bw = config["reference_genome"]+".{width}bps.GC.bw",
                GC_bw_deeptools = "analysis/deeptools/{width}bps.GC.bw",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "get_GC",
    conda:
                "envs/bedgraphtobigwig.yaml",
    shell:
                """
                # make chromsizes
                cut -f1,2 {input.ref_fai} > {output.chrom_sizes}

                # create a BED file with windows of wished width
                bedtools makewindows \
                	-g {output.chrom_sizes} \
                	-w {params.width} \
                	> {output.width_bed}

                # compute GC with bedtools
                bedtools nuc \
	                   -fi {input.ref} \
	                   -bed {output.width_bed} |\
                # pipe into gawk to reformat as BDG
                gawk -v w={params.width} 'BEGIN{{FS="\t"; OFS="\t"}} \
                    {{if (FNR>1) {{print $1,$2,$3,$5}}}}' \
                    > {output.GC_bdg}

                sort -k1,1 -k2,2n {output.GC_bdg} > {output.GC_bdg_sorted}

                # convert BED/BDG to BigWig format
                bedGraphToBigWig {output.GC_bdg_sorted} {output.chrom_sizes} {output.GC_bw}
                cp {output.GC_bw} {output.GC_bw_deeptools}
                """

### Motif Analysis -------------------------------------------------------------

rule homer:
    input:
                bed = ancient("analysis/ecc_caller/{condition}.merged.bed"),
                motif = ancient("bin/my_motifs.motif"),
    params:
                ref = config["reference_genome"],
                size = "given", # 'given' uses to the full sequence provided
    output:
                homer_bed = "analysis/homer/{condition}/{condition}.homer.bed",
                dir = directory("analysis/homer/{condition}"),
                html = "analysis/homer/{condition}/knownResults.html",
    log:
                "logs/homer/{condition}.homer.log",
    conda:
                "envs/homer.yaml",
    resources:
                threads =   20,
                nodes =     1,
                mem_gb =    64,
                name =      "homer",
    shell:
                """
                [ ! -d {output.dir} ] && mkdir -p {output.dir}

                awk -v OFS='\t' '{{print $1, $2, $3, "ecc_"NR}}' {input.bed} > {output.homer_bed}

                findMotifsGenome.pl \
                {output.homer_bed} \
                {params.ref} \
                {output.dir} \
                -size {params.size} \
                -mknown {input.motif} \
                -p {resources.threads} \
                2>> {log}
                """

rule extract_ARS_details:
    input:
                ref = config["reference_genome"],
                bed = "analysis/ecc_caller/{condition}.merged.bed",
    params:
                conf_tsv = "analysis/ecc_caller/{condition}_*.ecccaller_output.renamed.details.conf.tsv",
    output:
                fa = "analysis/ecc_caller/{condition}.merged.fa",
                ARS_fa = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.fa",
                ARS_bed = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.bed",
                ARS_conf = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf",
    log:
                "logs/extract_ARS_details/{condition}.extract_ARS_details.log",
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "extract_ARS_details",
    conda:
                "envs/bedtools.yaml",
    shell:
                """
                bedtools getfasta -fi {input.ref} -bed {input.bed} -fo {output.fa} 2> {log}
                echo "bedtools getfasta COMPLETE." > {log}

                # use grep to extract .fa seqs with ARS motif and write to BED
                grep -B 1 -e '[AT]TTTA[TT][AG]TTT[AT]' {output.fa} > {output.ARS_fa}
                echo "ARS-motif.fa extraction COMPLETE." > {log}
                grep '>' {output.ARS_fa} |\
                awk -v FS=':|-' -v OFS='\t' '{{gsub(/^>/, ""); print $1, $2, $3}}' > {output.ARS_bed}
                echo "ARS-motif.bed extraction COMPLETE." > {log}

                # now get the ecc_caller details from whichever technical rep called
                # this exact eccDNA. these files contain the number of split reads.
                grep -f {output.ARS_bed} {params.conf_tsv} |  awk -v FS=':' -v OFS='\t' '{{gsub(/\.ecccaller_output\.renamed.details\.conf\.tsv/, ""); print $1, $2}}' > {output.ARS_conf}
                echo "ARS-motif.conf extraction COMPLETE." > {log}
                """

rule extract_random_eccs:
    input:
                ARS_conf = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf",
                ecc_conf = expand("analysis/ecc_caller/{units.sample}.ecccaller_output.renamed.details.conf.tsv", units=units.itertuples()),
    params:
                conf_tsv = "analysis/ecc_caller/{condition}_*.ecccaller_output.renamed.details.conf.tsv",
                out = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random",
    output:
                cat_tmp = temp("analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random.conf.cat_tmp"),
                random1 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random1",
                random2 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random2",
                random3 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random3",
                random4 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random4",
                random5 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random5",
                random6 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random6",
                random7 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random7",
                random8 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random8",
                random9 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random9",
                random10 = "analysis/ecc_caller/{condition}.merged.fa.ARS-motif.conf.random10",
    log:
                "logs/extract_random_eccs.{condition}.log"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "extract_random_regions",
    shell:
                """
                # make tmp file with all detected circles by ecc_caller
                cat {params.conf_tsv} > {output.cat_tmp}

                # get number of random circles to sample
                WC=`wc -l {input.ARS_conf}`
                N_LINES=`echo $WC | cut -d' ' -f1`

                for i in {{1..10}}; do
                    sort -R {output.cat_tmp} | tail -n ${{N_LINES}} > {params.out}${{i}} 2> {log}
                done
                """

rule extract_background_regions:
    input:
                ecc_bed = "analysis/ecc_caller/{condition}.merged.bed",
                chromsizes = config["reference_genome"]+".chromsizes",
    output:
                shuf_bed = "analysis/ecc_caller/{condition}.merged.shuf.bed",
    log:
                "logs/extract_background_regions.{condition}.log"
    resources:
                threads =   1,
                nodes =     1,
                mem_gb =    64,
                name =      "extract_background_regions",
    conda:
                "envs/bedtools.yaml",
    shell:
                """
                bedtools shuffle -i {input.ecc_bed} -g {input.chromsizes} -seed 1 > {output.shuf_bed} 2> {log}
                """
