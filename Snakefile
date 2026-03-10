"""
A Snakemake workflow to process ChIP-Seq data
"""

import json
import os
import re
import datetime
import sys

##############################################################################
# Load required variables
##############################################################################

configfile: "./config.yaml"

DIR_IN = config.get("input_dir")
DIR_OUT = config.get("output_dir")
try:
   os.mkdir(DIR_OUT)
except:
   pass

# Extract sample names
pattern = r'^(.+)_R[12]\.fastq.gz$'
FASTQS = [f for f in os.listdir(DIR_IN) if re.match(pattern,f)]
SAMPLES = set([re.match(pattern,f).group(1) for f in FASTQS])

# Load genome
bwa2_index = config.get("bwa2_index")
effective_genome_size = config.get("effective_genome_size")
blacklist = config.get("blacklist")

##############################################################################
# Generate rules outputs
##############################################################################

# Save config settings for each run
DIR_LOGS = os.path.join(DIR_OUT, "logs")
v = datetime.datetime.now()
run_date = v.strftime("%Y.%m.%d")
CONFIG = [os.path.join(DIR_LOGS, "config_" + run_date + ".json")]

# Quality control
DIR_QC = os.path.join(DIR_OUT, "qc")
try:
   os.mkdir(DIR_QC)
except:
   pass

MULTI_QC = [os.path.join(DIR_QC, "multiqc_report.html")]
FAST_QC = expand(
   os.path.join(DIR_QC, "{sample}_R{read}_fastqc.html"),
   sample=SAMPLES,
   read=[1,2])

# Trimming
TRIMGALORE = expand(
   os.path.join(DIR_OUT, "trimfq", "{sample}_R{read}_val_{read}.fq.gz"),
   sample=SAMPLES,
   read=[1,2])

# Alignment
DIR_ALIGN = os.path.join(DIR_OUT, "bwamem2")
try:
   os.mkdir(DIR_ALIGN)
except:
   pass

BWAMEM2 = expand(
   os.path.join(DIR_ALIGN, "{sample}_align.sam"),
   sample=SAMPLES)

# SAM Filtering and Sorting
BAMFILES = expand(
   os.path.join(DIR_ALIGN, "{sample}_filt.sort.bam"),
   sample=SAMPLES)

# Peak calling
CHIPS = [i for i in SAMPLES if "input" not in i]
MACS3 = expand(
   os.path.join(DIR_OUT, "macs3", "{chip}_peaks.broadPeak"),
   chip=CHIPS)

# Binning and fold-change
DIR_BW = os.path.join(DIR_OUT, "bigwigs")
try:
   os.mkdir(DIR_BW)
except:
   pass

binsizes=config.get("binsizes")

BAMCOV = expand(
   os.path.join(DIR_BW, "{sample}_bs{bs}.bigwig"),
   sample=SAMPLES,
   bs=binsizes)

BAMCOM = expand(
   os.path.join(DIR_BW, "{chip}_bs{bs}.log2fc.bigwig"),
   chip=CHIPS,
   bs=binsizes)

ALL_BIGWIGS = BAMCOV+BAMCOM
ALL_BIGWIGS = [b.split("/")[-1] for b in ALL_BIGWIGS]
ALL_BEDGRAPHS = [b.replace(".bigwig", ".bedgraph") for b in ALL_BIGWIGS]

BWTOBG = expand(
   os.path.join(DIR_OUT, "bedgraphs", "{bedgraphs}"),
   bedgraphs=ALL_BEDGRAPHS)

# Rule all input
FINAL = \
   FAST_QC + MACS3 + BAMCOV + BAMCOM + BWTOBG

##############################################################################
# RULE ALL
##############################################################################

rule all:
   input:
      FINAL

rule clean:
   params:
      dir_align=DIR_ALIGN,
      dir_trim=TRIMGALORE
   shell:
      '''
      #  Delete all alignment files except filt.sort.bam and filt.sort.bam.bai
      find {params.dir_align} -type f ! -name "*.filt.sort.bam*" -delete

      # Delete trimmed fastq files
      rm -r {params.dir_trim}
      '''

##############################################################################
# FASTQC and trimming
##############################################################################

rule fastqc:
   input:
      os.path.join(DIR_IN, "{sample}_R{read}.fastq.gz")
   output:
      os.path.join(DIR_QC, "{sample}_R{read}_fastqc.html")
   params:
      outdir=DIR_QC
   shell:
      '''
      fastqc -o {params.outdir} {input}
      '''

rule trimgalore:
   input:
      read1=os.path.join(DIR_IN, "{sample}_R1.fastq.gz"),
      read2=os.path.join(DIR_IN, "{sample}_R2.fastq.gz")
   output:
      trim1=os.path.join(DIR_OUT, "trimfq", "{sample}_R1_val_1.fq.gz"),
      trim2=os.path.join(DIR_OUT, "trimfq", "{sample}_R2_val_2.fq.gz")
   threads: 4
   resources:
      mem_mb=32000,
      runtime=480
   params:
      outdir=os.path.join(DIR_OUT, "trimfq")
   shell:
      '''
      trim_galore --paired -j 4 -o {params.outdir} {input.read1} {input.read2}
      '''

##############################################################################
# Alignment and filtering
##############################################################################

rule bwamem2:
   input:
      trim1=os.path.join(DIR_OUT, "trimfq", "{sample}_R1_val_1.fq.gz"),
      trim2=os.path.join(DIR_OUT, "trimfq", "{sample}_R2_val_2.fq.gz")
   output:
      os.path.join(DIR_ALIGN, "{sample}_align.sam")
   threads: 24
   resources:
      mem_mb=96000,
      runtime=480
   params:
      bwa2_index=bwa2_index
   shell:
      '''
      bwa-mem2 mem -t {threads} \
         {params.bwa2_index} {input.trim1} {input.trim2} > {output}
      '''

rule sort_markdupes:
   input:
      os.path.join(DIR_ALIGN, "{sample}_align.sam")
   output:
      sorted=os.path.join(DIR_ALIGN, "{sample}_qsort.bam"),
      marked=os.path.join(DIR_ALIGN, "{sample}_markdupe.sort.bam"),
      metrics=os.path.join(DIR_QC, "{sample}_markdup.metrics.txt")
   threads: 16
   resources:
      mem_mb=64000,
      runtime=480
   shell:
      '''
      samtools sort -@ {threads} -N {input} | \
      samtools fixmate -m - {output.sorted}
      samtools sort -@ {threads} {output.sorted} | \
      samtools markdup -@ {threads} -s -f {output.metrics} - {output.marked}
      '''

rule filter_bam:
   input:
       os.path.join(DIR_ALIGN, "{sample}_markdupe.sort.bam")
   output:
      bam=os.path.join(DIR_ALIGN, "{sample}_filt.sort.bam"),
      bai=os.path.join(DIR_ALIGN, "{sample}_filt.sort.bam.bai"),
   params:
      tempbam=os.path.join(DIR_ALIGN, "{sample}.temp.bam"),
      filtered=os.path.join(DIR_ALIGN, "{sample}.mapq.sub10.bam"),
      blacklist=blacklist
   threads: 16
   resources:
      mem_mb=64000,
      runtime=480
   shell:      
      '''
      samtools view -b -F 1796 -q 0 -@ {threads} -U {params.filtered} -o {params.tempbam} {input}
      bedtools intersect -v -abam {params.tempbam} -b {params.blacklist} > {output.bam}
      samtools index {output.bam}
      '''

##############################################################################
# Peak calling and binning
##############################################################################

rule macs3:
   input:
      control=os.path.join(DIR_ALIGN, "{sample}_input_filt.sort.bam"),
      chip=os.path.join(DIR_ALIGN, "{sample}_{chip}_filt.sort.bam")
   output:
      os.path.join(DIR_OUT, "macs3", "{sample}_{chip}_peaks.broadPeak")
   params:
      name="{sample}_{chip}",
      outdir=os.path.join(DIR_OUT, "macs3"),
      g=effective_genome_size,
      min_length=config.get("min_length"),
      max_gap=config.get("max_gap")
   wildcard_constraints:
      chip="(?!input).*"
   resources:
      mem_mb=64000,
      runtime=240
   shell:
      '''
      macs3 callpeak \
         --broad \
         --bdg \
         -t {input.chip} \
         -c {input.control} \
         -f BAMPE \
         -n {params.name} \
         -g {params.g} \
         --min-length {params.min_length} \
         --max-gap {params.max_gap} \
         --outdir {params.outdir} \
         -q 0.01 \
         --keep-dup auto 
      '''

rule bamcoverage:
   input:
      os.path.join(DIR_ALIGN, "{sample}_filt.sort.bam")
   output:
      os.path.join(DIR_BW, "{sample}_bs{bs}.bigwig")
   wildcard_constraints:
      bs="|".join(map(str, binsizes))
   params:
      binsize=lambda wildcards: wildcards.bs
   threads: 16
   resources:
      mem_mb=64000,
      runtime=480
   shell:
      '''
      bamCoverage -b {input} \
                  -o {output} \
                  -of bigwig \
                  -bs {params.binsize} \
                  -p max
      '''

rule bamcompare:
   input:
      control=os.path.join(DIR_ALIGN, "{sample}_input_filt.sort.bam"),
      chip=os.path.join(DIR_ALIGN, "{sample}_{chip}_filt.sort.bam")
   output:
      os.path.join(DIR_BW, "{sample}_{chip}_bs{bs}.log2fc.bigwig")
   wildcard_constraints:
      chip="(?!input).*",
      bs="|".join(map(str, binsizes))
   params:
      binsize=lambda wildcards: wildcards.bs
   threads: 16
   resources:
      mem_mb=64000,
      runtime=480
   shell:
      '''
      bamCompare -b1 {input.chip} \
                 -b2 {input.control} \
                 -o {output} \
                 -of bigwig \
                 -bs {params.binsize} \
                 -p max
      '''

rule bigwigtobedgraph:
   input:
      os.path.join(DIR_BW, "{fname}.bigwig")
   output:
      os.path.join(DIR_OUT, "bedgraphs", "{fname}.bedgraph")
   resources:
      mem_mb=32000,
      runtime=480
   shell:
      '''
      bigWigToBedGraph {input} {output}
      '''
