#!/usr/bin/env snakemake
samples = ["m54263_190710_094916.subreads"]
threads = 20

rule all:
    input:
        # Mapping
        expand("results/cupcake/mapping/{sample}.hq.fasta", sample=samples),
        "results/cupcake/genome.mmi",
        expand("results/cupcake/mapping/{sample}.sam", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.sam", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.bam", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.bam.bai", sample=samples),
        expand("results/cupcake/mapping/{sample}.sorted.stats", sample=samples),
       
rule get_hq_fasta:
    input:
        fsa = "results/isoseq/polished/{sample}.hq.fasta.gz"
    output:
        fsa = "results/cupcake/mapping/{sample}.hq.fasta"
    shell:
        """
        gzip -d -c {input.fsa} > {output.fsa}
        """

rule minimap2_index:
    input:
        fsa = "genome/GCF_001952655.1_M_albus_1.0_genomic.fa"
    output:
        mmi = "results/cupcake/genome.mmi"
    threads:
        threads
    shell:
        """
        minimap2 -x splice -d {output.mmi} {input.fsa}
        """

rule minimap2_align:
    input:
        mmi = "results/cupcake/genome.mmi",
        fsa = "results/cupcake/mapping/{sample}.hq.fasta"
    output:
        sam = "results/cupcake/mapping/{sample}.sam"
    log:
        log = "results/cupcake/mapping/{sample}.log"
    threads:
        threads
    shell:
        """
        minimap2 -ax splice -t {threads} --secondary=no -C5 -O6,24 -B4 {input.mmi} {input.fsa} > {output.sam} 2> {log}
        """

rule sort_sam:
    input:
        sam = "results/cupcake/mapping/{sample}.sam"
    output:
        sam = "results/cupcake/mapping/{sample}.sorted.sam"
    shell:
        """
        cat {input.sam} | grep '^@' > {output.sam}
        cat {input.sam} | grep -v '^@' | sort -k3,3 -k4,4n >> {output.sam}
        """

rule sam_to_bam:
    input:
        sam = "{prefix}.sam"
    output:
        tmp = temp("{prefix}.tmp.bam"),
        bam = "{prefix}.bam"
    shell:
        """
        samtools view -b {input.sam} > {output.tmp}
        samtools sort {output.tmp} > {output.bam}
        """

rule bam_index:
    input:
        bam = "{prefix}.bam"
    output:
        bai = "{prefix}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule bam_stats:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input.bam} > {output.txt}
        """
