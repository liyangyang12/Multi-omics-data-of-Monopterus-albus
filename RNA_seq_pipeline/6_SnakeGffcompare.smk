#!/usr/bin/env snakemake
names = [
    "stringtie_taco",
]

rule all:
    input:
        expand("results/gffcompare/{name}.tracking", name=names),

rule gffcompare:
    input:
        ref = "/DATA/lyy/sex_reversal_project/1_ngs/genome/GCF_001952655.1_M_albus_1.0_genomic.gtf",
        gtf = "results/collect/{name}.gtf"
    output:
        out1 = "results/gffcompare/{name}.loci",
        out2 = "results/gffcompare/{name}.tracking",
        out3 = "results/gffcompare/{name}.annotated.gtf",
        out4 = "results/gffcompare/{name}.stats",
    log:
        log = "results/gffcompare/{name}.log"
    params:
        prefix = "results/gffcompare/{name}"
    shell:
        """
        gffcompare -r {input.ref} -R -o {params.prefix} {input.gtf} &> {log}
        """
