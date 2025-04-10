#!/usr/bin/env snakemake
names = ["cupcake"]
threads = 8

rule all:
    input:
        expand("results/collect/{name}.gtf.gz", name=names),
        expand("results/collect/{name}.isoform.txt", name=names),
        expand("results/collect/{name}.sorted.gtf.gz", name=names),
        expand("results/sqanti3/{name}", name=names),
        expand("results/gffcompare/{name}.tracking", name=names),

def get_input_gtf(name):
    if name == "cupcake":
        return "results/cupcake/collapsed/m54263_190710_094916.subreads/out.collapsed.filtered.sorted.gtf.gz"
    else:
        print(name)
        assert False

rule collect:
    input:
        gff = lambda wildcards: get_input_gtf(wildcards.name)
    output:
        gtf = "results/collect/{name}.gtf.gz",
        tbi = "results/collect/{name}.gtf.gz.tbi"
    shell:
        """
        cp {input.gff} {output.gtf}
        tabix -p gff {output.gtf}
        """

rule isoform_count:
    input:
        gtf = "results/collect/{name}.gtf.gz"
    output:
        txt = "results/collect/{name}.isoform.txt"
    shell:
        """
        zcat {input.gtf} | awk '$3=="transcript"' | sed -E 's/^.*gene_id "(\S*)";.*$/\\1/g' \
            | sort | uniq -c | awk '{{print $2"\\t"$1}}' > {output.txt}
        """

rule sqanti3:
    input:
        gtf = "results/collect/{name}.gtf",
        ref = "data/genome/annotation.gtf",
        fsa = "data/genome/genome.fasta"
    output:
        directory("results/sqanti3/{name}")
    log:
        log = "results/sqanti3/{name}.log"
    params:
        prefix = "results/sqanti3/{name}"
    threads:
        threads
    shell:
        """
        set +u; source activate SQANTI3.env
        ~/software/SQANTI3/sqanti3_qc.py -d {output} -t {threads} \
            --gtf {input.gtf} {input.ref} {input.fsa} &> {log}
        conda deactivate
        """

rule gffcompare:
    input:
        ref = "data/genome/annotation.gtf",
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
