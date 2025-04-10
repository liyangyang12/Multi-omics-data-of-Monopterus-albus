#!/usr/bin/env snakemake
names = ["ncbi", "ngs", "tgs"]
pairs = []
for name1 in names:
    for name2 in names:
        if name1 == name2:
            continue
        pairs.append("%s_vs_%s" % (name1, name2))
print(pairs)

indir = "results/merged/sources"
outdir = "results/sqanti"
rule all:
    input:
        expand(outdir + "/{pair}", pair=pairs),


rule sqanti3:
    input:
        gtf1 = indir + "/{name1}.gtf.gz",
        gtf2 = indir + "/{name2}.gtf.gz",
        fsa = "genome/GCF_001952655.1_M_albus_1.0_genomic.fa"
    output:
        directory(outdir + "/{name1}_vs_{name2}")
    log:
        log = outdir + "/{name1}_vs_{name2}.log"
    params:
        prefix = outdir + "/{name1}_vs_{name2}"
    threads:
        8
    shell:
        """
        mkdir {output}
        zcat {input.gtf1} > {output}/reference.gtf
        zcat {input.gtf2} > {output}/query.gtf
        set +u; source activate SQANTI3.env
        software/SQANTI3-4.0/sqanti3_qc.py -d {output} -t {threads} \
            --gtf {output}/query.gtf {output}/reference.gtf {input.fsa} &> {log}
        conda deactivate
        """
        

