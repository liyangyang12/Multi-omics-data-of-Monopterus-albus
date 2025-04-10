#!/usr/bin/env snakemake
names = ["GASR", "NCBI"]
outdir = "results/function"

rule all:
    input:
        expand(outdir + "/isoforms/{name}.fa", name=names),
        expand(outdir + "/eggnog/{name}.emapper.hits", name=names),
        "results/function/classify/kegg.txt",
        "results/function/gocount/tama_merge.modified_GO.txt",
        "results/function/gocount/Malbus_GO.txt",


rule get_isoform_protein:
    input:
        gtf = "genome/{name}.bed",
        fa = "genome/GCF_001952655.1_M_albus_1.0_genomic.fna"
    output:
        fa = outdir + "/isoforms/{name}.fa"
    shell:
        """
        gffread --in-bed -S {input.gtf} -g {input.fa} -y {output.fa}
        """


rule EggNOG:
    input:
        fa = rules.get_isoform_protein.output.fa
    output:
        ann = outdir + "/eggnog/{name}.emapper.annotations",
        hit = outdir + "/eggnog/{name}.emapper.hits",
        ort = outdir + "/eggnog/{name}.emapper.seed_orthologs"
    log:
        outdir + "/eggnog/{name}.emapper.log"
    params:
        prefix = outdir + "/eggnog/{name}"
    threads:
        40
    shell:
        """
        emapper.py -m diamond -i {input.fa} --output {params.prefix} -d euk --usemem --cpu {threads} &> {log}
        """

rule function_classify:
    input:
        "results/function/eggnog/GASR.emapper.annotations"
    output:
        allgene =  "results/function/classify/alggene.txt",
        GO =  "results/function/classify/GO.txt",
        kegg = "results/function/classify/kegg.txt"
    shell:
        """
        python scripts/eggnog_parser.py --eggnog {input} --resultall {output.allgene} --resultGO {output.GO} --resultKEGG {output.kegg}
        """

rule GOGASR:
    input:
        "results/function/eggnog/GASR.emapper.annotations"
    output:
        GO =  "results/function/gocount/GASR_GO.txt",
    shell:
        """
        python scripts/GO.py --eggnog {input}  --resultGO {output.GO} 
        """


rule GONCBI:
    input:
        anno = "results/function/eggnog/NVCI.emapper.annotations",
        refid = "results/function/refid/ref_id.txt"
    output:
        GO =  "results/function/gocount/NCBI_GO.txt",
    shell:
        """
        python scripts/GO2.py --eggnog {input.anno} --refid {input.refid}  --resultGO {output.GO} 
        """