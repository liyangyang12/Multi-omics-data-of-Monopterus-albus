configfile: "config.yaml"
samples = config["samples"]
indir = "clean"
outdir = "results/prepare"
reads = ["1", "2"]

rule all:
    input:
        expand(outdir + "/fastqc/{sample}_{read}.clean_fastqc.html", sample=samples, read=reads),

rule fastqc:
    input:
        indir + "/{name}.clean.fq.gz"
    output:
        outdir + "/fastqc/{name}.clean_fastqc.html"
    log:
        outdir + "/fastqc/{name}.log"
    params:
        odir = outdir + "/fastqc"
    shell:
        """
        fastqc -o {params.odir} {input} &> {log}
        """
