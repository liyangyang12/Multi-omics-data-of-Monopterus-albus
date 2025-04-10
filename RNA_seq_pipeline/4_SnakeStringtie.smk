configfile: "config.yaml"
samples = config["samples"]
indir = "results/mapping/filtered"
outdir = "results/stringtie"

rule all:
    input:
        expand(outdir + "/assembly/{sample}.gtf", sample=samples),
        outdir + "/taco/outdir",
        outdir + "/taco/stringtie_taco.sorted.gtf.gz", 


rule stringtie_assembly:
    input:
        bam = indir + "/{sample}.bam",
        bai = indir + "/{sample}.bam.bai",
    output:
        gtf = outdir + "/assembly/{sample}.gtf",
    log:
        log = outdir + "/assembly/{sample}.log"
    threads:
        10
    shell:
        """
        stringtie {input.bam} -p {threads} -o {output.gtf} -l {wildcards.sample} &> {log}
        """

rule taco_merge:
    input:
        gtfs = expand(rules.stringtie_assembly.output.gtf, sample=samples)
    output:
        txt = outdir + "/taco/filelist.txt",
        out = directory(outdir + "/taco/outdir")
    log:
        log = outdir + "/taco/outdir.log"
    threads:
        8
    shell:
        """
        for f in {input.gtfs}; do echo $f; done > {output.txt}
        taco_run -p {threads} -o {output.out} \
            --gtf-expr-attr FPKM \
            --filter-min-length 200 \
            --filter-min-expr 1.0 \
            --isoform-frac 1 {output.txt} &> {log}
        """

rule taco_gtf_index:
    input:
        outdir + "/taco/outdir"
    output:
        gtf = outdir + "/taco/stringtie_taco.sorted.gtf.gz",
        tbi = outdir + "/taco/stringtie_taco.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -header -i {input}/assembly.gtf | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

