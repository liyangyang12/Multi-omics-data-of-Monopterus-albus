threads = 20
indir = "/DATA/lyy/malbus/huwei_tgs/F19FTSCCKF1214_CYPvdhE_HuangShan"
outdir = "results/isoseq"
samples=["m54263_190710_094916.subreads"]

rule all:
    input:
        outdir + "/genome.isoseq.mmi",
        expand(outdir + "/ccs/{sample}.bam", sample=samples),
        expand(outdir + "/demux/{sample}.primer_5p--primer_3p.bam", sample=samples),
        expand(outdir + "/flnc/{sample}.bam", sample=samples),
        expand(outdir + "/clustered/{sample}.bam", sample=samples),
        expand(outdir + "/polished/{sample}.bam", sample=samples),
        expand(outdir + "/aligned/{sample}.bam", sample=samples),
        expand(outdir + "/collapsed/{sample}.gff", sample=samples),


rule ccs:
    input:
        indir + "/{sample}.bam",
    output:
        outdir + "/ccs/{sample}.bam"
    log:
        outdir + "/ccs/{sample}.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/ccs --noPolish --minPasses 1 -j {threads} {input} {output} &> {log}
        """

rule lima:
    input:
        bam = "results/isoseq/ccs/{sample}.bam",
        primer = "/DATA/lyy/malbus/workdir/isoseq3/primers.fasta"
    output:
        outdir + "/demux/{sample}.primer_5p--primer_3p.bam"
    log:
        outdir + "/demux/{sample}.log"
    params:
        bam = outdir + "/demux/{sample}.bam"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/lima --isoseq -j {threads} {input} {params.bam} &> {log}
        """

rule refine:
    input:
        bam = rules.lima.output,
        primer = "/DATA/lyy/malbus/workdir/isoseq3/primers.fasta"
    output:
        outdir + "/flnc/{sample}.bam"
    log:
        outdir + "/flnc/{sample}.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/isoseq3 refine -j {threads} {input} {output} &> {log}
        """

rule cluster:
    input:
        rules.refine.output
    output:
        outdir + "/clustered/{sample}.bam"
    log:
        outdir + "/clustered/{sample}.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/isoseq3 cluster --verbose -j {threads} {input} {output} &> {log}
        """

rule polish:
    input:
        bam1 = rules.cluster.output,
        bam2 = indir + "/{sample}.bam",
        pbi2 = indir + "/{sample}.bam.pbi",
    output:
        outdir + "/polished/{sample}.bam"
    log:
        outdir + "/polished/{sample}.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/isoseq3 polish -j {threads} {input.bam1} {input.bam2} {output} &> {log}
        """

rule mmindex_isoseq:
    input:
        "genome/GCF_001952655.1_M_albus_1.0_genomic.fa"
    output:
        outdir + "/genome.isoseq.mmi"
    log:
        outdir + "/genome.isoseq.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/pbmm2 index -j {threads} --preset ISOSEQ {input} {output} &> {log}
        """

rule align:
    input:
        rules.mmindex_isoseq.output,
        rules.polish.output
    output:
        outdir + "/aligned/{sample}.bam"
    log:
        outdir + "/aligned/{sample}.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/pbmm2 align -j {threads} --preset ISOSEQ --sort {input} {output} &> {log}
        """

rule collapse:
    input:
        bam = rules.align.output,
        ccs = rules.ccs.output
    output:
        outdir + "/collapsed/{sample}.gff"
    log:
        outdir + "/collapsed/{sample}.log"
    threads:
        threads
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/isoseq3 collapse -j {threads} {input.bam} {input.ccs} {output}
        """


rule pbindex:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.pbi"
    shell:
        """
        /DATA/lyy/smrtlink/smrtcmds/bin/pbindex {input}
        """

