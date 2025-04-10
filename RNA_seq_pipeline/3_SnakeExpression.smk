configfile: "config.yaml"
samples = config["samples"]
threads = 8
indir = "results/mapping/filtered"
outdir = "results/expression"

REPS = {"O_I":["results/mapping/filtered/O_I1.bam"],"O_IV":["results/mapping/filtered/O_IV1.bam","results/mapping/filtered/O_IV2.bam","results/mapping/filtered/O_IV3.bam"],"O_V":["results/mapping/filtered/O_V1.bam","results/mapping/filtered/O_V2.bam","results/mapping/filtered/O_V3.bam"],"OT_I":["results/mapping/filtered/OT_I1.bam","results/mapping/filtered/OT_I2.bam"],"OT_II":["results/mapping/filtered/OT_II1.bam","results/mapping/filtered/OT_II2.bam"],"OT_III":["results/mapping/filtered/OT_III1.bam","results/mapping/filtered/OT_III2.bam","results/mapping/filtered/OT_III3.bam"],"T_I":["results/mapping/filtered/T_I1.bam","results/mapping/filtered/T_I2.bam","results/mapping/filtered/T_I3.bam"],"T_III":["results/mapping/filtered/T_III1.bam","results/mapping/filtered/T_III2.bam"],"T_IV":["results/mapping/filtered/T_IV1.bam","results/mapping/filtered/T_IV2.bam","results/mapping/filtered/T_IV3.bam"]}
SAMPLE = ["O_I","O_IV","O_V","OT_I","OT_II","OT_III","T_I","T_III","T_IV"]
diff = []
for i in range(0,len(SAMPLE)-1):
    for j in range(i+1,len(SAMPLE)):
        diff.append(SAMPLE[i] + '_vs_' + SAMPLE[j])

print(diff)
number = {}
for n in SAMPLE:
    number[n] = len(REPS[n])

rule all:
    input:
        outdir + "/featureCount/counts.txt",
        expand(outdir + "/featureCountpair/{diffname}.txt",diffname = diff),
        expand(outdir + "/deseq2/{diffname}.csv",diffname = diff),
        expand(outdir + "/deseq2/table/{diffname}.flag",diffname = diff)

rule featureCount:
    input:
        bam = expand(indir + "/{sample}.bam",sample=samples),
        gtf = "/beegfs/zhoulab/lyy/sex_reversal_project/1_ngs/genome/GCF_001952655.1_M_albus_1.0_genomic.noMT.gff"
    output:
        txt = outdir + "/featureCount/counts.txt"
    log:
        outdir + "/featureCount/count.log"
    threads:
        threads
    shell:
        """
        featureCounts -T {threads} -p -C -B -t exon -g gene -a {input.gtf} -o {output} {input.bam} >{log} 2>&1
        """

rule featurecountpair:
    input:
        bam1 = lambda wildcards: REPS[wildcards.sample1],
        bam2  =lambda wildcards: REPS[wildcards.sample2],
        gtf = "genome/GCF_001952655.1_M_albus_1.0_genomic.noMT.gff",
    output:
        txt = outdir + "/featureCountpair/{sample1}_vs_{sample2}.txt"
    log:
        outdir + "/featureCountpair/{sample1}_vs_{sample2}.log"
    threads:
        threads
    shell:
        """
        featureCounts -T {threads} -p -C -B -t exon -g gene -a {input.gtf} -o {output} {input.bam1} {input.bam2}>{log} 2>&1
        """   

rule deseqdiff:
    input:
        txt = "results/expression/featureCountpair/{sample1}_vs_{sample2}.txt",
    output:
        outdir + "/deseq2/{sample1}_vs_{sample2}_matrix.txt",
        outdir + "/deseq2/{sample1}_vs_{sample2}.csv",
    params:
        length1 = lambda wildcards: number[wildcards.sample1],
        length2 = lambda wildcards: number[wildcards.sample2],
    shell:
        """
        tail -n +2 {input.txt} | awk '{{$2=null;$3=null;$4=null;$5=null;$6=null;print $0}}' > {output[0]}
        Rscript scripts/diffExp.R {output[0]} {output[1]} {params.length1} {params.length2}
        """

rule filterDEG:
    input:
        outdir + "/deseq2/{sample1}_vs_{sample2}.csv"
    output:
        outdir + "/deseq2/{sample1}_vs_{sample2}_filter.txt"
    shell:
        """
        python scripts/filterDEG.py --deg {input} --result {output}
        """

rule table:
    input:
        outdir + "/deseq2/{sample1}_vs_{sample2}_filter.txt"
    output:
        flag = outdir + "/deseq2/table/{sample1}_vs_{sample2}.flag",
    shell:
        """
        touch {output.flag}
        echo {wildcards.sample1}_vs_{wildcards.sample2} "`cat {input} |grep "Up-reg" |wc -l`" Up-reg >>{output.flag}
        echo {wildcards.sample1}_vs_{wildcards.sample2} "`cat {input} |grep "Down-reg" |wc -l`" Down-reg >>{output.flag}
        """
