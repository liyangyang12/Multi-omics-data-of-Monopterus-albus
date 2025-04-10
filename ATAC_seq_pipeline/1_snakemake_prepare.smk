SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep2","MO_rep3","ISE_rep1","ISE_rep2","ISE_rep3","ISM_rep1","ISM_rep2","ISM_rep3","ISL_rep1","ISL_rep2","ISL_rep3","ES_rep1","ES_rep2","ES_rep3","MS_rep1","MS_rep2","MS_rep3","LS_rep1","LS_rep2","LS_rep3","ES_rep4","ES_rep5","MS_rep4","MS_rep5","MS_rep6","MS_rep7","MS_rep8","MS_rep9"]

rule all:
        input:
            expand("results/trim/{sample}_1_val_1.fq.gz",sample=SAMPLE),
            "results/index/malbus.genome.1.bt2",
            expand("results/bowtie2/{sample}.sort.bam",sample=SAMPLE),
            expand("results/bowtie2/{sample}.sort.bam.bai",sample=SAMPLE),
            expand("results/filter/{sample}.tmp2.bam",sample=SAMPLE),
            expand("results/filter/{sample}.final.bam",sample=SAMPLE),
            expand("results/filter/{sample}.rmdup.bam.bai",sample=SAMPLE),  
            expand("results/filter/{sample}.filter.bam",sample=SAMPLE),
            expand("results/filter/{sample}.filter.bam.bai",sample=SAMPLE),
            expand("results/bw/{sample}.bw",sample=SAMPLE),


rule trim:
    input:
        r1 = "data/{sample}_1.fq.gz",
        r2 = "data/{sample}_2.fq.gz"
    output:
        r1 = "results/trim/{sample}_1_val_1.fq.gz",
        r2 = "results/trim/{sample}_2_val_2.fq.gz"
    shell:
        """
        trim_galore -q 25 --phred33 --length 45 -e 0.1 --paired -o results/trim  {input.r1} {input.r2}
        """
rule index:
    input:
        fa = "genome/GCF_001952655.1_M_albus_1.0_genomic.fna"
    output:
        "results/index/malbus.genome.1.bt2"
    threads: 60
    shell:
        """
        bowtie2-build {input.fa} results/index/malbus.genome
        """

rule mapping:
    input:
        r1 = "results/trim/{sample}_1_val_1.fq.gz",
        r2 = "results/trim/{sample}_2_val_2.fq.gz"
    output:
        bam =  "results/bowtie2/{sample}.sort.bam",
        index = "results/bowtie2/{sample}.sort.bam.bai",
        log = "results/bowtie2/{sample}.log",
    params:
        prefiex = "results/bowtie2/{sample}",
    shell:
        """
        bowtie2 -k 5 -X 2000 --mm --threads 5 -x results/index/malbus.genome -1 {input.r1} -2 {input.r2} -S {params.prefiex}.tmp.sam 1>{output.log} 2>&1
        samtools view -bS {params.prefiex}.tmp.sam > {params.prefiex}.tmp.bam
        samtools sort {params.prefiex}.tmp.bam -o {output.bam}
        samtools index {output.bam}
        rm {params.prefiex}.tmp.sam {params.prefiex}.tmp.bam
        """


rule filterbam1_1:
    input:
        bam = "results/bowtie2/{sample}.sort.bam",
    output:
        bam2 = "results/filter/{sample}.tmp2.bam",
    params:
        prefiex = "results/filter/{sample}",
    shell:
        """
        samtools view  -@ 20 -F 524 -f 2 -u {input.bam} > {params.prefiex}.tmp1.bam
        samtools view -h -@ 20 {params.prefiex}.tmp1.bam | scripts/assign_multimappers.py -k 5 --paired-end | samtools sort -n -O bam -@ 20 -o - > {output.bam2}
        rm {params.prefiex}.tmp1.bam
        """

rule filterbam1_2:
    input:
        bam = "results/filter/{sample}.tmp2.bam",
    output:
        bam3 = "results/filter/{sample}.final.bam",
    params:
        prefiex = "results/filter/{sample}",
    shell:
        """
        samtools fixmate -@ 20 -r {input.bam} {params.prefiex}.fixmate.bam
        samtools view -h -@ 20 -F 1804 -f 2 {params.prefiex}.fixmate.bam | samtools sort -O bam -@ 20 -o {output.bam3}
        rm {params.prefiex}.fixmate.bam
        """

rule filterbam2:
    input:
        bam = "results/filter/{sample}.final.bam",
    output:
        log = "results/filter/{sample}.log",
        bai1="results/filter/{sample}.rmdup.bam.bai",
        bam2 = "results/filter/{sample}.filter.bam",
        bai2 = "results/filter/{sample}.filter.bam.bai",
    params:
        prefiex = "results/filter/{sample}",
    shell:
        """
        picard MarkDuplicates -Xmx60g I={input.bam} O={params.prefiex}.rmdup.bam M={output.log} REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT
        samtools index {params.prefiex}.rmdup.bam
        samtools view -h -@ 60 -F 1804 -f 2 {params.prefiex}.rmdup.bam | grep -v "NC_003192.1" | samtools sort -O bam -@ 60 -o - > {output.bam2}
        samtools index {output.bam2}
        rm {params.prefiex}.rmdup.bam 
        """
rule bw:
    input:
        bam = "results/filter/{sample}.filter.bam",
    output:
        bw = "results/bw/{sample}.bw"
    shell:
        """
        bamCoverage --bam {input.bam} -o {output.bw} --binSize 10 --normalizeUsing RPKM --extendReads
        """

rule stat:
    input:
        "results/filter/{sample}.filter.bam"
    output:
        "results/filter/{sample}.stat"
    shell:
        """
        samtools flagstat -@ 5 {input} > {output}
        """
