SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep3","ISE_rep1","ISE_rep2","ISM_rep1","ISM_rep3","ISL_rep1","ISL_rep2","ES_rep3","ES_rep4","MS_rep4","MS_rep6","LS_rep1","LS_rep3"]

rule all:
        input:
            "results/computematrix/results_rep1.pdf",
            expand("results/ATACQC/{sample}.pdf",sample = SAMPLE)
            expand("results/fragment/{sample}_fragment.txt",sample = SAMPLE)



rule atacqc:
    input:
        "results/filter/{sample}.filter.bam"
    output:
        "results/ATACQC/{sample}.pdf"
    shell:
        """
        Rscript scripts/atacqc.R {input} {output}
        """

rule fragment:
    input:
        "../results/filter/{sample}.filter.bam"
    output:
        "results/fragment/{sample}_fragment.txt"
    shell:
        """
        samtools  view {input}| awk -F'\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print $1"\t"abs($9)}}' | sort | uniq | cut -f2 > {output}
        """


rule computematrix:
    input:
        bw =expand("../results/bw/{sample}.bw",sample=SAMPLE),
        bed = "../genome/Malbus.bed12"
    output:
        gz = "results/computematrix/results_rep1.gz",
        tab = "results/computematrix/results_rep1.tab"
    shell:
        """
        computeMatrix scale-regions -p 40 -R {input.bed}  -S {input.bw} -b 2000 -a 2000 --regionBodyLength 5000 --skipZeros -o {output.gz} --outFileNameMatrix {output.tab}
        """

rule plot2:
    input:
        "results/computematrix/results_rep1.gz"
    output:
        "results/computematrix/results_rep1.pdf"
    shell:
        """
        plotProfile -m {input} -out {output} --plotHeight 10 --plotWidth 12 --plotFileFormat pdf --perGroup --colors "#a23124" "#cb3228" "#e86c44" "#f4aa62" "#ecc264" "#a9d3e2" "#72a6c8" "#4470aa" "#1d2e62"
        """