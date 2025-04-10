SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep2","MO_rep3","ISE_rep1","ISE_rep2","ISE_rep3","ISM_rep1","ISM_rep2","ISM_rep3","ISL_rep1","ISL_rep2","ISL_rep3","ES_rep1","ES_rep2","ES_rep3","MS_rep1","MS_rep2","MS_rep3","LS_rep1","LS_rep2","LS_rep3","ES_rep4","ES_rep5","MS_rep4","MS_rep5","MS_rep6","MS_rep7","MS_rep8","MS_rep9"]
SAMPLE2=["VO","MO","ISE","ISM","ISL","ES","MS","LS"]

rule all:
        input:
            expand("results/bed/{sample}_shift.bed",sample=SAMPLE),
            expand("results/peak/{sample}/{sample}_peaks.xls",sample=SAMPLE),
            expand("results/peaksort/{sample}_peaks.sort.narrowPeak",sample=SAMPLE),
            expand("results/idr1_2/{sample}.txt",sample=SAMPLE2),
            expand("results/idr2_3/{sample}.txt",sample=SAMPLE2),
            expand("results/idr1_3/{sample}.txt",sample=SAMPLE2),
           
           

rule bamtobed:
    input:
        "results/filter/{sample}.filter.bam"
    output:
        "results/bed/{sample}.bed"
    shell:
        """
        bedtools bamtobed -i {input}  > {output}
        """

rule shift:
    input:
        "results/bed/{sample}.bed"
    output:
        "results/bed/{sample}_shift.bed"
    shell:
        """
        cat {input} | awk -F $'\t' 'BEGIN {{OFS = FS}}{{ if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} {{print $0}}}}' >| {output}
        """

rule peak:
    input:
        "results/bed/{sample}_shift.bed"
    output:
        "results/peak/{sample}/{sample}_peaks.xls",
        "results/peak/{sample}/{sample}_peaks.narrowPeak"
    params:
        prefix ="results/peak/{sample}"
    shell:
        """
        macs2 callpeak -t {input} -f BED -n {wildcards.sample} -g 684144148 -q 0.01 --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits --outdir {params.prefix}
        """

rule peak_sort:
    input:
        "results/peak/{sample}/{sample}_peaks.narrowPeak"
    output:
        "results/peaksort/{sample}_peaks.sort.narrowPeak"
    shell:
        """
        sort -k8,8nr {input} > {output}
        """

rule idr1_2:
    input:
        peak1 = "results/peaksort/{sample}_rep1_peaks.sort.narrowPeak",
        peak2 = "results/peaksort/{sample}_rep2_peaks.sort.narrowPeak",
    output:
        "results/idr1_2/{sample}.txt"
    log:
        "results/idr1_2/{sample}.log"
    shell:
        """
        idr --samples {input.peak1} {input.peak2} --input-file-type narrowPeak --output-file {output} --rank p.value --plot --use-best-multisummit-IDR --log-output-file {log}
        """

rule idr2_3:
    input:
        peak1 = "results/peaksort/{sample}_rep2_peaks.sort.narrowPeak",
        peak2 = "results/peaksort/{sample}_rep3_peaks.sort.narrowPeak",
    output:
        "results/idr2_3/{sample}.txt"
    log:
        "results/idr2_3/{sample}.log"
    shell:
        """
        idr --samples {input.peak1} {input.peak2} --input-file-type narrowPeak --output-file {output} --rank p.value --plot --use-best-multisummit-IDR --log-output-file {log}
        """

rule idr1_3:
    input:
        peak1 = "results/peaksort/{sample}_rep1_peaks.sort.narrowPeak",
        peak2 = "results/peaksort/{sample}_rep3_peaks.sort.narrowPeak",
    output:
        "results/idr1_3/{sample}.txt"
    log:
        "results/idr1_3/{sample}.log"
    shell:
        """
        idr --samples {input.peak1} {input.peak2} --input-file-type narrowPeak --output-file {output} --rank p.value --plot --use-best-multisummit-IDR --log-output-file {log}
        """




