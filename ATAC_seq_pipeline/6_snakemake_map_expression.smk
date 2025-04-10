SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
SAMPLE3 = ["closing","opening"]
diff = []
for i in range(0,len(SAMPLE2)-1):
    for j in range(i+1,len(SAMPLE2)):
        diff.append(SAMPLE2[i] + '_vs_' + SAMPLE2[j])
diff2 = []
for i in range(0,len(SAMPLE2)-1):
    diff2.append(SAMPLE2[i] + '_vs_' + SAMPLE2[i+1])
SAMPLE1 = ["O_I","O_IV","O_V","OT_I","OT_II","OT_III","T_I","T_III","T_IV"]
diff1 = []
for i in range(0,len(SAMPLE1)-1):
    for j in range(i+1,len(SAMPLE1)):
        diff1.append(SAMPLE1[i] + '_vs_' + SAMPLE1[j])
diff3=[]
for i in range(0,len(diff1)):
    diff3.append(diff1[i] + '_' + diff[i])

rule all:
    input:
        "results/identify/closing.txt", 
        "results/identify/opening.txt",
        expand("results/epi/{diffname}_dar_dmr.txt",diffname=diff),
        expand("results/closet/{diffname}_closet.txt",diffname=diff),
        expand("results/deseq2/{diffname2}_filter.txt",diffname2=diff1),
        expand("results/map_expression/{diffname3}.txt",diffname3=diff3)


rule get_bed:
    input:
        dar = expand("DAR/DARfilter/{diffname}_dar_filter.txt",diffname=diff),
    output:
        closing = "results/identify/closing.txt",
        opening = "results/identify/opening.txt",
    shell:
        """
        cat {input.dar} | grep "closing" | cut -f 1,2,3|sort|uniq > {output.closing}
        cat {input.dar} | grep "opening" | cut -f 1,2,3|sort|uniq > {output.opening}
        """

rule get_epi:
    input:
        dar = "DAR/DARfilter/{diffname}_dar_filter.txt",
    output:
        dar = "DAR/DARfilter/{diffname}_dar_filter_cut.txt",
        epi = "results/epi/{diffname}_dar.txt"
    shell:
        """
        cat {input.dar} | cut -f 1,2,3,12 > {output.dar}
        cat {output.dar}  > {output.epi}
        """

rule sort_DAR:
    input:
        DAMR =   "results/epi/{diffname}_dar.txt"
    output:
        result =  "results/epi/{diffname}_dar_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n {input.DAMR} > {output.result}
        """

rule closet:
    input:
        promotor = "results/promoter/malbus_promoter_sort_filter.bed",
        DAMR = "results/epi/{diffname}_dar_sort.txt"
    output:
        result = "results/closet/{diffname}_closet.txt"
    shell:
        """
        bedtools closest -d -a {input.DAMR} -b {input.promotor} > {output.result}
        """


rule filterDEG:
    input:
        "results/expression/deseq2/{sample1}_vs_{sample2}.csv"
    output:
        "results/deseq2/{sample1}_vs_{sample2}_filter.txt"
    shell:
        """
        python scripts/filterDEG.py --deg {input} --result {output}
        """


rule map_expression:
    input:
        result = "results/closet/{sample3}_vs_{sample4}_closet.txt",
        DEG = "results/deseq2/{sample1}_vs_{sample2}_filter.txt"
    output:
        mapresult = "results/map_expression/{sample1}_vs_{sample2}_{sample3}_vs_{sample4}.txt"
    shell:
        """
        python scripts/map_expression.py --closet {input.result} --deg {input.DEG} --sample1 {wildcards.sample3} --sample2 {wildcards.sample4} --result {output.mapresult}
        """
