SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
diff2 = []
for i in range(0,len(SAMPLE2)-1):
    for j in range(i+1,len(SAMPLE2)):
        diff2.append(SAMPLE2[i] + '_vs_' + SAMPLE2[j])

rule all:
    input:
        expand("results/DARtable/{diffname}.flag",diffname=diff2)

rule filterDMR:
    input:
        "results/diffbind/{sample1}_vs_{sample2}.txt"
    output:
        "results/DARfilter/{sample1}_vs_{sample2}_dar_filter.txt"
    shell:
        """
        python ../scripts/filterDAR.py --dar {input} --result {output}
        """

rule table:
    input:
        "results/DARfilter/{sample1}_vs_{sample2}_dar_filter.txt"
    output:
        flag = "results/DARtable/{sample1}_vs_{sample2}.flag",
    shell:
        """
        echo {wildcards.sample1}_vs_{wildcards.sample2} "`cat {input} |grep "closing" |wc -l`" closing >>{output.flag}
        echo {wildcards.sample1}_vs_{wildcards.sample2} "`cat {input} |grep "opening" |wc -l`" opening >>{output.flag}
        """
