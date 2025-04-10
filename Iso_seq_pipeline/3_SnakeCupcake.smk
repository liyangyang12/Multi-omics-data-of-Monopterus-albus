#!/usr/bin/env snakemake
samples = ["m54263_190710_094916.subreads"]
threads = 20

rule all:
    input:
        expand("results/cupcake/collapsed/{sample}/out.collapsed.gff", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.sorted.gtf.gz", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.abundance.txt", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.filtered.gff", sample=samples),
        expand("results/cupcake/collapsed/{sample}/out.collapsed.filtered.sorted.gtf.gz", sample=samples),
        expand("results/cupcake/fusion/{sample}/out.gff", sample=samples),
        expand("results/cupcake/fusion.sqanti3/{sample}/out_classification.txt", sample=samples),
        expand("results/cupcake/fusion/{sample}/out.annotated.txt", sample=samples),

rule collapse:
    input:
        fsa = "results/cupcake/mapping/{sample}.hq.fasta", # 序列名称里面包含有FLNC的数量
        sam = "results/cupcake/mapping/{sample}.sorted.sam"
    output:
        out1 = "results/cupcake/collapsed/{sample}/out.collapsed.gff",
        out2 = "results/cupcake/collapsed/{sample}/out.collapsed.gff.unfuzzy",
        out3 = "results/cupcake/collapsed/{sample}/out.collapsed.group.txt",
        out4 = "results/cupcake/collapsed/{sample}/out.collapsed.group.txt.unfuzzy",
        out5 = "results/cupcake/collapsed/{sample}/out.ignored_ids.txt",
        out6 = "results/cupcake/collapsed/{sample}/out.collapsed.rep.fa",
    params:
        out = "results/cupcake/collapsed/{sample}/out"
    log:
        log = "results/cupcake/collapsed/{sample}/out.log"
    shell:
        """
        collapse_isoforms_by_sam.py --input {input.fsa} -s {input.sam} \
            --dun-merge-5-shorter -o {params.out} &> {log}
        """

rule compress_collapsed_gff:
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.gff"
    output:
        gtf = "results/cupcake/collapsed/{sample}/out.collapsed.sorted.gtf.gz",
        tbi = "results/cupcake/collapsed/{sample}/out.collapsed.sorted.gtf.gz.tbi",
    shell:
        """
        bedtools sort -i {input.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule get_abundance_post_collapse: # 获取每个isoform上的FLNC的数量
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.gff",
        csv = "results/isoseq/polished/{sample}.cluster_report.csv"
    output:
        out1 = "results/cupcake/collapsed/{sample}/out.collapsed.abundance.txt",
        out2 = "results/cupcake/collapsed/{sample}/out.collapsed.read_stat.txt"
    params:
        prefix = "results/cupcake/collapsed/{sample}/out.collapsed"
    shell:
        """
        get_abundance_post_collapse.py {params.prefix} {input.csv} &> /dev/null
        """

# Filter away 5' degraded isoforms
rule filter_away_subset:
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.gff",
        abd = "results/cupcake/collapsed/{sample}/out.collapsed.abundance.txt"
    output:
        out1 = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.gff",
        out2 = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.rep.fa",
        out3 = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.abundance.txt"
    params:
        prefix = "results/cupcake/collapsed/{sample}/out.collapsed"
    shell:
        """
        filter_away_subset.py {params.prefix}
        """
    
rule compress_filtered_gff:
    input:
        gff = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.gff"
    output:
        gtf = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.sorted.gtf.gz",
        tbi = "results/cupcake/collapsed/{sample}/out.collapsed.filtered.sorted.gtf.gz.tbi"
    shell:
        """
        bedtools sort -i {input.gff} | bgzip -c > {output.gtf}
        tabix -p gff {output.gtf}
        """

# Fusion gene
rule fusion_finder:
    input:
        fsa = "results/cupcake/mapping/{sample}.hq.fasta",
        sam = "results/cupcake/mapping/{sample}.sorted.sam",
        csv = "results/isoseq/polished/{sample}.cluster_report.csv"
    output:
        out1 = "results/cupcake/fusion/{sample}/out.abundance.txt",
        out2 = "results/cupcake/fusion/{sample}/out.abundance.txt.bak",
        out3 = "results/cupcake/fusion/{sample}/out.gff",
        out4 = "results/cupcake/fusion/{sample}/out.group.txt",
        out5 = "results/cupcake/fusion/{sample}/out.read_stat.txt",
        out6 = "results/cupcake/fusion/{sample}/out.rep.fa",
    params:
        prefix = "results/cupcake/fusion/{sample}/out"
    log:
        log = "results/cupcake/fusion/{sample}/out.fusion.log"
    shell:
        """
        fusion_finder.py --input {input.fsa} -s {input.sam} --cluster_report {input.csv} \
            -o {params.prefix}  &> {log}
        mv {output.out1} {output.out2}
        grep -v '#' {output.out2} > {output.out1}
        """

rule sqanti3_qc_for_fusion:
    input:
        gff = "results/cupcake/fusion/{sample}/out.gff",
        gtf = "genome/GCF_001952655.1_M_albus_1.0_genomic.gtf",
        fsa = "genome/GCF_001952655.1_M_albus_1.0_genomic.fa",
    output:
        out1 = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
        out2 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.fasta",
        out3 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.genePred",
        out4 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.gtf",
        out5 = "results/cupcake/fusion.sqanti3/{sample}/out_corrected.gtf.cds.gff",
        out6 = "results/cupcake/fusion.sqanti3/{sample}/out_junctions.txt",
        out7 = "results/cupcake/fusion.sqanti3/{sample}/out.params.txt",
        out8 = "results/cupcake/fusion.sqanti3/{sample}/out_SQANTI3_report.html",
        out9 = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred"
    log:
        log = "results/cupcake/fusion.sqanti3/{sample}.log"
    params:
        prefix = "results/cupcake/fusion.sqanti3/{sample}"
    shell:
        """
        set +u; source activate SQANTI3.env
        /DATA/lyy/software/SQANTI3-4.0/sqanti3_qc.py --gtf {input.gff} {input.gtf} {input.fsa} \
            --is_fusion -d {params.prefix} &> {log}
        conda deactivate
        """

rule fusion_collate_info:
    input:
        in1 = "results/cupcake/fusion.sqanti3/{sample}/out_classification.txt",
        #in2 = "results/cupcake/fusion.sqanti3/{sample}/refAnnotation_out.genePred",
        in2 = "genome/GCF_001952655.1_M_albus_1.0_genomic.gtf",
        fsa = "genome/GCF_001952655.1_M_albus_1.0_genomic.fa"
    output:
        out1 = "results/cupcake/fusion/{sample}/out.annotated.txt",
        out2 = "results/cupcake/fusion/{sample}/out.annotated_ignored.txt"
    params:
        prefix = "results/cupcake/fusion/{sample}/out"
    shell:
        """
        fusion_collate_info.py {params.prefix} {input.in1} {input.in2} --genome {input.fsa}
        """

