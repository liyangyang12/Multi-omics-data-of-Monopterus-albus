configfile: "config.yaml"
samples = config["samples"]
threads = config["threads"]
indir = "clean"
outdir = "results/mapping"

rule all:
    input:
        # Mapping
        outdir + "/star_index",
        expand(outdir + "/mapped.1st/{sample}", sample=samples),
        expand(outdir + "/mapped.2nd/{sample}", sample=samples),
        # Filter
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/filtered/{sample}.bam.bai", sample=samples),
        expand(outdir + "/filtered/{sample}.infer.txt", sample=samples),
        expand(outdir + "/filtered/{sample}.stats", sample=samples),
        expand(outdir + "/bedgraph/{sample}Signal.UniqueMultiple.str1.out.bg", sample=samples),
        expand(outdir + "/bw/{sample}.bw",sample=samples),
        expand(outdir + "/bedgraphsort_unnorm/{sample}.bg",sample=samples)        

rule removing_adapter:
    input:
        fq1 = indir + "/{sample}_1.clean.fq.gz",
        fq2 = indir + "/{sample}_2.clean.fq.gz"
    output:
        fq1 = "results/trimmed/{sample}_1.clean_val_1.fq.gz",
        fq2 = "results/trimmed/{sample}_2.clean_val_2.fq.gz",
    shell:
        """
        trim_galore -q 20 --stringency 3 --length 20 -e 0.1 --paired {input.fq1} {input.fq2} --gzip -o results/trimmed
        """

rule star_index:
    input:
        fas = "genome/GCF_001952655.1_M_albus_1.0_genomic.fna",
        gtf = "genome/GCF_001952655.1_M_albus_1.0_genomic.gff"
    output:
        directory(outdir + "/star_index")
    log:
        outdir + "/star_index.log"
    threads:
        threads
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --sjdbGTFfile {input.gtf} \
            --sjdbGTFtagExonParentTranscript Parent \
            --genomeSAindexNbases 13 \
            --genomeFastaFiles {input.fas} &> {log}
        """

rule star_mapping_1st:
    input:
        fq1 = "results/trimmed/{sample}_1.clean_val_1.fq.gz",
        fq2 = "results/trimmed/{sample}_2.clean_val_2.fq.gz",
        idx = rules.star_index.output
    output:
        directory(outdir + "/mapped.1st/{sample}")
    log:
        outdir + "/mapped.1st/{sample}.log"
    threads:
        threads
    params:
        prefix = outdir + "/mapped.1st/{sample}",
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads \
            --alignEndsType EndToEnd \
            --outSAMtype None \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule star_index2pass:
    input:
        fas = "genome/GCF_001952655.1_M_albus_1.0_genomic.fna",
        gtf = "genome/GCF_001952655.1_M_albus_1.0_genomic.gff"
    output:
        directory(outdir + "/star_index2pass")
    log:
        outdir + "/star_index2pass.log"
    threads:
        threads
    params:
        tabs = " ".join([outdir + "/mapped.1st/%s/SJ.out.tab" % sample for sample in samples])
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --sjdbGTFfile {input.gtf} \
            --sjdbGTFtagExonParentTranscript Parent \
            --genomeSAindexNbases 13 \
            --sjdbFileChrStartEnd {params.tabs} \
            --genomeFastaFiles {input.fas} &> {log}
        """

rule star_mapping_2nd:
    input:
        fq1 = "results/trimmed/{sample}_1.clean_val_1.fq.gz",
        fq2 = "results/trimmed/{sample}_2.clean_val_2.fq.gz", 
        idx = outdir + "/star_index2pass",
    output:
        directory(outdir + "/mapped.2nd/{sample}")
    log:
        outdir + "/mapped.2nd/{sample}.log"
    threads:
        threads
    params:
        prefix = outdir + "/mapped.2nd/{sample}",
        tabs = " ".join([outdir + "/mapped.1st/%s/SJ.out.tab" % sample for sample in samples])
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads \
            --alignEndsType EndToEnd \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS \
            --limitBAMsortRAM 10000000000 \
            --limitSjdbInsertNsj 2000000 \
            --quantMode TranscriptomeSAM \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --sjdbFileChrStartEnd {params.tabs} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule filter_alignment:
    input:
        rules.star_mapping_2nd.output,
    output:
        bam = outdir + "/filtered/{sample}.bam"
    shell:
        """
        bamtools filter -in {input}/Aligned.sortedByCoord.out.bam -out {output} \
            -tag NH:1 -isProperPair true -isPrimaryAlignment true
        """

rule infer_experiment:
    input:
        bam = rules.filter_alignment.output,
        bed = "genome/Malbus.bed12"
    output:
        txt = outdir + "/filtered/{sample}.infer.txt"
    log:
        outdir + "/filtered/{sample}.infer.log"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        """

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        """
        bamtools index -in {input}
        """

rule bam_stat:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """

rule bedgraph:
    input:
        bam = outdir + "/filtered/{sample}.bam"
    output:
        bedgraph= outdir + "/bedgraph/{sample}Signal.UniqueMultiple.str1.out.bg"
    params:
        prefix = outdir + "/bedgraph/{sample}",
    shell:
        """
        STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile {input.bam} \
        --outWigType bedGraph --outFileNamePrefix {params.prefix} \
        --outWigNorm RPM --outWigStrand Unstranded
        """

rule bedSort:
    input:
        bedgraph= outdir + "/bedgraph/{sample}Signal.UniqueMultiple.str1.out.bg"
    output:
        bedsort = outdir + "/bedgraphsort/{sample}.bg"
    shell:
        """
        bedSort {input.bedgraph} {output.bedsort}
        """

rule bw:
    input:
        bedsort = outdir + "/bedgraphsort/{sample}.bg",
        si = "/beegfs/zhoulab/lyy/sex_reversal_project/1_ngs/genome/chrsize.txt"
    output:
        bw = outdir + "/bw/{sample}.bw"
    shell:
        """
        bedGraphToBigWig {input.bedsort} {input.si} {output.bw}
        """







