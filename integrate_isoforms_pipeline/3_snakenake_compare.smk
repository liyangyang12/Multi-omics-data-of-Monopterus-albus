rule all:
    input:
        "results/gffcompare_v2/compare.annotated.gtf",
        "results/3UTR/gasr_3_UTR.bed"


rule gffcompare:
    input:
        ref = "genome/ncbi.gtf",
        gtf = "results/merged2/gasr.gtf",        
    output:
        out3 = "results/gffcompare_v2/compare.annotated.gtf",
    log:
        log = "results/gffcompare_v2/compare.log"
    params:
        prefix = "results/gffcompare_v2/compare"
    shell:
        """
        gffcompare -r {input.ref} -R -o {params.prefix} {input.gtf} &> {log}
        """

rule extract_3UTR:
    input:
        refid = "results/merged2/gasr_refid.txt",
        bed = "results/merged2/tama_merge.modified.bed"
    output:
        "results/3UTR/gasr_3_UTR.bed"
    shell:
        """
        python software/dapars/src/DaPars_Extract_Anno_modified.py -b {input.bed} -s {input.refid} -o {output}
        """