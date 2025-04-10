names = ["ncbi", "ngs", "tgs"]
outdir = "results/merged"
rule all:
    input:
        expand(outdir + "/sources/{name}.gtf.gz", name=names),
        expand(outdir + "/sources/{name}.bed", name=names),
        outdir + "/tama_filelist.txt",
        outdir + "/tama_merge.modified.sorted.bed.gz",
        "results/gffcompare/compare.annotated.bed",

rule gtf_ncbi:
    input:
        gtf = "genome/GCF_001952655.1_M_albus_1.0_genomic.gtf.gz"
    output:
        gtf = outdir + "/sources/ncbi.gtf.gz"
    shell:
        """
        cp {input.gtf} {output.gtf}
        """

rule gtf_ngs:
    input:
        gtf = "results/collect/stringtie_taco.gtf.gz"
    output:
        gtf = outdir + "/sources/ngs.gtf.gz"
    shell:
        """
        cp {input.gtf} {output.gtf}
        """

rule gtf_tgs:
    input:
        gtf = "results/collect/cupcake.gtf.gz"
    output:
        gtf = outdir + "/sources/tgs.gtf.gz"
    shell:
        """
        cp {input.gtf} {output.gtf}
        """

rule gtf_to_bed:
    input:
        gtf = "{prefix}.gtf.gz"
    output:
        gtf = "{prefix}.gtf",
        gp = "{prefix}.genepred",
        bed1 = "{prefix}.bed",
    shell:
        """
        gunzip -c {input.gtf} > {output.gtf}
        gtfToGenePred {output.gtf} {output.gp}
        genePredToBed {output.gp} {output.bed1}
        """

rule tama_merge:
    input:
        bed1 = outdir + "/sources/modify/ncbi.bed",
        bed2 = outdir + "/sources/modify/ngs.bed",
        bed3 = outdir + "/sources/modify/tgs.bed"
    output:
        cfg = outdir + "/tama_filelist.txt",
    log:
        log = outdir + "/tama_merge.log"
    params:
        prefix = outdir + "/tama_merge"
    shell:
        """
        echo -e "{input.bed1}\tno_cap\t1,1,2\tncbi" >> {output.cfg}
        echo -e "{input.bed2}\tno_cap\t1,1,2\tngs" >> {output.cfg}
        echo -e "{input.bed3}\tno_cap\t1,2,1\ttgs" >> {output.cfg}
        #set +u; source activate py27
        #python /DATA/lyy/software/tama-master/tama_merge.py -f {output.cfg} -cds ncbi -s ncbi -p {params.prefix} &> {log}
        #conda deactivate
        """

rule modify_name_of_tama_merge_bed:
    input:
        bed = outdir + "/tama_merge.bed"
    output:
        bed = outdir + "/tama_merge.modified.bed",
    shell:
        """
        ./scripts/modify_name_of_tama_merge_bed.py {input.bed} {output.bed}
        """

rule compress_bed:
    input:
        bed = outdir + "/tama_merge.modified.bed"
    output:
        bed = outdir + "/tama_merge.modified.sorted.bed.gz",
        tbi = outdir + "/tama_merge.modified.sorted.bed.gz.tbi"
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule bedtogtf:
    input:
        bed = outdir + "/tama_merge.modified.bed",
    output:
        gp = outdir + "/tama_merge.modified.genepred",
        gtf = outdir + "/tama_merge.modified.gtf",
    shell:
        """
        bedToGenePred {input.bed} {output.gp}
        genePredToGtf file {output.gp} {output.gtf}
        """

rule gffcompare:
    input:
        ref = "genome/GCF_001952655.1_M_albus_1.0_genomic.gtf",
        gtf = outdir + "/tama_merge.modified.gtf",
        
    output:
        out3 = "results/gffcompare/compare.annotated.gtf",
    log:
        log = "results/gffcompare/compare.log"
    params:
        prefix = "results/gffcompare/compare"
    shell:
        """
        gffcompare -r {input.ref} -R -o {params.prefix} {input.gtf} &> {log}
        """

rule gtftobed:
    input:
        gtf="results/gffcompare/compare.annotated.gtf",
    output:
        gp = "results/gffcompare/compare.annotated.gp",
        bed = "results/gffcompare/compare.annotated.bed"
    shell:
        """ 
        gtfToGenePred  {input.gtf} {output.gp}
        genePredToBed {output.gp} {output.bed} 
        """

