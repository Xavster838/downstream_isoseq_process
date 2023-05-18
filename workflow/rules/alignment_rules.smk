#these are the rules to align isoseq CCS libraries (fastq or unmapped bam) to the T2T and the associated NHP reference.
rule fastq_from_bam:
    input:
        read = get_flnc
    output:
        fastq = temp("tmp/{SMP}_{SPRPOP}.fastq.gz"),
    resources:
        mem_mb=8000
    conda:
        "../envs/alignment.yml"
    log:
        "logs/{SMP}_{SPRPOP}_fastq_from_bam.log"
    threads: 2
    shell:'''
touch {output.fastq}
'''
#    run:
#        # chekc if bam
#        if(input["read"][-4:]==".bam"):
#            shell("bedtools bamtofastq -i {input.read} -fq /dev/stdout | gzip > {output.fastq}")
#        elif(input["read"][-6:]==".fastq" or input["read"][-3:]==".fq"):
#            shell("cp {input.read} {output.fastq}")
#            shell("gzip {output.fastq}")
#        elif( input["read"][-9:]==".fastq.gz" or input["read"][-6:]==".fq.gz" ):
#            shell("ln -s {input.read} {output.fastq}")
#        else:
#            raise Exception(f"Input Error : iso_seq supported formats are fastqs and bams. File passed : {input.read}")
