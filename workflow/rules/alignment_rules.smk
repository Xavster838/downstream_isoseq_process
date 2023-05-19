#these are the rules to align isoseq CCS libraries (fastq or unmapped bam) to the T2T and the associated NHP reference.
rule fastq_from_bam:
    input:
        read = get_flnc
    output:
        fastq = temp("tmp/{SMP}_{SPRPOP}.fastq.gz"),
    resources:
        mem_mb=8000
    # conda:
    #     "../envs/alignment.yml"
    log:
        "logs/{SMP}_{SPRPOP}_fastq_from_bam.log"
    threads: 2
    run:
       # chekc if bam
       if(input["read"][-4:]==".bam"):
           shell("bedtools bamtofastq -i {input.read} -fq /dev/stdout | gzip > {output.fastq}")
       elif(input["read"][-6:]==".fastq" or input["read"][-3:]==".fq"):
           shell("cp {input.read} {output.fastq}")
           shell("gzip {output.fastq}")
       elif( input["read"][-9:]==".fastq.gz" or input["read"][-6:]==".fq.gz" ):
           shell("ln -s {input.read} {output.fastq}")
       else:
           raise Exception(f"Input Error : iso_seq supported formats are fastqs and bams. File passed : {input.read}")


rule filter_fastq:
    input:
        fastq = rules.fastq_from_bam.output.fastq
    output:
        fastq = temp("tmp/iso_fastqs/{SMP}_{SPRPOP}_FILTERED.fastq.gz")
    resources:
        mem_mb = 4000
    threads : 4
    conda:
        "../envs/alignment.yml"
    params:
        min_len = config['min_read_length'] ,
        min_qual = config['min_read_quality']
    shell:'''
        seqkit seq -j {threads} --min-len {params.min_len} --min-qual {params.min_qual} {input.fastq} -o - | gzip -c > {output.fastq}
'''

fracIDs=list(range(20))
rule split_fastq:
    input:
        fastq = rules.filter_fastq.output.fastq #"iso_fastas/{species}/{read}.fasta",
    output:
        fastq = temp(expand("tmp/iso_fastqs/{{SMP}}_{{SPRPOP}}_FILTERED_{frac}.fastq", frac=fracIDs)),
    resources:
        mem_mb=8000
    threads:1
    run:
        import gzip
        from Bio import SeqIO
        outs = output["fastq"]
        with gzip.open( input.fastq , "rt") as handle:
            recs = list( SeqIO.parse(handle, "fastq") )
            nrecs = len(recs)
            nouts = len(outs)

            for idx, out in enumerate(outs):
                start = int( idx*nrecs/nouts)
                end = int(  (idx + 1)*nrecs/nouts)
                print(start, end, nrecs)
                SeqIO.write(recs[start:end], out, "fastq")

MMCMD = "minimap2 -ax splice --sam-hit-only --secondary=yes -p 0.8 --eqx -K 2G"

Hsa_ref = config['Hsa_ref']
t2t_ref = config['T2T_ref']

rule mm_index_t2t:
    input:
        fasta = t2t_ref
    output:
        mmi = "mmdb/{t2t_version}_ref.mmi",
    resources:
        mem_mb=8000
    conda:
        "../envs/alignment.yml"
    log:
        "logs/{t2t_version}_index.log"
    wildcard_constraints:
        t2t_version = Path(config['T2T_ref']).stem
    threads: 4
    shell:"""
{MMCMD} -t {threads} -d {output.mmi} {input.fasta}
"""

rule mm_index_hg38:
    input:
        fasta = Hsa_ref
    output:
        mmi = "mmdb/hg38_ref.mmi",
    resources:
        mem_mb=8000
    conda:
        "../envs/alignment.yml"
    log:
        "logs/hg38_index.log"
    threads:4
    shell:"""
{MMCMD} -t {threads} -d {output.mmi} {input.fasta}
"""


#get alignment to reference (NHP not T2T or hg38)
ref_link_list = [ get_species_ref_link(cur_row) for i, cur_row in manifest.iterrows() ]
nhp_ref_dict = { (r['species'], r['ref_name'] ) : r['ref_path'] for r in ref_link_list}

def get_species_sample_ref_path(wc):
    '''given species sample combo, return reference fasta to index.'''
    return nhp_ref_dict[wc['species'], wc['ref_name']]

rule mm_index_nhp_ref:
    input:
        fasta = get_species_sample_ref_path,
        human_mmi = rules.mm_index_t2t.output.mmi
    output:
        mmi = "mmdb/{species}/{ref_name}_ref.mmi" # get_species_ref_path 
    wildcard_constraints:
        ref_name = "|".join( [Path(x).stem for x in manifest['nhp_ref']] )
    resources:
        mem_mb = 10000
    threads:4
    shell:"""
if [[ {input.fasta} -ef {Hsa_ref} ]]; then 
    ln -s {input.human_mmi} {output.mmi}
else
    {MMCMD} -t {threads} -d {output.mmi} {input.fasta}
fi
"""