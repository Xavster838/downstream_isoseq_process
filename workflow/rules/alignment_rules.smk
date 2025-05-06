#these are the rules to align isoseq CCS libraries (fastq or unmapped bam) to the T2T and the associated NHP reference.
rule fastq_from_bam:
    '''get fastq from flnc bam file... or T2T bam file'''
    input:
        read = get_flnc
    output:
        fastq = temp("tmp/{SMP}_{SPRPOP}.fastq"),
    resources:
        mem_mb=8000,
        runtime_hrs=2
    # conda:
    #     "../envs/alignment.yml"
    log:
        "logs/{SMP}_{SPRPOP}_fastq_from_bam.log"
    threads: 2
    run:
       # chekc if bam
       if(input["read"][-4:]==".bam"):
           shell("samtools fastq -@ {threads} {input.read} > {output.fastq} 2> {log}")
       elif(input["read"][-4:]==".sam"):
           shell("samtools fastq {input.read} > {output.fastq} 2> {log} ")
       elif(input["read"][-5:]==".cram"):
           cram_ref = get_cram_ref(input['read'])
           if(cram_ref == -1):
              raise Exception(f"Input Error : could not identify reference of cram file : {input.read}")
           else:
              shell("samtools fastq -T {cram_ref} {input.read} > {output.fastq}")
       elif(input["read"][-6:]==".fastq" or input["read"][-3:]==".fq"):
           shell("cp {input.read} {output.fastq}")
           #shell("gzip {output.fastq}")
       elif( input["read"][-9:]==".fastq.gz" or input["read"][-6:]==".fq.gz" ):
           shell("gunzip -c {input.read} > {output.fastq}")
       else:
           raise Exception(f"Input Error : iso_seq supported formats are fastqs, bams, and crams. File passed : {input.read}")


rule filter_fastq:
    input:
        fastq = rules.fastq_from_bam.output.fastq
    output:
        fastq = temp("tmp/iso_fastqs/{SMP}_{SPRPOP}_FILTERED.fastq.gz")
    resources:
        mem_mb = 8000,
        runtime_hrs=4
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
        mem_mb=16000,
        runtime_hrs=4
    threads:2
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


Hsa_ref = config['Hsa_ref']
t2t_ref = config['T2T_ref']
EXTRA_MMCMD = config['extra_minimap_flags'] if 'extra_minimap_flags' in config.keys() else "" #incase want to further specify minimap commands with splice size and other flags

MMCMD = f"minimap2 -ax splice --sam-hit-only --secondary=yes -p 0.5 --eqx -K 2G {EXTRA_MMCMD}"




rule mm_index_t2t:
    input:
        fasta = t2t_ref
    output:
        mmi = "mmdb/{t2t_version}_ref.mmi",
    resources:
        mem_mb=8000,
        runtime_hrs=10
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
        mem_mb=8000,
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    log:
        "logs/hg38_index.log"
    threads:4
    shell:"""
{MMCMD} -t {threads} -d {output.mmi} {input.fasta}
"""


#get alignment to reference (NHP not T2T or hg38)
ref_link_list = [ get_species_ref_link(cur_row) for i, cur_row in manifest_df.iterrows() ]

def get_manifest_ref_path(wc):
    '''for mm_index_nhp_ref: get original assembly path'''
    return get_nhp_ref(manifest_df.loc[wc.SMP,:])

rule mm_index_nhp_ref:
    input:
        fasta = get_sample_reference 
    output:
        mmi = "mmdb/{SMP}/{SPRPOP}_{ref_name}_ref.mmi" # get_species_ref_path 
    wildcard_constraints:
        ref_name = "|".join( [Path(x).stem for x in manifest_df['reference']] )
    resources:
        mem_mb = 10000,
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    threads:4
    shell:"""
    {MMCMD} -t {threads} -d {output.mmi} {input.fasta}
"""

#map isoseq to nhp reference
rule map_mm:
    input:
        fasta = "tmp/iso_fastqs/{SMP}_{SPRPOP}_FILTERED_{frac}.fastq",  #rules.split_fastq.output.fastq , #"iso_fastqs/{species}/{sample}/{species}_{sample}_{frac}.fastq" ,  #"iso_fastas/{species}/{read}_{frac}.fasta",
        mmi = get_species_ref_path 
    output:
        bam = temp("tmp/alignments/{SMP}/{ref_name}/{SMP}_{SPRPOP}_FILTERED_{frac}_{ref_name}.mm.bam"),
    benchmark:
        "benchmarks/{SMP}_{SPRPOP}_FILTERED_{frac}_{ref_name}.mm.bam.bench",
    wildcard_constraints:
        ref_name = "|".join( [Path(x).stem for x in manifest_df['reference']] )
    resources:
        mem_mb= 4000 , #lambda wildcards, attempt: 3 + 2 * attempt,
        mem_sw=lambda wildcards, attempt: 3 + 0 * attempt, # the mmi index is ~ 7Gb for a hg38
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    threads: 8
    shell:"""
{MMCMD} -t {threads} \
    {input.mmi} {input.fasta} | \
    samtools view -F 2052 -b - | \
    samtools sort -T tmp/{wildcards.SPRPOP}_{wildcards.SMP}_{wildcards.frac}_{wildcards.ref_name} -m {resources.mem_mb}M - > {output.bam}
"""

rule mergeBams:
    input:
        bams = get_nhp_ref_input #expand("tmp/alignments/{{ref_name}}/{{SMP}}_{{SPRPOP}}_FILTERED_{frac}_{{ref_name}}.mm.bam" , frac = fracIDs)
    output:
        bam= "alignments/{SMP}/{ref_name}/{SMP}_{SPRPOP}_FILTERED_{ref_name}.mm.bam", #"../aln_2_species_own_ref_clr/{species}/{read}.bam",
        bai= "alignments/{SMP}/{ref_name}/{SMP}_{SPRPOP}_FILTERED_{ref_name}.mm.bam.bai" #"../aln_2_species_own_ref_clr/{species}/{read}.bam.bai",
    resources:
        mem_mb = 8000, #lambda wildcards, attempt: 8 + 8 * attempt,
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    wildcard_constraints:
        ref_name = "|".join( [Path(x).stem for x in manifest_df['reference']] )
    threads: 4
    shell:"""
        if [ $(echo {input.bams} | wc -w)  -eq 1 ]; then
            ln -s {input.bams} {output.bam}
        else
            samtools merge -@ {threads} {output.bam} {input.bams}
        fi
        samtools index -@ {threads} {output.bam}
    """
#     shell:"""
# samtools merge -@ {threads} {output.bam} {input.bams}
# samtools index -@ {threads} {output.bam}
# """


#map isoseq to T2T
rule map_mm_t2t:
    input:
        fasta = "tmp/iso_fastqs/{SMP}_{SPRPOP}_FILTERED_{frac}.fastq",  #rules.split_fastq.output.fastq , #"iso_fastqs/{species}/{sample}/{species}_{sample}_{frac}.fastq" ,  #"iso_fastas/{species}/{read}_{frac}.fasta",
        mmi = "mmdb/{t2t_version}_ref.mmi" 
    output:
        bam = temp("tmp/alignments/{SMP}/t2t/{SMP}_{SPRPOP}_FILTERED_{frac}_{t2t_version}.mm.bam"),
    benchmark:
        "benchmarks/{SMP}_{SPRPOP}_FILTERED_{frac}_{t2t_version}.mm.bam.bench",
    wildcard_constraints:
        t2t_version = Path(config['T2T_ref']).stem
    resources:
        mem_mb= 4000 , #lambda wildcards, attempt: 3 + 2 * attempt,
        mem_sw=lambda wildcards, attempt: 3 + 0 * attempt, # the mmi index is ~ 7Gb for a hg38
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    threads: 8
    shell:"""
{MMCMD} -t {threads} \
    {input.mmi} {input.fasta} | \
    samtools view -F 2052 -b - | \
    samtools sort -T tmp/{wildcards.SPRPOP}_{wildcards.SMP}_{wildcards.frac}_{wildcards.t2t_version} -m {resources.mem_mb}M - > {output.bam}
"""

rule mergeBams_t2t:
    input:
        bams = get_t2t_input 
    output:
        bam= "alignments/{SMP}/t2t/{SMP}_{SPRPOP}_FILTERED_{t2t_version}.mm.bam", #"../aln_2_species_own_ref_clr/{species}/{read}.bam",
        bai= "alignments/{SMP}/t2t/{SMP}_{SPRPOP}_FILTERED_{t2t_version}.mm.bam.bai" #"../aln_2_species_own_ref_clr/{species}/{read}.bam.bai",
    resources:
        mem_mb = 8000, #lambda wildcards, attempt: 8 + 8 * attempt,
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    wildcard_constraints:
        t2t_version = Path(config['T2T_ref']).stem
    threads: 4
    shell:"""
        if [ $(echo {input.bams} | wc -w)  -eq 1 ]; then
            ln -s {input.bams} {output.bam}
        else
            samtools merge -@ {threads} {output.bam} {input.bams}
        fi
        samtools index -@ {threads} {output.bam}
    """


#map isoseq to T2T
rule map_mm_hg38:
    input:
        fasta = "tmp/iso_fastqs/{SMP}_{SPRPOP}_FILTERED_{frac}.fastq",  #rules.split_fastq.output.fastq , #"iso_fastqs/{species}/{sample}/{species}_{sample}_{frac}.fastq" ,  #"iso_fastas/{species}/{read}_{frac}.fasta",
        mmi = "mmdb/hg38_ref.mmi" 
    output:
        bam = temp("tmp/alignments/{SMP}/hg38/{SMP}_{SPRPOP}_FILTERED_{frac}_hg38.mm.bam"),
    benchmark:
        "benchmarks/{SMP}_{SPRPOP}_FILTERED_{frac}_hg38.mm.bam.bench",
    resources:
        mem_mb= 4000 , #lambda wildcards, attempt: 3 + 2 * attempt,
        mem_sw=lambda wildcards, attempt: 3 + 0 * attempt, # the mmi index is ~ 7Gb for a hg38
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    threads: 8
    shell:"""
{MMCMD} -t {threads} \
    {input.mmi} {input.fasta} | \
    samtools view -F 2052 -b - | \
    samtools sort -T tmp/{wildcards.SPRPOP}_{wildcards.SMP}_{wildcards.frac}_hg38 -m {resources.mem_mb}M - > {output.bam}
"""

rule mergeBams_hg38:
    input:
        bams = get_hg38_input 
    output:
        bam= "alignments/{SMP}/hg38/{SMP}_{SPRPOP}_FILTERED_hg38.mm.bam", #"../aln_2_species_own_ref_clr/{species}/{read}.bam",
        bai= "alignments/{SMP}/hg38/{SMP}_{SPRPOP}_FILTERED_hg38.mm.bam.bai" #"../aln_2_species_own_ref_clr/{species}/{read}.bam.bai",
    resources:
        mem_mb = 8000, #lambda wildcards, attempt: 8 + 8 * attempt,
        runtime_hrs=10
    conda:
        "../envs/alignment.yml"
    threads: 4
    shell:"""
        if [ $(echo {input.bams} | wc -w)  -eq 1 ]; then
            ln -s {input.bams} {output.bam}
        else
            samtools merge -@ {threads} {output.bam} {input.bams}
        fi
        samtools index -@ {threads} {output.bam}
    """
