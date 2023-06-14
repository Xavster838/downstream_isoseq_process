##NOT locus SPECIFIC
rule get_alignment_stats:
    input:
        bam = "alignments/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam"
    output:
        stats = "alignments/{SMP}/{ref1}/stats/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam.stats.tbl"
    resources:
        mem_mb = 4000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
rb stats {input.bam} > {output.stats}    
"""

rule collapse_to_isoforms:
    input:
        bam = "alignments/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam"
    output:
        gff = "alignments/{SMP}/{ref1}/collapsed_gff/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam.collapsed.gff",
        read_iso_tbl = "alignments/{SMP}/{ref1}/collapsed_gff/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam.collapsed.read_stat.txt",
        abundance_tbl = "alignments/{SMP}/{ref1}/collapsed_gff/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam.collapsed.abundance.txt"
    log:
        "logs/isoform_collapse/{SMP}_{SPRPOP}_FILTERED_{ref1}_{ref2}.collapse.log"
    resources:
        mem_mb = 4000
    threads : 4
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
isoseq3 collapse -j {threads} --log-file {log} {input.bam} {output.gff}
"""

# ####    locus specific
rule annotate_reference_locus:
    '''given a sequence, map and identify regions on the reference corresponding to that sequence'''
    input:
        ref = get_sample_reference,
        loc_seq = get_loc_path,
    output:
        temp_mapping_psl = temp( "tmp/ref_mappings/{loc_name}/{ref1}/{ref2}__{loc_name}.psl" ),
        mapping_bed = "reference_annotations/{loc_name}/{ref1}/{ref2}__{loc_name}_mappings.bed"
    resources:
        mem_mb = 4000
    threads: 4
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] ),
    shell:"""
blat -t=dna -q=dna -minScore=100 -maxIntron=500 -minMatch=3 {input.ref} {input.loc_seq} {output.temp_mapping_psl}
tail -n +6 {output.temp_mapping_psl} | cut -f14,16,17,10,9 | \
    awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $3,$4,$5,"{wildcards.loc_name}" , ".",$1}}' | bedtools sort | awk 'BEGIN {{FS="\\t"; OFS="\\t"}}{{$4=$4"_"NR; print $0}}' > {output.mapping_bed}
"""

rule subset_alignment_bam_by_locus:
    '''subset alignment bam to reads that overlap annotated regions by annotate_reference_locus'''
    input:
        bam = "alignments/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam",
        ref_loc_bed = "reference_annotations/{loc_name}/{ref1}/{ref2}__{loc_name}_mappings.bed"
    output:
        bam = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam",
        bai = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam.bai"
    resources: 
        mem_mb = 8000
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:"""
bedtools intersect -abam {input.bam} -b {input.ref_loc_bed} | samtools sort > {output.bam} 
samtools index {output.bam}
"""

rule get_locus_alignment_stats:
    '''get alignment stats for reads subset by rule subset_alignemnt_bam_by_locus'''
    input:
        locus_bam = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam",
    output:
        stats = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam.stats.tbl"
    resources:
        mem_mb = 8000
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:"""
    rb stats {input.locus_bam} > {output.stats}
"""


# rule visualize_locus_dot_plots:
#     '''visualize alignment stats of a locus from rule get_locus_alignment_stats'''

# rule visualize_locus_alignments:
#     '''visualize actual reads aligned to a particular locus'''


# rule visualize_locus_isoform_support:
#     '''plot collapsed isoforms and dense alignments beneath'''

rule merge_locus_gff_info:
    '''add intersection of collapsed gff with reference annotation information'''
    input:
        gff = rules.collapse_to_isoforms.output.gff,
        locus_bed = rules.annotate_reference_locus.output.mapping_bed
    output:
        locus_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed.gff"
    shell:'''
bedtools intersect -a {input.gff} -b {input.locus_bed} -wb | awk 'BEGIN{{FS="\\t"; OFS="\\t"}}{{$9=$9 " paralog="$13";" ;print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' > {output.locus_gff} 
'''

# rule pull_isoform_genomic_sequence:

# rule pull_isoform_intronic_sequence:

# rule pull_isoform_genomic_mRNA_sequence:

# rule get_isoform_ORF:

# rule get_isoform_aa_sequence:
