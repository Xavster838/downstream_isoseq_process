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
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] )
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
rb stats {input.bam} > {output.stats}    
"""

rule get_sample_ref_stats:
    input:
        stats = get_all_ref_stats

rule get_t2t_stats:
    input:
        stats = get_all_t2t_stats

rule get_hg38_stats:
    input:
        stats = get_all_hg38_stats

# rule collapse_to_isoforms:


# ####    locus specific
# rule annotate_reference_locus:
#     '''given a sequence, map and identify regions on the reference corresponding to that sequence'''

# rule subset_alignment_bam_by_locus:
#     '''subset alignment bam to reads that overlap annotated regions by annotate_reference_locus'''

# rule get_locus_alignment_stats:
#     '''get alignment stats for reads subset by rule subset_alignemnt_bam_by_locus'''

# rule visualize_locus_dot_plots:
#     '''visualize alignment stats of a locus from rule get_locus_alignment_stats'''

# rule visualize_locus_alignments:
#     '''visualize actual reads aligned to a particular locus'''


# rule visualize_locus_isoform_support:
#     '''plot collapsed isoforms and dense alignments beneath'''

# rule pull_isoform_genomic_sequence:

# rule pull_isoform_intronic_sequence:

# rule pull_isoform_genomic_mRNA_sequence:

# rule get_isoform_ORF:

# rule get_isoform_aa_sequence:
