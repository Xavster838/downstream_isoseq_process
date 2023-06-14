
ref_link_list = [ get_species_ref_link(cur_row) for i, cur_row in manifest_df.iterrows() ]
nhp_ref_dict = { r['ref_name'] : r['ref_path'] for r in ref_link_list}

def get_species_sample_ref_path(wc):
    '''given species sample combo, return reference fasta to index.'''
    return nhp_ref_dict[wc['ref2']]

rule pull_paralog_genomic_sequence:
    '''given locus annotations from annotation_rule. pull full genomic sequence'''
    input:
        ref  = get_species_sample_ref_path,
        locus_bed = "reference_annotations/{loc_name}/{ref1}/{ref2}__{loc_name}_mappings.bed"
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}_genomic_sequence.fa" ,
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}_genomic_sequence.fa.fai"
    resources:
        mem_mb = 4000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    bedtools getfasta -fi {input.ref} -bed {input.locus_bed} -name -fo {output.fa}
    samtools faidx {output.fa}
"""



#rule pull_isoform_genomic_mRNA_sequence:
# rule pull_isoform_intronic_sequence:
# rule pull_isoform_genomic_mRNA_sequence:
# rule get_isoform_ORF:
# rule get_isoform_aa_sequence: