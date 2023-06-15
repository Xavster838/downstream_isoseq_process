
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
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_genomic_sequence.fa" ,
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_genomic_sequence.fa.fai"
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

rule fold_ref:
    '''make sure references are folded to 80 characters or less per line. AGAT does not work with single line fastas.'''
    input:
        ref  = get_species_sample_ref_path,
    output:
        tmp_folded_ref = temp("tmp/{ref2}__folded.fa"),
        fai = temp("tmp/{ref2}__folded.fa.fai"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    fold -w 80 {input.ref} > {output.tmp_folded_ref}
    samtools faidx {output.tmp_folded_ref}
"""

rule pull_isoform_intronic_sequence:
    '''given subset gff file. get intronic sequence for each selected isoform.'''
    input:
        ref  = rules.fold_ref.output.tmp_folded_ref,  #get_species_sample_ref_path,
        isoform_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons_topIsoforms.gff"
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_intronic_sequence.fa" ,
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_intronic_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    agat_sp_extract_sequences.pl --gff {input.isoform_gff} --fasta {input.ref} -t intron --merge --output {output.fa}
    sed -i 's/[^>]*>\([^ ]*\) \(.*\)/>\\1/' {output.fa} #get rid of extranious info for later running ORFfinder
    samtools faidx {output.fa}
"""


rule pull_isoform_genomic_mRNA_sequence:
    '''given subset gff file. get all exons for each selected isoform.'''
    input:
        ref  = rules.fold_ref.output.tmp_folded_ref,  #get_species_sample_ref_path,
        isoform_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons_topIsoforms.gff"
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_exon_sequence.fa" ,
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_exon_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    agat_sp_extract_sequences.pl --gff {input.isoform_gff} --fasta {input.ref} -t exon --merge --output {output.fa}
    sed -i 's/[^>]*>\([^ ]*\) \(.*\)/>\\1/' {output.fa} #get rid of extranious info for later running ORFfinder
    samtools faidx {output.fa}
"""

rule get_isoform_ORF:
    '''given CDS file, predict ORFs.'''
    input:
        isoform_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons_topIsoforms.gff"
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_ORF_sequence.fa" ,
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_ORF_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    orfipy {input.fa} --dna {output.fa} --min 100 --max 10000 --start ATG
    samtools faidx {output.fa}
"""

# rule get_isoform_aa_sequence:
