
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
    agat_sp_extract_sequences.pl --gff {input.isoform_gff} --fasta {input.ref} -t intron --keep_attributes --merge --output {output.fa}
    sed -i -n '/^>/ s/.*paralog=\([^ ]*\).*transcript_id=\([^ ]*\).*/>{wildcards.SMP}__\\1__\\2/p; /^>/! p' {output.fa} #fix names
    samtools faidx {output.fa}
"""

rule pull_all_isoform_genomic_mRNA_sequence:
    '''given the {loc.name} gff (not top isoform subset gff), get all genomic mRNA'''
    input:
        ref  = rules.fold_ref.output.tmp_folded_ref,  #get_species_sample_ref_path,
        gff = rules.merge_locus_gff_info.output.locus_gff
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_exon_sequence.fa" 
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_exon_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    agat_sp_extract_sequences.pl --gff {input.gff} --fasta {input.ref} -t exon --merge --keep_attributes --output {output.fa}
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
    agat_sp_extract_sequences.pl --gff {input.isoform_gff} --fasta {input.ref} -t exon --merge --keep_attributes --output {output.fa}
    sed -i 's/[^>]*>\([^ ]*\) \(.*\)/>\\1/' {output.fa} #get rid of extranious info for later running ORFfinder
    samtools faidx {output.fa}
"""

rule get_all_isoform_ORF_and_AA:
    '''given the {loc.name} gff (not top isoform subset gff), get all the ORFs and AA'''
    input:
        mRNA_fa = rules.pull_all_isoform_genomic_mRNA_sequence.output.fa 
    output:
        orf_fa = temp( "tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_ORF_sequence.fa" ) ,
        orf_fai = temp( "tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_ORF_sequence.fa.fai" ) ,
        aa_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_aa_sequence.fa"),
        aa_fai = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_aa_sequence.fa.fai"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    orfipy {input.mRNA_fa} --dna $( basename "{output.orf_fa}" ) --pep $( basename "{output.aa_fa}") --outdir $( dirname "{output.orf_fa}" ) --min 100 --max 10000 --start ATG
    samtools faidx {output.orf_fa}
    samtools faidx {output.aa_fa}
'''


rule get_isoform_ORF_and_AA:
    '''given CDS file, predict ORFs.'''
    input:
        mRNA_fa = rules.pull_isoform_genomic_mRNA_sequence.output.fa 
    output:
        orf_fa = temp( "tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_ORF_sequence.fa" ) ,
        aa_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_aa_sequence.fa"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    orfipy {input.mRNA_fa} --dna $( basename "{output.orf_fa}" ) --pep $( basename "{output.aa_fa}") --outdir $( dirname "{output.orf_fa}" ) --min 100 --max 10000 --start ATG
'''

rule fix_ORF_AA_names:
    input:
        orf_fa = rules.get_isoform_ORF_and_AA.output.orf_fa ,
        aa_fa = rules.get_isoform_ORF_and_AA.output.aa_fa ,
        tbl = rules.get_top_paralog_isoforms.output.tbl
    output:
        orf_fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_ORF_sequence.fa" ,
        orf_fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_ORF_sequence.fa.fai",
        aa_fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_aa_sequence.fa",
        aa_fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_aa_sequence.fa.fai",
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    script: "../scripts/process_orf_aa_fastas.py"

# rule get_isoform_aa_sequence:
