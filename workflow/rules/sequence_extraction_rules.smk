
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

rule add_introns_gff:
    '''Add introns to loc_name gff file'''
    input:
        gff = rules.merge_locus_gff_info.output.locus_gff,
    output:
        temp_intron_gff = temp("tmp/alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_all_collapsed_withIntrons.gff"),
        intron_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_all_collapsed_withIntrons.gff"
    resources:
        mem_mb = 8000
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:'''
agat_sp_add_introns.pl --gff {input.gff} --out {output.temp_intron_gff}
bedtools sort -i {output.temp_intron_gff} > {output.intron_gff}
'''

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

rule pull_all_isoform_intron_sequence:
    '''given the {loc.name} gff (not top isoform subset gff), get all genomic mRNA'''
    input:
        ref  = rules.fold_ref.output.tmp_folded_ref,  #get_species_sample_ref_path,
        gff = rules.add_introns_gff.output.intron_gff
    output:
        fa = "tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_intron_sequence.fa", 
        fai = "tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_intron_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
    agat_sp_extract_sequences.pl --gff {input.gff} --fasta {input.ref} -t intron --merge --keep_attributes --output {output.fa}
    sed -i 's/[^>]*>\([^ ]*\) \(.*\)/>\\1/' {output.fa} #get rid of extranious info for later running ORFfinder
    samtools faidx {output.fa}
"""

rule pull_all_isoform_genomic_mRNA_sequence:
    '''given the {loc.name} gff (not top isoform subset gff), get all genomic mRNA'''
    input:
        ref  = rules.fold_ref.output.tmp_folded_ref,  #get_species_sample_ref_path,
        gff = rules.merge_locus_gff_info.output.locus_gff
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_exon_sequence.fa", 
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
    sed -i 's/[^>]*>\([^ ]*\) \(.*\)/>\\1/' {output.orf_fa} #get rid of extranious info for later running ORFfinder
    sed -i 's/[^>]*>\([^ ]*\) \(.*\)/>\\1/' {output.aa_fa} #get rid of extranious info for later running ORFfinder
    samtools faidx {output.orf_fa}
    samtools faidx {output.aa_fa}
'''

rule get_longest_paralog_isoform_orfs_aa_list:
    '''given AA sequence from get_all_isoform_ORF_and_AA, determine the longest reading frame for each paralog and output just those sequences'''
    input:
        aa_fai = rules.get_all_isoform_ORF_and_AA.output.aa_fai
    output:
        isoform_list = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoforms.lst"),
        tmp_sorted_fai = temp("sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_all_aa_sequence_SORTED.fa.fai"),
        tmp_fai = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_tmp.fai")
    resources:
        mem_mb = 8000
    threads : 2
    log: "logs/{loc_name}__{SMP}__{SPRPOP}__{ref1}__{ref2}__get_longest_paralog_isoform_orfs_aa_list.log"
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    #sort by longest isoform
    sort -nr -k 2,2 {input.aa_fai} > {output.tmp_sorted_fai} 2> {log}
    #get top (process each paralog separately)
    mapfile -t paralog_array < <(cut -f 1 {output.tmp_sorted_fai} | sort | sed "s/\.[0-9]\+_ORF\.[0-9]\+//g" | sort | uniq) #get array of paralog names
    # Print the array elements
    for cur_paralog in "${{paralog_array[@]}}"; do
        grep "${{cur_paralog}}\." {output.tmp_sorted_fai} | sort -nr -k 2,2 > {output.tmp_fai} 
        top_paralog=$( head -n 1 {output.tmp_fai} | cut -f 1 )
        printf "${{top_paralog}}\n" >> {output.isoform_list} 2> {log}
    done
'''

rule pull_longest_isoform_introns_mRNAs:
    '''given list from rule get_longest_paralog_isoform_orfs_aa_list, pull introns and genomic mRNAs sequences.'''
    input:
        lst = rules.get_longest_paralog_isoform_orfs_aa_list.output.isoform_list,
        intron_fa = rules.pull_all_isoform_intron_sequence.output.fa,
        mRNA_fa = rules.pull_all_isoform_genomic_mRNA_sequence.output.fa ,
    output:
        tmp_isoform_tbl = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoforms_no_ORF.lst"),
        intron_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_intron_sequence.fa"),
        mRNA_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_genomic_mRNA_sequence.fa"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    sed 's/_.*//g' {input.lst} > {output.tmp_isoform_tbl}
    seqtk subseq {input.intron_fa} {output.tmp_isoform_tbl} > {output.intron_fa}
    seqtk subseq {input.mRNA_fa} {output.tmp_isoform_tbl} > {output.mRNA_fa}
''' 

rule pull_longest_supported_isoform_introns_mRNAs:
    '''given list from annotation rule get_long_supported_isoforms, pull introns and genomic mRNAs sequences.'''
    input:
        lst = rules.get_long_supported_isoforms.output.keep_isos_lst,
        intron_fa = rules.pull_all_isoform_intron_sequence.output.fa,
        mRNA_fa = rules.pull_all_isoform_genomic_mRNA_sequence.output.fa ,
    output:
        tmp_isoform_tbl = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms.lst"),
        intron_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_intron_sequence.fa"),
        mRNA_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_genomic_mRNA_sequence.fa"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    sed 's/_.*//g' {input.lst} > {output.tmp_isoform_tbl}
    seqtk subseq {input.intron_fa} {output.tmp_isoform_tbl} > {output.intron_fa}
    seqtk subseq {input.mRNA_fa} {output.tmp_isoform_tbl} > {output.mRNA_fa}
''' 

rule pull_longest_paralog_isofrom_ORFs_AAs:
    '''given list from rule get_longest_paralog_isoform_orfs_aa_list, pull ORF, and AA sequences.'''
    input:
        lst = rules.get_longest_paralog_isoform_orfs_aa_list.output.isoform_list,
        orf_fa = rules.get_all_isoform_ORF_and_AA.output.orf_fa ,
        aa_fa = rules.get_all_isoform_ORF_and_AA.output.aa_fa,
    output:
        orf_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_ORF_sequence.fa") ,
        aa_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_aa_sequence.fa"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    seqtk subseq {input.orf_fa} {input.lst} > {output.orf_fa}
    seqtk subseq {input.aa_fa} {input.lst} > {output.aa_fa}
''' 

rule get_longest_ORF_per_isoform:
    '''given sorted aa fa, output a list of names of the longest ORFs for each isoform annotated in rule get_all_isoform_ORF_and_AA '''
    input:
        sorted_aa_fai = get_longest_paralog_isoform_orfs_aa_list.output.tmp_sorted_fai, #to get only longest reading frames for each
    output:
        tmp_fai = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_tmp_longest_iso_aa.fai"),
        longest_iso_ORF_lst = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__longest_aa.lst")
    resources:
        mem_mb = 8000
    threads: 2
    log: "logs/{loc_name}__{SMP}__{SPRPOP}__{ref1}__{ref2}__get_longest_ORF_per_isoform.log"
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
    shell:'''
    #get top (process each paralog separately)
    mapfile -t paralog_array < <(cut -f 1 {input.sorted_aa_fai} | sort | sed "s/_ORF\.[0-9]\+//g" | sort | uniq) #get array of isoform names
    # for each isoform, sort and take top and add to list
    for cur_paralog in "${{paralog_array[@]}}"; do
        grep "${{cur_paralog}}\." {input.sorted_aa_fai} | sort -nr -k 2,2 > {output.tmp_fai} 
        top_paralog=$( head -n 1 {output.tmp_fai} | cut -f 1 )
        printf "${{top_paralog}}\n" >> {output.isoform_list} 2> {log}
    done
'''

rule pull_longest_supported_paralog_isofrom_ORFs_AAs:
    '''given list from rule get_long_supported_isoforms, pull ORF, and AA sequences.'''
    input:
        lst = rules.get_long_supported_isoforms.output.keep_isos_lst,
        longest_iso_aa_lst = rules.get_longest_ORF_per_isoform.output.longest_iso_ORF_lst,
        orf_fa = rules.get_all_isoform_ORF_and_AA.output.orf_fa ,
        aa_fa = rules.get_all_isoform_ORF_and_AA.output.aa_fa,
    output:
        tmp_orf_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__all_ORF_sequence.fa")
        tmp_aa_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__all_aa_sequence.fa")
        tmp_longest_aa_lst = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__longest_aa.lst"),
        orf_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_ORF_sequence.fa") ,
        aa_fa = temp("tmp/sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_aa_sequence.fa"),
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:'''
    seqtk subseq {input.orf_fa} {input.longest_iso_aa_lst} | sed '/^>/ s/_.*//' > {output.tmp_orf_fa}
    seqtk subseq {input.aa_fa} {input.longest_iso_aa_lst} | sed '/^>/ s/_.*//' > {output.tmp_aa_fa}
    seqtk subseq {output.tmp_orf_fa} {input.lst} > {output.orf_fa}
    seqtk subseq {output.tmp_aa_fa} {input.lst} > {output.aa_fa}
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

rule fix_longest_orf_isoform_orf_names:
    '''fix names for longest orf isoform names for ORF sequence fastas'''
    input: 
        fa = rules.pull_longest_paralog_isofrom_ORFs_AAs.output.orf_fa,
        gff = rules.add_introns_gff.output.intron_gff
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_ORF_sequence.fa",
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_ORF_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    script: "../scripts/rename_sequence_fa_reads.py"

rule fix_longest_orf_isoform_AA_names:
    '''fix names for longest orf isoform names for AA sequence fastas'''
    input: 
        fa = rules.pull_longest_paralog_isofrom_ORFs_AAs.output.aa_fa,
        gff = rules.add_introns_gff.output.intron_gff
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_aa_sequence.fa",
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_aa_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    script: "../scripts/rename_sequence_fa_reads.py"

rule fix_longest_orf_isoform_intron_names:
    '''fix names for longest orf isoform names for AA sequence fastas'''
    input: 
        fa = rules.pull_longest_isoform_introns_mRNAs.output.intron_fa,
        gff = rules.add_introns_gff.output.intron_gff
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_intron_sequence.fa",
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_intron_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    script: "../scripts/rename_sequence_fa_reads.py"

rule fix_longest_orf_isoform_mRNA_names:
    '''fix names for longest orf isoform names for AA sequence fastas'''
    input: 
        fa = rules.pull_longest_isoform_introns_mRNAs.output.mRNA_fa,
        gff = rules.add_introns_gff.output.intron_gff
    output:
        fa = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_genomic_mRNA_sequence.fa",
        fai = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_genomic_mRNA_sequence.fa.fai"
    resources:
        mem_mb = 8000
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    script: "../scripts/rename_sequence_fa_reads.py"

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
