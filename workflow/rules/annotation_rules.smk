##NOT locus SPECIFIC
rule get_alignment_stats:
    input:
        bam = "alignments/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam"
    output:
        stats = "alignments/{SMP}/{ref1}/stats/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam.stats.tbl"
    resources:
        mem_mb = 4000,
        runtime_hrs=4
    threads : 2
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    shell:"""
paste <(rb stats {input.bam}) <(samtools view {input.bam} | cut -f 2 | awk 'BEGIN{{print "alignment_flag"}}{{ptin $0}}') > {output.stats}
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
        mem_mb = 4000,
        runtime_hrs=4
    threads : 4
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ), #dealing with fact that t2t has two different reference names
        SPRPOP="|".join(manifest_df['superpop']),
        SMP = "|".join(manifest_df["sample"])
    shell:"""
isoseq3 collapse -j {threads} --log-file {log} {input.bam} {output.gff}
sed -i 's/gene_id /gene_id=/g' {output.gff} 
sed -i 's/transcript_id /transcript_id=/g' {output.gff} #add equals signs to attributes
sed -i 's/"//g' {output.gff} #get rid of double quotes
"""

# ####    locus specific
rule annotate_reference_locus:
    '''given a sequence, map and identify regions on the reference corresponding to that sequence. Output bed'''
    input:
        ref = get_sample_reference,
        loc_seq = get_loc_path,
    output:
        loc_bed = "reference_annotations/{loc_name}/{SMP}__{SPRPOP}/{ref1}/{ref2}__{loc_name}_mappings.bed"
        #temp_mapping_psl = temp( "tmp/ref_mappings/{loc_name}/{SMP}/{ref1}/{ref2}__{loc_name}.psl" ),
        #mapping_bed = "reference_annotations/{loc_name}/{SMP}/{ref1}/{ref2}__{loc_name}_mappings.bed"
    resources:
        mem_mb = 4000,
        runtime_hrs=4
    threads: 4
    params:
        min_aln_size = config['min_map_len']
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] ),
    shell:"""
minimap2 -x asm20 --secondary=yes -p 0.8 --eqx -t {threads} -r 500 -K 100M {input.ref} {input.loc_seq} | 
	awk -v min_size={params.min_aln_size} 'BEGIN{{OFS="\\t"}} ($9 - $8 >= min_size) {{print $6, $8, $9, $1, $10, $5}}' > {output.loc_bed}
"""

rule annotate_canonical_ref_mRNA:
    '''given a canonical mRNA, map and identify regions on the reference corresponding to that mRNA. Output: bed12'''
    input:
        ref = get_sample_reference,
        ref_loc_bed = rules.annotate_reference_locus.output.loc_bed,
        can_mRNA = get_can_mRNA_path,
    output:
        temp_bam = temp( "tmp/ref_mappings/{loc_name}/{SMP}__{SPRPOP}/{ref1}/{ref2}__{loc_name}_canonical_mRNA_mappings.bam" ),
        bed12 = "reference_annotations/{loc_name}/{SMP}__{SPRPOP}/{ref1}/{ref2}__{loc_name}_canonical_mRNA_mappings.bed12"
    resources:
        mem_mb = 4000,
        runtime_hrs=4
    threads: 4
    conda:
        "../envs/annotation.yml"
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] ),
    shell:"""
minimap2 -a -k19 -w5 --splice -g 2k  -A1 -B2 -O2,32 \
    -E1,0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=no -G50k \
    --secondary=yes -N 100 \
    {input.ref} {input.can_mRNA} | \
    samtools view -F 2052 -b - | \
    samtools sort -@ 4 - > {output.temp_bam}

bedtools intersect -a {output.temp_bam} -b {input.ref_loc_bed} -wa -wb -bed | \
    awk -F '\t' 'BEGIN{{OFS=FS}}{{$4=$16; print}}' > {output.bed12}
"""


rule subset_alignment_bam_by_locus:
    '''subset alignment bam to reads that overlap annotated regions by annotate_reference_locus'''
    input:
        bam = "alignments/{SMP}/{ref1}/{SMP}_{SPRPOP}_FILTERED_{ref2}.mm.bam",
        ref_loc_bed = rules.annotate_reference_locus.output.loc_bed 
    output:
        bam = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam",
        bai = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam.bai"
    resources: 
        mem_mb = 8000,
        runtime_hrs=4
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
        mem_mb = 8000,
        runtime_hrs=4
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:"""
    paste <(rb stats {input.locus_bam}) <(samtools view {input.locus_bam} | cut -f 2 | awk 'BEGIN{{print "alignment_flag"}}{{print $0}}') > {output.stats}
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
        locus_bed = rules.annotate_reference_locus.output.loc_bed
    output:
        locus_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed.gff"
    resources:
        mem_mb = 8000,
        runtime_hrs=4
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:'''
bedtools intersect -header -a {input.gff} -b {input.locus_bed} -wb | awk 'BEGIN{{FS="\\t"; OFS="\\t"}}{{$9=$9 " paralog="$13";" ;print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' > {output.locus_gff} 
'''

rule get_top_paralog_isoforms:
    '''find the most abundantly expressed isoforms for a given paralog defined by the locus annotations. Write to a table'''
    input:
        locus_gff = rules.merge_locus_gff_info.output.locus_gff,
        abundance_tbl = rules.collapse_to_isoforms.output.abundance_tbl
    output:
        tbl = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_top_paralog_isoforms.tbl"
    resources:
        mem_mb = 8000,
        runtime_hrs=4
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    script: "../scripts/get_top_paralog_isoforms.py"

rule subset_gff_top_isoforms:
    '''given locus gff, add introns for later extracting sequence'''
    input:
        intron_gff = rules.merge_locus_gff_info.output.locus_gff,
        isoform_tbl = rules.get_top_paralog_isoforms.output.tbl
    output:
        subset_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_topIsoforms.gff"
    resources:
        mem_mb = 8000,
        runtime_hrs=4
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:'''
top_isoforms=( $(cut -f 2 {input.isoform_tbl} | tail -n +2) )

if [ ${{#top_isoforms[@]}} -gt 0 ]; then
    for isoform in "${{top_isoforms[@]}}"; do
        if grep "${{isoform}};" {input.intron_gff} >> {output.subset_gff}; then
            echo "match_found"
        else
            echo "no match found for ${{isoform}} in {wildcards.SMP}"
        fi
    done
else
    echo "no isoforms in {input.isoform_tbl}"
fi

touch {output.subset_gff}
'''

rule add_introns_locus_gff:
    '''given locus gff, add introns for later extracting sequence'''
    input:
        gff = rules.subset_gff_top_isoforms.output.subset_gff,
    output:
        temp_intron_gff = temp("tmp/alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons.gff"),
        intron_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons_topIsoforms.gff"
    resources:
        mem_mb = 8000,
        runtime_hrs=4
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

rule add_introns_to_isoform_gff:
    '''given locus gff, add introns for later extracting sequence'''
    input:
        gff = rules.merge_locus_gff_info.output.locus_gff #rules.collapse_to_isoforms.output.gff,
    output:
        temp_intron_gff = temp("tmp/alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons.gff"),
        intron_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons.gff"
    resources:
        mem_mb = 8000,
        runtime_hrs=4
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

rule get_long_supported_isoforms:
    '''find the isoforms that retain the majority of exons in canonical mRNA. Return tbl'''
    input:
        canonical_bed12 = rules.annotate_canonical_ref_mRNA.output.bed12 ,#canonical mRNA bed12
        abundance_tbl = rules.collapse_to_isoforms.output.abundance_tbl ,
        ref_gene_mappings_bed = rules.annotate_reference_locus.output.loc_bed,
        isoform_gff = rules.add_introns_to_isoform_gff.output.intron_gff,
    output:
        keep_isos_tbl = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_long_supported_isoforms.tbl",
        keep_isos_lst = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_long_supported_isoforms.lst"
    resources:
        mem_mb = 8000,
        runtime_hrs=4
    threads: 2
    params:
        flank_tolerance = 100 ,#keep isoforms if their alignments are at least within 100bp of canonical mRNA mapping.
        min_abundance = 3 #only keep isoforms that have at least min_abundance reads in them
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    script: "../scripts/get_long_supported_isoforms.py"
