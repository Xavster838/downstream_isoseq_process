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
sed -i 's/gene_id /gene_id=/g' {output.gff} 
sed -i 's/transcript_id /transcript_id=/g' {output.gff} #add equals signs to attributes
sed -i 's/"//g' {output.gff} #get rid of double quotes
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
    awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $3,$4,$5,"{wildcards.loc_name}" , ".",$1}}' | bedtools sort | awk 'BEGIN {{FS="\\t"; OFS="\\t"}}{{$4=$4"_"NR; print $0}}' | \
	awk -v min_len={config[min_map_len]} '{{ if (($3 - $2) > min_len) print }}' > {output.mapping_bed}
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
        mem_mb = 8000
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
        mem_mb = 8000
    threads: 2
    wildcard_constraints:
        loc_name = "|".join( list(config["ref_map_loci"].keys() ) ),
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join( [get_nhp_ref_name( ref_path ) for ref_path in manifest_df["reference"] ] )
    conda:
        "../envs/annotation.yml"
    shell:'''
top_isoforms=( $(cut -f 2 {input.isoform_tbl} | tail -n +2) )
for isoform in "${{top_isoforms[@]}}"; do
    if grep "${{isoform}};" {input.intron_gff} >> {output.subset_gff}; then
        echo "match_found"
    else
        echo "no match found for ${{isoform}} in {wildcards.SMP}"
    fi
done
'''

rule add_introns_locus_gff:
    '''given locus gff, add introns for later extracting sequence'''
    input:
        gff = rules.subset_gff_top_isoforms.output.subset_gff,
    output:
        temp_intron_gff = temp("tmp/alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons.gff"),
        intron_gff = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_collapsed_withIntrons_topIsoforms.gff"
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
