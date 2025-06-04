rule get_iso_dot_plots:
    input:
        tbl = "alignments/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_mappings.bam.stats.tbl" , #rules.get_locus_alignment_stats.output.stats,
        rgns = rules.annotate_reference_locus.output.loc_bed
    output:
        plt = "plots/{loc_name}/{SMP}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_isoseq_dotplots.pdf"
    resources:
        mem_mb = 4000
    threads : 2
    conda:
        "../envs/plotting.yml"
    wildcard_constraints:
        ref1 = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    params:
        cur_nhp = "{SPRPOP}",
        cur_sample = "{SMP}"
    script: "../scripts/plot_dot_plots.R"
