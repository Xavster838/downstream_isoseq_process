rule get_iso_dot_plots:
    input:
    output:
    resources:
        mem_mb = 4000
    threads : 2
    conda:
        "../envs/plotting.yml"
    wildcard_constraints:
        ref = "|".join(["hg38", "t2t"] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) ,
        ref2 = "|".join(["hg38", Path(config['T2T_ref']).stem ] + [get_nhp_ref_name(x) for x in manifest_df["reference"]] ) #dealing with fact that t2t has two different reference names
    script: ../scripts/plot_dot_plots.R