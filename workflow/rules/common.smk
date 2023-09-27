#common functions I need for all snakemake
def get_flnc(wc):
  '''given species sample info, find flnc file from manifest.'''
  return manifest_df.loc[wc.SMP].isoseq_flnc

#functions to process and get information about NHP reference.
def get_nhp_ref(manifest_row):
    '''given sample species combo, get nhp_ref path'''
    return manifest_row.reference #need to use list indexing for sample variable because sample is a function of dataframes.

def get_nhp_ref_name( ref_path):
    '''given sample reference path from get_nhp_ref, return name to use for naming files.'''
    return Path(ref_path).stem

def get_species_ref_link(manifest_row):
    '''get species, ref_name, and path and return dictionary with this information.'''
    ref_path = get_nhp_ref(manifest_row)
    ref_name = get_nhp_ref_name(ref_path) 
    return {'superpop' : manifest_row["superpop"], 'ref_name' : ref_name , 'ref_path' : ref_path}


def get_species_ref_path(wc):
    '''given sample species combo, return all nhp_ref mmi_index to use as reference'''
    ref_mmi_string = "mmdb/{superpop}_{ref_name}_ref.mmi"
    ref_dict = get_species_ref_link(manifest_df.loc[ wc['SMP'] ])
    ref_out_name = ref_mmi_string.format(superpop = ref_dict['superpop'], ref_name = ref_dict['ref_name'])
    return ref_out_name


def get_all_ref_alignments(wc):
    '''get all sample reference defined alignment mergeBams outputs for files in manifest'''
    output_string = "alignments/{SMP}/{ref_name}/{SMP}_{SPRPOP}_FILTERED_{ref_name}.mm.bam"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"],ref_name = Path(cur_row["reference"]).stem  ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_t2t_alignments(wc):
    '''get all t2t mergeBams ouputs for files in manifest'''
    t2t_version = Path(config['T2T_ref']).stem
    output_string = "alignments/{SMP}/t2t/{SMP}_{SPRPOP}_FILTERED_{t2t_version}.mm.bam"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"] ,t2t_version = t2t_version ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_hg38_alignments(wc):
    '''get all t2t mergeBams ouputs for files in manifest'''
    output_string = "alignments/{SMP}/hg38/{SMP}_{SPRPOP}_FILTERED_hg38.mm.bam"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"]) for i, cur_row in manifest_df.iterrows()]
    return out_paths


def get_col_bam(cur_sample , alignment_type = "sample_reference"):
    '''given a merge rule, check if manifest.df already has an alignment bam for the given alignment or return False if not.'''        
    if alignment_type == 'sample_reference':
        aln_col = "isoseq_ref_aligned_bam"
    elif alignment_type == "T2T":
        aln_col = "isoseq_T2T_aligned_bam"
    elif alignment_type == "hg38":
        aln_col = "isoseq_hg38_aligned_bam"
    else:
        raise ValueError(f"Invalid alignment_type argument: can only be either: ['sample_reference', 'T2T', 'hg38']: passed: {alignment_type} ")
    cur_aln_path = manifest_df.loc[cur_sample][aln_col]
    if isinstance( cur_aln_path, str ):
        if(os.path.exists(cur_aln_path)):
            if( os.path.splitext(cur_aln_path)[1] == ".bam" ):
                return cur_aln_path
            else:
                print(f"file in {aln_col} column for sample {cur_sample} in manifest not of type bam... generating alignments.")
        else:
            print(f"file path in {aln_col} column for sample {cur_sample} in manifest not found... generating alignments.")
    else:
        print(f"not seeing file passed for sample {cur_sample} in column {aln_col} for manifest... generating alignments")
    return False

def get_nhp_ref_input( wc ):
    '''figure out if previous sample reference bam exists and return that or expanded sub bams of fracIDs.'''
    cur_sample = wc.SMP    
    prev_aln = get_col_bam(cur_sample, alignment_type = "sample_reference")
    if not prev_aln:
        return expand("tmp/alignments/{{ref_name}}/{{SMP}}_{{SPRPOP}}_FILTERED_{frac}_{{ref_name}}.mm.bam" , frac = fracIDs)
    return prev_aln


def get_t2t_input(wc):
    '''figure out if previous T2T aligned bam exists and return that or expanded sub bams from splitting by fracIDs'''
    cur_sample = wc.SMP
    prev_aln = get_col_bam(cur_sample, alignment_type = "T2T")
    return prev_aln or expand("tmp/alignments/{{SMP}}/t2t/{{SMP}}_{{SPRPOP}}_FILTERED_{frac}_{{t2t_version}}.mm.bam" , frac = fracIDs)

def get_hg38_input(wc):
    '''figure out if previous T2T aligned bam exists and return that or expanded sub bams from splitting by fracIDs'''
    cur_sample = wc.SMP
    prev_aln = get_col_bam(cur_sample, alignment_type = "hg38")
    return prev_aln or expand("tmp/alignments/{{SMP}}/hg38/{{SMP}}_{{SPRPOP}}_FILTERED_{frac}_hg38.mm.bam" , frac = fracIDs)


def get_all_ref_stats(wc):
    '''get all sample reference defined alignment stat outputs for rule get_sample_ref_stats'''
    output_string = "alignments/{SMP}/{ref_name}/stats/{SMP}_{SPRPOP}_FILTERED_{ref_name}.mm.bam.stats.tbl"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"],ref_name = Path(cur_row["reference"]).stem  ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_t2t_stats(wc):
    '''get all t2t mergeBams ouputs for files in manifest'''
    t2t_version = Path(config['T2T_ref']).stem
    output_string = "alignments/{SMP}/t2t/stats/{SMP}_{SPRPOP}_FILTERED_{t2t_version}.mm.bam.stats.tbl"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"] ,t2t_version = t2t_version ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_hg38_stats(wc):
    '''get all hg38 mergeBams ouputs for files in manifest'''
    output_string = "alignments/{SMP}/hg38/stats/{SMP}_{SPRPOP}_FILTERED_hg38.mm.bam.stats.tbl"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"]) for i, cur_row in manifest_df.iterrows()]
    return out_paths

#isoseq3 collapse
def get_all_ref_collapse(wc):
    '''get all sample reference defined alignment stat outputs for rule get_sample_ref_stats'''
    output_string = "alignments/{SMP}/{ref_name}/collapsed_gff/{SMP}_{SPRPOP}_FILTERED_{ref_name}.mm.bam.collapsed.gff"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"],ref_name = Path(cur_row["reference"]).stem  ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_t2t_collapse(wc):
    '''get all t2t mergeBams ouputs for files in manifest'''
    t2t_version = Path(config['T2T_ref']).stem
    output_string = "alignments/{SMP}/t2t/collapsed_gff/{SMP}_{SPRPOP}_FILTERED_{t2t_version}.mm.bam.collapsed.gff"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"] ,t2t_version = t2t_version ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_hg38_collapse(wc):
    '''get all hg38 mergeBams ouputs for files in manifest'''
    output_string = "alignments/{SMP}/hg38/collapsed_gff/{SMP}_{SPRPOP}_FILTERED_hg38.mm.bam.collapsed.gff"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"]) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_loc_path(wc):
    '''given a locus name, identify the path in the config variable ref_map_loci that contains the sequence to that locus'''
    if "ref_map_loci" in config.keys():
        if wc.loc_name in config["ref_map_loci"].keys():
            return config["ref_map_loci"][wc.loc_name]
    return None

def get_can_mRNA_path(wc):
    '''given a locus name, identify the path in the config variable ref_canonical_mRNA that contains the reference mRNA sequence to map. return path'''
    if "ref_canonical_mRNA" in config.keys():
        if wc.loc_name in config["ref_canonical_mRNA"].keys():
            return config["ref_canonical_mRNA"][wc.loc_name]
    return None

def get_sample_reference(wc):
    '''given a reference name, like CHM13 or Jim_h1, return the path to that reference identified by one of the samples in the manifest'''
    return [ref_path for ref_path in manifest_df["reference"] if get_nhp_ref_name( ref_path ) == wc.ref2 ][0]

def get_all_intron_fas(wc):
    '''given manifest, get all intron fastas'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_intronic_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_exon_fas(wc):
    '''given manifest, get all exon (mRNA) fastas'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_exon_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_ORF_fas(wc):
    '''given manifest, get all ORF and AA fastas'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_ORF_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas
    

def get_all_dot_plots(wc):
    '''given manifest, get all isoseq dotplots for choosing functional paralogs'''
    output_string = "plots/{loc_name}/{ref1}/{SMP}__{SPRPOP}__{ref2}__{loc_name}_isoseq_dotplots.pdf"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_longest_paralog_isoform_orf(wc):
    '''given manifest: get all isoseq paralog longest ORF isoforms'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_ORF_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_longest_paralog_isoform_aa(wc):
    '''given manifest: get all isoseq paralog longest AA isoforms'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_aa_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_longest_paralog_isoform_introns(wc):
    '''given manifest: get all isoseq paralog longest intron sequences'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_intron_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_longest_paralog_isoform_mRNA(wc):
    '''given manifest: get all isoseq paralog longest intron sequences'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}_longest_paralog_isoform_genomic_mRNA_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

# def get_all_longest_supported_isoform_intron(wc):
#     '''given manifest: get all isoseq paralog longest supported intron sequences'''

def get_all_longest_supported_isoform_mRNA(wc):
    '''given manifest: get all isoseq paralog longest supported mRNA sequences'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_genomic_mRNA_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_longest_supported_isoform_orf(wc):
    '''given manifest: get all isoseq paralog longest supported orf sequences'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_ORF_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas

def get_all_longest_supported_isoform_aa(wc):
    '''given manifest: get all isoseq paralog longest supported aa sequences'''
    output_string = "sequence/{loc_name}/{SMP}/{ref1}/{SMP}_{SPRPOP}_{ref2}__{loc_name}__long_supported_isoforms_aa_sequence.fa"
    out_fastas = []
    for cur_loc in list(config["ref_map_loci"].keys()):
        print(cur_loc)
        out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"], ref1 = Path(cur_row["reference"]).stem, ref2 = Path(cur_row["reference"]).stem, loc_name = cur_loc) for i, cur_row in manifest_df.iterrows()]
        print(out_paths)
        out_fastas = out_fastas + out_paths
    return out_fastas