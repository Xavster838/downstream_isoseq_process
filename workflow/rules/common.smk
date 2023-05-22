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
    output_string = "alignments/{ref_name}/{SMP}_{SPRPOP}_FILTERED_{ref_name}.mm.bam"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"],ref_name = Path(cur_row["reference"]).stem  ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_t2t_alignments(wc):
    '''get all t2t mergeBams ouputs for files in manifest'''
    t2t_version = Path(config['T2T_ref']).stem
    output_string = "alignments/t2t/{SMP}_{SPRPOP}_FILTERED_{t2t_version}.mm.bam"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"] ,t2t_version = t2t_version ) for i, cur_row in manifest_df.iterrows()]
    return out_paths

def get_all_hg38_alignments(wc):
    '''get all t2t mergeBams ouputs for files in manifest'''
    output_string = "alignments/hg38/{SMP}_{SPRPOP}_FILTERED_hg38.mm.bam"
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"]) for i, cur_row in manifest_df.iterrows()]
    return out_paths

  
def get_col_bam(cur_sample , alignment_type = "sample_reference"):
    '''given a merge rule, check if manifest.df already has an alignment bam for the given alignment or return False if not.'''        
    if alignment_type == 'sample_reference':
        aln_col = "isoseq_ref_aligned_bam"
        out_bam_str = get_nhp_ref_name( get_nhp_ref(manifest_df.loc[cur_sample]))
    elif alignment_type == "T2T":
        aln_col = "isoseq_T2T_aligned_bam"
        out_ref_str = Path(config['T2T_ref']).stem
    elif alignment_type == "hg38":
        aln_col = "isoseq_hg38_aligned_bam"
        out_ref_str = "hg38"
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
    ref_name = Path(manifest_df.loc[cur_sample]["reference"]).stem
    out_name = f"alignments/{ref_name}/{cur_sample}_{ref_name}_FILTERED_{ref_name}.mm.bam"
    prev_aln = get_col_bam(cur_sample, alignment_type = "sample_reference")
    return prev_aln or expand("tmp/alignments/{{ref_name}}/{{SMP}}_{{SPRPOP}}_FILTERED_{frac}_{{ref_name}}.mm.bam" , frac = fracIDs)
