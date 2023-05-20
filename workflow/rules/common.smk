#common functions I need for all snakemake
def get_flnc(wc):
  '''given species sample info, find flnc file from manifest.'''
  return manifest_df.loc[wc.SMP].isoseq_flnc

#functions to process and get information about NHP reference.
def get_nhp_ref(manifest_row):
    '''given sample species combo, get nhp_ref path'''
    return manifest_row.reference #need to use list indexing for sample variable because sample is a function of dataframes.

def get_species_ref_link(manifest_row):
    '''get species, ref_name, and path and return dictionary with this information.'''
    ref_path = get_nhp_ref(manifest_row)
    ref_name = Path(ref_path).stem
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
    out_paths = [ output_string.format(SMP = cur_row["sample"], SPRPOP = cur_row["superpop"],ref_name = Path(cur_row["nhp_ref"]).stem  ) for i, cur_row in manifest_df.iterrows()]
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