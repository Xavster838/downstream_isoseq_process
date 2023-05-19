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
