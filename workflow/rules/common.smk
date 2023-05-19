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
