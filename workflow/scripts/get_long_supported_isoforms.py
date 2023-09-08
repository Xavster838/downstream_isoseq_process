import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


canonical_bed12_path = snakemake.input['canonical_bed12'] 
abundance_tbl_path = snakemake.input['abundance_tbl']
ref_gene_mappings_bed_path = snakemake.input['ref_gene_mappings_bed']
isoform_gff_path = snakemake.input['isoform_gff ']

#step 1: filter abundance tbl 
MIN_ISO_COUNT = snakemake.params["min_abundance"]
abund_tbl = pd.read_csv(abundance_tbl_path, comment = "#", sep = "\t")
abund_tbl = abund_tbl.loc[abund_tbl["count_fl"] > MIN_ISO_COUNT,:]


#step 2 filter TBC gff with canonical gff
#step 2a. load bed12
bed_12_cols = ["contig", "start", "end", "aln", ".", "strand", "_start2", "_end2", "rgb", "n_exons", "exon_size", "exon_start"]
bed12 = pd.read_csv(ref_mRNA_bed12_path, sep = "\t", usecols=list(range(12)), names = bed_12_cols)
bed12 = bed12.set_index("aln", drop = False)

#load gff file
trans_df = pd.DataFrame(columns = ['id', 'start', 'end', 'strand', 'paralog', 'exon_start', 'exon_size'], )
trans_df = trans_df.set_index("id", drop = False)
with open(isoform_gff_path, 'r') as file:
    for line in file:
        # Skip comment lines starting with '#'
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        attributes = columns[8].strip().split(';')
        attribute_dict = {key_value.split('=')[0]: key_value.split('=')[1] for key_value in attributes}
        if columns[2] == "transcript":
            if "nbis" in attribute_dict["ID"]: #not sure why, but gff has some extra annotations with an nbis- name in ID...
                continue
            new_transcript = pd.Series({"id" : attribute_dict["transcript_id"] , "start" : int(columns[3]), "end" : int(columns[4]), 
                                        "strand" : columns[6], "paralog" : attribute_dict["paralog"],
                                        "exon_start" : [], "exon_size" : []})
            trans_df = pd.concat([trans_df, new_transcript.to_frame().T.set_index("id", drop = False)])

#add exons : note : can't put in same for loop because sometimes run into exons for transcripts that haven't been found yet
with open(isoform_gff_path, 'r') as file:
    for line in file:
        # Skip comment lines starting with '#'
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        attributes = columns[8].strip().split(';')
        attribute_dict = {key_value.split('=')[0]: key_value.split('=')[1] for key_value in attributes}
        if columns[2] == "exon":
            assert attribute_dict["transcript_id"] in trans_df.index, f"not finding {attribute_dict['transcript_id']} in trans_df index"
            cur_start = int(columns[3])
            cur_size = int(columns[4]) - cur_start
            cur_strand = columns[6]
            trans_df.loc[attribute_dict["transcript_id"], "exon_start"].append(cur_start)
            trans_df.loc[attribute_dict["transcript_id"], "exon_size"].append(cur_size)

#filter on length
FLANK_TOLERANCE = snakemake.params["flank_tolerance"]
length_filter_isos = pd.DataFrame(columns = ['id', 'start', 'end', 'strand', 'paralog', 'exon_start', 'exon_size'], )
for transcript_id, row in trans_df.iterrows():
    if(row.paralog not in bed12.index):
        missed += 1
        continue
    cur_aln = bed12.loc[row.paralog]   
    if(row.start <= cur_aln.start + 100 and row.end >= cur_aln.end - 100):
        length_filter_isos = pd.concat([length_filter_isos, row.to_frame().T] )
#only take isoforms that fall in filtering of abundance_tbl
#now take al these that also have enough copies
keep_isos = pd.DataFrame(columns = ['id', 'start', 'end', 'strand', 'paralog', 'exon_start', 'exon_size', 'abundance'], )
for transcript_id, row in length_filter_isos.iterrows():
    if(transcript_id in list(abund_tbl['pbid'])):
        row['abundance'] = int(abund_tbl.loc[abund_tbl["pbid"]== transcript_id, "count_fl"])
        #keep largest transcript
        if(row.paralog in list(keep_isos['paralog']) ):
            compare_row = keep_isos.loc[ keep_isos['paralog'] == row.paralog]
            if(sum(row.exon_size) > sum(compare_row.exon_size[0] ) ):
                keep_isos = keep_isos.drop([compare_row.id[0]])
                keep_isos = pd.concat([keep_isos, row.to_frame().T] )
        else:
            keep_isos = pd.concat([keep_isos, row.to_frame().T] )

#write to tbl
keep_isos.to_csv(snakemake.output["keep_isos_tbl"], sep = "\t", header=True, index = False )
#write list
with open(snakemake.output["keep_isos_lst"], "w") as file:
    for element in my_list:
        file.write(element + '\n')