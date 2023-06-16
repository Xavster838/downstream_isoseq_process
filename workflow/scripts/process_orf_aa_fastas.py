#select largest reading frames and fix names of reading frames.
import pandas as pd
from Bio import SeqIO
import os

isoform_paralog_df = pd.read_csv(snakemake.input['tbl'], sep = "\t")

def add_info(cur_record):
    '''add isoform, ORF, and ORF length as features for a given SeqRecord'''
    cur_record.isoform, cur_record.orf = cur_record.description.strip().split(" ")[0].split("_")
    cur_record.length = int([x for x in cur_record.description.strip().split(" ") if "length" in x][0].split(":")[1])
    cur_record.paralog = isoform_paralog_df.loc[isoform_paralog_df['top_isoform']==cur_record.isoform, 'paralog'].iloc[0]
    return cur_record

def process_fasta(in_fasta, out_fasta):
    '''same steps for processing orf fasta and AA fasta: select largest isoform, fix name, print to out_fasta'''
    #load fasta
    records = list(SeqIO.parse(in_fasta, "fasta"))
    records = list(map(add_info, records)) #add orf, isoform, length to each record
    #get largest reading frame for each isoform
    #now go through records and get largest reading frame for an isoform
    largest_orf_isoforms = {}
    for cur_record in records:
        if(cur_record.isoform in list(largest_orf_isoforms.keys())):
            if cur_record.length <= largest_orf_isoforms[cur_record.isoform]['length']:
                continue
        largest_orf_isoforms[cur_record.isoform] = {'orf':cur_record.orf, 'length' : cur_record.length,
                                                'paralog' : cur_record.paralog, 'sequence' : str(cur_record.seq) }
    #write to fasta
    with open(out_fasta, 'w') as new_file:
        for cur_isoform, info_dict in largest_orf_isoforms.items():
            out_name = f"{snakemake.wildcards['SMP']}__{info_dict['paralog']}__{cur_isoform}"
            out_seq = info_dict['sequence']
            new_file.write(f'>{out_name}\n')
            new_file.write(f'{out_seq}\n')

#process ORF fa
process_fasta(snakemake.input['orf_fa'], snakemake.output['orf_fa'])
# Run samtools
exit_code = os.system(f"samtools faidx {snakemake.output['orf_fa']}")
# Check if the command was successful
if exit_code == 0:
    print("orf_fa processed and indexed.")
else:
    print(f"Error. Unable to faidx {snakemake.output['orf_fa']}")

#process AA fa
process_fasta(snakemake.input['aa_fa'], snakemake.output['aa_fa'])
# Run the bash command
exit_code = os.system(f"samtools faidx {snakemake.output['aa_fa']}")
# Check if the command was successful
if exit_code == 0:
    print("aa_fa processed and indexed.")
else:
    print(f"Error. Unable to faidx {snakemake.output['aa_fa']}")
