#select largest reading frames and fix names of reading frames.
from Bio import SeqIO
import os

#step 1 : process gff
def read_gff_lines(gff_path):
    '''read in gff lines and skip comment lines'''
    gff_lines = []
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                gff_lines.append(line.rstrip('\n'))
    return gff_lines

gff_lines = read_gff_lines(snakemake.input['gff'])
if(len(gff_lines) == 0):
    print(f"{snakemake.input['gff']} is empty. Opening empty file for {snakemake.output['fa']} and fai files.")
    with open(snakemake.output['fa'], 'w'):
        pass
    with open(snakemake.output['fai'], 'w'):
        pass
    break

#link transcript dict to paralog
par_iso_dict = {}
for line in gff_lines:
    columns = line.strip().split('\t')
    col_9= columns[8].split(";")
    col_9_dict = { x.split("=")[0]:x.split("=")[1] for x in col_9 }
    par_iso_dict[col_9_dict["transcript_id"]] = col_9_dict["paralog"]

#step 2 : process fasta
records = list(SeqIO.parse(fa, "fasta"))
#get new name
for record in records:
    transcript_id = record.id.split('_')[0] #ORFs and AA have an ending _ORF_10 that I need to remove
    paralog = par_iso_dict[transcript_id] # 
    new_name = f"{snakemake.wildcards['SMP']}__{paralog}__{transcript_id}"
    record.id = new_names
    record.description = ""
#write file
SeqIO.write(records, snakemake.output['fa'] , "fasta")

#step 3: index new fasta
os.system(f"samtools faidx {snakemake.output['fa']}")

print("Done.")

