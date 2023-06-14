import pandas as pd

gff_file = snakemake.input['locus_gff']

#get paralog isoforms
paralog_isoforms = {}
with open(gff_file, 'r') as file:
    # Open the output GFF file for writing
    for line in file:
        # Skip comment lines starting with '#'
        if line.startswith('#'):
            continue
        # Split the line by tab to separate the columns
        columns = line.strip().split('\t')
        # Split the 9th column by semicolon to separate the key-value pairs
        attributes = columns[8].strip().split(';')
        # Iterate through the key-value pairs
        for attribute in attributes:
            if attribute.strip() == "":
                continue
            # Split each key-value pair by '=' to get the key and value
            key, value = attribute.strip().split('=')
            # Check if the attribute key and value match your desired criteria
            if key == 'paralog':
                # Write the line to the output file
                cur_paralog = value
            elif key == 'transcript_id':
                cur_isoform = value
        if cur_paralog not in list(paralog_isoforms.keys()):
            paralog_isoforms[cur_paralog] = set([cur_isoform])
        else:
            paralog_isoforms[cur_paralog].add(cur_isoform)

#get abundance of isoforms
abundance_file = snakemake.input.abundance_tbl
abundance_df = pd.read_csv(abundance_file, comment="#", sep = "\t")
abundance_dict = dict(zip(abundance_df['pbid'], abundance_df['count_fl']))

#get top isoforms for each paralog
out_tbl = pd.DataFrame()
for cur_paralog, isoforms in paralog_isoforms.items():
    max_value_key = None
    max_value = float('-inf')
    for key in list(isoforms):
        assert key in abundance_dict , f"Error: not finding {key} in abundance dictionary"
        if abundance_dict[key] > max_value:
            max_value_key = key
            max_value = abundance_dict[key]
    out_ser = pd.Series({"paralog" : cur_paralog , "top_isoform" : max_value_key , "n_reads" : max_value })
    out_tbl = pd.concat( [out_tbl, out_ser.to_frame().T ])

#write to output
out_tbl.to_csv(path_or_buf=snakemake.output.tbl ,sep = "\t", header=True, index = False )