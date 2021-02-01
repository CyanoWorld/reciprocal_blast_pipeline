import pandas as pd

frec_prot=pd.read_table(snakemake.input['rec_res'])
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns=["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                  "stitle"]

fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])


for i in range(0, len(fw_res), 1):
    taxids = fw_res.iat[i, 7]
    scientific_names = fw_res.iat[i, 8]
    common_names = fw_res.iat[i, 9]

    if (len(taxids.split(";")) > 3):
        taxid_string = ';'.join(taxids.split(";")[0:3]) + "..."
        scientific_names_string = ';'.join(scientific_names.split(";")[0:3]) + "..."
        common_names_string = ';'.join(common_names.split(";")[0:3]) + "..."
        fw_res.iloc[i, 7] = taxid_string
        fw_res.iloc[i, 8] = scientific_names_string
        fw_res.iloc[i, 9] = common_names_string

rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
result_data = rec_prot.merge(fw_res,how='inner', on=['sacc','qseqid'])
result_data = result_data.drop_duplicates('sacc', keep='first')

pd.set_option('colheader_justify', 'left')
html_string = '''
<html>
  <head>
    <title>BLAST Result Table</title>
    <style>
    .mystyle {{
        font-size: 11pt; 
        font-family: Arial;
        border-collapse: collapse; 
        border: 1px solid silver;

    }}

    .mystyle td, th {{
        padding: 5px;
        max-width: 200px;
    }}

    .mystyle tr:nth-child(even) {{
        background: #E0E0E0;
    }}

    .mystyle tr:hover {{
        background: silver;
        cursor: pointer;
    }}
    </style>

  </head>
  <body>

    {table}
  </body>
</html>
'''

# OUTPUT AN HTML FILE
# with open('fw_results.html', 'w') as f:
#    f.write(html_string.format(table=fw_res.to_html(classes='mystyle')))

# with open('bw_results.html', 'w') as f:
#    f.write(html_string.format(table=bw_res.to_html(classes='mystyle')))

with open(snakemake.output['rec_html'], 'w') as f:
    f.write(html_string.format(table=result_data.to_html(classes='mystyle')))