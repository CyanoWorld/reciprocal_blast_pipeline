#this script writes the RBHs identified by the reciprocal BLAST pipeline into an html table
import pandas as pd
rec_prot=pd.read_table(snakemake.input['rec_res'])
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns = ["qseqid","sseqid","pident","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
rec_prot = rec_prot.rename(columns={"forward_genome_id": "sseqid"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
result_data = rec_prot.merge(fw_res,how='inner', on=['sseqid','qseqid'])
result_data = result_data.drop_duplicates('sseqid', keep='first')

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

#with open('fw_results.html', 'w') as f:
#    f.write(html_string.format(table=fw_res.to_html(classes='mystyle')))

#with open('bw_results.html', 'w') as f:
#    f.write(html_string.format(table=bw_res.to_html(classes='mystyle')))

with open('reciprocal_results.html', 'w') as f:
    f.write(html_string.format(table=result_data.to_html(classes='mystyle')))