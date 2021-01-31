import pandas as pd
fw_res = pd.read_table(snakemake.input['fw_res'],header=None)
#bw_res = pd.read_table("backward_blast_output.csv",header=None)
fw_res.columns = ["qseqid","sseqid","pident","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
#bw_res.columns = ["qseqid","sseqid","pident","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

rec_prot = pd.read_table(snakemake.input['rec_res'])
rec_prot['qseqid'] = rec_prot['backward_genome_id']
rec_prot['sseqid'] = rec_prot['forward_genome_id']
result_data = rec_prot.merge(fw_res, on=['sseqid','qseqid'])

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