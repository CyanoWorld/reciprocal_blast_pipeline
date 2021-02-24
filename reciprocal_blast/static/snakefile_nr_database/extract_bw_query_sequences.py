#backward query preparation
#this scripts reads the forward BLAST output table in order to write the subject accession identifier into a file that is used by the blastdbcmd program
import pandas as pd
fw_results = pd.read_table(snakemake.input['fw_res'], header=None)
out = open(snakemake.output['gi_list'],"w+")
for gi in fw_results[6][:]:
    out.write(gi+"\n")
out.close()
