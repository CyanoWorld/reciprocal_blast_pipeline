#this script is used by the snakemake pipeline in order to produce result images
#currently it produces a subplot with statistical informations, based on the forward BLAST output of RBHs
import pandas as pd
import matplotlib.pyplot as plt, mpld3

rec_prot=pd.read_table(snakemake.input['rec_res'])
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns = ["qseqid","sseqid","pident","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
rec_prot = rec_prot.rename(columns={"forward_genome_id": "sseqid"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
result_data = rec_prot.merge(fw_res,how='inner', on=['sseqid','qseqid'])
result_data = result_data.drop_duplicates('sseqid', keep='first')

fig, axs = plt.subplots(3,figsize=(10,8))
axs[0].hist(result_data['pident'])
axs[0].set_title("percentage identity")
axs[1].hist(result_data['evalue'],color='r')
axs[1].set_title("evalue")
axs[2].hist(result_data['bitscore'],color='g')
axs[2].set_title("bitscore")
#fig.text(0.04, 0.5, 'sequences', va='center', rotation='vertical')
plt.tight_layout()
plt.savefig(snakemake.params['statistics'])
plt.savefig(snakemake.output[0])
mpld3.save_html(fig,snakemake.output[1])