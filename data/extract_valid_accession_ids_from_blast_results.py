import pandas as pd
import numpy as np
import itertools as it


def get_seq_match_dict_and_flat_list(df):
    #extract protein identifier for matches: key=protein identifier for forward_input_sequences
    #and values=protein identifier for matches 
    seq_matches_dict = {}
    for i in df[0].unique():
        seq_matches_dict[i] = list(set(df[df[0] == i][1]))

    return seq_matches_dict

def extract_reciprocal_best_hits_and_return_protein_ids(seq_matches_backward_dict,seq_matches_forward_dict):
    result_set = []
    for forward_key in seq_matches_forward_dict.keys():
        for forward_value in seq_matches_forward_dict[forward_key]:
            if forward_value in seq_matches_backward_dict.keys():
                if forward_key in seq_matches_backward_dict[forward_value]:
                    print("[+] human: ",forward_value,"[+] mouse: ", seq_matches_backward_dict[forward_value])
                    result_set.append([forward_value,seq_matches_backward_dict[forward_value]])
    return result_set
    
forward_df = pd.read_csv(snakemake.input['fw_result'],delimiter=',',header=None)
backward_df = pd.read_csv(snakemake.input['bw_result'],delimiter=',',header=None)

seq_matches_forward_dict = get_seq_match_dict_and_flat_list(forward_df)
seq_matches_backward_dict = get_seq_match_dict_and_flat_list(backward_df)

best_hits = extract_reciprocal_best_hits_and_return_protein_ids(seq_matches_backward_dict,seq_matches_forward_dict)

out = open(snakemake.output[0],'w')
out.write('forward_genome_id\tbackward_genome_id\n')
for prot_id_pair in best_hits:
	out.write(str(prot_id_pair[0])+'\t'+str(prot_id_pair[1])+'\n')
out.close()
