import pandas as pd
import itertools as it

forward_df = pd.read_csv('forward_blast_output.csv',delimiter=',',header=None)

def get_seq_match_dict_and_flat_list(df):
    #extract protein identifier for matches: key=protein identifier for forward_input_sequences
    #and values=protein identifier for matches 
    seq_matches_dict = {}
    for i in df[0].unique():
        seq_matches_dict[i] = list(set(df[df[0] == i][1]))

    #flat list of protein identifier from matches (database genome) to the query sequences
    matches = [prot_id for prot_id in seq_matches_dict.values()]
    matches = list(it.chain.from_iterable(matches))
    return seq_matches_dict,matches

#checks if line contains sequence information (fasta header)
def true_if_line_is_fasta_header(first_line_symbol):
    check = True if first_line_symbol == '>' else False
    return check

#checks if next sequence is a match of the forward blast
def true_if_identifier_is_in_forward_matches(prot_id,forward_matches):
    check = True if prot_id in forward_matches else False
    return check

#check if current line is still part of sequence that should be written into the input for backward blast
#this is true as long as the iteration over the forward_blast_db is on a matched sequence and false if the iteration
#is on a sequence that has not matched with the query sequences
def line_is_part_of_matched_sequence(prot_id,forward_matches):
    check = True if true_if_identifier_is_in_forward_matches(prot_id,forward_matches) else False
    return check

#function for preparing input of the backward blast
def write_file_for_backward_blast_based_on_matches_of_forward_blast(path_to_forward_database, forward_matches):
    backward_input = open('./query_sequences/input_backward.faa', 'w')
    with open(path_to_forward_database) as forward_db_genome:
        
        write_line_to_output = False
        for line in forward_db_genome:
            if true_if_line_is_fasta_header(line[0]):
                prot_id = line[1:].split(' ')[0]
                write_line_to_output = line_is_part_of_matched_sequence(prot_id,forward_matches)
                if write_line_to_output:
                    print("[+] writing {} sequence into input_backward.faa".format(prot_id))
            if write_line_to_output:
                backward_input.write(line)
                
    backward_input.close()
    
seq_matches_forward_dict,forward_matches = get_seq_match_dict_and_flat_list(forward_df)
write_file_for_backward_blast_based_on_matches_of_forward_blast('./forward_genome/database/fw_prot_db.faa',forward_matches)
