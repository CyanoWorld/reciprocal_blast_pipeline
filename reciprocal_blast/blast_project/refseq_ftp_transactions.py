from gzip import decompress
from os import remove, chdir, getcwd
from os.path import isfile
from wget import download
import pandas as pd

#downloads the current refseq assembly file into an specified directory
def download_current_assembly_summary_into_specific_directory(directory):
    refseq_url = "ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    try:
        #/blast/reciprocal_blast
        current_working_directory = getcwd()
        chdir(directory)
        downloaded_file = download(refseq_url)
        if isfile(downloaded_file) == False:
            raise FileExistsError
        chdir(current_working_directory)
    except Exception as e:
        raise FileExistsError("[-] Couldn't download assembly_summary_refseq.txt file ...\n\tException: {}".format(e))

    return downloaded_file

#deletes the downloaded assembly summary file - returns True if successfull else raise an error
def delete_downloaded_assembly_summary(directory):
    if(isfile(directory+'assembly_summary_refseq.txt')):
        try:
            remove(directory+'assembly_summary_refseq.txt')
        except Exception as e:
            raise OSError("[-] Couldn't delete assembly_summary_refseq.txt file ...\n\tException: {}".format(e))
    return True

def refseq_file_exists(directory):
    return isfile(directory+'assembly_summary_refseq.txt')

#completeness_level = 'Chromosome', 'Scaffold', 'Complete Genome', 'Contig'
def read_current_assembly_summary_with_pandas(summary_file_path,completeness_level):

    try:
        refseq_table = pd.read_table(summary_file_path, skiprows=[0, 1], header=None)
        header = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid",
                  "species_taxid", "organism_name", "infraspecific_name", "isolate", "version_status", "assembly_level",
                  "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm",
                  "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material"]
        refseq_table.columns = header
    except Exception as e:
        raise ValueError("[-] Exception during pandas parsing of assembly_summary_refseq.txt file ...\n\tException: {}".format(e))

    try:
        summary_dict = {}
        accession = []
        organism_name = []
        for index, row in refseq_table[refseq_table['assembly_level'] == 'Complete Genome'].iterrows():
            protein_genome = row['ftp_path'].split('/')[-1:][0]
            protein_genome = row['ftp_path'] + '/' + str(protein_genome) + '_protein.faa.gz'
            accession.append(row['assembly_accession'])
            organism_name.append(str(row['assembly_accession']) + " " + str(row['organism_name']))
            summary_dict[row['assembly_accession']] = [protein_genome, row['taxid'], row['species_taxid'],
                                                       row['organism_name']]
        html_input_list = tuple(zip(accession, organism_name))
    except Exception as e:
        raise ValueError("[-] Exception during creation of dictionary for assembly_summary_refseq.txt file ...\n\tException: {}".format(e))
    return summary_dict, html_input_list