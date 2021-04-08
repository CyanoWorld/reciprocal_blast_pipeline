import gzip
from os import remove, chdir, getcwd, waitpid
from os.path import isfile, isdir
from wget import download
import json
import pandas as pd
import re
import psutil
from .services import upload_file, create_and_save_refseq_database_model, create_refseq_genome_directory, \
    write_pandas_table_to_project_dir, create_and_save_refseq_transaction,get_refseq_database_title


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

#TODO: extend function for reading just specific genome files that can serve as input for the form ?
#TODO: add parsing of taxid files
#completeness_level = 'Chromosome', 'Scaffold', 'Complete Genome', 'Contig'
def read_current_assembly_summary_with_pandas(summary_file_path, refseq_level_checklist):
    #function for changing the ftp_header in the pandas table
    def set_protein_assembly_file(ftp_path):
        protein_genome = ftp_path.split('/')[-1:][0]
        protein_genome = ftp_path + '/' + str(protein_genome) + '_protein.faa.gz'
        return protein_genome

    #init parsing refseq table with pandas
    try:
        refseq_table = pd.read_table(summary_file_path+'assembly_summary_refseq.txt', skiprows=[0, 1], header=None)
        header = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid",
                  "species_taxid", "organism_name", "infraspecific_name", "isolate", "version_status", "assembly_level",
                  "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm",
                  "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material"]
        refseq_table.columns = header
    except Exception as e:
        raise ValueError("[-] Exception during pandas parsing of assembly_summary_refseq.txt file ...\n\tException: {}".format(e))

    #extract necessary data fields: assembly number, names, taxids and the correct ftp_filepath for downloading with gzip
    try:
        refseq_table = refseq_table[['assembly_accession', 'organism_name', 'taxid', 'species_taxid','assembly_level', 'ftp_path']]
        refseq_table['ftp_path'] = refseq_table['ftp_path'].apply(lambda row: set_protein_assembly_file(row))

        pandas_genome_level_dataframes = []
        for genome_level in refseq_level_checklist:
            pandas_genome_level_dataframes.append(refseq_table[refseq_table['assembly_level'] == genome_level])

        desired_refseq_genomes_dataframe = pd.concat(pandas_genome_level_dataframes)

        #tuple list for dropdown menu, not implemented yet
        #html_input_list = tuple(zip(refseq_table['assembly_accession'], refseq_table['organism_name']))
    except Exception as e:
        raise ValueError("[-] Exception during extraction of smaller dataframe from refseq_table dataframe ...\n\tException: {}".format(e))
    return [desired_refseq_genomes_dataframe]

def transform_data_table_to_json_dict(df):
    json_records = df.reset_index().to_json(orient='records')
    data = []
    data = json.loads(json_records)
    return data


# filter refseq table with taxids optained by the get_species_taxids.sh script
def read_taxonomy_table(filepath):
    if (isfile(filepath) == False):
        raise Exception("[-] There is no taxonomy file called: {}".format(filepath))
    taxonomy_file = pd.read_table(filepath, header=None)
    # species_taxid and taxid should normally be interchangeable, the species_taxid may inherit more informations
    # to current strain (have a look at the README description of the refseq summary file)
    taxonomy_file.columns = ['species_taxid']
    return taxonomy_file


def filter_table_by_taxonomy(refseq_table, taxonomy_table):
    # species_taxid
    return refseq_table.merge(taxonomy_table, how='inner', on=['species_taxid'])


# returns the amount of filtered genomes
def reduction_amount(refseq_table, filtered_table):
    reduct = len(refseq_table) - len(filtered_table)
    if (reduct <= 0):
        raise Exception("[-] After filtering, there is no data remaining!")
    return reduct


# downloading genomes with wget
def download_genome_from_ftp_path(ftp_path):
    try:
        return download(ftp_path)
    except:
        return None


# function with two responsibilities in order to save memory usage
def extract_downloaded_file_and_write_taxid_file(genome_downloaded_file, taxid):
    # decompression to fasta file
    try:
        output = open(genome_downloaded_file + ".decompressed.faa", "w")
        with gzip.open(genome_downloaded_file, "rb") as bytes_out:
            bytes_from_file = bytes_out.read()
            line = bytes_from_file.decode("utf-8")
            output.write(line)
            output.close()
        # remove gz file
        remove(genome_downloaded_file)
    except Exception as e:
        output.close()
        raise Exception("[-] Exception during decompressing: {}".format(e))

    # creation of taxmap file | taxid \t acc_id
    try:
        accession_id_pattern = re.compile('>(\S*)')
        taxmap = open('acc_taxid_map.table', 'a')
        for acc_id in re.findall(accession_id_pattern, line):
            taxmap.write(str(taxid) + "\t" + str(acc_id) + "\n")
        taxmap.close()
        # print("\t[*] Done writing taxonomic informations into taxmap file ...")
        return True 
    except Exception as e:
        taxmap.close()
        raise Exception("[-] Exception during writing taxmap file: {}".format(e))


def download_assemblies(filtered_table):
    try:
        for genome_url, taxid in zip(filtered_table['ftp_path'], filtered_table['species_taxid']):
            genome_file = download_genome_from_ftp_path(genome_url)
            #need additional import (import time)
            #time.sleep(1)
            # genome_file = wget.download(genome_url)
            #print("[+] downloaded genome: {}\n[+] taxonomic node: {}".format(genome_url, taxid))
            if (genome_file):
                extract_downloaded_file_and_write_taxid_file(genome_file, taxid)
    except Exception as e:
        #print("[-] ERROR during download_assemblies")
        raise Exception("[-] Error during downloading assemblies with Exception: {}".format(e))

#function that is triggered by POST submit in download_refseq_genomes.html
def refseq_download_project(refseq_form):
    #check if a taxid file is given in order to limit the refseq download process by taxonomy
    if refseq_form.cleaned_data.get('taxid_file', False):

        print("USING NEW FILE")
        #upload taxonomic information file
        taxid_file = refseq_form.cleaned_data['taxid_file']
        taxid_file_path = 'media/' + 'databases/' + 'refseq_databases/' + 'taxonomic_node_files/' + taxid_file.name
        upload_file(taxid_file,taxid_file_path)

        #get data from refseq_summary_file and limit by assembly_level (=refseq_levels)
        refseq_table = read_current_assembly_summary_with_pandas("./static/refseq_summary_file/",
                                                                 refseq_form.cleaned_data['refseq_levels'])[0]
        #limit file by taxonomy
        taxonomy_table = read_taxonomy_table(taxid_file_path)
        filtered_table = filter_table_by_taxonomy(refseq_table, taxonomy_table)

        if len(filtered_table) == 0:
            raise Exception("[-] The database doesnt contain any entries, pls apply an other filter method!")

        #create new refseq database model and save it into the database
        new_refseq_model = create_and_save_refseq_database_model(database_description=refseq_form.cleaned_data['database_name'],
                                                                 assembly_levels=refseq_form.cleaned_data['refseq_levels'],
                                                                 assembly_entries=len(filtered_table),
                                                                 attached_taxonomic_file=taxid_file_path)
        #create directory in media/databases/refseq_databases
        refseq_database_table_path = create_refseq_genome_directory(new_refseq_model.id)

        write_pandas_table_to_project_dir(refseq_database_table_path,
                                          filtered_table,
                                          refseq_form.cleaned_data['database_name'])

    elif refseq_form.cleaned_data.get('taxid_uploaded_file',False):
        print("USING UPLOADED FILE")

        taxid_file = refseq_form.cleaned_data['taxid_uploaded_file']
        taxid_file_path = 'media/' + 'databases/' + 'refseq_databases/' + 'taxonomic_node_files/' + taxid_file

        refseq_table = read_current_assembly_summary_with_pandas("./static/refseq_summary_file/",
                                                                 refseq_form.cleaned_data['refseq_levels'])[0]
        taxonomy_table = read_taxonomy_table(taxid_file_path)
        filtered_table = filter_table_by_taxonomy(refseq_table, taxonomy_table)

        if len(filtered_table) == 0:
            raise Exception("[-] The database doesnt contain any entries, pls apply an other filter method!")

            # create new refseq database model and save it into the database
        new_refseq_model = create_and_save_refseq_database_model(
            database_description=refseq_form.cleaned_data['database_name'],
            assembly_levels=refseq_form.cleaned_data['refseq_levels'],
            assembly_entries=len(filtered_table),
            attached_taxonomic_file=taxid_file_path)
        # create directory in media/databases/refseq_databases
        refseq_database_table_path = create_refseq_genome_directory(new_refseq_model.id)

        write_pandas_table_to_project_dir(refseq_database_table_path,
                                          filtered_table,
                                          refseq_form.cleaned_data['database_name'])
    else:
        print("NO TAXONOMIC FILTERING!")
        filtered_table = read_current_assembly_summary_with_pandas("./static/refseq_summary_file/",
                                                                 refseq_form.cleaned_data['refseq_levels'])[0]
        new_refseq_model = create_and_save_refseq_database_model(
            database_description=refseq_form.cleaned_data['database_name'],
            assembly_levels=refseq_form.cleaned_data['refseq_levels'],
            assembly_entries=len(filtered_table))

        refseq_database_table_path = create_refseq_genome_directory(new_refseq_model.id)

        write_pandas_table_to_project_dir(refseq_database_table_path,
                                          filtered_table,
                                          refseq_form.cleaned_data['database_name'])

    print("[+] Created Table")
    print("[*] \t length filtered table: {}".format(len(filtered_table)))


def trigger_database_download(database_id):
    title = get_refseq_database_title(database_id)
    title = title.replace(' ', '_').upper()
    table = 'media/databases/refseq_databases/'+str(database_id)+'/'+title
    working_dir = 'media/databases/refseq_databases/'+str(database_id)
    static_file_path = 'static/python_scripts/download_refseq_genomes.py'
    if(isfile(static_file_path) and isdir(working_dir) and isfile(table)):
        try:
            pid = psutil.Popen(['python','static/python_scripts/download_refseq_genomes.py','--ifile',table,'--outputdir',working_dir],stdout=open(working_dir+'/'+'download.log','w'))
            waitpid(pid.pid,0)
            #create_and_save_refseq_transaction(pid.pid,database_id,'DOWNLOAD DATABASE ID {} REFSEQ TRANSACTION'.format(database_id))
        except Exception as e:
            raise Exception("[-] ERROR during database download in trigger_database_download. Exception: {}".format(e))