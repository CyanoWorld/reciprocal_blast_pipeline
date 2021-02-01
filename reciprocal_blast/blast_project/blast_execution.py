from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences, \
    TaxNodesForBackwardDatabase, TaxNodesForForwardDatabase
from django.shortcuts import get_object_or_404
from django.db import IntegrityError
from shutil import copy
from os.path import isfile, isdir
from os import remove, chmod, listdir, getcwd, system
from datetime import datetime
from pathlib import Path
import subprocess

def write_nr_snakefile(project_id):
    try:
        project = get_object_or_404(BlastProject, pk=project_id)
        query_sequences = get_object_or_404(QuerySequences, associated_project=project)
        forward_db_organisms = TaxNodesForForwardDatabase.objects.filter(associated_project=project_id)
        backward_db_organisms = TaxNodesForBackwardDatabase.objects.filter(associated_project=project_id)
        fw_settings = get_object_or_404(ForwardBlastSettings, associated_project=project)
        bw_settings = get_object_or_404(BackwardBlastSettings, associated_project=project)
    except Exception as e:
        raise IntegrityError("[-] Couldn't get NR project specific data out of the database with exception: {}".format(e))
    try:
        if isdir('media/'+str(project_id)+'/'):
            #copy_nr_snakefile(project_id)
            write_nr_snakemake_configuration_and_taxid_file(project_id,query_sequences, forward_db_organisms, backward_db_organisms, fw_settings, bw_settings)
        else:
            raise ValueError("[-] Couldn't write snakemake config for NR project, there is no project folder.")
    except Exception as e:
        raise ValueError("[-] Couldn't write snakemake config for NR project with exception: {}".format(e))

def write_nr_snakemake_configuration_and_taxid_file(project_id, query_sequences, fw_taxids, bw_taxids, fw_settings, bw_settings):
    try:
        snakemake_config = open("media/" + str(project_id) + '/snakemake_config', 'w')
        snakemake_config.write("fw_query: \'query_sequences/"+query_sequences.query_file_name+"\'\n")
        snakemake_config.write("fw_word_size: " + str(fw_settings.word_size) + "\n")
        snakemake_config.write("fw_num_alignments: " + str(fw_settings.num_alignments) + "\n")
        snakemake_config.write("fw_evalue: " + str(fw_settings.e_value) + "\n")
        snakemake_config.write("fw_num_descriptions: " + str(fw_settings.num_descriptions) + '\n')
        snakemake_config.write("fw_max_target_seqs: " + str(fw_settings.max_target_seqs) + '\n')
        snakemake_config.write("fw_max_hsps: " + str(fw_settings.max_hsps) + '\n')
        snakemake_config.write("fw_num_threads: " + str(fw_settings.num_threads) + '\n')

        snakemake_config.write("bw_word_size: " + str(bw_settings.word_size) + "\n")
        snakemake_config.write("bw_num_alignments: " + str(bw_settings.num_alignments) + "\n")
        snakemake_config.write("bw_evalue: " + str(bw_settings.e_value) + "\n")
        snakemake_config.write("bw_num_descriptions: " + str(bw_settings.num_descriptions) + '\n')
        snakemake_config.write("bw_max_target_seqs: " + str(bw_settings.max_target_seqs) + '\n')
        snakemake_config.write("bw_max_hsps: " + str(bw_settings.max_hsps) + '\n')
        snakemake_config.write("bw_num_threads: " + str(bw_settings.num_threads) + '\n')

        snakemake_config.write("project_id: "+str(project_id)+"\n")

        snakemake_config.close()

        write_taxid_file(project_id,"taxonomic_nodes_forward.txt",fw_taxids)
        write_taxid_file(project_id,"taxonomic_nodes_backward.txt",bw_taxids)
    except Exception as e:
        raise ValueError("[-] Error creating project associated configuration files with Exception: {}".format(e))

def write_taxid_file(project_id,filename,taxonomic_nodes):
    try:
        taxid_file = open("media/" + str(project_id) + "/"+filename, "w")
        for organism in taxonomic_nodes:
            taxid_file.write(str(organism.taxonomic_node) + "\n")
        taxid_file.close()
    except Exception as e:
        raise Exception("[-] Coudn't write taxonomic node file with exception: {}".format(e))



def copy_nr_snakefile(project_id):
    try:
        copy('static/snakefile_nr_database/Snakefile', 'media/' + str(project_id))
        chmod("./media/" + str(project_id) + "/Snakefile", 777)

    except Exception as e:
        raise ValueError("[-] Failed copying snakefile to project dir with exception: {}".format(e))

#TODO one function for copying
def write_snakefile(project_id):
    #prepare project relevant data
    try:
        project = get_object_or_404(BlastProject, pk=project_id)
        dbtype = "prot" if project.search_strategy == "blastp" else "nucl"
        forward_genome = get_object_or_404(Genomes,associated_project=project, reciprocal_type='forward')
        backward_genome = get_object_or_404(Genomes,associated_project=project, reciprocal_type='backward')
        query_sequences = get_object_or_404(QuerySequences, associated_project=project)
        fw_settings = get_object_or_404(ForwardBlastSettings, associated_project=project)
        bw_settings = get_object_or_404(BackwardBlastSettings, associated_project=project)
    except Exception as e:
        raise IntegrityError("[-] Couldn't get project specific data out of the database with exception: {}".format(e))

    #check if media dir of project exists if not return Http404
    try:
        if isdir('media/' + str(project_id) + '/'):
            #write project specific snakemake configuration and grant permissions
            write_project_associated_snakefile_configuration_or_throw_valueerror(backward_genome, bw_settings, dbtype, forward_genome,
                                                                                 fw_settings, project, project_id, query_sequences)
            chmod("./media/" + str(project_id) + "/snakemake_config", 777)
        else:
            raise ValueError("[-] Couldn't write snakemake configuration, there is no project folder.")
    except Exception as e:
        raise ValueError("[-] Couldn't write snakemake configuration with exception: {}".format(e))

def write_project_associated_snakefile_configuration_or_throw_valueerror(backward_genome, bw_settings, dbtype, forward_genome,
                                                                         fw_settings, project, project_id, query_sequences):
    try:
        forward_genome_path = getcwd() + '\\media\\databases\\'+forward_genome.genome_name
        backward_genome_path = getcwd() + '\\media\\databases\\' + backward_genome.genome_name
        snakemake_config = open("media/" + str(project_id) + '/snakemake_config','w')
        snakemake_config.write("fw_db: \'"+forward_genome_path+"\'\n")
        snakemake_config.write("bw_db: \'"+backward_genome_path+"\'\n")
        snakemake_config.write("fw_query: \'query_sequences/"+query_sequences.query_file_name+"\'\n")
        snakemake_config.write("bw_query: \'query_sequences/input_backward.faa\'\n")

        snakemake_config.write("fw_word_size: "+str(fw_settings.word_size)+"\n")
        snakemake_config.write("fw_num_alignments: "+str(fw_settings.num_alignments)+"\n")
        snakemake_config.write("fw_evalue: "+str(fw_settings.e_value)+"\n")
        snakemake_config.write("fw_num_descriptions: "+str(fw_settings.num_descriptions)+'\n')
        snakemake_config.write("fw_max_target_seqs: "+str(fw_settings.max_target_seqs)+'\n')
        snakemake_config.write("fw_max_hsps: "+str(fw_settings.max_hsps)+'\n')
        snakemake_config.write("fw_num_threads: "+str(fw_settings.num_threads)+'\n')


        snakemake_config.write("bw_word_size: "+str(bw_settings.word_size)+"\n")
        snakemake_config.write("bw_num_alignments: "+str(bw_settings.num_alignments)+"\n")
        snakemake_config.write("bw_evalue: "+str(bw_settings.e_value)+"\n")
        snakemake_config.write("bw_num_descriptions: " + str(bw_settings.num_descriptions) + '\n')
        snakemake_config.write("bw_max_target_seqs: " + str(bw_settings.max_target_seqs) + '\n')
        snakemake_config.write("bw_max_hsps: " + str(bw_settings.max_hsps) + '\n')
        snakemake_config.write("bw_num_threads: " + str(bw_settings.num_threads) + '\n')

        snakemake_config.write("dbtype: "+dbtype+"\n")
        snakemake_config.write("search_type: "+str(project.search_strategy)+"\n")
        snakemake_config.write("project_id: "+str(project_id)+"\n")

        snakemake_config.close()
    except Exception as e:
        raise ValueError("[-] Failed writing snakemake_config to project dir with exception: {}".format(e))


def snakefile_exists(project_id):
    switch = True if isfile("media/" + str(project_id) + "/snakemake_config") else False
    return switch

#
#TODO: render dictionary based on snakemake config file for output in pipeline_dashboard
def view_builded_snakefile(project_id,nr_or_upload):
    #snakemake_config_file = 'media/' + str(project_id) + '/snakemake_config'
    if nr_or_upload == 'nr':
        snakefile_path = 'static/snakefile_nr_database/Snakefile'
    else:
        snakefile_path = 'static/snakefile_genome_upload/Snakefile'
    try:
        snakefile = open(snakefile_path,'r')
        lines = snakefile.readlines()
        snakefile.close()

        content = {}
        for line in lines[2:]:
            if 'rule' in line:
                rule = line
                content[rule] = []
            if line != rule:
                content[rule].append(line)

        return content
    except Exception as e:
        raise ValueError('[-] Error during reading Snakefile at: {} with Exception: {}'.format("media/" + str(project_id) + "/Snakefile", e))

def read_snakemake_logs(project_id):
    try:
        files = []
        time_modified = []
        logfile_path = 'media/'+str(project_id)+'/.snakemake/log'
        for file in listdir(logfile_path):
            files.append(file)
            current_file = Path('media/'+str(project_id)+'/.snakemake/log/' + file)
            timestamp = datetime.fromtimestamp(current_file.stat().st_mtime)
            time_modified.append(timestamp)
        latest_file_name = files[time_modified.index(max(time_modified))]

        content={}
        content['logfile'] = []
        latest_file = open('media/'+str(project_id)+'/.snakemake/log/'+latest_file_name)
        latest_percentage = 0
        for line in latest_file.readlines():
            content['logfile'].append(line)
            if 'steps' and '%' in line:
                latest_percentage =  line.split(' ')[4].split('(')[1].split(')')[0]
        latest_file.close()
        return content, latest_percentage
    except Exception as e:
        return [{'no_logs':"no logs are available, exception: ".format(e)}]


#don't output something in snakemake's invoked scripts
#TODO add error handling
def exec_snakemake(project_id):
    project = get_object_or_404(BlastProject, pk=project_id)
    if project.using_nr_database == True:
        snakemake_working_dir = 'media/'+str(project_id)+'/'
        snakemake_config_file = 'media/'+str(project_id)+'/snakemake_config'
        #print("[+] execute snakemake ...")
        #snakemake --snakefile 'C:\\Users\\lujeb\\Documents\\github_projects\\django_blast_project\\reciprocal_blast\\static\\snakefile_nr_database\\Snakefile' --configfile 'snakemake_config' --cores 2
        #subprocess.Popen(['snakemake','--dry-run','--cores','2','--wms-monitor','http://127.0.0.1:5000'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=snakemake_working_dir)
        subprocess.Popen(['snakemake','--snakefile','C:\\Users\\lujeb\\Documents\\github_projects\\django_blast_project\\reciprocal_blast\\static\\snakefile_nr_database\\Snakefile','--wms-monitor','http://127.0.0.1:5000','--cores','2','--configfile',snakemake_config_file,'--directory',snakemake_working_dir], shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    else:
        snakemake_working_dir = 'media/'+str(project_id)+'/'
        snakemake_config_file = 'media/'+str(project_id)+'/' + 'snakemake_config'
        # print("[+] execute snakemake ...")
        # subprocess.Popen(['snakemake','--dry-run','--cores','2','--wms-monitor','http://127.0.0.1:5000'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=snakemake_working_dir)
        #print("[+] WORKING DIR: "+getcwd())
        snakefile_dir = getcwd()+'\\static\\snakefile_genome_upload\\Snakefile'
        #print('\n[+] SNAKEFILE DIR: '+snakefile_dir,'\n[+] SNAKEMAKE WORKING DIR: '+snakemake_working_dir,'\n[+] SNAKEMAKE CONFIG DIR: '+snakemake_config_file)
        #snakemake --snakefile F:\programming\github_projects\bachelor_project_github\reciprocal_blast_pipeline\reciprocal_blast\static\snakefile_genome_upload\Snakefile --cores 2 --configfile media/1/snakemake_config --directory media/1/ --dry-run
        print('\n[+] COMMAND: snakemake --snakefile {} --cores 2 --configfile {} --directory {}\n'.format(snakefile_dir,snakemake_config_file,snakemake_working_dir))
        subprocess.Popen(['snakemake', '--snakefile', snakefile_dir, '--cores', '2', '--configfile',snakemake_config_file, '--directory', snakemake_working_dir], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
