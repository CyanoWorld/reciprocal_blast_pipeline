from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences
from django.shortcuts import get_object_or_404
from django.http import Http404
from shutil import copy
from os.path import isfile, isdir
from os import remove, chmod
import subprocess

def write_snakefile(project_id):
    #prepare project relevant data
    project = get_object_or_404(BlastProject, pk=project_id)
    dbtype = "prot" if project.search_strategy == "blastp" else "nucl"
    forward_genome = get_object_or_404(Genomes,associated_project=project, reciprocal_type='forward')
    backward_genome = get_object_or_404(Genomes,associated_project=project, reciprocal_type='backward')
    query_sequences = get_object_or_404(QuerySequences, associated_project=project)
    fw_settings = get_object_or_404(ForwardBlastSettings, associated_project=project)
    bw_settings = get_object_or_404(BackwardBlastSettings, associated_project=project)

    #check if media dir of project exists if not return Http404
    if isdir('media/' + str(project_id) + '/'):
        #edit and copy relevant python scripts from static folder to project folder

        edit_process_fw_blast_results_with_real_fw_genome_file(
            'static/snakefile_templates/process_fw_blast_results.py', 'media/' + str(project_id) + '/process_fw_blast_results.py',
            3,'path_to_fw_genome='+'\"./forward_genome/database/'+forward_genome.genome_name+"\"\n")
        #ensure execution rights for snakemake
        chmod("./media/" + str(project_id) + "/process_fw_blast_results.py", 777)
        copy_post_processing_script_accession_ids(project_id)
        chmod("./media/" + str(project_id) + "/extract_valid_accession_ids_from_blast_results.py", 777)

        #write project associated snakefile
        write_project_associated_snakefile_or_throw_valueerror(backward_genome, bw_settings, dbtype, forward_genome,
                                                               fw_settings, project, project_id, query_sequences)
    else:
        return Http404()


def write_project_associated_snakefile_or_throw_valueerror(backward_genome, bw_settings, dbtype, forward_genome,
                                                           fw_settings, project, project_id, query_sequences):
    try:
        snakefile = open("./media/" + str(project_id) + "/Snakefile", "w")
        snakefile.write('rule all:\n')
        snakefile.write('\tinput:\n')
        snakefile.write('\t\t\"reciprocal_best_hits_protein_ids.txt\"\n\n')

        snakefile.write("rule make_forward_database:\n")
        snakefile.write("\tinput: \"" + "forward_genome/" + forward_genome.genome_name + "\"\n")
        snakefile.write('\toutput: \"forward_genome/database/' + forward_genome.genome_name + '\"\n')
        snakefile.write(
            "\tshell: \"makeblastdb -in {{input}} -dbtype {} -out {{output}} && mv {{input}} forward_genome/database/\"\n\n".format(
                dbtype))

        snakefile.write("rule forward_blast:\n")
        snakefile.write(
            "\tinput: query=\"query_sequences/" + query_sequences.query_file_name + "\", db=\"forward_genome/database/" + forward_genome.genome_name + "\"\n")
        snakefile.write("\toutput: \"forward_blast_output.csv\"\n")
        snakefile.write(
            "\tshell: \"" + project.search_strategy + " -db " + "{input.db}" + " -query {input.query} " + " -outfmt 10 " + "-word_size "
            + str(fw_settings.word_size) + " -evalue " + str(fw_settings.e_value) + " -num_alignments " + str(
                fw_settings.num_alignments) + " -out {output} " + "\"\n\n")

        snakefile.write("rule get_backward_sequences:\n")
        snakefile.write("\tinput: \"forward_blast_output.csv\"\n")
        snakefile.write("\toutput: \"query_sequences/input_backward.faa\"\n")
        snakefile.write("\tscript: \n")
        snakefile.write("\t\t\"process_fw_blast_results.py\"\n\n")

        snakefile.write("rule make_backward_database:\n")
        snakefile.write("\tinput: \"" + "backward_genome/" + backward_genome.genome_name + "\"\n")
        snakefile.write('\toutput: \"backward_genome/database/' + backward_genome.genome_name + '\"\n')
        snakefile.write(
            "\tshell: \"makeblastdb -in {{input}} -dbtype {} -out {{output}} && mv {{input}} backward_genome/database/\"\n\n".format(
                dbtype))

        snakefile.write("rule backward_blast:\n")
        snakefile.write(
            "\tinput: query=\"query_sequences/input_backward.faa\", db=\"backward_genome/database/" + backward_genome.genome_name + "\"\n")
        snakefile.write("\toutput: \"backward_blast_output.csv\"\n")
        snakefile.write(
            "\tshell: \"" + project.search_strategy + " -db " + "{input.db}" + " -query {input.query} " + " -outfmt 10 " + "-word_size "
            + str(bw_settings.word_size) + " -evalue " + str(bw_settings.e_value) + " -num_alignments " + str(
                bw_settings.num_alignments) + " -out {output} " + "\"\n\n")

        snakefile.write("rule get_reciprocal_best_hits_proteins_ids:\n")
        snakefile.write("\tinput: \"backward_blast_output.csv\", \"forward_blast_output.csv\"\n")
        snakefile.write("\toutput: \"reciprocal_best_hits_protein_ids.txt\"\n")
        snakefile.write("\tscript: \n")
        snakefile.write("\t\t\"extract_valid_accession_ids_from_blast_results.py\"\n\n")

        snakefile.close()
        chmod("./media/" + str(project_id) + "/Snakefile",777)
    except Exception as e:
        raise ValueError("[-] Failed writing snakefile to disc with exception: {}".format(e))


#copy relevant snakemake file into project directory
def copy_post_processing_script_accession_ids(project_id):
    try:
        copy('static/snakefile_templates/extract_valid_accession_ids_from_blast_results.py', 'media/' + str(project_id))
    except Exception as e:
        raise ValueError(
            '[-] Couldnt copy: extract_valid_accession_ids_from_blast_results.py with exception: {}'.format(e))


def snakefile_exists(project_id):
    switch = True if isfile("media/" + str(project_id) + "/Snakefile") else False
    return switch

#render dictionary based ib snakefile for output in pipeline_dashboard
def view_builded_snakefile(project_id):
    try:
        snakefile = open("media/" + str(project_id) + "/Snakefile")
        lines = snakefile.readlines()
        snakefile.close()

        content = {}
        for line in lines:
            if 'rule' in line:
                rule = line
                content[rule] = []
            if line != rule:
                content[rule].append(line)

        return content
    except Exception as e:
        raise ValueError('[-] Error during reading Snakefile at: {} with Exception: {}'.format("media/" + str(project_id) + "/Snakefile", e))

def delete_snakefile(project_id):
    if snakefile_exists(project_id):
        remove("media/" + str(project_id) + "/Snakefile")
        remove("media/" + str(project_id) + "/extract_valid_accession_ids_from_blast_results.py")
        remove("media/" + str(project_id) + "/process_fw_blast_results.py")
    else:
        raise ValueError('[-] Failed deleting snakemake files: {}'.format("media/" + str(project_id) + "/Snakefile"))

#this function is used for edititng the process_fw_blast_results.py for the snakefile
#it will write the path to the forward genome int line 4
def edit_process_fw_blast_results_with_real_fw_genome_file(file_name,path_to_new_file, line_num, text):
    try:
        scriptfile = open(file_name, 'r')
        lines = scriptfile.readlines()
        scriptfile.close()
        lines[line_num] = text
        out = open(path_to_new_file, 'w')
        out.writelines(lines)
        out.close()
    except Exception as e:
        raise ValueError("[-] Error during editing {} with exception: {}".format(file_name,e))


def exec_snakemake(project_id):
    snakemake_working_dir = 'media/'+str(project_id)+'/'
    #chdir(snakemake_working_dir)
    subprocess.Popen(['snakemake','--cores','2'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=snakemake_working_dir)
    #system('snakemake --cores 2')
    #chdir('../..')
    #subprocess.Popen(['snakemake', '--dry-run'], cwd=snakemake_working_dir, shell=True)
    #subprocess.Popen(['blastn','-h'], cwd=snakemake_working_dir, shell=True)