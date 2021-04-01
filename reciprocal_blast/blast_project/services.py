'''
Content: functions for file uploading / loading procedures and database transactions

Use this file for additional database transactions or file uploading / loading procedures that are triggered in the views

additionally two functions are present for reading html results produced by snakemake
one function that checks if the pipeline has finished or not, this function is used in the detail view in order to load html results
'''
from django.db import IntegrityError
from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences, TaxNodesForForwardDatabase, TaxNodesForBackwardDatabase, RefseqGenome
from os import walk, mkdir
from os.path import isfile,  isdir
from shutil import rmtree

#provides the possibility to delete files if a user has accidently moved files into that directory
def delete_files_without_projects():
    try:
        all_project_ids = []
        for project in BlastProject.objects.all():
            all_project_ids.append(project.id)
        #list all projects in media folder except the database folder
        projects_in_media = next(walk('media/'))[1]
        projects_in_media.remove('databases')

        for folder in projects_in_media:
            if int(folder) not in all_project_ids:
                print("[+] removing folder {} ...".format(folder))
                rmtree('media/'+str(folder))
    except Exception as e:
        raise ValueError('[-] unable to remove project directory with exception: {}'.format(e))

#This won't work for some files on windows systems
def delete_project_files_by_project_id(project_id):
    try:
        if isdir('media/' + str(project_id) + '/'):
            rmtree('media/'+str(project_id)+'/')
        if isdir('static/result_images/'+str(project_id)):
            rmtree('static/result_images/'+str(project_id))
    except Exception as e:
        raise ValueError('[-] unable to remove project directory: {} with exception: {}'.format(project_id,e))

#this function is executed after deleting a project, thus it will pretend that uploaded genomes gets deleted
#if they are not connected to any project after deleting a project
def check_if_genomes_should_be_resaved(project_id):
    try:
        # if the associated genome database just has this project as its use, the genome database won't appear in the database anymore
        # so save the genome into the database without a project association
        genomes = Genomes.objects.filter(associated_project=project_id)
        # len(genomes) == 2 because there are always two genomes associated to genome upload projects
        if len(genomes) == 2:
            genome_first = Genomes(genome_name=genomes[0].genome_name, path_to_file=genomes[0].path_to_file)
            genome_second = Genomes(genome_name=genomes[1].genome_name, path_to_file=genomes[1].path_to_file)
            genome_first.save()
            genome_second.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't resave genome databases during deletion of the project with Exception: {}".format(e))

def save_project_from_form_or_raise_exception(new_title, new_strategy, user):
    try:
        project = BlastProject(project_title=new_title, search_strategy=new_strategy,project_username=user)
        project.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save project in database with exception: {}".format(e))
    return project

def save_nr_project_from_form_or_raise_exception(new_title, user):
    try:
        project = BlastProject(project_title=new_title, search_strategy='blastp', project_username=user, using_nr_database=True)
        project.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save project in database with exception: {}".format(e))
    return project

def save_query_file_in_db(query_sequences, project):
    try:
        uploaded_file_url_queries = 'media/'+str(project.id)+'/'+query_sequences.name
        new_query_sequences = QuerySequences(associated_project=project, query_file_name=query_sequences.name,
                                             path_to_query_file=uploaded_file_url_queries)
        new_query_sequences.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save query file in database with exception: {}".format(e))

def save_new_genome_file_in_db(genome_file):
    try:
        uploaded_file_url = 'media/databases/'+genome_file.name
        new_genome_database = Genomes(genome_name=genome_file.name,path_to_file=uploaded_file_url)
        new_genome_database.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save database files with exception: {}".format(e))


#function will be executed in the following order
def create_project_dir(project):
    try:
        mkdir('static/result_images/'+str(project.id))
        mkdir('media/' + str(project.id))
        #mkdir('media/'+str(project.id)+'/forward_genome')
        #mkdir('media/'+str(project.id)+'/backward_genome')
        mkdir('media/' + str(project.id) + '/query_sequences')
    except Exception as e:
        raise IntegrityError('[-] A very specific bad thing happened during creation of your project folder: {}'.format(e))

def create_nr_project_dir(project):
    try:
        mkdir('static/result_images/' + str(project.id))
        mkdir('media/' + str(project.id))
        mkdir('media/' + str(project.id) + '/query_sequences')
    except Exception as e:
        raise IntegrityError('[-] A very specific bad thing happened during creation of your project folder: {}'.format(e))

def upload_file(project_file,destination):
    try:
        with open(destination,'wb+') as dest:
            for chunk in project_file.chunks():
                dest.write(chunk)
    except Exception as e:
        raise IntegrityError('[-] A very specific bad thing happened during file upload of: {} Exception: {}'.format(project_file.name,e))


def save_refseq_genome_in_db(database_description,assembly_levels):
    pass

def save_genomes_and_query_in_db(query_sequences, forward_genome_name, backward_genome_name, project):
    try:
        uploaded_file_url_forward = 'media/databases/'+forward_genome_name
        uploaded_file_url_backward = 'media/databases/'+backward_genome_name
        uploaded_file_url_queries = 'media/'+str(project.id)+'/'+query_sequences.name
        new_forward_genome = Genomes(associated_project=project, genome_name=forward_genome_name,
                                 reciprocal_type='forward', path_to_file=uploaded_file_url_forward)
        new_backward_genome = Genomes(associated_project=project, genome_name=backward_genome_name,
                                  reciprocal_type='backward', path_to_file=uploaded_file_url_backward)
        new_query_sequences = QuerySequences(associated_project=project,query_file_name=query_sequences.name,
                                             path_to_query_file=uploaded_file_url_queries)
        new_forward_genome.save()
        new_backward_genome.save()
        new_query_sequences.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save database and query files with exception: {}".format(e))



def save_forward_settings_from_form_or_raise_exception(project,settings_form_forward):
    try:
        settings_fw = ForwardBlastSettings(associated_project=project, e_value=settings_form_forward['fw_e_value'],
                                           word_size=settings_form_forward['fw_word_size'], num_alignments=settings_form_forward['fw_num_alignments'],
                                            num_descriptions=settings_form_forward['fw_num_descriptions'],num_threads=settings_form_forward['fw_num_threads'])
        settings_fw.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save forward BLAST settings with exception: {}".format(e))

def save_backward_settings_from_form_or_raise_exception(project,settings_form_backward):
    try:
        settings_bw = BackwardBlastSettings(associated_project=project, e_value=settings_form_backward['bw_e_value'],
                                            word_size=settings_form_backward['bw_word_size'], num_alignments=settings_form_backward['bw_num_alignments'],
                                            num_descriptions=settings_form_backward['bw_num_descriptions'],num_threads=settings_form_backward['bw_num_threads'],
                                           max_hsps=settings_form_backward['bw_max_hsps'])
        settings_bw.save()
    except Exception as e:
        raise IntegrityError("[-] Couldn't save backward BLAST settings with exception: {}".format(e))

#TODO add function that resets the execution boolean in order to allow a re-execution of the pipeline if an error occured
def set_executed_on_true_and_save_project(current_project):
    try:
        current_project.pipeline_executed = True
        current_project.save()
    except Exception as e:
        raise IntegrityError('[-] Could not set pipeline_executed on true with Exception: {}'.format(e))

#taxids are a list of integers produced in the respective forms.py class
def validate_fw_taxids_and_save_into_database(project, user_email, taxids):
    for id in taxids:
            tax_fw_node = TaxNodesForForwardDatabase(associated_project=project, taxonomic_node=id)
            if tax_fw_node.if_valid_save_organism_name(user_email):
                tax_fw_node.save()
            else:
                raise IntegrityError("[-] Couldn't save taxonomic node {} for FW database into project database! There is no organism with such a name!".format(id))

#taxids are a list of integers produced in the respective forms.py class
def validate_bw_taxids_and_save_into_database(project, user_email, taxids):
    for id in taxids:
            tax_bw_node = TaxNodesForBackwardDatabase(associated_project=project, taxonomic_node=id)
            if tax_bw_node.if_valid_save_organism_name(user_email):
                tax_bw_node.save()
            else:
                raise IntegrityError("[-] Couldn't save taxonomic node {} for BW database into project database! There is no organism with such a name!".format(id))

#checks if the pipeline finishes (if result html table is present)
def snakemake_project_finished(project_id):
    if isfile('media/'+str(project_id)+'/reciprocal_results.html'):
        return True
    else:
        return False

#loads the reciprocal results table that is written with one of the last rules in the snakefiles
def get_html_results(project_id,filename):
    try:
        with open("media/"+str(project_id)+"/"+filename) as res:
            data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("[-] Couldn't read file {} with Exception: {}".format(e))

#loads the html graphs that are written with one of the last rules in the snakefiles
def load_html_graph(project_id,filename):
    try:
        html_file = open('media/'+str(project_id)+'/'+filename)
        html_lines = ''
        for line in html_file.readlines():
            html_lines += line
        html_file.close()
        return html_lines
    except Exception as e:
        raise FileNotFoundError("[-] Couldn't create result graph with Exception: {}".format(e))