from django.core.files.storage import FileSystemStorage


from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences
from os import listdir, walk, mkdir
from os.path import isfile, join, isdir
from shutil import rmtree

#TODO add exceptions and error handling
#TODO deprecated do deleting while deleting project
def delete_files_without_projects():
    all_project_ids = []
    for project in BlastProject.objects.all():
        all_project_ids.append(project.id)
    projects_in_media = next(walk('media/'))[1]
    print(projects_in_media)
    print(all_project_ids)
    for folder in projects_in_media:
        if int(folder) not in all_project_ids:
            print("[+] removing folder {} ...".format(folder))
            rmtree('media/'+str(folder))

#
def delete_genomfiles_by_project_id(project_id):
    if isdir('media/' + str(project_id) + '/'):
        rmtree('media/'+str(project_id)+'/')

#TODO return error pages
def save_project_from_form_or_raise_exception(new_title, new_strategy, user):
    project = BlastProject(project_title=new_title, search_strategy=new_strategy,project_username=user)
    try:
        project.save()
        return project
    except Exception as e:
        raise ValueError('A very specific bad thing happened during project saving: {}'.format(e))

#function will be executed in the following order
def create_project_dir(project):
    try:
        mkdir('media/' + str(project.id))
        mkdir('media/'+str(project.id)+'/forward_genome')
        mkdir('media/'+str(project.id)+'/backward_genome')
        mkdir('media/' + str(project.id) + '/query_sequences')
    except Exception as e:
        raise ValueError('A very specific bad thing happened during file upload: {}'.format(e))

def upload_forward_genome_file(sequence_file, project):
    try:
        with open('media/'+str(project.id)+"/forward_genome/"+sequence_file.name, 'wb+') as destination:
            for chunk in sequence_file.chunks():
                destination.write(chunk)
    except Exception as e:
        raise ValueError('A very specific bad thing happened during file upload: {}'.format(e))

def upload_backward_genome_file(sequence_file, project):
    try:
        with open('media/'+str(project.id)+"/backward_genome/"+sequence_file.name, 'wb+') as destination:
            for chunk in sequence_file.chunks():
                destination.write(chunk)
    except Exception as e:
        raise ValueError('A very specific bad thing happened during file upload: {}'.format(e))

def upload_query_sequences_file(sequence_file, project):
    try:
        with open('media/'+str(project.id)+"/query_sequences/"+sequence_file.name, 'wb+') as destination:
            for chunk in sequence_file.chunks():
                destination.write(chunk)
    except Exception as e:
        raise ValueError('A very specific bad thing happened during file upload: {}'.format(e))


def save_genomes_and_query_in_db(query_sequences, forward_genome, backward_genome, project):
    try:
        uploaded_file_url_forward = 'media/'+str(project.id)+forward_genome.name
        uploaded_file_url_backward = 'media/'+str(project.id)+backward_genome.name
        uploaded_file_url_queries = 'media/'+str(project.id)+query_sequences.name
        new_forward_genome = Genomes(associated_project=project, genome_name=forward_genome.name,
                                 reciprocal_type='forward', path_to_file=uploaded_file_url_forward)
        new_backward_genome = Genomes(associated_project=project, genome_name=backward_genome.name,
                                  reciprocal_type='backward', path_to_file=uploaded_file_url_backward)
        new_query_sequences = QuerySequences(associated_project=project,query_file_name=query_sequences.name,
                                             path_to_query_file=uploaded_file_url_queries)
        new_forward_genome.save()
        new_backward_genome.save()
        new_query_sequences.save()
    except Exception as e:
        raise ValueError('A very specific bad thing happened during file upload: {}'.format(e))


def save_forward_settings_from_form_or_raise_exception(project,settings_form_forward):
    try:
        settings_fw = ForwardBlastSettings(associated_project=project, e_value=settings_form_forward['fw_e_value'],
                                           word_size=settings_form_forward['fw_word_size'], num_alignments=settings_form_forward['fw_num_alignments'])
        settings_fw.save()
    except Exception as e:
        raise ValueError('A very specific bad thing happened during settings saving: {}'.format(e))

def save_backward_settings_from_form_or_raise_exception(project,settings_form_backward):
    try:
        settings_bw = BackwardBlastSettings(associated_project=project, e_value=settings_form_backward['bw_e_value'],
                                            word_size=settings_form_backward['bw_word_size'], num_alignments=settings_form_backward['bw_num_alignments'])
        settings_bw.save()
    except Exception as e:
        raise ValueError('A very specific bad thing happened during settings saving: {}'.format(e))