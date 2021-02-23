'''
Content: anchor point for project creation. Functions are used by the views.py function project_creation.
External functions are used for creating the project directories and snakemake files, this is done in services.py and blast_execution.py.
'''

from django.db import IntegrityError, transaction
from .models import Genomes
from .blast_execution import prepare_data_and_write_genome_upload_snakemake_configurations, prepare_data_and_write_nr_snakemake_configurations_and_taxid_files
from .services import save_genomes_and_query_in_db, save_project_from_form_or_raise_exception, \
    upload_file, create_project_dir, \
    save_forward_settings_from_form_or_raise_exception, save_backward_settings_from_form_or_raise_exception, \
    save_nr_project_from_form_or_raise_exception, save_query_file_in_db, \
    create_nr_project_dir, validate_fw_taxids_and_save_into_database, validate_bw_taxids_and_save_into_database

def create_project_with_nr_database(request, project_creation_form,settings_form_forward,settings_form_backward):
    try:
        # ensures that everything is correctly saved into the database, if an error occurres saving would not be transmitted
        with transaction.atomic():
            new_title = project_creation_form.cleaned_data['project_title']
            taxonomic_nodes_bw = project_creation_form.cleaned_data['taxid_bw']
            taxonomic_nodes_fw = project_creation_form.cleaned_data['taxid_fw']
            query_sequences = request.FILES['query_sequence_file']

            project = save_nr_project_from_form_or_raise_exception(new_title, request.user)

            #validation is previously handled in the forms.py class
            validate_fw_taxids_and_save_into_database(project, request.user.email, taxonomic_nodes_fw)
            validate_bw_taxids_and_save_into_database(project, request.user.email, taxonomic_nodes_bw)

            save_query_file_in_db(query_sequences, project)

            save_forward_settings_from_form_or_raise_exception(project, settings_form_forward.cleaned_data)
            save_backward_settings_from_form_or_raise_exception(project,
                                                                settings_form_backward.cleaned_data)

            create_nr_project_dir(project)
            upload_file(query_sequences,
                        'media/' + str(project.id) + '/' + 'query_sequences' + '/' + query_sequences.name)

            #this function of the blast_execution.py file triggers a function that translates high order taxids into species ids for limiting BLAST searches
            prepare_data_and_write_nr_snakemake_configurations_and_taxid_files(project.id)
    except Exception as e:
        raise IntegrityError("[-] Couldn't perform project creation due to following exception: {}".format(e))

def create_project_with_uploaded_files(request, project_creation_form,settings_form_forward,settings_form_backward):
    try:
        # ensures that everything is correctly saved into the database, if an error occurres saving would not be transmitted
        with transaction.atomic():
            new_title = project_creation_form.cleaned_data['project_title']
            new_strategy = project_creation_form.cleaned_data['search_strategy']
            project = save_project_from_form_or_raise_exception(new_title, new_strategy, request.user)

            query_sequences = request.FILES['query_sequence_file']

            save_forward_settings_from_form_or_raise_exception(project, settings_form_forward.cleaned_data)
            save_backward_settings_from_form_or_raise_exception(project, settings_form_backward.cleaned_data)

            #these conditions check wether a previously genome database should be used or if a file gets uploaded
            #four conditions:
            # available = previously uploaded
            # fw and bw not available
            # fw available and bw not
            # bw available and fw not
            # bw and fw available
            if  project_creation_form.cleaned_data['forward_genome_file'] == None and project_creation_form.cleaned_data['backward_genome_file'] == None:
                forward_genome = project_creation_form.cleaned_data['forward_genome_uploaded_file']
                forward_genome_data = Genomes.objects.filter(genome_name=forward_genome).order_by('id').first()
                backward_genome = project_creation_form.cleaned_data['backward_genome_uploaded_file']
                backward_genome_data = Genomes.objects.filter(genome_name=backward_genome).order_by('id').first()
                save_genomes_and_query_in_db(query_sequences, forward_genome_data.genome_name, backward_genome_data.genome_name,
                                             project)

            elif project_creation_form.cleaned_data['forward_genome_file'] == None and project_creation_form.cleaned_data['backward_genome_file'] != None:
                forward_genome = project_creation_form.cleaned_data['forward_genome_uploaded_file']
                forward_genome_data = Genomes.objects.filter(genome_name=forward_genome).order_by('id').first()
                backward_genome_data = request.FILES['backward_genome_file']
                save_genomes_and_query_in_db(query_sequences, forward_genome_data.genome_name, backward_genome_data.name,
                                             project)
                upload_file(backward_genome_data, 'media/' + 'databases/' + backward_genome_data.name)


            elif project_creation_form.cleaned_data['forward_genome_file'] != None and project_creation_form.cleaned_data['backward_genome_file'] == None:
                backward_genome = project_creation_form.cleaned_data['backward_genome_uploaded_file']
                backward_genome_data = Genomes.objects.filter(genome_name=backward_genome).order_by('id').first()
                forward_genome_data = request.FILES['forward_genome_file']
                save_genomes_and_query_in_db(query_sequences, forward_genome_data.name, backward_genome_data.genome_name,
                                             project)
                upload_file(forward_genome_data, 'media/' + 'databases/' + forward_genome_data.name)


            elif project_creation_form.cleaned_data['forward_genome_file'] != None and project_creation_form.cleaned_data['backward_genome_file'] != None:
                forward_genome_data = request.FILES['forward_genome_file']
                backward_genome_data = request.FILES['backward_genome_file']
                save_genomes_and_query_in_db(query_sequences, forward_genome_data.name, backward_genome_data.name,
                                             project)
                upload_file(forward_genome_data, 'media/' + 'databases/' + forward_genome_data.name)
                upload_file(backward_genome_data, 'media/' + 'databases/' + backward_genome_data.name)
            else:
                raise IntegrityError("[-] Something with the database specifications is not correct")



            create_project_dir(project)
            upload_file(query_sequences,
                        'media/' + str(project.id) + '/' + 'query_sequences' + '/' + query_sequences.name)

            prepare_data_and_write_genome_upload_snakemake_configurations(project.id)
    except Exception as e:
        raise IntegrityError("[-] Couldn't perform project creation due to exception: {}".format(e))
