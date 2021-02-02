from django.shortcuts import render, get_object_or_404, redirect
from django.contrib import messages
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group
from django.db import transaction, IntegrityError

from django.http import HttpResponse
from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences, TaxNodesForForwardDatabase, TaxNodesForBackwardDatabase
from .services import delete_files_without_projects, \
    save_genomes_and_query_in_db, save_project_from_form_or_raise_exception, \
    delete_project_files_by_project_id, upload_file, create_project_dir, \
    save_forward_settings_from_form_or_raise_exception, save_backward_settings_from_form_or_raise_exception, \
    set_executed_on_true_and_save_project, save_nr_project_from_form_or_raise_exception, save_query_file_in_db, \
    create_nr_project_dir, validate_fw_taxids_and_save_into_database, snakemake_project_finished,\
    get_html_results, load_html_graph, validate_bw_taxids_and_save_into_database, save_new_genome_file_in_db, \
    check_if_genomes_should_be_resaved

from .biopython_functions import get_species_taxid, get_scientific_name_by_taxid

from .blast_execution import write_snakefile, view_builded_snakefile, snakefile_exists,\
    exec_snakemake, write_nr_snakefile, read_snakemake_logs

from .forms import BlastProjectForm, BlastProjectNrForm, AdvancedSettingsForm_Forward, AdvancedSettingsForm_Backward, CreateUserForm, SpeciesNameForm, BlastProjectUploadedForm, UploadDatabaseForm

from .decorators import unauthenticated_user, allowed_user
# TODO refactoring and annotating

@unauthenticated_user
def loginPage(request):

    if request.method == 'POST':
        username = request.POST.get('username')
        password = request.POST.get('password')
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('main')
        else:
            messages.info(request,'Username OR password is incorrect')
            return render(request, 'blast_project/login.html')
    return render(request,'blast_project/login.html')

def logoutUser(request):
    logout(request)
    return redirect('login')

@unauthenticated_user
def registrationPage(request):
    userForm = CreateUserForm()
    if request.method == 'POST':
        userForm = CreateUserForm(request.POST)
        if userForm.is_valid():
            try:
                userForm.save()
                username = userForm.cleaned_data.get('username')
                #group = Group.objects.get(name='customer')
                #user.groups.add(group)
                messages.success(request,'Account was created for '+ username)
                return redirect('login')
            except Exception as e:
                return failure_view(request,e)
    context = {'form': userForm, }
    return render(request,'blast_project/register.html',context)

#homepage: lists all projects associated to the user
#main.html : lists projects and links to their detail and pipeline dashboard pages
#in the navigation bar are all other necessary links to project creation and supplementary pages
#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def main(request):
    try:
        projects = BlastProject.objects.filter(project_username=request.user)
        #query set of project specific databases
        target_genomes = Genomes.objects.filter(associated_project__in=projects)
        #deletes folders that are not included as projects ids in the database
        delete_files_without_projects()
    except Exception as e:
        return failure_view(request,e)
    context = {
        'projects':projects,
        'genomes':target_genomes,
        }
    return render(request,'blast_project/main.html',context)

#this page is rendered when post forms are valid
def success_view(request):
    return render(request,'blast_project/success.html')

#if an exception occurres this page is rendered in order to evaluate the exceptions context
def failure_view(request,exception):
    context={'exception':exception}
    return render(request,'blast_project/failure.html', context)

#view for deleting specific projects and all associated data
#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def delete_project(request,project_id):
    project = get_object_or_404(BlastProject, pk=project_id)
    if request.method == 'POST':
        #deleting associated files
        try:
            #execute this function before deletion of project
            check_if_genomes_should_be_resaved(project_id)
            delete_project_files_by_project_id(project_id)
            
            project.delete()
        except Exception as e:
            return failure_view(request,e)
        return success_view(request)
    else:
        target_genomes = Genomes.objects.filter(associated_project=project_id)
        forward_settings = ForwardBlastSettings.objects.get(associated_project=project_id)
        backward_settings = BackwardBlastSettings.objects.get(associated_project=project_id)
        query_sequences = QuerySequences.objects.get(associated_project=project_id)
        context = {'project': project, 'genomes': target_genomes, 'forward_settings': forward_settings,
                   'backward_settings': backward_settings,'query_sequences':query_sequences }
        return render(request, 'blast_project/delete_project.html', context)

#project detail view; if snakemake finished this view also displays result graphs
@login_required(login_url='login')
def project_details(request, project_id):
    #associated data for both project types
    project = get_object_or_404(BlastProject,pk=project_id)
    forward_settings = ForwardBlastSettings.objects.get(associated_project=project_id)
    backward_settings = BackwardBlastSettings.objects.get(associated_project=project_id)
    query_sequences = QuerySequences.objects.get(associated_project=project_id)

    #project details for projects that use uploaded files as databases
    if project.using_nr_database == False:
        target_genomes = Genomes.objects.filter(associated_project=project_id)

        context = {'project': project, 'genomes': target_genomes, 'forward_settings': forward_settings,
                   'backward_settings': backward_settings, 'query_sequences': query_sequences}
        if snakemake_project_finished(project.id):
            context['html_results_upload'] = True
            try:
                context['statistics_graph'] = load_html_graph(project_id, 'statistics.html')
            except Exception as e:
                context['error_html_graph_reading'] = "Couldn't create result plots with exception: " + e

    #project details for projects that use the nr database as database
    else:
        forward_db_organisms = TaxNodesForForwardDatabase.objects.filter(associated_project=project_id)
        backward_db_organisms = TaxNodesForBackwardDatabase.objects.filter(associated_project=project_id)

        context = {'project': project, 'forward_settings': forward_settings,'backward_db_organisms':backward_db_organisms,
                   'backward_settings': backward_settings, 'query_sequences': query_sequences, 'forward_db_organisms':forward_db_organisms}

        #view results if the snakemake execution has finished; if the last output file was created
        if snakemake_project_finished(project.id):
            context['html_results'] = True
            try:
                context['organisms_orthologs_graph'] = load_html_graph(project_id,'organisms_orthologs_graph.html')
                context['amount_hits'] = load_html_graph(project_id,'amount_hits.html')
            except Exception as e:
                context['error_html_graph_reading'] = "Couldn't create result plots with exception: "+e

    return render(request, 'blast_project/project_details.html', context)

#this view is used when a button is triggered in project_details it displays the reciprocal results
def display_reciprocal_result_table(request,project_id):
    try:
        html_data = get_html_results(project_id,"reciprocal_results.html")
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request,e)

#view for executing and displaying the status of the genome upload blast project pipeline
#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def pipeline_dashboard(request,project_id):
    project = get_object_or_404(BlastProject,pk=project_id)
    context = {'project': project}

    #code for displaying the actual status of the snakemake execution
    if snakefile_exists(project_id):
        try:
            content = view_builded_snakefile(project_id, 'upload')
            log_content = read_snakemake_logs(project_id)
            if 'no_logs' not in log_content[0].keys():
                context['log'] = log_content[0]
                context['pipeline_percentage'] = log_content[1]
            else:
                context['no_logs'] = log_content[0]
            context['content'] = content

        except Exception as e:
            return failure_view(request, e)
    return render(request,'blast_project/pipeline_dashboard.html', context)

#view for executing and displaying the status of the nr blast project pipeline
@login_required(login_url='login')
def pipeline_nr_dashboard(request,project_id):
    project = get_object_or_404(BlastProject, pk=project_id)
    context = {'project': project}

    #code for displaying the actual status of the snakemake execution
    if snakefile_exists(project_id):
        try:
            content = view_builded_snakefile(project_id,'nr')
            log_content = read_snakemake_logs(project_id)
            if 'no_logs' not in log_content[0].keys():
                context['log'] = log_content[0]
                context['pipeline_percentage'] = log_content[1]
            else:
                context['no_logs'] = log_content[0]
            context['content'] = content

        except Exception as e:
            return failure_view(request,e)

    return render(request,'blast_project/pipeline_nr_dashboard.html', context)

def execute_nr_snakefile(request):
    project_id = request.GET['execute_nr_snakefile']

    current_project=get_object_or_404(BlastProject,pk=project_id)
    try:
        set_executed_on_true_and_save_project(current_project)
    except IntegrityError as e:
        return failure_view(request,e)

    exec_snakemake(project_id)
    return redirect('pipeline_nr_dashboard', project_id)

def execute_snakefile(request):
    project_id = request.GET['execute_snakefile']

    current_project=get_object_or_404(BlastProject,pk=project_id)
    try:
        set_executed_on_true_and_save_project(current_project)
    except Exception as e:
        return failure_view(request,e)

    exec_snakemake(project_id)
    return redirect('pipeline_dashboard', project_id)

#view for project creation based on uploaded genome files
#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def create_project(request):
    if request.method == 'POST':
        project_creation_form = BlastProjectForm(request.POST,request.FILES)
        settings_form_forward = AdvancedSettingsForm_Forward(request.POST)
        settings_form_backward = AdvancedSettingsForm_Backward(request.POST)
        if project_creation_form.is_valid() and settings_form_forward.is_valid() and settings_form_backward.is_valid():
            try:
                # ensures that everything is correctly saved into the database, if an error occurres saving would not be transmitted
                with transaction.atomic():
                    new_title = project_creation_form.cleaned_data['project_title']
                    new_strategy = project_creation_form.cleaned_data['search_strategy']
                    project = save_project_from_form_or_raise_exception(new_title, new_strategy, request.user)

                    forward_genome = request.FILES['forward_genome_file']
                    backward_genome = request.FILES['backward_genome_file']
                    query_sequences = request.FILES['query_sequence_file']

                    save_genomes_and_query_in_db(query_sequences, forward_genome.name, backward_genome.name, project)
                    save_forward_settings_from_form_or_raise_exception(project, settings_form_forward.cleaned_data)
                    save_backward_settings_from_form_or_raise_exception(project, settings_form_backward.cleaned_data)

                    create_project_dir(project)
                    upload_file(forward_genome,'media/'+'databases/'+forward_genome.name)
                    upload_file(backward_genome,'media/'+'databases/'+backward_genome.name)
                    upload_file(query_sequences,'media/'+str(project.id)+'/'+'query_sequences'+'/'+query_sequences.name)

                    write_snakefile(project.id)
                    return success_view(request)
            except IntegrityError as e:
                return failure_view(request,e)

    else:
        project_creation_form = BlastProjectForm()
        settings_form_forward = AdvancedSettingsForm_Forward()
        settings_form_backward = AdvancedSettingsForm_Backward()
    context = {'BlastProjectForm': project_creation_form,
               'usage info': 'Forward Genome: Genome of query sequences\nBackward Genome: Genome of sequences for comparison to query sequences',
               'AdvancedSettingsForm_Forward':settings_form_forward,'AdvancedSettingsForm_Backward':settings_form_backward}


    return render(request, 'blast_project/project_creation.html', context)

@login_required(login_url='login')
def create_uploaded_based_project(request):
    if request.method == 'POST':
        project_creation_form = BlastProjectUploadedForm(request.POST,request.FILES)
        settings_form_forward = AdvancedSettingsForm_Forward(request.POST)
        settings_form_backward = AdvancedSettingsForm_Backward(request.POST)
        if project_creation_form.is_valid() and settings_form_forward.is_valid() and settings_form_backward.is_valid():
            try:
                # ensures that everything is correctly saved into the database, if an error occurres saving would not be transmitted
                with transaction.atomic():
                    new_title = project_creation_form.cleaned_data['project_title']
                    new_strategy = project_creation_form.cleaned_data['search_strategy']
                    project = save_project_from_form_or_raise_exception(new_title, new_strategy, request.user)

                    forward_genome = project_creation_form.cleaned_data['forward_genome_file']
                    backward_genome = project_creation_form.cleaned_data['backward_genome_file']

                    #just the genome/database name is required in order to reuse this object
                    forward_genome_data = Genomes.objects.filter(genome_name=forward_genome).order_by('id').first()
                    backward_genome_data = Genomes.objects.filter(genome_name=backward_genome).order_by('id').first()

                    query_sequences = request.FILES['query_sequence_file']

                    create_project_dir(project)
                    upload_file(query_sequences,'media/'+str(project.id)+'/'+'query_sequences'+'/'+query_sequences.name)
                    save_genomes_and_query_in_db(query_sequences, forward_genome_data.genome_name, backward_genome_data.genome_name, project)

                    save_forward_settings_from_form_or_raise_exception(project,settings_form_forward.cleaned_data)
                    save_backward_settings_from_form_or_raise_exception(project,settings_form_backward.cleaned_data)
                    write_snakefile(project.id)
                    return success_view(request)
            except IntegrityError as e:
                return failure_view(request,e)

    else:
        project_creation_form = BlastProjectUploadedForm()
        settings_form_forward = AdvancedSettingsForm_Forward()
        settings_form_backward = AdvancedSettingsForm_Backward()
    context = {'BlastProjectUploadedForm': project_creation_form,
               'usage info': 'Forward Genome: Genome of query sequences\nBackward Genome: Genome of sequences for comparison to query sequences',
               'AdvancedSettingsForm_Forward':settings_form_forward,'AdvancedSettingsForm_Backward':settings_form_backward}


    return render(request, 'blast_project/project_creation_based_on_uploaded_files.html', context)

@login_required(login_url='login')
def create_nr_based_project(request):
    if request.method == 'POST':

        project_creation_form = BlastProjectNrForm(request.POST, request.FILES)
        settings_form_forward = AdvancedSettingsForm_Forward(request.POST)
        settings_form_backward = AdvancedSettingsForm_Backward(request.POST)
        if project_creation_form.is_valid() and settings_form_forward.is_valid() and settings_form_backward.is_valid():
            try:
                #ensures that everything is correctly saved into the database, if an error occurres saving would not be transmitted
                with transaction.atomic():
                    new_title = project_creation_form.cleaned_data['project_title']
                    taxonomic_nodes_bw = project_creation_form.cleaned_data['taxid_bw']
                    taxonomic_nodes_fw = project_creation_form.cleaned_data['taxid_fw']
                    query_sequences = request.FILES['query_sequence_file']

                    project = save_nr_project_from_form_or_raise_exception(new_title, request.user)

                    validate_fw_taxids_and_save_into_database(project, request.user.email, taxonomic_nodes_fw)
                    validate_bw_taxids_and_save_into_database(project, request.user.email, taxonomic_nodes_bw)

                    save_query_file_in_db(query_sequences, project)


                    save_forward_settings_from_form_or_raise_exception(project, settings_form_forward.cleaned_data)
                    save_backward_settings_from_form_or_raise_exception(project,
                                                                        settings_form_backward.cleaned_data)

                    create_nr_project_dir(project)
                    upload_file(query_sequences,'media/'+str(project.id)+'/'+'query_sequences'+'/'+query_sequences.name)

                    write_nr_snakefile(project.id)
                    return success_view(request)

            except IntegrityError as e:
                return failure_view(request,e)

    else:
        project_creation_form = BlastProjectNrForm()
        settings_form_forward = AdvancedSettingsForm_Forward()
        settings_form_backward = AdvancedSettingsForm_Backward()

    context = {'BlastProjectNrForm': project_creation_form,
               'AdvancedSettingsForm_Forward': settings_form_forward,
               'AdvancedSettingsForm_Backward': settings_form_backward}

    return render(request, 'blast_project/project_creation_based_on_nr_database.html', context)

@login_required(login_url='login')
def species_taxid(request):
    user_mail = request.user.email
    context={}
    if request.method == 'POST':
        taxid_form = SpeciesNameForm(request.POST)
        if taxid_form.is_valid():
            organism_name = taxid_form.cleaned_data.get('organism_name')
            try:
                taxid_and_record_info = get_species_taxid(user_mail,organism_name)
                context['taxid']=taxid_and_record_info[0]
                context['translation']=taxid_and_record_info[1]
            except Exception as e:
                return failure_view(request,e)
    else:
        taxid_form = SpeciesNameForm()
    context['SpeciesNameForm'] = taxid_form
    return render(request,'blast_project/species_taxid.html',context)

def upload_databases(request):
    context={}
    if request.method == 'POST':
        upload_database_form = UploadDatabaseForm(request.POST, request.FILES)
        if upload_database_form.is_valid():
            try:
                new_database_file = upload_database_form.cleaned_data['genome_file']
                upload_file(new_database_file,'media/'+'databases/'+new_database_file.name)
                save_new_genome_file_in_db(new_database_file)
                return success_view(request)
            except IntegrityError as e:
                return failure_view(request,e)
    else:
        upload_database_form = UploadDatabaseForm()
    context['UploadDatabaseForm'] = upload_database_form
    return render(request,'blast_project/upload_databases.html',context)