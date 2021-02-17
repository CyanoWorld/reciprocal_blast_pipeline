from django.contrib import messages
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.db import IntegrityError
from django.http import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect

from .biopython_functions import get_species_taxid
from .blast_execution import view_builded_snakefile, snakemake_config_exists, \
    exec_snakemake, read_snakemake_logs
from .decorators import unauthenticated_user
from .forms import BlastProjectForm, BlastProjectNrForm, AdvancedSettingsForm_Forward, AdvancedSettingsForm_Backward, \
    CreateUserForm, SpeciesNameForm, UploadDatabaseForm
#BlastProjectUploadedForm
from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences, \
    TaxNodesForForwardDatabase, TaxNodesForBackwardDatabase
from .project_creation import create_project_with_uploaded_files, create_project_with_nr_database
from .services import delete_files_without_projects, \
    delete_project_files_by_project_id, upload_file, set_executed_on_true_and_save_project, snakemake_project_finished, \
    get_html_results, load_html_graph, save_new_genome_file_in_db, \
    check_if_genomes_should_be_resaved


# TODO refactoring and annotating
#General informations concerning views:
#request.method == 'GET' view



#This view loads the project_creation.html template which in turn includes the upload_genomes_form.html as well as
#the nr_database_form.html. It then
@login_required(login_url='login')
def project_creation(request):
    if request.method == 'POST':

        settings_form_forward = AdvancedSettingsForm_Forward(request.POST)
        settings_form_backward = AdvancedSettingsForm_Backward(request.POST)

        try:
            if request.POST['project_type'] == 'upload':
                project_creation_nr_from = BlastProjectNrForm()

                project_creation_form = BlastProjectForm(request.POST, request.FILES)

                if project_creation_form.is_valid() and settings_form_forward.is_valid() and settings_form_backward.is_valid():
                    create_project_with_uploaded_files(request, project_creation_form,settings_form_forward,settings_form_backward)
                    return success_view(request)

            if request.POST['project_type'] == 'nr':
                project_creation_form = BlastProjectForm()

                project_creation_nr_from = BlastProjectNrForm(request.POST, request.FILES)
                if project_creation_nr_from.is_valid() and settings_form_forward.is_valid() and settings_form_backward.is_valid():
                    create_project_with_nr_database(request, project_creation_nr_from,settings_form_forward,settings_form_backward)
                    return success_view(request)

        except Exception as e:
            return failure_view(request,e)


    else:
        project_creation_form = BlastProjectForm()
        project_creation_nr_from = BlastProjectNrForm()
        #project_creation_uploaded_form = BlastProjectUploadedForm()
        settings_form_forward = AdvancedSettingsForm_Forward()
        settings_form_backward = AdvancedSettingsForm_Backward()

    #user_email just for nr project creation
    context = {'BlastProjectForm': project_creation_form,
               'BlastProjectNrForm': project_creation_nr_from,
               #'BlastProjectUploadedForm': project_creation_uploaded_form,
           'AdvancedSettingsForm_Forward': settings_form_forward,
           'AdvancedSettingsForm_Backward': settings_form_backward,
               'user_email':request.user.email}

    return render(request,'blast_project/project_creation.html',context)

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
@login_required(login_url='login')
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
#the view is divided into a view for the nr associated projects and the uploaded genome database projects
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

        #view results if the snakemake execution has finished; if the last output file was created; currently 'reciprocal_results.html'
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

        #view results if the snakemake execution has finished; if the last output file was created; currently 'reciprocal_results.html'
        if snakemake_project_finished(project.id):
            context['html_results'] = True
            try:
                context['organisms_orthologs_graph'] = load_html_graph(project_id,'organisms_orthologs_graph.html')
                context['amount_hits'] = load_html_graph(project_id,'amount_hits.html')
            except Exception as e:
                context['error_html_graph_reading'] = "Couldn't create result plots with exception: "+e

    return render(request, 'blast_project/project_details.html', context)

#this view is used when a button is triggered in project_details it displays the reciprocal results
@login_required(login_url='login')
def display_reciprocal_result_table(request,project_id):
    try:
        html_data = get_html_results(project_id,"reciprocal_results.html")
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request,e)

#view for executing and displaying the status of the genome upload project pipeline
#the view reads the snakefile and the last snakemake log file
#the view loads the pipeline_dashboard.html template in which buttons for executing and monitoring are included
#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def pipeline_dashboard(request,project_id):
    project = get_object_or_404(BlastProject,pk=project_id)
    context = {'project': project}

    #code for displaying the actual status of the snakemake execution
    if snakemake_config_exists(project_id):
        try:
            content = view_builded_snakefile(project_id, 'upload')
            log_content = read_snakemake_logs(project_id)
            #logs are available
            if 'no_logs' not in log_content[0].keys():
                context['log'] = log_content[0]
                context['pipeline_percentage'] = log_content[1]
            #there are no logs available, yet
            else:
                context['no_logs'] = log_content[0]
            context['content'] = content

        except Exception as e:
            return failure_view(request, e)
    return render(request,'blast_project/pipeline_dashboard.html', context)

#this view is similar but slightly different to the pipeline_dashboard view, maybe a combination would be appropriate
#view for executing and displaying the status of the nr blast project pipeline
#the view reads the snakefile and the last snakemake log file
#the view loads the pipeline_nr_dashboard.html template in which buttons for executing and monitoring are included
@login_required(login_url='login')
def pipeline_nr_dashboard(request,project_id):
    project = get_object_or_404(BlastProject, pk=project_id)
    context = {'project': project}

    #code for displaying the actual status of the snakemake execution
    if snakemake_config_exists(project_id):
        try:
            content = view_builded_snakefile(project_id,'nr')
            log_content = read_snakemake_logs(project_id)
            #logs are available
            if 'no_logs' not in log_content[0].keys():
                context['log'] = log_content[0]
                context['pipeline_percentage'] = log_content[1]
            #there are no logs available, yet
            else:
                context['no_logs'] = log_content[0]
            context['content'] = content

        except Exception as e:
            return failure_view(request,e)

    return render(request,'blast_project/pipeline_nr_dashboard.html', context)

#this view is triggered if the user pressed the EXECUTE NR SNAKEMAKE button in the pipeline_nr_dashboard.html template
#it
@login_required(login_url='login')
def execute_nr_snakefile(request):
    project_id = request.GET['execute_nr_snakefile']

    current_project=get_object_or_404(BlastProject,pk=project_id)
    try:
        set_executed_on_true_and_save_project(current_project)
    except IntegrityError as e:
        return failure_view(request,e)

    exec_snakemake(project_id)
    return redirect('pipeline_nr_dashboard', project_id)

@login_required(login_url='login')
def execute_snakefile(request):
    project_id = request.GET['execute_snakefile']

    current_project=get_object_or_404(BlastProject,pk=project_id)
    try:
        set_executed_on_true_and_save_project(current_project)
    except Exception as e:
        return failure_view(request,e)

    exec_snakemake(project_id)
    return redirect('pipeline_dashboard', project_id)

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

@login_required(login_url='login')
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