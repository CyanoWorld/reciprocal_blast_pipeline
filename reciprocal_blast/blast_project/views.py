from django.shortcuts import render, get_object_or_404, redirect
from django.contrib import messages
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group

from django.http import Http404
from .models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences
from .services import delete_files_without_projects, \
    save_genomes_and_query_in_db, save_project_from_form_or_raise_exception, \
    delete_genomfiles_by_project_id, upload_forward_genome_file, upload_backward_genome_file, upload_query_sequences_file, create_project_dir, \
    save_forward_settings_from_form_or_raise_exception, save_backward_settings_from_form_or_raise_exception
from .blast_execution import write_snakefile, view_builded_snakefile, snakefile_exists, delete_snakefile, exec_snakemake
from .forms import BlastProjectForm, AdvancedSettingsForm_Forward, AdvancedSettingsForm_Backward, CreateUserForm

from .decorators import unauthenticated_user, allowed_user
# TODO refactoring

# Create your views here.
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

#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def main(request):
    #print("[+] {} | {] | {}".format(BASE_DIR,MEDIA_URL,MEDIA_ROOT))
    try:
        projects = BlastProject.objects.filter(project_username=request.user)
        target_genomes = Genomes.objects.filter(associated_project__in=projects)
        delete_files_without_projects()
    except Exception as e:
        return failure_view(request,e)
    context = {
        'projects':projects,
        'genomes':target_genomes,
        }
    return render(request,'blast_project/main.html',context)

def success_view(request):
    return render(request,'blast_project/success.html')

def failure_view(request,exception):
    context={'exception':exception}
    return render(request,'blast_project/failure.html', context)

#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def delete_project(request,project_id):
    project = get_object_or_404(BlastProject, pk=project_id)
    if request.method == 'POST':
        #deleting associated files
        delete_genomfiles_by_project_id(project_id)
        project.delete()
        return success_view(request)
    else:
        target_genomes = Genomes.objects.filter(associated_project=project_id)
        forward_settings = ForwardBlastSettings.objects.get(associated_project=project_id)
        backward_settings = BackwardBlastSettings.objects.get(associated_project=project_id)
        query_sequences = QuerySequences.objects.get(associated_project=project_id)
        context = {'project': project, 'genomes': target_genomes, 'forward_settings': forward_settings,
                   'backward_settings': backward_settings,'query_sequences':query_sequences }
        return render(request, 'blast_project/delete_project.html', context)

@login_required(login_url='login')
def project_details(request, project_id):
    project = get_object_or_404(BlastProject,pk=project_id)
    target_genomes = Genomes.objects.filter(associated_project=project_id)
    forward_settings = ForwardBlastSettings.objects.get(associated_project=project_id)
    backward_settings = BackwardBlastSettings.objects.get(associated_project=project_id)
    query_sequences = QuerySequences.objects.get(associated_project=project_id)
    context = {'project': project, 'genomes': target_genomes, 'forward_settings': forward_settings,
               'backward_settings': backward_settings, 'query_sequences': query_sequences}
    return render(request, 'blast_project/project_details.html', context)

#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def pipeline_dashboard(request,project_id):
    project = get_object_or_404(BlastProject,pk=project_id)
    context = {'project': project}
    if snakefile_exists(project_id):
        try:
            content = view_builded_snakefile(project_id)
            context['content'] = content
        except Exception as e:
            return failure_view(request,e)
    return render(request,'blast_project/pipeline_dashboard.html', context)

#TODO GET request for deleting, this breaks the REST Framework definition...
#@allowed_user(['admin','customer'])
@login_required(login_url='login')
def remove_snakefile(request):
    if request.method == 'GET':
        project_id = request.GET['remove_snakefile']
        try:
            delete_snakefile(project_id)
        except Exception as e:
            return failure_view(request,e)

        return redirect('pipeline_dashboard',project_id)
    else:
        return Http404()

#@allowed_user(['admin','customer'])
#placeholder for a realy exiting function
@login_required(login_url='login')
def build_snakefile(request):
    if request.method == 'POST':
        project_id = request.POST['build_snakefile']
        #write project associated snakefile --> load project and all associated objects and get access to projects media folder
        try:
            write_snakefile(project_id)
        except Exception as e:
            return failure_view(request,e)
        return redirect('pipeline_dashboard',project_id)
    else:
        return Http404()

def execute_snakefile(request):
    project_id = request.GET['execute_snakefile']
    exec_snakemake(project_id)
    return redirect('pipeline_dashboard', project_id)

#@allowed_user(['admin','customer'])
#TODO add not valid features ...
@login_required(login_url='login')
def create_project(request):
    if request.method == 'POST':
        project_creation_form = BlastProjectForm(request.POST,request.FILES)
        settings_form_forward = AdvancedSettingsForm_Forward(request.POST)
        settings_form_backward = AdvancedSettingsForm_Backward(request.POST)
        if project_creation_form.is_valid():
            new_title = project_creation_form.cleaned_data['project_title']
            new_strategy = project_creation_form.cleaned_data['search_strategy']
            project = save_project_from_form_or_raise_exception(new_title, new_strategy, request.user)

            try:
                forward_genome = request.FILES['forward_genome_file']
                backward_genome = request.FILES['backward_genome_file']
                query_sequences = request.FILES['query_sequence_file']

                create_project_dir(project)
                upload_forward_genome_file(forward_genome, project)
                upload_backward_genome_file(backward_genome, project)
                upload_query_sequences_file(query_sequences, project)
                save_genomes_and_query_in_db(query_sequences, forward_genome, backward_genome, project)
            except Exception as e:
                return failure_view(request,e)

            if settings_form_forward.is_valid() and settings_form_backward.is_valid():
                save_forward_settings_from_form_or_raise_exception(project,settings_form_forward.cleaned_data)
                save_backward_settings_from_form_or_raise_exception(project,settings_form_backward.cleaned_data)
                return success_view(request)

    else:
        project_creation_form = BlastProjectForm()
        settings_form_forward = AdvancedSettingsForm_Forward()
        settings_form_backward = AdvancedSettingsForm_Backward()
    context = {'BlastProjectForm': project_creation_form,
               'usage info': 'Forward Genome: Genome of query sequences\nBackward Genome: Genome of sequences for comparison to query sequences',
               'AdvancedSettingsForm_Forward':settings_form_forward,'AdvancedSettingsForm_Backward':settings_form_backward}


    return render(request, 'blast_project/project_creation.html', context)
