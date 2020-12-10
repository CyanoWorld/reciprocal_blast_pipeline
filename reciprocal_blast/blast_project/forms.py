from django import forms
from .models import BlastProject, Genomes
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

class BlastProjectForm(forms.Form):
    project_title = forms.CharField(label="Project title",
                                    error_messages={'required':"A project title is required for saving project metadata into the database"})
    search_strategy = forms.ChoiceField(choices=(('blastp','blastp'),('blastn','blastn')),label="Search strategy",
                                        error_messages={'required':"Please specify a search strategy otherwise you won't be able to execute a BLAST run .."})
    forward_genome_file = forms.FileField(error_messages={'required':"Upload a forward genome file fasta file, this is the genome for the first BLAST (not your query sequence genome)"})
    backward_genome_file = forms.FileField(error_messages={'required':"Upload a backward genome file fasta file, this is the genome for the second BLAST (from your query sequences)"})
    query_sequence_file = forms.FileField(error_messages={'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

class AdvancedSettingsForm_Forward(forms.Form):
    fw_e_value = forms.DecimalField(label="FW E-Value", required=False, initial=0.001)
    fw_word_size = forms.IntegerField(label="FW Word Size", required=False, initial=3)
    fw_num_alignments = forms.IntegerField(label="FW Number of possible alignment outputs", required=False, initial=10000)

class AdvancedSettingsForm_Backward(forms.Form):
    bw_e_value = forms.DecimalField(label="BW E-Value", required=False, initial=0.001)
    bw_word_size = forms.IntegerField(label="BW Word Size", required=False,initial=3)
    bw_num_alignments = forms.IntegerField(label="BW Number of possible alignment outputs",required=False, initial=1)

class CreateUserForm(UserCreationForm):
    class Meta:
        model = User
        fields = ['username','email','password1','password2']
