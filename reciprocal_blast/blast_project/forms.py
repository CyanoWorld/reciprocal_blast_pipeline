from django import forms
from .models import BlastProject, Genomes
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from .biopython_functions import get_species_taxid, get_species_taxid_without_email

def get_genomes_tuple():
    try:
        genomes = []
        for genome_name in Genomes.objects.all().values('genome_name'):
            if genome_name['genome_name'] not in genomes:
                genomes.append(genome_name['genome_name'])
        genomes = tuple(zip(genomes,genomes))
        return genomes
    except:
        return (('no genomes available','no genomes available'))

class BlastProjectForm(forms.Form):
    project_title = forms.CharField(label="Project title",
                                    error_messages={'required':"A project title is required for saving project metadata into the database"})
    search_strategy = forms.ChoiceField(choices=(('blastp','blastp'),('blastn','blastn')),label="Search strategy",
                                        error_messages={'required':"Please specify a search strategy otherwise you won't be able to execute a BLAST run .."})
    forward_genome_file = forms.FileField(error_messages={'required':"Upload a forward genome file fasta, this is the database for the first BLAST (not your query sequence genome)"},required=True)
    backward_genome_file = forms.FileField(error_messages={'required':"Upload a backward genome file fasta, this is the database for the second BLAST (from your query sequences)"},required=True)

    forward_genome_uploaded_file = forms.TypedChoiceField(choices=get_genomes_tuple(),required=False, error_messages={'required':'Please specify an uploaded genome database'})
    backward_genome_uploaded_file = forms.TypedChoiceField(choices=get_genomes_tuple(),required=False, error_messages={'required':'Please specify an uploaded genome database'})
    #forward_genome_uploaded_file = forms.ModelChoiceField(queryset=Genomes.objects.all().values('genome_name'),required=False, error_messages={'required':'Please specify an uploaded genome database'})
    #backward_genome_uploaded_file = forms.ModelChoiceField(queryset=Genomes.objects.all().values('genome_name'),required=False, error_messages={'required':'Please specify an uploaded genome database'})

    query_sequence_file = forms.FileField(error_messages={'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})


    def __init__(self,data=None,*args,**kwargs):
        super(BlastProjectForm,self).__init__(data,*args,**kwargs)
        if data:
            print(self.fields['backward_genome_uploaded_file'])
            if data.get('backward_genome_file', None) == '':
                self.fields['backward_genome_uploaded_file'].required = True
                self.fields['backward_genome_file'].required = False

            if data.get('forward_genome_file', None) == '':
                self.fields['forward_genome_uploaded_file'].required = True
                self.fields['forward_genome_file'].required = False

    def clean_query_sequence_file(self):
        query_file = self.cleaned_data['query_sequence_file']
        if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
            raise ValidationError("[-] Please upload only fasta files!")
        else:
            return query_file

    def clean_forward_genome_file(self):
        genome = self.cleaned_data['forward_genome_file']
        if genome != None:
            if genome.name.endswith('.faa') != True and genome.name.endswith('.fasta') != True:
                raise ValidationError("[-] Please upload only fasta files!")

            for genome_name in Genomes.objects.all().values('genome_name'):
                if genome_name['genome_name'] == genome.name:
                    raise ValidationError(
                        '[-] A genome database with that name already exist, please specify another file! Or try to use previously uploaded databases!')
        return genome

    def clean_backward_genome_file(self):
        genome = self.cleaned_data['backward_genome_file']
        if genome != None:
            if genome.name.endswith('.faa') != True and genome.name.endswith('.fasta') != True:
                raise ValidationError("[-] Please upload only fasta files!")

            for genome_name in Genomes.objects.all().values('genome_name'):
                if genome_name['genome_name'] == genome.name:
                    raise ValidationError(
                        '[-] A genome database with that name already exist, please specify another file! Or try to use previously uploaded databases!')
        return genome

'''
class BlastProjectUploadedForm(forms.Form):
    project_title = forms.CharField(label="Project title",
                                    error_messages={
                                        'required': "A project title is required for saving project metadata into the database"})
    search_strategy = forms.ChoiceField(choices=(('blastp', 'blastp'), ('blastn', 'blastn')), label="Search strategy",
                                        error_messages={
                                            'required': "Please specify a search strategy otherwise you won't be able to execute a BLAST run .."})

    #error_messages={'required': "Use a database file, this is the database for the second BLAST (from your query sequences)"}
    forward_genome_file = forms.ChoiceField(choices=get_genomes_tuple())
    backward_genome_file = forms.ChoiceField(choices=get_genomes_tuple())
    query_sequence_file = forms.FileField(error_messages={'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

    #this is just possible because choice fields are previously 'cleaned'
    #it is not a good practice to do that
    def clean_backward_genome_file(self):
        bw_genome = self.cleaned_data['backward_genome_file']
        fw_genome = self.cleaned_data['forward_genome_file']
        if bw_genome == fw_genome:
            raise ValidationError('[-] Do not use the same databases for the forward and backward blast.')
        else:
            return bw_genome
'''

#TODO: use the user_email in order to check the scientific names
class BlastProjectNrForm(forms.Form):
    project_title = forms.CharField(label="Project title",
                                    error_messages={
                                        'required': "A project title is required for saving project metadata into the database"})
    taxid_bw = forms.CharField(required=True ,label='Scientific Names (conversion to Taxonomic Nodes) for Backward BLAST', error_messages={'required':"Specify a Scientific Name for your backward BLAST - use a comma separated list - names will be converted to taxids that will be written to a file which will serve as the -taxidlist parameter of your backward BLAST"})
    query_sequence_file = forms.FileField(error_messages={'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})
    taxid_fw = forms.CharField(required=False, label = 'Scientific Names (conversion to Taxonomic Nodes) for Forward BLAST',
                               error_messages={'required':"Specify a Scientific Name for your forward BLAST - use a comma seperated list"})
    user_email = forms.CharField(required=True)

    def clean_taxid_fw(self):
        taxids_fw = self.cleaned_data['taxid_fw']
        taxids = []
        if taxids_fw != '':
            for name in taxids_fw.split(','):
                try:
                    taxid = get_species_taxid_without_email(name)
                    taxids.append(int(taxid))
                except Exception as e:
                    raise ValidationError('[-] Use a comma separated list of scientific names. An error occured during parsing of the scientific names. Exception: {}'.format(e))
        return taxids

    def clean_taxid_bw(self):
        taxids_bw = self.cleaned_data['taxid_bw']
        if len(taxids_bw.split(',')) > 1:
            raise ValidationError(
                '[-] Specify at least one scientific name!')
        try:
            taxid = get_species_taxid_without_email(taxids_bw)
        except Exception as e:
            raise ValidationError('[-] An error occured during parsing of the scientific names. Exception: {}'.format(e))
        return [taxid]

    def clean_query_sequence_file(self):
        query_file = self.cleaned_data['query_sequence_file']
        if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
            raise ValidationError("[-] Please upload only fasta files!")
        else:
            return query_file

class AdvancedSettingsForm_Forward(forms.Form):
    fw_e_value = forms.DecimalField(label="FW E-Value", required=False, initial=0.001)
    fw_word_size = forms.IntegerField(label="FW Word Size", required=False, initial=3)
    fw_num_alignments = forms.IntegerField(label="FW Number of possible alignment outputs", required=False, initial=10000)
    fw_num_descriptions = forms.IntegerField(label="FW Number of possible alignment description outputs", required=False, initial=500)
    fw_num_threads = forms.IntegerField(label="Number of threads used for executing this BLAST search", required=False, initial=1)
    #fw_max_hsps = forms.IntegerField(label="FW Number of high scoring pairs to output", required=False, initial=0)
    #fw_max_target_seqs = forms.IntegerField(label="FW Number of maximal alignment outputs", required=False, initial=500)

class AdvancedSettingsForm_Backward(forms.Form):
    bw_e_value = forms.DecimalField(label="BW E-Value", required=False, initial=0.001)
    bw_word_size = forms.IntegerField(label="BW Word Size", required=False,initial=3)
    bw_num_alignments = forms.IntegerField(label="BW Number of possible alignment outputs",required=False, initial=1)
    bw_num_descriptions = forms.IntegerField(label="FW Number of possible alignment description outputs",
                                              required=False, initial=1)
    bw_num_threads = forms.IntegerField(label="Number of threads used for executing this BLAST search", required=False,
                                         initial=1)
    bw_max_hsps = forms.IntegerField(label="FW Number of high scoring pairs to output", required=False, initial=1)
    #bw_max_target_seqs = forms.IntegerField(label="FW Number of maximal alignment outputs", required=False,initial=1)

class CreateUserForm(UserCreationForm):
    class Meta:
        model = User
        fields = ['username','email','password1','password2']

class SpeciesNameForm(forms.Form):
    organism_name = forms.CharField(label="Organism Name", initial='Cyanobacteria')

class UploadDatabaseForm(forms.Form):
    genome_file = forms.FileField(error_messages={'required':"Upload a genome fasta file for database preparement"})
    def clean_genome_file(self):
        genome = self.cleaned_data['genome_file']
        for genome_name in Genomes.objects.all().values('genome_name'):
            if genome_name['genome_name'] == genome.name:
                raise ValidationError('[-] A genome database with that name already exist, please specify another file! Or try to use previously uploaded databases!')
        return genome