from django.db import models
from django.contrib.auth.models import User
from .biopython_functions import get_scientific_name_by_taxid

class BlastProject(models.Model):
    project_title = models.CharField(max_length=200, blank=False, unique=True)
    search_strategy = models.CharField(max_length=20, choices=[('blastp', 'blastp'), ('blastn', 'blastn')], default='blastp')
    project_username = models.ForeignKey(User,on_delete=models.CASCADE)
    pipeline_executed = models.BooleanField(default=False)
    using_nr_database = models.BooleanField(default=False)
    def __str__(self):
        return "{}".format(self.project_title)

class TaxNodesForForwardDatabase(models.Model):
    associated_project = models.ForeignKey(BlastProject,on_delete=models.CASCADE)
    taxonomic_node = models.IntegerField()
    organism_name = models.CharField(max_length=200)
    #at this point the taxnode should be considered as true as filtering for wrong taxnodes is done in the form
    def if_valid_save_organism_name(self,user_email):
        try:
            scientific_name = get_scientific_name_by_taxid(user_email,self.taxonomic_node)
            self.organism_name = scientific_name
        except Exception:
            return False
        return True

    def __str__(self):
        return "forward taxonomic node: {} of project {}".format(self.organism_name, self.associated_project.id)


class TaxNodesForBackwardDatabase(models.Model):
    associated_project = models.ForeignKey(BlastProject,on_delete=models.CASCADE)
    taxonomic_node = models.IntegerField()
    organism_name = models.CharField(max_length=200)

    def if_valid_save_organism_name(self,user_email):
        try:
            scientific_name = get_scientific_name_by_taxid(user_email,self.taxonomic_node)
            self.organism_name = scientific_name
        except Exception:
            return False
        return True

    def __str__(self):
        return "backward taxonomic node: {} of project {}".format(self.organism_name, self.associated_project.id)

class QuerySequences(models.Model):
    associated_project = models.ForeignKey(BlastProject, on_delete=models.CASCADE)
    query_file_name = models.CharField(max_length=200, blank=False)
    path_to_query_file = models.CharField(max_length=300, blank=False)

    def __str__(self):
        return "query sequence: {} of project {}".format(self.query_file_name, self.associated_project.id)


class Genomes(models.Model):
    associated_project = models.ForeignKey(BlastProject, on_delete=models.CASCADE,blank=True, null=True)
    reciprocal_type = models.CharField(max_length=200, choices=[('forward', 'forward'), ('backward', 'backward')], blank=True, null=True)
    genome_name = models.CharField(max_length=200, blank=False)
    path_to_file = models.CharField(max_length=300, blank=False)

    def __str__(self):
        return "genome: {}".format(self.genome_name)


class ForwardBlastSettings(models.Model):
    associated_project = models.ForeignKey(BlastProject, on_delete=models.CASCADE)
    e_value = models.DecimalField(max_digits=30, decimal_places=15, default=0.0001)
    word_size = models.IntegerField(default=3)
    num_alignments = models.IntegerField(default=10000)

    def __str__(self):
        return "fw settings of {}".format(self.associated_project.__str__())


class BackwardBlastSettings(models.Model):
    associated_project = models.ForeignKey(BlastProject, on_delete=models.CASCADE)
    e_value = models.DecimalField(max_digits=30, decimal_places=15, default=0.0001)
    word_size = models.IntegerField(default=3)
    num_alignments = models.IntegerField(default=1)

    def __str__(self):
        return "bw settings of {}".format(self.associated_project.__str__())