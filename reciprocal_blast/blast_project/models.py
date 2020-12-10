from django.db import models
from django.contrib.auth.models import User


class BlastProject(models.Model):
    project_title = models.CharField(max_length=200, blank=False, unique=True)
    search_strategy = models.CharField(max_length=20, choices=[('blastp', 'blastp'), ('blastn', 'blastn')], default='blastp')
    project_username = models.ForeignKey(User,on_delete=models.CASCADE)

    def __str__(self):
        return "{}".format(self.project_title)


class QuerySequences(models.Model):
    associated_project = models.ForeignKey(BlastProject, on_delete=models.CASCADE)
    query_file_name = models.CharField(max_length=200, blank=False)
    path_to_query_file = models.CharField(max_length=300, blank=False)

    def __str__(self):
        return "{}".format(self.query_file_name)


class Genomes(models.Model):
    associated_project = models.ForeignKey(BlastProject, on_delete=models.CASCADE)
    reciprocal_type = models.CharField(max_length=200, choices=[('forward', 'forward'), ('backward', 'backward')], blank=False)
    genome_name = models.CharField(max_length=200, blank=False)
    path_to_file = models.CharField(max_length=300, blank=False)

    def __str__(self):
        return "{}".format(self.genome_name)


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