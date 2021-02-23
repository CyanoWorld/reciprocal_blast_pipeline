'''
For further information see Django developer guidelines
'''
from django.contrib import admin
from .models import BlastProject, Genomes
# Register your models here.
admin.site.register(BlastProject)
admin.site.register(Genomes)