from django.test import TestCase, Client
from django.urls import reverse
from blast_project.models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences
from django.contrib.auth.models import User

class TestViews(TestCase):
    def setUp(self):
        self.user = User.objects.create(username='testuser', password='12345', is_active=True, is_staff=True, is_superuser=True)
        self.user.set_password('hello')
        self.user.save()
        self.client = Client()
        self.client.login(username='testuser', password='hello')
        self.project = BlastProject.objects.create(project_title='test case',search_strategy='blastn',project_username=self.user)
        self.project.save()

        self.fw_genome = Genomes.objects.create(associated_project=self.project,reciprocal_type='forward',genome_name="forward_genome.faa",path_to_file="great_filepath/")
        self.bw_genome = Genomes.objects.create(associated_project=self.project,reciprocal_type='backward',genome_name="backward_genome.faa",path_to_file="great_filepath/")

        self.fw_genome.save()
        self.bw_genome.save()

        self.query_sequences = QuerySequences.objects.create(associated_project=self.project,query_file_name="best_query_file.faa",path_to_query_file="best_filepath/")
        self.query_sequences.save()

        self.forward_settings = ForwardBlastSettings.objects.create(associated_project=self.project)
        self.forward_settings.save()

        self.backward_settings = BackwardBlastSettings.objects.create(associated_project=self.project)
        self.backward_settings.save()


        self.main_url = reverse('main')
        self.download_refseq_genomes_url = reverse('refseq_genome_download')
        self.project_details_url = reverse('project_details',args='2')
        #


    def test_main_GET(self):
        #self.client.login(username='User1', password='checkthisgreat546!')
        response = self.client.get(self.main_url)
        self.assertEquals(response.status_code, 200)
        self.assertTemplateUsed(response,'blast_project/main.html')

    def test_project_details_GET(self):
        #self.client.login(username='User1', password='checkthisgreat546!')
        response = self.client.get(self.project_details_url)
        self.assertEquals(response.status_code, 200)
        self.assertTemplateUsed(response,'blast_project/project_details.html')


