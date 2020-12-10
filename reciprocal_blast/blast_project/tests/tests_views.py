from django.test import TestCase, Client
from django.urls import reverse
from blast_project.models import BlastProject, Genomes

class TestViews(TestCase):

    def setUp(self):
        self.client = Client()
        self.main_url = reverse('main')
        self.project_details_url = reverse('project_details',args='2')
        BlastProject.objects.create(project_title='test case',search_strategy='blastn')

    def test_main_GET(self):
        response = self.client.get(self.main_url)
        self.assertEquals(response.status_code, 200)
        self.assertTemplateUsed(response,'blast_project/main.html')

    def test_project_details_GET(self):
        response = self.client.get(self.project_details_url)
        self.assertEquals(response.status_code, 200)
        self.assertTemplateUsed(response,'blast_project/project_details.html')
