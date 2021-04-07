from django.test import TestCase, Client
from django.urls import reverse
from blast_project.models import BlastProject, Genomes, ForwardBlastSettings, BackwardBlastSettings, QuerySequences
from django.contrib.auth.models import User

class TestRefseqViews(TestCase):
    def setUp(self):
        self.user = User.objects.create(username='testuser', password='12345', is_active=True, is_staff=True,
                                        is_superuser=True)
        self.user.set_password('hello')
        self.user.save()
        self.client = Client()
        self.client.login(username='testuser', password='hello')
        self.download_refseq_genomes_url = reverse('refseq_genome_download')

    def test_download_refseq_genomes_GET(self):
        response = self.client.get(self.download_refseq_genomes_url)
        self.assertEqual(response.status_code,200)
        self.assertTemplateUsed(response,'blast_project/download_refseq_genomes.html')

    def test_download_refseq_genomes_POST(self):
        response = self.client.post( self.download_refseq_genomes_url,
                                    data={'refseq_levels':['Complete Genome'],'database_description':'test'})
        self.assertEqual(response.status_code,200)



