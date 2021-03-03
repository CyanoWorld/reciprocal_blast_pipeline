from django.test import SimpleTestCase

from django.test import TestCase, Client
from os.path import isfile
from os import remove
from blast_project.refseq_ftp_transactions import download_current_assembly_summary_into_specific_directory, delete_downloaded_assembly_summary

class TestRefseqFtpTransactions(TestCase):
    def test_download(self):
        downloaded_file = download_current_assembly_summary_into_specific_directory('./blast_project/tests/')
        self.assertEquals(downloaded_file,'assembly_summary_refseq.txt')
        self.assertTrue(isfile('./blast_project/tests/assembly_summary_refseq.txt'))
        self.assertTrue(delete_downloaded_assembly_summary('./blast_project/tests/'))
