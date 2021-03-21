from django.test import SimpleTestCase
import subprocess
from django.test import TestCase, Client
from os.path import isfile
from os import remove
from blast_project.models import RefSeqTransaction
from time import sleep
from blast_project.refseq_ftp_transactions import download_current_assembly_summary_into_specific_directory, delete_downloaded_assembly_summary, read_current_assembly_summary_with_pandas

class TestRefseqFtpTransactions(TestCase):

    def test_download(self):

        if(isfile('./blast_project/tests/assembly_summary_refseq.txt') == False):
            downloaded_file = download_current_assembly_summary_into_specific_directory('./blast_project/tests/')
            self.assertEquals(downloaded_file,'assembly_summary_refseq.txt')
        else:
            self.assertTrue(isfile('./blast_project/tests/assembly_summary_refseq.txt'))

        #self.assertTrue(delete_downloaded_assembly_summary('./blast_project/tests/'))

    # completeness_level = 'Chromosome', 'Scaffold', 'Complete Genome', 'Contig'
    def test_parsing_refseq_file(self):
        refseq_summary = read_current_assembly_summary_with_pandas('./blast_project/tests/assembly_summary_refseq.txt',['Complete Genome','Chromosome'])
        self.assertGreater(len(refseq_summary[0]['assembly_accession'].keys()), 1000)

    def test_get_process_id(self):
        pid = subprocess.Popen(['python','./blast_project/tests/static/long_process.py'])
        #TODO capture cmdline output in proc
        self.assertIsNotNone(pid.pid)
        #sleep(10.5)

