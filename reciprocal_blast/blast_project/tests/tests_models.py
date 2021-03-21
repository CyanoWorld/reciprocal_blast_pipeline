from blast_project.models import RefSeqTransaction
from django.test import TestCase
class RefSeqDatabaseTestCase(TestCase):
    def setUp(self):
        RefSeqTransaction.objects.create(pid=2200,process_title="test case")

    def test_access_refseqTransaction(self):
        transaction = RefSeqTransaction.objects.get(pid=2200)
        print(transaction.process_title)
        self.assertEqual("2200",transaction.return_pid())