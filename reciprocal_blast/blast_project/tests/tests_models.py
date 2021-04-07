from blast_project.models import RefSeqTransaction, RefseqGenome, RefseqGenomeAssemblyLevels
from django.test import TestCase
class RefSeqDatabaseTestCase(TestCase):

    def setUp(self):
        RefSeqTransaction.objects.create(pid=2200,process_title="test case")
        assembly_1 = RefseqGenomeAssemblyLevels.objects.create(assembly_level='Complete Genome')
        assembly_2 = RefseqGenomeAssemblyLevels.objects.create(assembly_level='Chromosome')
        self.refseq_genome = RefseqGenome.objects.create(database_description="Test database complete genome chromosome",assembly_entries=2000)
        self.refseq_genome.assembly_levels.add(assembly_1)
        self.refseq_genome.assembly_levels.add(assembly_2)

    def test_access_refseqTransaction(self):
        transaction = RefSeqTransaction.objects.get(pid=2200)
        print(transaction.process_title)
        self.assertEqual("2200",transaction.return_pid())

    def test_refseqGenome(self):
        print("[****]")
        print("ID {}".format(self.refseq_genome.id))
        self.assertEqual(self.refseq_genome.database_description,'Test database complete genome chromosome')