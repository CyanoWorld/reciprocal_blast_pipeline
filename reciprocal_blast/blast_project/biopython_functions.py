from Bio import Entrez
from django.db import IntegrityError


def get_species_taxid(user_email,scientific_name):
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList'][0]
        translation = record['QueryTranslation']
    except Exception as e:
        raise IntegrityError("[-] There is no taxonomic node defined by your specified scientific name: {} Exception: {}".format(scientific_name,e))
    return taxid, translation

#try to fetch the user_email in the form
def get_species_taxid_without_email(scientific_name):
    try:
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList'][0]
    except Exception as e:
        raise IntegrityError("[-] There is no taxonomic node defined by your specified scientific name: {} Exception: {}".format(scientific_name,e))
    return taxid

def get_scientific_name_by_taxid(user_email,taxid):
    try:
        Entrez.email = user_email
        handle = Entrez.efetch('taxonomy', id=taxid, rettype='xml')
        record = Entrez.read(handle)
        scientific_name = record[0]['ScientificName']
    except Exception as e:
        raise IntegrityError("[-] There is no organism with your specified taxonomic node: {} Exception: {}".format(id,e))
    return scientific_name

