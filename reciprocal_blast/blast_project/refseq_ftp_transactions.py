from gzip import decompress
from os import remove, chdir, getcwd
from os.path import isfile
from wget import download
import pandas as pd

#downloads the current refseq assembly file into an specified directory
def download_current_assembly_summary_into_specific_directory(directory):
    refseq_url = "ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    try:
        current_working_directory = getcwd()
        chdir(directory)
        downloaded_file = download(refseq_url)
        if isfile(downloaded_file) == False:
            raise FileExistsError
        chdir(current_working_directory)
    except Exception as e:
        raise FileExistsError("[-] Couldn't download assembly_summary_refseq.txt file ...")

    return downloaded_file


def read_current_assembly_summary_with_pandas(summary_file):
    pass