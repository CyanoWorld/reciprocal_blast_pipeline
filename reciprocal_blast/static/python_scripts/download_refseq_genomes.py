#!/usr/bin/python

import sys, getopt
import subprocess
import pandas as pd
import wget
import os
import gzip
from os.path import isfile, isdir
import time
#re.compile and re.findall
import re
#different package for downloading files from ncbi's ftp server
#import signal
from Bio import Entrez

#downloading genomes with wget 
def download_genome_from_ftp_path(ftp_path):
    try:
        return wget.download(ftp_path)
    except:
        return None

    
def makeblastdatabase(database_description):
    print("[*] EXECUTING makeblastdb")
    blast_db_cmd = subprocess.Popen(['makeblastdb','-in',database_description+'.database.faa','-dbtype','prot','-taxid_map','acc_taxid_map.table','-parse_seqids'],stdout=open('makeblastdb.log','w'))
    os.waitpid(blast_db_cmd.pid,0)
    print("[*] END DATABASE CREATION, finished makeblastdb, logs are written into makeblastdb.log")
    
def concatenate_assemblies(database_description):
    print("[*] STARTING CONCATENATION")
    fasta_files = [file for file in os.listdir() if 'decompressed.faa' in file]
    database_file = open(database_description,'w')
    errorlist=[]
    for ffile in fasta_files:
        try:
            cat_cmd = subprocess.Popen(['cat',str(ffile)],stdout=database_file)
            os.waitpid(cat_cmd.pid,0)
            os.remove(str(ffile))
        except Exception as e:
            print("[-] Exception during concatenation of downloaded files: {}".format(e))
            errorlist.append(e)
    database_file.close()
    print("[*] ENDED CONCATENATION")
    return errorlist

#function with two responsibilities in order to save memory usage
def extract_downloaded_file_and_write_taxid_file(genome_downloaded_file,taxid):
    
    #decompression to fasta file
    try:
        output = open(genome_downloaded_file+".decompressed.faa","w")
        with gzip.open(genome_downloaded_file,"rb") as bytes_out:
            bytes_from_file = bytes_out.read()
            line = bytes_from_file.decode("utf-8")
            output.write(line)
            output.close()
        #remove gz file
        os.remove(genome_downloaded_file)
    except Exception as e:
        output.close()
        raise Exception("[-] Exception during decompressing: {}".format(e))
    
    #creation of taxmap file | acc_id /t taxid
    try:
        accession_id_pattern = re.compile('>(\S*)')
        taxmap = open('acc_taxid_map.table','a')
        for acc_id in re.findall(accession_id_pattern,line):
            taxmap.write(str(acc_id)+"\t"+str(taxid)+"\n")
        taxmap.close()
        #print("\t[*] Done writing taxonomic informations into taxmap file ...")
        return True
    
    except Exception as e:
        taxmap.close()
        raise Exception("[-] Exception during writing taxmap file: {}".format(e))
        
def download_assemblies(filtered_table,errorlist):
    
    try:
        print("[*] STARTING DOWNLOAD")
        genomes_number = len(filtered_table)
        percentage_per_run = 100.0/genomes_number
        percentage = 0
        for genome_url,taxid in zip(filtered_table['ftp_path'],filtered_table['species_taxid']):
            genome_file = download_genome_from_ftp_path(genome_url)
            time.sleep(1)
            #genome_file = wget.download(genome_url)
            percentage += percentage_per_run
            print(percentage)
            if(genome_file):
                extract_downloaded_file_and_write_taxid_file(genome_file,taxid)
            else:
                errorlist.append(genome_url)
                print("[ERROR] {}".format(genome_url))
                
        print("[*] ENDED DOWNLOAD")
        return errorlist
    except Exception as e:
        raise Exception("[-] Error during downloading assemblies with Exception: {}".format(e))
    
    

def read_refseq_table(refseq_table):
    try:
        return pd.read_csv(refseq_table,index_col=0)
    except Exception as e:
        raise Exception("[-] Couldnt read pandas table: {} with Exception: {}".format(refseq_table,e))

def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hio:",["ifile=","outputdir="])
    except getopt.GetoptError:
        print('test.py -i <inputfile>')
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfilepath>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o","--outputdir"):
            outputdir = arg
        else:
            raise getopt.GetoptError("[-] ERROR: Please specify the input file with --ifile and the outputdirectory with --outputdir")
    
    maindir = os.getcwd()
    print(str(maindir))
    print(os.listdir())
    print(isdir(outputdir))
    print(isfile(inputfile))
    os.chdir(outputdir)
    #process_pid = os.getpid()
    #process = psutil.Process(process_pid)
    try:
        errorlist=[]
        refseq_table=read_refseq_table(inputfile.split('/')[-1])
        errorlist=download_assemblies(refseq_table,errorlist)
     
        second_error_list=concatenate_assemblies(inputfile.split('/')[-1]+".database.faa")
        for error in second_error_list:
            errorlist.append(error)
            
        if(len(errorlist)!=0):
            print("[-] ERRORS OCCURED: take a look at the errors.txt file")
            errorfile = open('errors.txt','w')
            for i in errorlist:
                errorfile.write(str(i)+"\n")
            errorfile.close()
        makeblastdatabase(inputfile.split('/')[-1])
        os.chdir(maindir)
        sys.exit(0)
    except Exception as e:
        os.chdir(maindir)
        print("[-] ERROR: Exception: {}".format(e))
            
        sys.exit(-1)
    #return 0

if __name__ == "__main__":
    main(sys.argv[1:])