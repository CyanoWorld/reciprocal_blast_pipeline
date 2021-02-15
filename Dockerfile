#AUTHOR Lukas Becker
# miniconda - django - snakemake - NCBI BLAST+ 2.11.0+ - Dockerfile
#

#base image; maybe choose another image
FROM ubuntu:focal

# Download and install required software
RUN apt-get update -y && apt-get upgrade -y && apt-get install curl -y && apt-get install wget bzip2 -y
RUN mkdir /opt/blast
# set working directory
WORKDIR /blast

# Update path environment variable
ENV PATH /blast/miniconda3/bin:$PATH

# Download and install anaconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#-b Batch mode with no PATH modifications to ~/.bashrc
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /blast/miniconda3
RUN rm Miniconda3-latest-Linux-x86_64.sh

# Download & install BLAST
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz | tar -zxvpf- 

# Update path environment variable
ENV PATH /blast/ncbi-blast-2.11.0+/bin:$PATH

# Update miniconda
RUN conda update conda
RUN conda update --all

# Add bioconda as channel
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Download and install snakemake via the bioconda channel
RUN conda install -c bioconda snakemake

# Set-up django and postgres packages
COPY requirements.txt /blast/
RUN pip install -r requirements.txt 
COPY /reciprocal_blast /blast/

#create default BLASTDB environment variable
ENV BLASTDB /blast/media/databases
# Delete not required packages etc..
RUN apt-get autoremove --purge --yes && apt-get clean && rm -rf /var/lib/apt/lists/*
# Optional commands e.g. initiating scripts
CMD ["bash"]
