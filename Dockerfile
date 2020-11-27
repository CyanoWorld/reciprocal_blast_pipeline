#
# miniconda snakemake - NCBI BLAST+ 2.11.0+ Dockerfile
#

FROM ubuntu:focal

MAINTAINER Lukas Becker

# Download and install required software
RUN apt-get update -y && apt-get upgrade -y && apt-get install curl -y && apt-get install wget bzip2 -y
RUN mkdir /opt/blast
# set working directory
WORKDIR /blast

ENV PATH /blast/miniconda3/bin:$PATH

# Download and install anaconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#-b Batch mode with no PATH modifications to ~/.bashrc
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /blast/miniconda3
RUN rm Miniconda3-latest-Linux-x86_64.sh

# Download & install BLAST
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz | tar -zxvpf- 

# update environment variable
ENV PATH /blast/ncbi-blast-2.11.0+/bin:$PATH

#update miniconda
RUN conda update conda
RUN conda update --all
#add bioconda as channel
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

#installing snakemake
RUN conda install -c bioconda snakemake

# folder for volume binding
RUN mkdir /data

# Copying genomes and data for BLAST analysis
#COPY data .

# Optional commands e.g. initiating scripts 
CMD ["bash"]