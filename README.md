# reciprocal_blast_pipeline

# Attention
- currently the input files for the genomes and forward_input.faa are too big 

## Things to do next week: TODO
- write a backward blast result processing script --> integrate available scripts
- integrate BioPython or another software package to process the accession IDs in order to get taxon specifications --> integrate available scripts
- test available pipeline with less time

## snakemake, miniconda and BLAST in a docker container
Build the docker image with:
`docker build --tag snakeblast:1.0` this will take a while. The image has a size of 2.38 GB. 

Start a container and bind a volume to it:
`docker run -v path_on_your_local_machine/desired_directory_for_volume_sharing:/blast/data -i -d --name snakeblast snakeblast:1.0`

Start a shell inside the container and start working with snakemake and BLAST:
`docker exec -it snakeblast snakeblast:1.0 /bin/bash`

Execute snakemake where X are your desired cores to use:
`cd /data && snakemake --cores X `

This will probably take some time. Once snakemake has finished the "all" job, there will be two .csv output files in the /data directory, one is called backward_blast_results.csv which is currently the result file for this first approach towards a full reciprocal pipeline. 

Stop the container:
`docker stop snakeblast`

Restart the container:
`docker start snakeblast`