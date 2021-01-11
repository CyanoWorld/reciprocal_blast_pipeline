# reciprocal_blast_pipeline

## Things to do next week: TODO
- refactoring of the webinterface
- big updates in relevant branch
- webinterface in relevant branch
- integrate the whole nr database into the webinterface - tried to use the update_blastdb.pl perl script but unfortunatly this was no solution
- refactoring of the genomes model: uploaded genomes could be used multiple times! 

## snakemake, miniconda, BLAST, django and postgres in a docker container
Edit the docker-compose.yml file with right file paths to volumes. Start the network container with:
`docker-compose up`. Keep in mind, if you delete the containers with `docker-compose down` all data in the volume folders will be lost.

## snakemake, miniconda and BLAST in a docker container
Build the docker image with:
`docker build --tag snakeblast:1.0 .` this will take a while. The image has a size of 2.38 GB. 

Start a container and bind a volume to it:
`docker run -v path_on_your_local_machine/desired_directory_for_volume_sharing:/blast/data -i -d --name snakeblast snakeblast:1.0`

Start a shell inside the container and start working with snakemake and BLAST:
`docker exec -it snakeblast /bin/bash`

Execute snakemake where X are your desired cores to use:
`cd /data && snakemake --cores X `

This will probably take some time. Once snakemake has finished the "all" job, there will be two .csv output files in the /data directory, one is called backward_blast_results.csv which is currently the result file for this first approach towards a full reciprocal pipeline. 

Stop the container:
`docker stop snakeblast`

Restart the container:
`docker start snakeblast`