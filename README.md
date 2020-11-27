# reciprocal_blast_pipeline

## snakemake, miniconda and BLAST in a docker container
Build the docker image with:
`docker build --tag snakeblast:1.0` this will take a while. The image has a size of 2.38 GB. 

Start a container and bind a volume to it:
`docker run -v path_on_your_local_machine/desired_directory_for_volume_sharing:/blast/data -i -d --name snakeblast snakeblast:1.0`

Start a shell inside the container and start working with snakemake and BLAST:
`docker exec -it snakeblast snakeblast:1.0`

Stop the container:
`docker stop snakeblast`

Re-Start the container:
`docker start snakeblast`