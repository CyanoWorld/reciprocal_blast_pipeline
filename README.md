# reciprocal_blast_pipeline
## Installation
In order to start developing and using this application you need to download the application from the git-repository. There are two possibilities for installing the application either you install the package distribution platform anaconda or the container virtualization platform docker. Installation via docker involves alteration and editing of the docker-compose.yaml file with subsequent execution of `docker-compose up`. Particularly one would have to change the path of the volumes into: 

`path_to_project_directory/data:/blast/media` and `path_to_project_directory/tmp:/blast/tmp`. Once this is done you can start the container-network with the command: `docker-compose up`. This will download and install the desired images and software requirements defined in the Dockerfile and requirements.txt file. Additionally, the application gets copied into the docker container.
For installing the application within a conda environment you need the anaconda distribution platform and the biosnakedjango-\{system\}.yml file for your specific system. You can install the required packages in a virtual environment by starting a shell and submitting the command: `conda env create -f biosnakedjango-{system}.yml`
In order to activate the environment use the command: `conda activate biosnakedjango`. If you want to use the PostgreSQL database you need the additional python package psycopg2 and an active database server on your local machine, that runs on the port 5432. Otherwise, you can change the database to SQLite3 by editing the DATABASE section in the settings.py file, or by replacing everything in the DATABASE section with:

```
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}
```


Afterwards, change your active working directory to the reciprocal_blast project root folder and type following commands: `python manage.py makemigration` and `python manage.py migrate`. These commands will directly create the needed database tables. In the same working directory you can start the server with the command `python manage.py runserver 8080`. Snakemake is executed with the --wms-monitor http://127.0.0.1:5000 flag which tells snakemake to communicate with the monitoring service panoptes (https://github.com/panoptes-organization). In order to allow the correct snakemake execution and communication with panoptes, open a new shell in the project directory, activate the biosnakedjango environment and start the panoptes server by submitting the command `panoptes`. The last installation step for the conda setup involves downloading and installing the BLAST command line application and databases, as well as the E-Direct tool. For installation of those executables follow the instructions that are given in the BLAST and E-Direct user manuals. If you want to develop on windows, you cannot install the E-Direct tool. You need to out-comment the function

## Developer Guidelines

The recommended installation for development is to use the conda environment files. If you have installed the application via docker, make sure to share the reciprocal_blast directory with the container. Otherwise you would need to copy all changed files into the container. If you share the application volume with the web container you might run into errors with the command argument. Thus, you can delete the argument and submit the commands for model migrations and server start directly in the container. In order to do this you need to start the container network with the command docker-compose up. Open a new shell and type `docker exec -it reciprocal_blast_pipeline_web_1 /bin/bash`.  
Load the project folder reciprocal_blast into your development IDE. One potential IDE that has plugins for snakemake support is PyCharm. The next step involves changes in the database section of the reciprocal_blast/reciprocal_blast/settings.py file. Replace the current postgres database settings with the settings for the SQLite database. Those settings are out commented in the relevant part of the settings file.

Keep in mind that django is a MVT framework, which has its own markup language for designing templates and an own class based database monitoring. In addition to models and templates the core hook for execution of background tasks such as the execution of snakemake is the view.py file. Functions in this file will process the individual url specific HTTP GET and POST methods. This aspect is dealt with in more detail in Section DJANGO and on the [django project website](https://www.djangoproject.com/)). 

HTML template files are rendered with Bootstrap 4, which is an HTML, CSS and JavaScript framework primarily for designing websites. A good overview of all design elements and options Bootstrap offers can be found on the [Bootstrap project page](https://getbootstrap.com/docs/5.0/getting-started/introduction/).

## User Guidelines

The recommended installation procedure is the Docker setup. Especially windows users should install the application with docker, thus automated deletion of project files will not work correctly on windows platforms and installation of NCBI E-Direct tools is only possible on Mac or Linux operating systems. In order to enable the use of the preformatted nr database you have to download and decompress the database files from the [NCBI FTP server ](ftp.ncbi.nlm.nih.gov/blast/db/). Additionally you have to download the taxonomy database (taxdb) in order to get taxonomic information in your blast output files. You can either download the files manually or use the update_blastdb.pl perl script that is provided by the BLAST executables. The size of the database is, after decompressing, roughly around 253 Gigabyte. Store the database in an appropriate directory and create a path variable called BLASTDB that points to the directory. If you choose the docker-compose.yaml file for the installation procedure you have to add two things to the container specifications: the database folder as a volume and an environment variable that points to this volume. This can be done in the docker-compose.yaml file with the volume and environment parameters. The application starts on port 8080, the postgres database starts with port 5432. 

### Web interface operations

Once you have started the application you can visit the browser interface by typing and submitting localhost:8080/blast_project/ into your browser search bar. You need to register and create an account in order to use the application. After successful login you get redirected to the homepage of this application, which serves as a dashboard for your created projects. If you press the navigation bar button you can view a bunch of subsites, e.g. the project creation or the species taxonomic node checkup. Utilization of most subpages are self explaining. 


| URL           | Supported HTTP methods| Name in urls.py | Short Description  |
| ------------- |:------:| --------------:| -----------------------------------|
| localhost:8080/blast_project| GET| main |Homepage and project dashboard, gives an overview about current projects and serves as an anchor point for reaching the project specific subsites.|
| /register      | GET, POST     |   register | Account registration.|
| /login | GET, POST      |    login | Login for application access.|
|/logout|GET|logout|Logut.|



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
