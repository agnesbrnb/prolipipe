# Prolific

## Description

The prolific pipeline is about genome annotation and metabolic pathway reconstruction. From .fasta genomes, it enables to find if strains of bacteria are theorically able to produce metabolites.

## Installation

Because of the amount of calculations, we advise to create a genouest account and install all the requirements on your cluster account.  
We advise you to create a conda environment with bioperl 5.22 (not above) to be able to use prokka.  
You may have to install python and the followng packages : os, optparse, supervenn, matplotlib.pyplot, pandas, argparse, numpy.

#### Prokka
You can install prokka many ways, as described in the [prokka github page](https://github.com/tseemann/prokka)

#### Pathway tools and mpwt
To install and use pathway tools and mpwt with singularity, you can follow the tutorial of Metage2Metabo [here](https://metage2metabo.readthedocs.io/en/latest/install.html#installation-with-singularity-e-g-on-a-cluster)


## Usage

To begin, you can download genomes using [ncbi genome download](https://github.com/kblin/ncbi-genome-download) or use your own ones.  
The genomes must be ordonned like the genome folder in toy_example.  

  


#### rename.py

Use rename.py to rename practically your genome folders and .fasta files :  
`rename.py -d [path to the directory of the folders to be renamed]`

#### pipeline_yael.py

`pipeline_yael.py --option [option]`
	
options:  

	--file		Directory of the data read
	--padmet_ref	Path to the padmet_ref needed for the module mpwt
	--output	Path to the folder where you want to put the results in
	--ptsc		Path to the scratch folder (on the cluster)
	--ptsi		Name of the sigularity image
	-e 		Launch extract module to decompress all .tar and prepare data for prokka module
	-t		Module to help you construct the taxon_id file that is need for Pathwaytools
	-p		Launch the prokka module that will create all .gbk with the .fasta for each folder organism
	-m		Launch mpwt module that will run the pathwaytools singularity and do multiple things like : create database for each organism, create a padmet file and wikipages with the padmet toolbox
	-q		quiet mode
	-k		Keep .faa files that can be need to use other annotation software like eggNOG-mapper
	-v		Activate verbose  

Before using the option -t, please make sure you have completed the all_taxons.tsv file with the information of your own genomes.

##### Examples

###### Launch prokka :  
`pipeline_yael.py -p --file /toy_example/genomes/ `

###### Create the id_taxon.tsv file for mpwt :  
`pipeline_yael.py -t --file /toy_example/genomes/ `

###### Launching mpwt :  
`pipeline_yael.py -m --output [name for output file] --ptsc [Path to the scratch folder] --ptsi [Name of the sigularity image] --padmet_ref [path to padmet ref]`


#### pathways.py

Before launching the script, make sure you have created a .txt file which contains a list of metacyc reaction corresponding to a pathway, as in the pathway_pyruvate.txt file.

`pathways.py -r [path and name of reaction file] -p [path and name of pathway file : reaction list] -o [path and name for the output .csv file]`

optional argument that creates a graph :

	-g		path and name of the output .png graph

