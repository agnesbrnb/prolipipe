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
To install and use pathway tools and mpwt, you can follow the tutorial of Metage2Metabo [here](https://metage2metabo.readthedocs.io/en/latest/install.html#installation-with-singularity-e-g-on-a-cluster)


## Usage

To begin, you can download genomes using [ncbi genome download](https://github.com/kblin/ncbi-genome-download) or use your own ones.  
The genomes must be ordonned like the genome folder in toy_example.  

You also will have to complete the all_taxons.tsv file with the information of your own genomes.  


#### rename.py

Use rename.py to rename practically your genome folders and .fasta files :  
rename.py -d [path to the directory of the folders to be renamed]

#### pipeline_yael.py

pipeline_yael.py --option [option]
	
options:  
-- file		:  Directory of the data read  
-- padmet_ref	 : Path to the padmet_ref needed for the module mpwt  
-- output	 : Path to the folder where you want to put the results in  
-- ptsc		:  Path to scratch folder (link to the singularity)  
-- ptsi		:  Path to the singularity  
-- id_wiki	:  The wiki id that will be use as reference for all the wikipages that are created by the padmet toolbox inside the mpwt module  
	
	-v		Activate verbose
	-e 		Launch extract module to decompress all .tar and prepare data for prokka module
	-t		Module to help you construct the taxon_id file that is need for Pathwaytools
	-q		quiet mode
	-k		Keep .faa files that can be need to use other annotation software like eggNOG-mapper
	-p		Launch the prokka module that will create all .gbk with the .fasta for each folder organism
	-m		Launch mpwt module that will run the pathwaytools singularity and do multiple things like : create database for each organism, create a padmet file and wikipages with the padmet toolbox	

